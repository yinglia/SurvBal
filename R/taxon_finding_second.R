#' STEP 2: select more taxon one by one
#'
#' Search for next taxon that is the most significant and see whether to add it to numerator or denominator
#'
#' @param Surv_obj the object containing information about the survival outcome
#' @param data the pre-processed data (produced from STEP 0)
#' @param covariates the covariates for adjustment
#' @param candidates the pool of remaining taxa 
#' @param tracker the index tracking how many taxa have been added to form the balance 
#' @param selected the taxa that have been selected so far and associated statistics 
#' @param model the model that the user wants to use 
#' @param dist if \code{model} is ``\code{parametric}'', than dist is the distribution the user wants to use, ``\code{weibull}'' is default
#' @param sequential_test whether sequential testing is used to stop the forward search
#' @param sequential_alpha the level of significance for the sequential testing
#'
#' @noRd

taxon_finding_second <- function(Surv_obj, data, covariates, candidates, tracker, selected, model, dist, sequential_test, sequential_alpha){

  len_cand = length(candidates)

  record_positive = matrix(ncol=2, nrow=len_cand) # keeps track of coefficient and p value for candidate with positive balance
  record_negative = matrix(ncol=2, nrow=len_cand) # keeps track of coefficient and p value for candidate with negative balance
  selected_positive = selected[, 1][which( selected[,'sign']== 1 )] # the selected candidates with a positive coefficient
  selected_negative = selected[, 1][which( selected[,'sign']== 0 )] # the selected candidates with a negative coefficient

  mean_positive = apply( data[, selected_positive, drop=F], 1, mean)
  mean_negative = apply( data[, selected_negative, drop=F], 1, mean)

  for (ii in 1:len_cand) {
    cur_cand = candidates[ii]

    # B_positive is given by the mean of the selected positive with current candidate
    # minus the mean of the selected negative
    B_positive = apply(data[, c(selected_positive, cur_cand), drop=F], 1, mean) - mean_negative
    if (is.null( covariates )){
      if (model == "coxph") {
        mod = coxph(Surv_obj ~ B_positive)  # get coefficient and p value for current B_positive
        record_positive[ii,] = summary(mod)$coef[, c(1, 5)] # records the coefficient and p value with positive balance
      } else if (model == "parametric"){
        if (is.na(dist)) {
          dist = "weibull"
        }
        mod = survreg(Surv_obj ~ B_positive, dist=dist)
        record_positive[ii, ] = summary(mod)$table[2, c(1, 4)] # records the coefficient and p value for current candidate
      }
    } else{
      d_step2 = data.frame(Surv_obj, B_positive, covariates)
      
      if (model == "coxph") {
        mod = coxph(Surv_obj ~ ., data=d_step2)  # get coefficient and p value for current B_positive
        record_positive[ii,] = summary(mod)$coef[1, c(1, 5)] # records the coefficient and p value with positive balance
      } else if (model == "parametric"){
        if (is.na(dist)) {
          dist = "weibull"
        }
        mod = survreg(Surv_obj ~ ., dist=dist, data=d_step2)
        record_positive[ii, ] = summary(mod)$table[2, c(1, 4)] # records the coefficient and p value for current candidate
      }
    }
    
    # B_negative is given by mean of the selected positive minus the mean of the selected negative
    # with the current candidate
    B_negative = mean_positive - apply(data[, c(selected_negative, cur_cand), drop=F], 1, mean)
    if (is.null( covariates )){
      if (model == "coxph") {
        mod = coxph(Surv_obj ~ B_negative)  # get coefficient and p value for current B_negative
        record_negative[ii, ] = summary(mod)$coef[, c(1, 5)]
      } else if (model == "parametric"){
        if (is.na(dist)) {
          dist = "weibull"
        }
        mod = survreg(Surv_obj ~ B_negative, dist=dist)
        record_negative[ii, ] = summary(mod)$table[2, c(1, 4)] # records the coefficient and p value for current candidate
      }
    } else{
      d_step2 = data.frame(Surv_obj, B_negative, covariates) 
      
      if (model == "coxph") {
        mod = coxph(Surv_obj ~., data=d_step2)  # get coefficient and p value for current B_negative
        record_negative[ii, ] = summary(mod)$coef[1, c(1, 5)]
      } else if (model == "parametric"){
        if (is.na(dist)) {
          dist = "weibull"
        }
        mod = survreg(Surv_obj ~., dist=dist, data=d_step2)
        record_negative[ii, ] = summary(mod)$table[2, c(1, 4)] # records the coefficient and p value for current candidate
      }
    }
  }

  # arrange the result: taxon_id, sign (num or denom), p_value
  ind_positive_selected = which.min(record_positive[, 2])

  result_positive = c(candidates[ind_positive_selected], # taxon_id
                      1,
                      p_value = record_positive[ind_positive_selected, 2]) # p_value

  # same as above but for negative
  ind_negative_selected = which.min(record_negative[, 2])

  result_negative = c(candidates[ind_negative_selected],
                      0,
                      p_value = record_negative[ind_negative_selected, 2])

  # selects the result with the smaller p value
  if (result_negative[3] < result_positive[3]) {
    result = result_negative
  } else {
    result = result_positive
  }
  
  # update the candidate pool and tracker
  candidates = setdiff(candidates, result[1]) # removes selected candidate from pool
  tracker = tracker + 1

  sig = NA
  # sequential test if it is chosen
  if (sequential_test == TRUE){
    
    B_original = mean_positive - mean_negative # original balance
    if (result[2] == 1){
      B_updated = apply(data[, c(selected_positive, result[1]), drop=F], 1, mean) - mean_negative
    } else if (result[2] == 0){
      B_updated = mean_positive - apply(data[, c(selected_negative, result[1]), drop=F], 1, mean)
    } # updated balance
  
    if (is.null(  covariates) ){
      if (model == "coxph") {
        mod_original = coxph(Surv_obj ~ B_original, x=T) # original cox model
        mod_updated = coxph(Surv_obj ~ B_updated, x=T) # updated cox model
        test = plrtest(mod_updated, mod_original, nested = FALSE, adjusted = "none") # could also adjust for model complexity, AIC/BIC, no need here as balance is one variable
        sig = ( test$pLRTAB < sequential_alpha ) # non-nested LRT, model 1 fits better than model 2 (pLRTA), Model 2 fits better than Model 1 (pLRTB), Model fits not equally close to true Model  
        # variance test, whether distinguishable (pOmega)
      } else if (model == "parametric"){
        if (is.na(dist)) {
          dist = "weibull"
        }
        data_toboot = data.frame(Surv_obj, B_original, B_updated)
        reps = boot(data=data_toboot, dist_boot=dist, statistic=AICdiff_function, R=1000)
        ci = boot.ci(reps, type="bca", conf=1-sequential_alpha)
        sig = ( ci$bca[4] * ci$bca[5] > 0 )
      }
    } else{
      if (model == "coxph") {
        d_original = data.frame(Surv_obj, B_original, covariates)
        mod_original = coxph(Surv_obj ~., x=T, data=d_original) # original cox model
        
        d_updated = data.frame(Surv_obj, B_updated, covariates)
        mod_updated = coxph(Surv_obj ~., x=T, data=d_updated) # updated cox model
        
        test = plrtest(mod_updated, mod_original, nested = FALSE, adjusted = "none") # could also adjust for model complexity, AIC/BIC, no need here as balance is one variable
        sig = ( test$pLRTAB < sequential_alpha ) # non-nested LRT, model 1 fits better than model 2 (pLRTA), Model 2 fits better than Model 1 (pLRTB), Model fits not equally close to true Model  
        # variance test, whether distinguishable (pOmega)
      } else if (model == "parametric"){
        if (is.na(dist)) {
          dist = "weibull"
        }
        data_toboot = data.frame(Surv_obj, B_original, B_updated, covariates)
        reps = boot(data=data_toboot, dist_boot=dist, statistic=AICdiff_function_adjusted, R=1000)
        ci = boot.ci(reps, type="bca", conf=1-sequential_alpha)
        sig = ( ci$bca[4] * ci$bca[5] > 0 )
      }
    }
    
  }
  
  return (list(result=result, candidates=candidates, tracker=tracker, sig=sig)) # return the result in the end, use list for multiple items

}

