#' STEP 1: initial ratio
#'
#' Find the first ratio of two taxon (difference of the log of the two taxon) to put in the model by
#' looping through the bacteria to find the bacterium that is the most statistically significant to survival
#'
#' @param Surv_obj the object containing information about the survival outcome
#' @param data the pre-processed data (produced from STEP 0)
#' @param covariates the covariates for adjustment
#' @param candidates the pool of remaining taxa 
#' @param tracker the index tracking how many taxa have been added to form the balance 
#' @param model the model that the user wants to use 
#' @param dist if \code{model} is ``\code{parametric}'', than dist is the distribution the user wants to use, ``\code{weibull}'' is default
#'
#' @import utils
#' 
#' @noRd

taxon_finding_first <- function(Surv_obj, data, covariates, candidates, tracker, model, dist) {

  len_cand = length(candidates)
  record = t(combn(candidates, 2)) # combinations

  record = cbind(record, t( apply(record, 1, function(z){
    B_ij = data[, z[1]] - data[, z[2]]

    if (is.null( covariates )){
      if (model == "coxph") {
        mod1 = coxph(Surv_obj ~ B_ij) # fit coxph for B_ij
        summary(mod1)$coef[, c(1, 5)] # records the coefficient and p value for current combination
      } else if (model == "parametric"){
        if (is.na(dist)) {
          dist = "weibull"
        }
        mod1 = survreg(Surv_obj ~ B_ij, dist=dist) # fit parametric for B_ij
        summary(mod1)$table[2, c(1, 4)] # records the coefficient and p value for current combination
      } else {
        stop("The only options of model specification are 'coxph' and 'parametric'.")
      }
    } else{
      d_step1 = data.frame(Surv_obj, B_ij, covariates)
      
      if (model == "coxph") {
        mod1 = coxph(Surv_obj ~., data=d_step1) # fit coxph for B_ij
        summary(mod1)$coef[1, c(1, 5)] # records the coefficient and p value for current combination
      } else if (model == "parametric"){
        if (is.na(dist)) {
          dist = "weibull"
        }
        mod1 = survreg(Surv_obj ~ ., dist=dist, data=d_step1) # fit parametric for B_ij
        summary(mod1)$table[2, c(1, 4)] # records the coefficient and p value for current combination
      } else {
        stop("The only options of model specification are 'coxph' and 'parametric'.")
      }
    }
  }) ) )

  # stores the result with three components
  result = matrix(ncol=3, nrow=2)
  ind_selected = which.min(record[, 4])

  result[1,1] = record[ind_selected, 1] # taxon_id i
  result[2,1] = record[ind_selected, 2] # taxon_id j
  result[1,3] = result[2,3] = record[ind_selected, 4] # p value
  if (record[ind_selected, 3] >= 0){ # the sign (num or denom)
    result[1,2] = 1; result[2,2] = 0
  } else{
    result[1,2] = 0; result[2,2] = 1
  }

  # update the candidate pool and tracker
  candidates = setdiff(candidates, result[, 1])
  tracker = tracker + 2

  return (list(result=result, candidates=candidates, tracker=tracker)) # return the result in the end

}
