#' STEP 2, help function
#'
#' Sequential test for parametric model using boostrap, the criterion is AIC 
#' (no need to adjust complexity since balance is one variable)
#' 
#' @param data the data to boot, containing Surv_obj, B_original, and B_updated
#' @param dist_boot the parametric distribution
#' @param indices placeholder, allowing boot to select sample
#'
#' @noRd

AICdiff_function <- function(data, dist_boot, indices) {
  
  d <- data[indices, ] #allows boot to select sample
  
  fit1 = survreg(Surv_obj ~ B_original, dist=dist_boot, data=d)
  fit2 = survreg(Surv_obj ~ B_updated, dist=dist_boot, data=d)
  
  return(AIC(fit2) - AIC(fit1))
  
}

#' Sequential test for parametric model using boostrap, when covariates are adjusted, the criterion is AIC 
#' (no need to adjust complexity since the total number of explanatory variables is the same)
#' 
#' @param data the data to boot, containing Surv_obj, B_original, and B_updated, also covariates
#' @param dist_boot the parametric distribution
#' @param indices placeholder, allowing boot to select sample
#'
#' @noRd

AICdiff_function_adjusted <- function(data, dist_boot, indices) {
  
  d <- data[indices, ] #allows boot to select sample
  d1 <- subset(d, select = -c(B_updated))
  d2 <- subset(d, select = -c(B_original))
    
  fit1 = survreg(Surv_obj ~ ., dist=dist_boot, data=d1)
  fit2 = survreg(Surv_obj ~ ., dist=dist_boot, data=d2)
  
  return(AIC(fit2) - AIC(fit1))
  
}
