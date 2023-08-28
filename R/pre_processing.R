#' STEP 0: pre_processing
#'
#' Process the given data by handling the inflated zeros, 
#' then computing the relative abundance and taking log transformation of the relative abundance
#'
#' @param data a \eqn{n \times p} matrix which is the raw taxon count table and needs to be processed, where rows are the samples, columns are the taxa
#' @param min_prevalence the minimum prevalence of taxa that could be considered in the selection procedure
#' @param mult_repl true or false, if the user wants to use Bayesian Multiplicative Replacement or not to handle inflated zeros, respectively
#'
#' @noRd

pre_processing <- function(data, min_prevalence, mult_repl) {

  # remove rare taxa
  killbugs = which(apply(data>0, 2, mean) < min_prevalence)
  if (length(killbugs) > 0) data = data[, -killbugs]

  # compute relative abundance by Bayesian Multiplicative Replacement or adding a constant then dividing sample library size
  if (mult_repl == TRUE){
    data_prepped = cmultRepl(data)
  } else{
    data = data + 0.5
    data_prepped = data / rowSums(data)
  }

  data_prepped = log(data_prepped) # take the log transformation

  return(data_prepped) # return the result in the end

}
