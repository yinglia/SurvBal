#' STEP 3: go over the selected to determine the final model
#'
#' Find model that gives us the smallest p val
#'
#' @param selected the candidate(s) that have been selected so far
#' @param criterion the criterion to select balance in the final model. ``\code{min_pvalue}'' or ``\code{min_decrement_pvalue}''
#' @param threshold the threshold of p-value decrement if \code{criterion} is ``\code{min_decrement_pvalue}''
#'
#' @noRd

select_model <- function(selected, criterion, threshold) {

  if (any(!is.na(selected[, 4])) & any(selected[, 4] == 0, na.rm=TRUE)) selected = selected[1:(which(selected[, 4] == 0)-1), ]
  
  if (criterion == "min_pvalue"){
    ind_selected_to = which.min(selected[, 3])
  } else if (criterion == "min_decrement_pvalue"){
    if (nrow(selected) == 2){
      ind_selected_to = 2
    } else{
      nrow_selected = nrow(selected)
      increment = c(NA, NA, (selected[3:nrow_selected, 3] - selected[2:(nrow_selected-1), 3]) / selected[2:(nrow_selected-1), 3])
      ind_selected_to = which(increment > -threshold)[1] - 1
      if (is.na(ind_selected_to)) ind_selected_to = 2
    }
  }

  ind_selected_to = max(2, ind_selected_to)
  result = selected[1:ind_selected_to, ]

  return (result)

}
