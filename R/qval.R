#' @title qval
#'
#' @description A function to retrieve a vector of q-values of each taxa
#'    for all variables that are included in the pseudo-value regression model.
#'
#' @param SOHPIEres An object called after running SOHPIE_DNA.
#'
#' @return A table that includes q-values for all predictor variables considered in the regresson.
#'
#' @examples
#' \donttest{
#' data(combinedamgut) # A complete data containing columns with taxa and clinical covariates.
#'
#' # Note: The line below will use a toy example with the first 30 out of 138 taxa.
#' OTUtab = combinedamgut[ , 8:37]
#' # Clinical/demographic covariates (phenotypic data):
#' # Note: All of these covariates will be included in the regression, so
#' # please make sure that phenodat includes the variables that will be analyzed only.
#' phenodat = combinedamgut[, 1:7] # first column is ID, so not using it.
#' # Obtain indices of each grouping factor
#' # In this example, a variable indicating the status of living with a dog was chosen (i.e. bin_dog).
#' # Accordingly, Groups A and B imply living without and with a dog, respectively.
#'  newindex_grpA = which(combinedamgut$bin_dog == 0)
#'  newindex_grpB = which(combinedamgut$bin_dog == 1)
#'
#' SOHPIEres <- SOHPIE_DNA(OTUdat = OTUtab, clindat = phenodat,
#' groupA = newindex_grpA, groupB = newindex_grpB)
#'
#' # Create an object to keep the table with q-values using qval() function.
#' qvaltab <- qval(SOHPIEres)
#' }
#' @export

qval <- function(SOHPIEres) {
  qvaltab <- SOHPIEres$q_values
  return(qvaltab = qvaltab)
}
