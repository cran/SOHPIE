#' @title coeff_specific_var
#'
#' @description A function to retrieve a vector of coefficient estimates of each taxa
#'    for one specific variable.
#'
#' @param coefftab A table that includes coefficient estimates for a specific variable.
#' @param varname Specify the name of the variable of interest.
#'
#' @return A vector of coefficient estimates for a single variable from the model.
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
#' groupA = newindex_grpA, groupB = newindex_grpB, c = 0.5)
#'
#' # coeff() function will return coefficient estimates only.
#' coefftab <- coeff(SOHPIEres)
#'
#' # coeff_specific_var() will return coefficient estimates of
#' # a single variable of interest.
#' coeff_specific_var(coefftab = coefftab, varname = "bin_dog")
#' }
#' @export

coeff_specific_var <- function(coefftab, varname) {
        coeffvec <- coefftab[varname]
        return(coeffvec = coeffvec)
}
