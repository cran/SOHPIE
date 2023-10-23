#' @title DCtaxa_tab
#'
#' @description A function to obtain a list consisting of taxa that are significantly differentially connected (DC)
#' between two biological or clinical conditions. These DC taxa are resulted from the pseudo-value regression method with
#' additional covariates. In addition, a user can extract the names of DC taxa only.
#'
#' @param qvaltab A table with adjusted p-values (or q-value in this package).
#' @param groupvar Specify the name of binary indicator variable.
#' @param alpha A level of significance (e.g. 0.05).
#'
#' @return q-values and names of significantly DC taxa (e.g. taxa name) based on SOHPIE_DNA function.
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
#' # Create an object to keep the table with q-values using qval() function.
#' qvaltab <- qval(SOHPIEres)
#'
#' # Please do NOT forget to provide the name of variable in DCtaxa_tab(groupvar = ).
#' DCtaxa_tab <- DCtaxa_tab(qvaltab = qvaltab, groupvar = "bin_dog", alpha = 0.3)
#' DCtaxa_tab
#' }
#' @export

DCtaxa_tab <- function(qvaltab, groupvar, alpha) {
  qvalvec = qvaltab[groupvar]
  DCtaxa_complete_tab = filter(qvalvec, qvalvec < alpha)
  DCtaxa_names_only = rownames(DCtaxa_complete_tab)
  return(list(DCtaxa_complete_tab = DCtaxa_complete_tab, DCtaxa_names_only = DCtaxa_names_only))
}
