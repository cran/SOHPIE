#' @title thetatildefun
#'
#' @description A function to calculate jackknife pseudo-values
#'
#' @param thetahatinput A network centrality calculated from association matrix
#'      for whole sample.
#' @param thetahatdropinput Network centralities calculated from re-estimated 
#'     association matrices for leave-one-out samples.
#' @param sizegroup Sample size for group.
#'
#' @return A jackknife pseudo-value will be returned.
#'   
#' @examples
#' \donttest{
#' # In this example, the subset of the American Gut Project data will be used.
#' data(combinedamgut) # A complete data containing columns with taxa and clinical covariates.
#'
#' # Note: The line below will use a toy example with the first 30 out of 138 taxa.
#' OTUtab = combinedamgut[ , 8:37]
#' 
#' # Obtain indices of each grouping factor
#' # In this example, a variable indicating the status of living
#' # with a dog was chosen (i.e. bin_dog).
#' # Accordingly, Groups A and B imply living without and with a dog, respectively.
#' newindex_grpA = which(combinedamgut$bin_dog == 0)
#' newindex_grpB = which(combinedamgut$bin_dog == 1)
#'
#' # Now, we estimate (and re-estimate) association matrices 
#' # for each group separately.
#' asso_matA = asso_mat(OTUdat=OTUtab, group=newindex_grpA) 
#' asso_matB = asso_mat(OTUdat=OTUtab, group=newindex_grpB) 
#' 
#' # Calculate the network centrality.
#' thetahat_grpA = thetahats(asso_matA$assomat) 
#' thetahat_grpB = thetahats(asso_matB$assomat)
#' 
#' # Obtain network centrality for the re-estimated association matrices.
#' thetahat_drop_grpA = sapply(asso_matA$reest.assomat, thetahats)
#' thetahat_drop_grpB = sapply(asso_matB$reest.assomat, thetahats)
#' 
#' # Sample sizes for each group.
#' n_A <- length(newindex_grpA) 
#' n_B <- length(newindex_grpB)
#' 
#' # Now calculate jackknife pseudo-values for each group.
#' thetatilde_grpA = thetatildefun(thetahat_grpA, thetahat_drop_grpA, n_A)
#' thetatilde_grpB = thetatildefun(thetahat_grpB, thetahat_drop_grpB, n_B)
#' }
#' @export


thetatildefun <- function(thetahatinput, thetahatdropinput, sizegroup) {
        thetatildeout = matrix(NA, ncol=length(thetahatinput), nrow=sizegroup)
        thetatildeout = sapply(1:nrow(thetahatdropinput), function(k) {
                sizegroup * thetahatinput[k] - (sizegroup - 1) * thetahatdropinput[k, ]
        })
        return(thetatildeout)
}
