#' @title pseudoreg.summary()
#'
#' @description A function to output summary results
#'      (p-values, q-values, and coefficient estimates)
#'      from the fitted pseudo-value regression.
#'
#' @param pseudo.reg.res A fitted pseudo-value regression using pseudoreg()
#' @param taxanames Names of taxa from the OTU table.
#'
#' @return A pseudo-value regression is fitted. Please use pseudoreg.summary()
#'        to output p-values, q-values, and coefficient estimates.
#'
#' @examples
#' \donttest{
#' # In this example, the subset of the American Gut Project data will be used.
#' data(combinedamgut) # A complete data containing columns with taxa and clinical covariates.
#'
#' # Note: The line below will use a toy example with the first 30 out of 138 taxa.
#' OTUtab = combinedamgut[ , 8:37]
#'
#' #Clinical/demographic covariates (phenotypic data):
#' # Note: All of these covariates will be included in the regression, so
#' # please make sure that phenodat includes the variables that will be analyzed only.
#' phenodat = combinedamgut[, 1:7] # first column is ID, so not using it.
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
#'
#' thetatilde = rbind(thetatilde_grpA, thetatilde_grpB)
#'
#' # Map the column names (taxa names)
#' colnames(thetatilde) = colnames(OTUtab)
#'
#' # Fit a pseudo-value regression using jackknife pseudovalues
#' # and phenotypic data. A reminder that the phenotypic data should
#' # contain a set of predictor variables to be fitted.
#' fitmod = pseudoreg(pseudoval=thetatilde, clindat=phenodat, c=0.5)
#'
#' # Extract summary results from the fitted model from fitmod object above.
#' summary.result = pseudoreg.summary(pseudo.reg.res=fitmod, taxanames=colnames(OTUtab))
#' }
#' @export
#' @import dplyr
#' @import fdrtool


pseudoreg.summary <- function(pseudo.reg.res, taxanames) {
        ### Obtain p-values (and coefficient estimates if interested) for each taxa:
        beta_hat = vector(mode = "list", length(pseudo.reg.res))
        p_values = vector(mode = "list", length(pseudo.reg.res))
        k = NULL
        for(k in 1:length(pseudo.reg.res)) {
                # Estimates for the beta coefficients:
                beta_hat[[k]] <- summary(pseudo.reg.res[[k]])$coef[-1, "Estimate"]
                p_values[[k]] <- summary(pseudo.reg.res[[k]])$coef[-1, "Pr(>|t|)"]

        }

        beta_hat = as.data.frame(bind_rows(beta_hat))  # Convert list into data.frame
        rownames(beta_hat) <- taxanames # Map the taxa names to the data.frame for betahats

        p_values = as.data.frame(bind_rows(p_values))  # Convert list into data.frame
        rownames(p_values) <- taxanames # Map the taxa names to the data.frame for p-values

        # Compute the q-values :
        suppressWarnings({q_values = as.data.frame(apply(p_values, 2, function(x) fdrtool(x, statistic = "pvalue", plot=FALSE, verbose = FALSE)$qval))})
        rownames(q_values) <- taxanames # Map the taxa names to the data.frame for q-values

        return(list(beta_hat = beta_hat, p_values = p_values, q_values = q_values))
}


