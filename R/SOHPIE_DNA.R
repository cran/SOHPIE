#' @title SOHPIE_DNA
#'
#' @description A pseudo-value regression approach for differential co-abundance
#'     network analysis that adjusts for additional covariates.
#'
#' @param OTUdat An OTU table with subjects in rows and taxa in columns.
#' @param clindat A subdata consisting of the clinical and demographic variables that
#'               the user wants to include in the regression.
#'               (e.g., binary group indicator for intervention vs. control, continuous age, ...)
#' @param groupA Indices of the subjects in the first category (e.g., not living with a dog; see example below with American Gut Project sample data)
#'               of binary group variable.
#' @param groupB Indices of the subjects in the second category (e.g., living with a dog; see example below with American Gut Project sample data)
#'               of binary group variable.
#' @param c  The choice of trimming proportion for the least trimmed estimator of robust regression.
#'           A value has to be between 0.5 and 1 as specified in ltsReg() function in robustbase package.
#'
#' @return A list containing three data frame objects returned from this SOHPIE_DNA main function. A user will see
#'         beta coefficients, p-values, and adjusted p-values (q-values)
#'         for each predictor variables that are included in the regression model.
#'
#' @references Ahn S, Datta S. Differential Co-Abundance Network Analyses for Microbiome Data Adjusted for Clinical Covariates Using Jackknife Pseudo-Values.
#'      ArXiv [Preprint]. 2023 Mar 23:arXiv:2303.13702v1. PMID: 36994149; PMCID: PMC10055480.
#'
#' @examples
#' \donttest{
#' # In this example, the subset of the American Gut Project data will be used.
#' data(combinedamgut) # A complete data containing columns with taxa and clinical covariates.
#'
#' # Note: The line below will use a toy example with the first 30 out of 138 taxa.
#' OTUtab = combinedamgut[ , 8:37]
#'  #Clinical/demographic covariates (phenotypic data):
#' # Note: All of these covariates will be included in the regression, so
#' # please make sure that phenodat includes the variables that will be analyzed only.
#' phenodat = combinedamgut[, 1:7] # first column is ID, so not using it.
#' # Obtain indices of each grouping factor
#' # In this example, a variable indicating the status of living
#' # with a dog was chosen (i.e. bin_dog).
#' # Accordingly, Groups A and B imply living without and with a dog, respectively.
#' newindex_grpA = which(combinedamgut$bin_dog == 0)
#' newindex_grpB = which(combinedamgut$bin_dog == 1)
#'
#' SOHPIEres <- SOHPIE_DNA(OTUdat = OTUtab, clindat = phenodat,
#' groupA = newindex_grpA, groupB = newindex_grpB, c = 0.5)
#' }
#' @export
#' @import robustbase
#' @import parallel
#' @import dplyr
#' @import gtools
#' @import fdrtool

SOHPIE_DNA <- function(OTUdat, clindat, groupA, groupB, c) {
      est_method = sparcc # estimate the association matrix/network via SparCC

      ########################################################################################
      #  STEP 1. Estimate an association matrix via SparCC from the microbiome data.
      ########################################################################################
      OTUtabA = OTUdat[groupA, ]
      OTUtabB = OTUdat[groupB, ]

      n_A <- length(groupA) # Sample size for Group A
      n_B <- length(groupB) # Sample size for Group B

      sparcc.matA = est_method(x = OTUtabA)$cor.w
      sparcc.matB = est_method(x = OTUtabB)$cor.w



      #############################################################################################
      #  STEP 2. Calculate \hat{\theta}_{k} by taking the marginal sum of the
      #          association matrix to obtain the degree centrality measure for each taxon.
      #############################################################################################
      thetahat_grpA = thetahats(sparcc.matA)
      thetahat_grpB = thetahats(sparcc.matB)



      #############################################################################################
      #  STEP 3. Re-estimate association matrix using the same OTU table without i-th subject.
      #          Then, calculate \hat{\theta}_{k(i)} from the re-estimated association matrix.
      #############################################################################################
      # Re-estimation part
      # Note: You can specify the core size by specifying mc.cores option within mclapply().
      sparcc.mat_drop_grpA <- mclapply(groupA, function(j) sparcc(OTUtabA[-j, ])$cor.w)
      sparcc.mat_drop_grpB <- mclapply(groupB, function(j) sparcc(OTUtabB[-j, ])$cor.w)
      # thetahat_{-i} for each taxa
      thetahat_drop_grpA <- sapply(sparcc.mat_drop_grpA, thetahats)
      thetahat_drop_grpB <- sapply(sparcc.mat_drop_grpB, thetahats)




      #############################################################################################
      #  STEP 4. Calculate jackknife pseudovalues (\tilde{\theta}_{ik} using \hat{\theta}_{k}
      #          and \hat{\theta}_{k(i)}.
      #############################################################################################
      thetatildefun <- function(thetahatinput, thetahatdropinput, sizegroup) {
        thetatildeout = matrix(NA, ncol=length(thetahatinput), nrow=sizegroup)
        thetatildeout = sapply(1:nrow(thetahatdropinput), function(k) {
          sizegroup * thetahatinput[k] - (sizegroup - 1) * thetahatdropinput[k, ]
        })
        return(thetatildeout)
      }

      thetatilde_grpA = thetatildefun(thetahat_grpA, thetahat_drop_grpA, n_A)
      thetatilde_grpB = thetatildefun(thetahat_grpB, thetahat_drop_grpB, n_B)
      thetatilde = rbind(thetatilde_grpA, thetatilde_grpB)
      colnames(thetatilde) = colnames(OTUdat) # Map the column names (taxa names)



      #############################################################################################
      #  STEP 5. Fit a robust regression model
      ##############################################################################################
      pseudo.reg.res <- lapply(1:ncol(thetatilde), function(i) {
        m <- thetatilde[, i]
        df <- data.frame(clindat,
                         m = m)
        fit <- ltsReg(m ~ ., data = df, mcd=FALSE, alpha=c) # include a set of covariates in this model
        return(fit)
      })

      ### Obtain p-values and coefficient estimates and  if interested) for each taxa:
      ### 10/20/2023: Standard Error (SE) is also available
      beta_hat = vector(mode = "list", ncol(thetatilde))
      p_values = vector(mode = "list", ncol(thetatilde))
      stderrs = vector(mode = "list", ncol(thetatilde))
      k = NULL
      for(k in 1:ncol(thetatilde)) {
        # Estimates for the beta coefficients:
        beta_hat[[k]] <- summary(pseudo.reg.res[[k]])$coef[-1, "Estimate"]
        p_values[[k]] <- summary(pseudo.reg.res[[k]])$coef[-1, "Pr(>|t|)"]
        stderrs[[k]] <- summary(pseudo.reg.res[[k]])$coef[-1, "Std. Error"]

      }

      beta_hat = as.data.frame(bind_rows(beta_hat))  # Convert list into data.frame
      rownames(beta_hat) <- colnames(OTUdat) # Map the taxa names to the data.frame for betahats

      p_values = as.data.frame(bind_rows(p_values))  # Convert list into data.frame
      rownames(p_values) <- colnames(OTUdat) # Map the taxa names to the data.frame for p-values

      stderrs = as.data.frame(bind_rows(stderrs))  # Convert list into data.frame
      rownames(stderrs) <- colnames(OTUdat) # Map the taxa names to the data.frame for stderr

      # Compute the q-values :
      suppressWarnings({q_values = as.data.frame(apply(p_values, 2, function(x) fdrtool(x, statistic = "pvalue", plot=FALSE, verbose = FALSE)$qval))})
      rownames(q_values) <- colnames(OTUdat) # Map the taxa names to the data.frame for q-values

      return(list(beta_hat = beta_hat, p_values = p_values, q_values = q_values, stderrs = stderrs))
}
