## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# library(robustbase) # To fit a robust regression.
# library(parallel) # To use mclapply() when reestimating the association matrix.
# library(dplyr)  # For the convenience of tabulating p-values, coefficients, and q-values.
# library(fdrtool) # For false discovery rate control.
# library(gtools) # To estimate an association matrix via SparCC.

## ----setup--------------------------------------------------------------------
library(SOHPIE)

## -----------------------------------------------------------------------------
set.seed(20050505)
data(combinedamgut) # A complete data containing columns with taxa and clinical covariates.

## -----------------------------------------------------------------------------
# Note: Again, we will use a toy example with the first 30 out of 138 taxa.
OTUtab = combinedamgut[ , 8:37]

# Clinical/demographic covariates (phenotypic data):
# Note: All of these covariates in phenodat below will be included in the regression 
#       when you use SOHPIE_DNA function later. Please make sure 
#       phenodat below include variables that will be analyzed only.
phenodat = combinedamgut[, 1:7] # first column is ID, so not using it.

## -----------------------------------------------------------------------------
# Obtain indices of each grouping factor.
# In this example, a variable indicating the status of living with a dog was chosen (i.e. bin_dog).
# Accordingly, Groups A and B imply living without and with a dog, respectively.
newindex_grpA = which(combinedamgut$bin_dog == 0)
newindex_grpB = which(combinedamgut$bin_dog == 1)

## -----------------------------------------------------------------------------
SOHPIEres <- SOHPIE_DNA(OTUdat = OTUtab, clindat = phenodat, 
                        groupA = newindex_grpA, groupB = newindex_grpB, c = 0.5)

## -----------------------------------------------------------------------------
# qval() function will get you a table with q-values.
qval(SOHPIEres)

## -----------------------------------------------------------------------------
# Create an object to keep the table with q-values.
qvaltab <- qval(SOHPIEres)
# Retrieve a vector of q-values for a single variable of interest.
qval_specific_var(qvaltab = qvaltab, varname = "bin_dog")

## -----------------------------------------------------------------------------
# Please do NOT forget to provide the name of variable in DCtaxa_tab(groupvar = )
# and the level of significance (0.3 in this example).
DCtaxa_tab <- DCtaxa_tab(qvaltab = qvaltab, groupvar = "bin_dog", alpha = 0.3)
DCtaxa_tab

## -----------------------------------------------------------------------------
data(combineddietswap)

OTUtab = combineddietswap[ , 5:ncol(combineddietswap)]
phenodat = combineddietswap[ , 1:4] # first column is ID, so not using it.


## -----------------------------------------------------------------------------
# Obtain indices for each groups 
# (i.e. African-Americans from Pittsburgh (AAM) vs. Africans from rural South Africa (AFR))
# at baseline (time = 1) and at 29-days (time = 6)

# Group A1 for AAM at baseline.
newindex_A1 = which(combineddietswap$timepoint == 1 & combineddietswap$nationality == "AAM")
# Group A6 for AAM at 29-days.
newindex_A6 = which(combineddietswap$timepoint == 6 & combineddietswap$nationality == "AAM")
# Group B1 for AFR at baseline.
newindex_B1 = which(combineddietswap$timepoint == 1 & combineddietswap$nationality == "AFR")
# Group A6 for AFR at 29-days.
newindex_B6 = which(combineddietswap$timepoint == 6 & combineddietswap$nationality == "AFR")

## -----------------------------------------------------------------------------
est_asso_matA1 = asso_mat(OTUdat=OTUtab, group=newindex_A1)
est_asso_matA6 = asso_mat(OTUdat=OTUtab, group=newindex_A6)
est_asso_matB1 = asso_mat(OTUdat=OTUtab, group=newindex_B1)
est_asso_matB6 = asso_mat(OTUdat=OTUtab, group=newindex_B6)

## For each group, we take the difference of estimated 
## association matrices between time points (29-days minus baseline).
asso_mat_diffA61 = est_asso_matA6$assomat - est_asso_matA1$assomat
asso_mat_diffB61 = est_asso_matB6$assomat - est_asso_matB1$assomat

## Similarly, for each group, we take the difference of re-estimated 
## association matrices between time points (29-days minus baseline).
asso_mat_drop_diffA61 = mapply('-', est_asso_matA6$reest.assomat,
                               est_asso_matA1$reest.assomat, SIMPLIFY=FALSE)
asso_mat_drop_diffB61 = mapply('-', est_asso_matB6$reest.assomat,
                               est_asso_matB1$reest.assomat, SIMPLIFY=FALSE)

## -----------------------------------------------------------------------------
## For changes in network centrality between association matrices estimated from the whole data.
thetahat_grpA = thetahats(asso_mat_diffA61)
thetahat_grpB = thetahats(asso_mat_diffB61)
## For changes in network centrality between association matrices 
## re-estimated from the leave-one-out sample.
thetahat_drop_grpA = sapply(asso_mat_drop_diffA61, thetahats)
thetahat_drop_grpB = sapply(asso_mat_drop_diffB61, thetahats)

## -----------------------------------------------------------------------------
# Sample sizes for each group.
n_A <- length(newindex_A1) 
n_B <- length(newindex_B1)

# Jackknife pseudo-values.
thetatilde_grpA = thetatildefun(thetahat_grpA, thetahat_drop_grpA, n_A)
thetatilde_grpB = thetatildefun(thetahat_grpB, thetahat_drop_grpB, n_B)

thetatilde = rbind(thetatilde_grpA, thetatilde_grpB)

## -----------------------------------------------------------------------------
fitmod = pseudoreg(pseudoval=thetatilde, clindat=phenodat, c=0.5)

## -----------------------------------------------------------------------------
summary.result = pseudoreg.summary(pseudo.reg.res=fitmod, 
                                   taxanames=colnames(OTUtab))

## In this study, the main grouping variable was nationality. 
# We can use qval_specific_var to see the q-values only.
# of nationality. Alternatively, we further use DCtaxa_tab to see
# the DC taxa with their p-values.
qval_specific_var(summary.result$q_values, varname = "nationalityAFR")
DCtaxa_tab(summary.result$q_values, groupvar = "nationalityAFR", alpha=0.05)


