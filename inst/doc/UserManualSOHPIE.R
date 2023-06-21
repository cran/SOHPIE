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
                        groupA = newindex_grpA, groupB = newindex_grpB)

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

