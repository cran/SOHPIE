% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stderrs.R
\name{stderrs}
\alias{stderrs}
\title{stderrs}
\usage{
stderrs(SOHPIEres)
}
\arguments{
\item{SOHPIEres}{An object called after running SOHPIE_DNA.}
}
\value{
A table that includes standard error of betahats
for all predictors regressed in the fitted model.
}
\description{
A function to retrieve a vector of standard error (stderr)
   of coefficient estimates (betahats) all predictor variables in the pseudo-value regression model.
}
\examples{
\donttest{
data(combinedamgut) # A complete data containing columns with taxa and clinical covariates.

# Note: The line below will use a toy example with the first 30 out of 138 taxa.
OTUtab = combinedamgut[ , 8:37]
# Clinical/demographic covariates (phenotypic data):
# Note: All of these covariates will be included in the regression, so
# please make sure that phenodat includes the variables that will be analyzed only.
phenodat = combinedamgut[, 1:7] # first column is ID, so not using it.
# Obtain indices of each grouping factor
# In this example, a variable indicating the status of living with a dog was chosen (i.e. bin_dog).
# Accordingly, Groups A and B imply living without and with a dog, respectively.
 newindex_grpA = which(combinedamgut$bin_dog == 0)
 newindex_grpB = which(combinedamgut$bin_dog == 1)

SOHPIEres <- SOHPIE_DNA(OTUdat = OTUtab, clindat = phenodat,
groupA = newindex_grpA, groupB = newindex_grpB, c = 0.5)

# stderrs() function will return standard error of betahats only.
stderrs(SOHPIEres)
}
}
