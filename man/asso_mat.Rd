% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asso_mat.R
\name{asso_mat}
\alias{asso_mat}
\title{asso_mat}
\usage{
asso_mat(OTUdat, group)
}
\arguments{
\item{OTUdat}{An OTU table with subjects in rows and taxa in columns.}

\item{group}{Indices of the subjects in a category of binary group variable.}
}
\value{
A list of an association matrix and reestimated association matrix
  is returned, which are estimated via SparCC.
}
\description{
A function to estimate an association matrix. 
  This function also includes re-estimation for leave-one-out sample.
}
\examples{
\donttest{
# In this example, the subset of the American Gut Project data will be used.
data(combinedamgut) # A complete data containing columns with taxa and clinical covariates.

# Note: The line below will use a toy example with the first 30 out of 138 taxa.
OTUtab = combinedamgut[ , 8:37]

# Obtain indices of each grouping factor
# In this example, a variable indicating the status of living
# with a dog was chosen (i.e. bin_dog).
# Accordingly, Groups A and B imply living without and with a dog, respectively.
newindex_grpA = which(combinedamgut$bin_dog == 0)
newindex_grpB = which(combinedamgut$bin_dog == 1)

# Now, we estimate (and re-estimate) association matrices 
# for each group separately.
asso_matA = asso_mat(OTUdat=OTUtab, group=newindex_grpA) 
asso_matB = asso_mat(OTUdat=OTUtab, group=newindex_grpB) 

}
}
