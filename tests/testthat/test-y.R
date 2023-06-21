# Load and process the data to run SOHPIE
data(combinedamgut)
OTUtab = as.data.frame(combinedamgut[ , 8:37])
phenodat = combinedamgut[ ,1:7]
newindex_grpA <- which(combinedamgut$bin_dog == 0)
newindex_grpB <- which(combinedamgut$bin_dog == 1)
SOHPIEres <- SOHPIE_DNA(OTUdat = OTUtab, clindat = phenodat, groupA = newindex_grpA, groupB = newindex_grpB)

# Create a function to check the sign of p-values
# and adjusted p-values:
sign_check <- function(pval) {
        if (sum(pval < 0) > 0) {
                return("ERROR: Negative p-values")
        } else {
                return("No negative p-values")
        }
}


# A function to check if results have any negative p-values or q-values:
test_that("Check any negative p-values after SOHPIE", {
        expect_equal(sign_check(SOHPIEres$p_values), "No negative p-values")
        expect_equal(sign_check(SOHPIEres$q_values), "No negative p-values")
})

# Use testthat to check if results return any missing or NA:
# Note: coefficient estimates and p-values
#       (and q-values) should not have any missing or NA.
test_that("Check the missingness after SOHPIE", {
        expect_false(any(is.na(SOHPIEres$beta_hat)))
        expect_false(any(is.na(SOHPIEres$p_values)))
        expect_false(any(is.na(SOHPIEres$q_values)))
})


