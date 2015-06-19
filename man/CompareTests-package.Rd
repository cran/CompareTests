\name{CompareTests-package}
\alias{CompareTests-package}
\docType{package}
\title{
Correct for Verification Bias in Diagnostic Accuracy & Agreement
}
\description{
A standard test is observed on all specimens.  We treat the second test (or sampled test) as being conducted on only a stratified sample of specimens.  Verification Bias is this situation when the specimens for doing the second (sampled) test is not under investigator control.  We treat the total sample as stratified two-phase sampling and use inverse probability weighting.  We estimate diagnostic accuracy (category-specific classification probabilities; for binary tests reduces to specificity and sensitivity) and agreement statistics (percent agreement, percent agreement by category, Kappa (unweighted), Kappa (quadratic weighted) and symmetry test (reduces to McNemar's test for binary tests)). 
}
\details{
\tabular{ll}{
Package: \tab CompareTests\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2015-06-19\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}

You have a dataframe with columns "stdtest" (no NAs allowed; all specimens with NA stdtest results are dropped), "sampledtest" (a gold standard which is NA for some specimens), sampling strata "strata1" "strata2" (values cannot be missing for any specimens). Correct for Verification Bias in the diagnostic and agreement statistics with CompareTests(stdtest,sampledtest,interaction(strata1,strata2),goldstd="sampledtest")


}
\author{
Hormuzd A. Katki and David W. Edelstein

Maintainer: Hormuzd Katki <katkih@mail.nih.gov>

}
\references{

Katki HA, Li Y, Edelstein DW, Castle PE.  Estimating the agreement and diagnostic accuracy of two diagnostic tests when one test is conducted on only a subsample of specimens.  Stat Med. 2012 Feb 28; 31(5): 10.1002/sim.4422.

}
\keyword{ package }
\examples{

# Get specimens dataset
data(specimens)

# Get diagnostic and agreement statistics if sampledtest is the gold standard
CompareTests(specimens$stdtest,specimens$sampledtest,specimens$stratum)

# Get diagnostic and agreement statistics if stdtest is the gold standard
CompareTests(specimens$stdtest,specimens$sampledtest,specimens$stratum,goldstd="stdtest")

# Get agreement statistics if neither test is a gold standard
CompareTests(specimens$stdtest,specimens$sampledtest,specimens$stratum,goldstd=FALSE)

}
