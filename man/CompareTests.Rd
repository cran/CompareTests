\name{CompareTests}
\alias{CompareTests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Correct for Verification Bias in Diagnostic Accuracy & Agreement
}
\description{
A standard test is observed on all specimens.  We treat the second test (or sampled test) as being conducted on only a stratified sample of specimens.  We treat the total sample as stratified two-phase sampling and use inverse probability weighting.  We estimate diagnostic accuracy (category-specific classification probabilities; for binary tests reduces to specificity and sensitivity) and agreement statistics (percent agreement, percent agreement by category, Kappa (unweighted), Kappa (quadratic weighted) and symmetry tests (reduces to McNemar's test for binary tests)).
}
\usage{
CompareTests(stdtest, sampledtest, strata = NA, goldstd = "sampledtest")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stdtest}{A vector of standard test results.  Any NA test results are dropped from the analysis entirely.
}
  \item{sampledtest}{A vector of test results observed only on a sample of specimens.  Test results with NA are assumed to no be observed for that specimen
}
  \item{strata}{The sampling stratum each specimen belongs to.  Set to NA if no sampling or simple random sampling.
}
  \item{goldstd}{For outputing diagnostic accuracy statistics, denote if "stdtest" or "sampledtest" is the gold standard.  If no gold standard, set to FALSE.
}
}

\value{

Outputs to screen the estimated contingency table of paired test results, agreement statistics, and diagnostic accuracy statistics.

Returns a list with the following components
\item{Cells }{Observed contingency tables of pair test results for each stratum}
\item{EstCohort }{Weighted contingency table of each pair of test results}
\item{Cellvars }{Variance of each weighted cell count}
\item{Cellcovars }{Variance-covariance matrix for each column of weighted cell counts}
\item{p0 }{Percent agreement}
\item{Varp0 }{Variance of percent agreement}
\item{AgrCat }{Percent agreement by each test category}
\item{VarAgrCat }{Variance of Percent agreement by each test category}
\item{uncondsymm }{Symmetry test test statistic}
\item{Margincovars }{covariance of each pair of margins}
\item{Kappa }{Kappa (unweighted)}
\item{Kappavar }{Variance of Kappa}
\item{iPV }{Each predictive value (for binary tests, NPV and PPV)}
\item{VarsiPV }{Variance of each predictive value (for binary tests, NPV and PPV)}
\item{iCSCP}{Each category-specific classification probability (for binary tests, specificity and sensitivity}
\item{VarsiCSCP}{Variance of each category-specific classification probability (for binary tests, specificity and sensitivity}
\item{WeightedKappa }{Kappa (quadratic weights)}
\item{varWeightedKappa }{Variance of quadratic-weighted Kappa}
}
\references{
Katki HA, Li Y, Edelstein DW, Castle PE.  Estimating the agreement and diagnostic accuracy of two diagnostic tests when one test is conducted on only a subsample of specimens. Stat Med. 2012 Feb 28; 31(5): 10.1002/sim.4422.
}
\author{
Hormuzd A. Katki and David W. Edelstein
}
\note{
Order the categories from least to most severe, for binary (-,+) or (0,1) to make sure that what is output as sensitivity is not the specificity, or that PPV is not reported as NPV.

If you have multiple variables to be crossed to represent the sampling strata, use interaction(), e.g. strata=interaction(strata1,strata2)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

##
# Stat Med Paper 2x2 Chlamydia testing verification bias example
# Note that p for symmetry test is 0.12 not 0.02 as reported in the Stat Med paper
##

# Convert 2x2 Chlamydia testing table to a dataframe for analysis
# Include NAs for the samples where CTDT test was not conducted (HC2 was conducted on all)
HC2stdtest <- c(rep(1,827),rep(0,4998))
stratum <- HC2stdtest
CTDTsampledtest <- c(           rep(1,800), # 1,1 cell
                                rep(0,27), # 1,0 cell HC2+ , CTDT-
                                rep(NA,827-800-27), # HC2+, and no CTDT test done
                                rep(1,6), # 0,1 cell: HC2-, CTDT+
                                rep(0,396),# 0,0 cell: HC2- and CTDT-
                                rep(NA,4998-6-396) # HC2-, no CTDT test done
)
chlamydia <- data.frame(stratum,HC2stdtest,CTDTsampledtest)					

# Analysis
temp <- CompareTests(chlamydia$HC2stdtest,
                     chlamydia$CTDTsampledtest,
                     chlamydia$stratum,
                     goldstd="sampledtest"
)


##
# Example analysis of fictitious data example
##
data(specimens)
temp <- CompareTests(specimens$stdtest,
                     specimens$sampledtest,
                     specimens$stratum,
                     goldstd="sampledtest")
## The output is
# The weighted contingency table:
#                             as.factor.stdtest.
# as.factor.sampledtest.     1       2      3      4
#                      1 47.88   7.158  3.322  0.000
#                      2 20.12 104.006 21.861  2.682
#                      3  0.00  10.836 97.494  8.823
#                      4  0.00   0.000  3.322 74.495
# 
# 
# Agreement Statistics
# 
# pct agree and 95% CI: 0.8057 ( 0.7438 0.8555 )
# pct agree by categories and 95% CI 
#      est   left  right
# 1 0.6101 0.4501 0.7494
# 2 0.6241 0.5315 0.7083
# 3 0.6693 0.5562 0.7658
# 4 0.8340 0.6340 0.9358
# Kappa and 95% CI: 0.734 ( 0.6509 0.8032 )
# Weighted Kappa (quadratic weights) and 95% CI: 0.8767 ( 0.7107 0.9536 )
# symmetry chi-square: 9.119 p= 0.167 
# 
# 
# 
# Diagnostic Accuracy statistics
# 
#        est   left  right
# 1PV 0.7041 0.5422 0.8271
# 2PV 0.8525 0.7362 0.9229
# 3PV 0.7738 0.6547 0.8605
# 4PV 0.8662 0.6928 0.9490
#          est   left  right
# 1CSCP 0.8204 0.6011 0.9327
# 2CSCP 0.6996 0.6169 0.7710
# 3CSCP 0.8322 0.7219 0.9046
# 4CSCP 0.9573 0.5605 0.9975
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
