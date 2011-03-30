CompareTests <- function(stdtest,sampledtest,strata=NA,goldstd="sampledtest")
{
## Purpose: Estimate % agreement, symmetry, test operating characteristics, and Kappa
## under two-phase sampling.  Also estimate sensitivity, specificity, PPV, NPV and their
## ordinal generalizations: the ith PV (iPV) and the ith category-specific classification
## probability (iCSCP).
##
## stdtest is observed on everyone, sampledtest is observed only on some
## fraction of those who get each stdtest result. Both tests have I >=2  categories
## 
## There can be many sampling strata, each stratum being a separate I by I matrix of cells.
## Set strata to NA if there is no sampling or if you have simple random sampling.
##
## For the iPVs, I presume that the sampledtest is the gold standard.  If instead stdtest 
## is the gold standard, then set "stdtest" in goldstd, or if none, set to FALSE
##
## Order the categories from least to most severe, for binary (-,+) or (0,1) to make sure
## that what is output as sensitivity is not the specificity, or that PPV is not reported
## as NPV.
## ----------------------------------------------------------------------
## Arguments: IxIxS array of cells plus IxS array of margins
## ----------------------------------------------------------------------
## Author: David Edelstein and Hormuzd Katki, Date: 27 Jul 2009, 12:24

##
# Take raw data on each unit make the IxIxS array of cells and IxS array of margins
# If no strata, cells are IxI and margins are an I element vector  
##
if (!any(is.na(strata))) {
  temp <- fulltable(as.factor(sampledtest),as.factor(stdtest),as.factor(strata))
  dimtemp <- dim(temp); I <- dim(temp)[1]-3; S <- dim(temp)[3]-3
  # For I ratings and S strata, temp is (I+3)x(I+3)
  # Get rid of the 3 NA, <NA>, and Sum row/cols/strata; cells is IxIxS
  cells <- temp[1:I,1:I,1:S] 
  # For each stratum (3rd index), Pick out the Sum row (last row) the p colsums; is 1xIxS
  margins <- temp[I+3,1:I,1:S] 
}
else { # note that if any strata have missing values, I ignore all the strata
  warning('Some sampling strata have missing values. Ignoring all strata.')
  temp <- fulltable(as.factor(sampledtest),as.factor(stdtest))
  dimtemp <- dim(temp); I <- dim(temp)[1]-3
  # For I ratings and S strata, temp is (I+3)x(I+3)
  # Get rid of the 3 NA, <NA>, and Sum row/cols/strata; cells is IxIx1
  cells <- array(dim=c(I,I,1)) # Ensure third index remains when only 1 stratum
  cells[,,1] <- temp[1:I,1:I,drop=FALSE] 
  # Pick out the Sum row (last row) the p colsums; is vector of I elements
  margins <- temp[I+3,1:I] 
}

##
# Calculate sampling weights and fractions for each column
# If a margin is zero, make the weight for that zero as well
##
weights <- margins/(colSums(cells)+1*(colSums(cells)==0))
if (dim(cells)[3]==1) # If no strata, need to transpose weights to keep as a row
  weights <- t(t(weights))
fractions <- 1/weights 


#Weight cells
weightedCells <- cells
for(j in 1:dim(cells)[2]){
	for(s in 1:dim(cells)[3]){
		weightedCells[,j,s] <- weightedCells[,j,s]*weights[j,s]
	}#end of for
}#end of for


#Create estimation of the cohort
estCohort <- cells[,,1]
for(i in 1:dim(cells)[1]){
	for(j in 1:dim(cells)[2]){
		estCohort[i,j] <- sum(weightedCells[i,j,])
	}#end of for
}#end of for


#calculate test 2 margin for estimated cohort
margins2 <- rowSums(estCohort)

#calculate test 1 margin for estimated cohort (this is fixed and known)
margins1 <- colSums(estCohort)


#Find the variance of each cell and generate hat p
hatp <- cells
vars <- estCohort
# If any zero cells, add 0.5 to all cells
if (any(estCohort==0)) {cells <- cells+0.5}
for(i in 1:dim(cells)[1]){
	for(j in 1:dim(cells)[2]){
		vars[i,j] <- 0
		for(s in 1: dim(cells)[3]){
			hatp[i,j,s] <- cells[i,j,s]/(sum(cells[,j,s])+1*(sum(cells[,j,s])==0))
			vars[i,j] <- vars[i,j] + sum(cells[,j,s])*hatp[i,j,s]*(1-hatp[i,j,s])*weights[j,s]^2
		}#end of for
	}#end of for
}#end of for


#Find the covariances
covars <- array(0, c(dim(cells)[1], dim(cells)[1], dim(cells)[2]))
for(i in 1:dim(cells)[1]){
	for(j in 1:dim(cells)[2]){
		for(k in 1:dim(cells)[1]){
			covars[i,k,j] <- 0
			for(s in 1: dim(cells)[3]){
				covars[i,k,j] <- covars[i,k,j] - sum(cells[,j,s])*hatp[i,j,s]*hatp[k,j,s]*weights[j,s]^2
			}#end of for
		}#end of for
	}#end of for
}#end of for
for(i in 1:dim(vars)[1]){
	for(j in 1:dim(vars)[2]){
		covars[i,i,j] <- vars[i,j]
	}#end of for
}#end of for

# If any zero cells, remove the 0.5 already added to all cells
if (any(estCohort==0)) {cells <- cells-0.5}


#Percent agreement
p0 <- sum(diag(estCohort))/sum(estCohort)
varP0 <- sum(diag(vars)/(sum(estCohort)^2))


#Percent agreement by category
percAgrByCat <- diag(estCohort)/(margins2+margins1-diag(estCohort))
varPercAgrByCat <- (diag(vars)+percAgrByCat^2*(apply(vars,1,sum)-diag(vars)))/(margins2+margins1-diag(estCohort))^2


# Naive conditional Symmetry testing (McNemar's test for 2x2) that ignores sampling
naivesymm <- 0
for(i in 1:dim(cells)[1]){
  for(j in 1:dim(cells)[1]){
    if (i>j) 
		naivesymm <- naivesymm + ((estCohort[i,j]-estCohort[j,i])^2)/(estCohort[i,j]+estCohort[j,i])
  }#end of for
}#end of for

# Conditional Symmetry test accounting for sampling
# Have to compute the proper null hypothesis for each off-diagonal pair
# This does it for all i,j.  For i=j, the null is zero, and the nulls sum to zero for each
# off diagonal pair.  By convention, I keep the positive null to compute the test.
nulls <- estCohort
for(i in 1:dim(cells)[1]){
	for(j in 1:dim(cells)[2]){
		nulls[i,j] <- 0
		for(s in 1: dim(cells)[3]){
			nulls[i,j] <- nulls[i,j] + (cells[i,j,s]+cells[j,i,s])*(fractions[i,s]-fractions[j,s])/(fractions[i,s]+fractions[j,s])
		}#end of for
	}#end of for
}#end of for
# Plug in the nulls to get the test
# since I use nulls[i,j]>0, then cells[j,i]>cells[i,j]
condsymm <- 0
# 2x2 table, 1 stratum
condsymm <- (cells[2,1,]-cells[1,2,]-nulls[1,2])^2/(cells[2,1,]+cells[1,2,])
# Fix this later for IxI tables with S strata
# for(i in 1:dim(cells)[1]){
#   for(j in 1:dim(cells)[1]){
#     if (i>j )
#       condsymm <- condsymm + (sum(cells[i,j,])-sum(cells[j,i,])-abs(nulls[i,j]))^2/(sum(cells[i,j,])+sum(cells[j,i,]))
#   }#end of for
# }#end of for

# Unconditional Symmetry test accounting for sampling (does not reduce to McNemar's test)
uncondsymm <- 0
for(i in 1:dim(cells)[1]){
  for(j in 1:dim(cells)[1]){
    if (i>j)
      uncondsymm <- uncondsymm + ((estCohort[i,j]-estCohort[j,i])^2)/(vars[i,j]+vars[j,i])
  }#end of for
}#end of for

# Conditional OR and its unconditional variance
condOR <- estCohort[2,1]/estCohort[1,2]
VlogcondOR <- estCohort[2,1]^-2 * vars[2,1] + estCohort[1,2]^-2 * vars[1,2]

#Kappa claculations
pe <- sum(margins1*margins2/(sum(estCohort)^2))
kappa <- (p0-pe)/(1-pe)
marginCovars <- array(0, c(dim(cells)[1], dim(cells)[1]))
for(i in 1:dim(cells)[1]){
	for(j in 1:dim(cells)[1]){
		marginCovars[i,j] <- sum(covars[i,j,])
	}#end of for
}#end of for
varPe <- 0
for(i in 1:dim(cells)[1]){
	for(j in 1:dim(cells)[2]){
		varPe <- varPe + margins1[i]*margins1[j]*marginCovars[i,j]
	}#end of for
}#end of for
varPe <- varPe/(sum(estCohort)^4)
covarP0Pe <- 0
for(i in 1:dim(cells)[1]){
	for(j in 1:dim(cells)[2]){
		covarP0Pe <- covarP0Pe + margins1[j]*covars[i,j,i]
	}#end of for
}#end of for
covarP0Pe <- covarP0Pe/(sum(estCohort)^3)
varKappa <- (varP0 + 2*covarP0Pe*(kappa-1) + varPe*(kappa-1)^2)/((1-pe)^2)


#Test operating characteristics
# PPV NPV
iPV <- diag(estCohort)/margins1
varsiPV <- diag(vars)/margins1^2
stderrlogitiPV <- sqrt(varsiPV)/iPV/(1-iPV)
CIiPV <- matrix(log(iPV/(1-iPV)),nrow=length(iPV),ncol=2)+stderrlogitiPV%*%t(c(-1.96,1.96))
CIiPV <- 1/(1+exp(-CIiPV))
# Sens Spec
iCSCP <- diag(estCohort)/margins2
varsiCSCP <-((1-iCSCP)^2*diag(vars) + (iCSCP^2)*(diag(marginCovars)-diag(vars)))/margins2^2
stderrlogitiCSCP <- sqrt(varsiCSCP)/iCSCP/(1-iCSCP)
CIiCSCP <- matrix(log(iCSCP/(1-iCSCP)),nrow=length(iCSCP),ncol=2)+stderrlogitiCSCP%*%t(c(-1.96,1.96))
CIiCSCP <- 1/(1+exp(-CIiCSCP))


# output
output <- list()
output <- append(output, list(cells))
output <- append(output, list(estCohort))
output <- append(output, list(vars))
output <- append(output, list(covars))
output <- append(output, list(p0))
output <- append(output, list(varP0))
output <- append(output, list(percAgrByCat))
output <- append(output, list(varPercAgrByCat))
#output <- append(output, list(naivesymm))
#output <- append(output, list(condsymm))
output <- append(output, list(uncondsymm))
#output <- append(output, list(condOR))
#output <- append(output, list(VlogcondOR))
output <- append(output, list(marginCovars))
output <- append(output, list(kappa))
output <- append(output, list(varKappa))
output <- append(output, list(iPV))
output <- append(output, list(varsiPV))
output <- append(output, list(iCSCP))
output <- append(output, list(varsiCSCP))
names(output) <- c("Cells", "EstCohort", "Cellvars", "Cellcovars", "p0", "Varp0", "AgrCat",
                   "VarAgrCat", "uncondsymm", "Margincovars",
                   "Kappa","Kappavar","iPV", "VarsiPV", "iCSCP", "VarsiCSCP")

# Print output
cat("The weighted contingency table:\n")
print(estCohort)
cat("\n\nAgreement Statistics\n\n")
cat("pct agree and 95% CI:",output$p0,
    "(",output$p0+1.96*c(-sqrt(output$Varp0),sqrt(output$Varp0)),")\n")
cat("pct agree by categories and 95% CI","\n")
leftci <- output$AgrCat-1.96*sqrt(output$VarAgrCat)
rightci <-output$AgrCat+1.96*sqrt(output$VarAgrCat)
agrcatstats <- cbind(output$AgrCat,leftci,rightci)
colnames(agrcatstats) <- c("est","left","right")
print(agrcatstats)
cat("Kappa and 95% CI:",output$Kappa,
    "(",output$Kappa+1.96*c(-sqrt(output$Kappavar),sqrt(output$Kappavar)),")\n")
cat("symmetry chi-square:",output$uncondsymm,"p=",1-pchisq(output$uncondsymm,6),"\n")
cat("\n")

if (goldstd != FALSE) {
   cat("\n\nDiagnostic Accuracy statistics\n\n")
	PV <- cbind(iPV,CIiPV)
	colnames(PV) <- c("est","left","right")
	sensspec <- cbind(iCSCP,CIiCSCP);
	colnames(sensspec) <- c("est","left","right")
	if (goldstd == "sampledtest") {
		if (length(iPV)==2) {
  			rownames(sensspec) <- c("Specificity","Sensitivity")
  			rownames(PV) <- c("NPV","PPV")
		}
		else {
  			rownames(PV) <- paste(1:dim(estCohort)[1],"PV",sep="")
  			rownames(sensspec) <- paste(1:dim(estCohort)[1],"CSCP",sep="")
  			# iCSCP: ith category-specific classification probability
		}
	}
	else {
		if (length(iPV)==2) {
  			rownames(sensspec) <- c("NPV","PPV")
  			rownames(PV) <- c("Specificity","Sensitivity")
		}
		else {
  			rownames(PV) <- paste(1:dim(estCohort)[1],"CSCP",sep="")
  			rownames(sensspec) <- paste(1:dim(estCohort)[1],"PV",sep="")
  			# iCSCP: ith category-specific classification probability
		}	
	}
print(PV);print(sensspec)
}

return(output)
}#end of program
