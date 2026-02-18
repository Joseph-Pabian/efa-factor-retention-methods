###########################################################################################################################################
# Important Note 1: You should go to the end of the program to input your correlation matrix as well
# as other specifications for your revised parallel analysis. 
#
# Important Note 2: The output will include many statements saying that the "maximum iteration exceeded".
# This is NOT an error. The factoring program is designed to conduct iterated principal axis factor analysis.
# We have chosen not to iterate the communalities so as a result the output includes the statement that
# the maximum iteration exceeded.
#
# General Information:
#
# This program is written in R code and conducts revised parallel analysis (R-PA). 
# The references for articles about R-PA are listed in the Publication section of this website. 
# In order to conduct this program, you need to have dowloaded the R software to your computer 
# (https://www.r-project.org/). In addition, you need to have installed the psych package  
# (http://personality-project.org/r/psych/) that is accessed to conduct the factor analyses.
# To install the psych package, type in "install.packages("psych")" in the console.
#
# In conducting an R-PA, you must input your correlation matrix. We refer to this matrix as the
# researcher's correlation matrix. The program inserts squared multiple correlations (SMCs) along 
# the diagonal of this correlation matrix in place of the 1s. We refer to the resulting matrix 
# as the researcher's reduced correlation matrix. Factors are extracted from the researcher's 
# reduced correlation matrix using principal axis factor analysis (with non-iterated communalities). 
# The eigenvalues and eigenvectors from the factoring of the reseacher's reduced correlation 
# matrix are saved.

# A defined number of parallel datasets (typically 100) are generated. Each dataset has the same number of 
# variables and the same sample size as involved in the computation of the researcher's 
# correlation matrix. Assuming normality, each parallel dataset is generated based on the factor 
# loadings of the principal axis factoring of the researcher's reduced correlation matrix. For details, see
# the relevant publications. A parallel correlation matrix is computed for each dataset. SMCs are inserted 
# along the diagonal of each parallel correlation matrix, yielding a parallel reduced correlation matrix.
# Factors are then extracted for each parallel reduced correlation matrix using principal axis 
# factor analysis (with non-iterated communalities). The eigenvalues from the factoring are saved
# for the parallel datasets. The program computes the 95%ile of the eigenvalues for each sequentially 
# extracted factor.
#
# For each sequential factor, the eigenvalue based on the researcher's reduced correlation matrix is compared 
# to the 95%ile of the eigenvalues based on the parallel reduced correlation matrices. Based on these
# results, the number of factors is estimated as well as other output.  To see a description of the output, see
# comments at the end of the program.
###########################################################################################################################################


###########################################################################################################################################
# General Orietation of the R-PA program:
# 
# This program contains two functions: "gen.pa.data" and "revised.pa". Within the R-PA program, 
# the "gen.pa.data" function is run first followed by the "revised.pa" function. 
#
# The "gen.pa.data" function generates the parallel datasets for the "revised.pa" function.
# Each parallel dataset contains N observations and P variables, which are set to be the same as 
# the researcher's dataset. 
#
# "revised.pa" is the main function for conducting revised parallel analysis.
# You (the researcher) need to supply your own correlation to matrix to run the "revised.pa" function.
# The supplied correlation matrix should have P rows and P columns. We provided a sample correlation matrix 
# at the end of this program file to illustrate how to conduct the R-PA. 
###########################################################################################################################################


###########################################################################################################################################
# Load the psych package.
# If the package is not installed, type in "install.packages("psych")" in the R console to install the package.
# We will use the "fa" function in the psych package to conduct principal axis factor analysis 
# for the "revised.pa" function.
###########################################################################################################################################

library(psych)


###########################################################################################################################################
# The "gen.pa.data" function generates parallel datasets for the "revised.pa" function. The parallel
# datasets contain N observations and P variables, which are set to be the same as the researcher's dataset. 
###########################################################################################################################################

gen.pa.data <- function(P, lambda, N, N.Factors){

	generated.data <- matrix(NA, nrow=N, ncol=P)
	
	factor.score <- replicate(N.Factors, rnorm(N, 0, 1))
	unique.score <- replicate(P, rnorm(N, 0, 1))

	unique.loading <- matrix(NA, nrow=P, ncol=1)
	unique.loading.sq <- matrix(NA, nrow=P, ncol=1)

	lambda.sq <- lambda*lambda

	for(p in 1:P){
		unique.loading.sq[p,1] <- 1-sum(lambda.sq[p, ])
		unique.loading[p,1] <- sqrt(unique.loading.sq[p,1])
		generated.data[ ,p] <- (factor.score %*% t(lambda))[ ,p] + unique.score[ ,p]*unique.loading[p,1]
	}

	generated.data

}


###########################################################################################################################################
# The "revised.pa" function conducts R-PA to estimate the number of factors underlying your variables.
# By running the "fa" command in this function each time, the screen will print
# "maximum iteration exceeded". This is because we set the maximum number of iterations 
# for convergence (max.iter) at 1 to make sure the squared multiple correlations (SMC=TRUE) 
# are used as the communality estimates. 
###########################################################################################################################################

revised.pa <- function(data, N, N.PD=100, percentile=.95) {

	# data is a correlation or covariance matrix or a raw data matrix supplied by users.
	# If a covariance matrix is supplied, it will be converted to a correlation matrix.
	# If raw data is supplied, the correlation matrix will be calcualted using pairwise deletion. 
	# N is the number of observations. Must be supplied by users.
	# N.PD is the number of parallel datasets for conducting the revised parallel analysis.
	# Default N.PD is 100.
	# percentile is the percentile to compare eigenvalues for observed dataset and generated
	# parallel datasets. Default percentile is .95. 

	P <- dim(data)[2] # Number of variables

	#############################################################################################
	# Conduct principal axis factor analysis and save the eigenvalues for the observed dataset
	#############################################################################################

	eigs.data <- fa(data, nfactors=1, rotate="none", SMC=TRUE, max.iter=1, fm="pa")$values

	# Define the max number of factors tested
	F.Max <- sum(eigs.data>0) 

	# Set up an array for saving eigenvalues for the parallel datasets
	eigs.para.array <- array(NA, c(N.PD, P, F.Max)) 

	# Set up a function for calcualting percentile rank
	perc.rank <- function(vec, x){length(vec[vec<=x])/length(vec)}

	#############################################################################################
	# Generate parallel datasets assuming k-1 factors underlying P variables
	# Compare eigenvalues for observed dataset and parallel datasets
	#############################################################################################

	# Start with testing the first factor
	T.NFactor=1 

	while(T.NFactor<=F.Max){ # Loop from the first factor to the max number of factors

		# Testing the first factor
		if(T.NFactor==1){

			# Generate parallel datasets assuming zero factor underlying data
			for(which.pd in 1:N.PD){
				para.data <- matrix(rnorm(N*P), nrow=N, ncol=P)
				eigs.para.array[which.pd, ,T.NFactor] <- fa(para.data, nfactors=1, rotate="none", SMC=TRUE, max.iter=1, fm="pa")$values
			}

			eigs.para <- as.numeric(quantile(eigs.para.array[ ,T.NFactor,T.NFactor], percentile))
			eigs.data.rank <- perc.rank(eigs.para.array[ ,T.NFactor,T.NFactor], eigs.data[T.NFactor])
			if(eigs.data[T.NFactor]>eigs.para[T.NFactor]){T.NFactor <- T.NFactor+1} else break

		}


		# Testing the T.NFactor-th factor
		if(T.NFactor>1 & T.NFactor<=F.Max){

			# Compute the factor loadings assuming (T.NFactor-1) factors underlying the observed data
			# These factor loadings will be used to generate parallel datasets
			lambda.data <- matrix(as.numeric(fa(data, nfactors=T.NFactor-1, rotate="none", SMC=TRUE, max.iter=1, fm="pa")$loadings), 
						nrow=P, ncol=T.NFactor-1, byrow=FALSE)

			# Generate parallel datasets assuming (T.NFactor-1) factor underlying data
			for(which.pd in 1:N.PD){
				para.data <- gen.pa.data(P, lambda=lambda.data, N, N.Factors=T.NFactor-1) 
				eigs.para.array[which.pd, ,T.NFactor] <- fa(para.data, nfactors=1, rotate="none", SMC=TRUE, max.iter=1, fm="pa")$values
			}

			eigs.para <- c(eigs.para, as.numeric(quantile(eigs.para.array[ ,T.NFactor,T.NFactor], percentile)))
			eigs.data.rank <- c(eigs.data.rank, perc.rank(eigs.para.array[ ,T.NFactor,T.NFactor], eigs.data[T.NFactor]))
			if(eigs.data[T.NFactor]>eigs.para[T.NFactor]){T.NFactor <- T.NFactor+1} else break

		}

	} # End the while loop
	
	n.factors = T.NFactor-1

	# Print out results
	cat("\n")
	cat("\n")
	cat("\n")
	cat("The estimated number of factors is ", n.factors, ".\n", sep="")
	cat("\n")
	cat(c("The eigenvalues extracted for each sequentially extracted factor based on the observed dataset are: \n", 
	      format(round(eigs.data,2), nsmall=2), "\n"), collapse=" ")
	cat("\n")
	cat(c(paste("The ", percentile*100, "%ile eigenvalues for each sequentially extracted factor based on the parallel datasets are: \n", sep=""), 
	      format(round(eigs.para,2), nsmall=2), "\n"), collapse=" ")
	cat("\n")
	cat(c("The percentile ranks of the eigenvalues based on the observed dataset compared to the eigenvalues based on the parallel datasets are: \n", 
	      paste(format(round(eigs.data.rank*100,0), nsmall=0), "%ile", sep=""), "\n"), collapse=" ")

} # End the R-PA function


###########################################################################################################################################
# Input values for the R-PA:
#
# We have included a researcher's correlation matrix to illustrate how to conduct a R-PA.  This sample correlation matrix was generated by 
# a factor model with three factors underlying nine indicators.  Substitute your own correlation matrix (researcher's correlation matrix) 
# for the one below to obtain the estimated number of factors for your dataset. Note that in reading in your correlation matrix, you must specify 
# "nrow" and "ncol", which is equal to the number of variables in you data. In the last statement, the sample size for your researcher's 
# correlation matrix is specified. For our illustration, the sample size was specified as 400 (N=400).  In place of 400, isert the sample size for your 
# researcher's correlation matrix.
#
# The program assumes no missing data.  There have been no manuscripts on how to conduct R-PA in the presence of missing data. 
#
# Currently the number of parallel datasets is set at 100 (N.PD=100) and the %ile at .95 (percentile=.95) as the default option.
# Thus these two arguments were not specified in the statement.  These are the values recommended in the articles on R-PA.  
# You can input a different number of parallel datasets and a different %ile as a cutoff for the eigenvalues for the parallel 
# datasets to override the default values.
#
# Output values for the R-PA:
# 
# R-PA returns the following output: 
# 1) the estimated number of factors; 
# 2) the eigenvalues for each sequentially extracted factor based on the researcher's reduced correlation matrix;
# 3) the 95%ile (or other pencentile values specified by the researcher) eigenvalues for each sequentially extracted factor based on
#    the parallel datasets. 
# 4) the percentile ranks of the eigenvalues based on the researcher's reduced correlation matrix compared to the eigenvalues based 
# on the parallel datasets
###########################################################################################################################################

