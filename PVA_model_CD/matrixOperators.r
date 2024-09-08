# Matrix operators for population models

## maximum lambda function
# Purpose: Calculates the maximum eigenvalue (lambda) of a Leslie matrix, which represents the dominant growth rate of the population.
# Usage: Used in your main script to assess the growth or decline of the population over a single time unit. It helps determine if the population is increasing, stable, or decreasing.
max.lambda <- function(x) Re((eigen(x)$values)[1]) ## where 'x' is a Leslie matrix

## Maximum r function
# Purpose: Calculates the natural logarithm of the dominant eigenvalue to provide the intrinsic rate of increase (r) of the population.
# Usage: Complements the max.lambda function by offering another perspective on population growth, useful in demographic analysis to express growth rate on a logarithmic scale.
max.r <- function(x) log(Re((eigen(x)$values)[1])) ## where 'x' is a Leslie matrix

## Stable stage distribution
# Purpose: Determines the stable stage distribution of a population from a Leslie matrix, indicating the long-term expected proportion of individuals in each age class.
# Usage: Essential for initializing the population vector (init.vec) in your script, which sets up the initial distribution of individuals across different age classes based on long-term expectations.
stable.stage.dist <- function(x) ((x %*% (Re((eigen(x)$vectors)[,1])))/(sum((x %*% (Re((eigen(x)$vectors)[,1]))))))[,1]

## Generation length function
# Purpose: Calculates the reproductive value (R0) for each age class within a Leslie matrix, reflecting the contribution of each age class to future generations.
# Usage: Helps in understanding the relative importance of different age classes in contributing to the population's future growth, crucial for conservation and management decisions.
# R.val Function Definition: Calculates the reproductive value (R0) for a given Leslie matrix.
R.val <- function(X,age.max) ## reproductive value (R0) where X = Leslie matrix; age.max = maximum age of females
{		
    # Transition Matrix (T): Extracts and modifies the Leslie matrix to focus only on survival transitions.
    #Purpose: This matrix represents the probabilities of surviving and transitioning from one age class to the next, excluding fertility.
    #Implementation: The function takes the first age.max rows and columns from the Leslie matrix X to focus on the transition probabilities, setting the first row to zero since it typically contains fertility rates, not survival probabilities.
		T <- X[1:age.max,1:age.max] # Copy the relevant portion of the Leslie matrix.
		T[1,1:(age.max)] <- 0 # Zero out the first row, typically containing fertility values, to isolate survival transitions.

		# Fertility Matrix (F): Isolates the fertility rates from the Leslie matrix.
		F <- X[1:age.max,1:age.max] # Copy the transition matrix setup for consistency.
		diag(F[2:age.max,1:(age.max-1)]) <- 0 # Zero out all but the first row to focus solely on fertility contributions.

		# Identity Matrix (I): Creates an identity matrix needed for matrix inversion.
		I <- matrix(data<-0,nrow<-age.max,ncol<-age.max)
		diag(I) <- 1  # Shortcut to create an identity matrix directly using diag function.

		# Fundamental Matrix (N.fund): Calculates the expected lifetime presence in each class.
		library(MASS)   # Ensure the MASS library is loaded for matrix operations.
		N.fund <- ginv(I-T) # Use the generalized inverse to calculate the fundamental matrix, which predicts lifetime transitions.

		# Reproductive Matrix (R): Combines fertility and survival to estimate overall reproductive output.
		R <- F %*% N.fund # Multiply the fertility matrix by the fundamental matrix to form the reproductive matrix.

		# Calculate R0: Finds the dominant eigenvalue of the reproductive matrix, representing the mean number of offspring produced per individual.
		# The reproductive value R0 calculated by this function is a fundamental demographic parameter that indicates whether the population is expected to grow, decline, or remain stable in the long term (without immigration or emigration). 
		# A value of R0 > 1 suggests population growth, R0 = 1 indicates stability, and	R0 < 1 signals potential decline.
		R0 <- Re((eigen(R)$values)[1]) # Extract the real part of the dominant eigenvalue to get R0, ensuring it is a real number.
		
		# Output the result: Print the calculated reproductive value to provide insight into the potential population growth.
		print("number of female offspring produced per female during its lifetime")
		print("_________________________________________________________________")
		print(R0)

}

## Mean generation time function
#Purpose: Computes the mean generation length using the Leslie matrix, which is the average time taken for one generation to replace itself with the next.
#Usage: Offers insights into the life expectancy and turnover rate of the population, providing a temporal scale for population dynamics analysis.
G.val <- function (X,age.max) ## where X is a Leslie Matrix
{	
		G <- (log(R.val(X,age.max)))/(log(Re((eigen(X)$values)[1])))
		print("mean generation time")
		print("____________________")
		print(G)
}
