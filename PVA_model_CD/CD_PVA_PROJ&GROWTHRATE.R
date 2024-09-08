############################################
## Common dolphin population projection model
############################################
## Population projection & growth rate 

source("CD_PVA_MATRIX.R")

#
##      Projection setup
#
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Initial Population Setup
start.pop <- total_population # Initialize the starting population size at total population size calculated previously individuals
# calculate the initial population vector based on the total population
init.vec.total.pop <- total_population*stable.stage.dist(popmat)

## Time Frame for Projection
yr.now <- 2024 #Define the current year for the start of the projection
yr.end <- 2100 # Define the end year for the projection, setting a long-term analysis

##Time Projection Setup
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Purpose: To establish the timeframe over which the population will be projected.
# Implementation: Calculates the total number of years (t) for the projection based on the difference between the start (yr.now) and end (yr.end) years, and creates a sequence of those years (yrs).
# Define the number of years for the projection based on the start and end years.
t <- (yr.end - yr.now) # Total years for projection
yrs <- seq(yr.now,yr.end,1) # Vector of years from start to end

## Density-Feedback Function for Population Regulation:
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Purpose: To model the effect of population density on fertility rates, simulating natural density-dependent mechanisms that regulate population growth.
# Implementation: Defines a carrying capacity (K) as a multiple of the starting population and creates a fertility reduction factor that decreases as the population size approaches this capacity. A linear model (lm) is used to fit these changes, which will adjust fertility rates in the projection.

K <- 3*start.pop # Define a carrying capacity as three times the starting population.

# Create a fertility reduction factor that scales with population size, aiming to simulate density-dependent fertility.
fert.min.mult <- 0.7 # Minimum multiplier for fertility reduction
i.K <- start.pop:K # Range from current population size to carrying capacity
# Calculate a multiplier for fertility that decreases as population size approaches carrying capacity.
fert.mult <- 1 - cumsum(rep((1 - fert.min.mult)/length(i.K), length(i.K)))
# Fit a linear model to describe how the fertility multiplier changes with population size.
i.K.fit <- lm(fert.mult ~ i.K)

plot(i.K, fert.mult, pch=10, xlab = "Population Size", ylab = "Fertility Multiplier", type="p", col = "blue")

# Add the linear model fit to the plot
abline(i.K.fit, col = "red")

# Display the summary of the linear model
summary(i.K.fit)


# STEP #7: Population Matrix + Projection iteration Setup + Storage
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Population matrix setup
# Purpose: To initialize the structure that will store population counts over time for simulation.
# Implementation: Creates a matrix (n.mat) where each column represents a year and each row an age class. The initial population distribution is set in the first column.
# Initialize a population matrix with zeros, where rows represent stages and columns represent time steps.
n.mat <- matrix(0,nrow=stages,ncol=(t+1))
n.mat[,1] <- init.vec.total.pop  # Set the initial population vector as the first column
popmat <- popmat.orig # Reset the population matrix to its original state for the simulation

## Stochastic Projection Iterations
# Purpose: To perform multiple simulations of population projections under varying conditions to capture uncertainty and variability in population dynamics.
# Implementation: Specifies the number of iterations (iter) for the stochastic projections and sets up for detailed simulation processes in the subsequent lines 
# Define the number of iterations for the stochastic projection to simulate variability.
iter <- 100000 # Number of stochastic iterations
itdiv <- iter/10000 # Interval for progress reporting

# Set up storage matrices and vectors for simulation results
# Create a matrix to store total population summaries for each iteration and each time step
# This matrix has 'iter' rows (one for each iteration) and 't+1' columns (one for each time point including the initial state)
n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1))

# Create a matrix to store the rate of population change ('r') for each year of each iteration
# This matrix has 'iter' rows and 't' columns (one for each year of the projection)
r.mat <- matrix(data = 0, nrow = iter, ncol = t)

# Initialize vectors to store the weighted mean age and total length (TL) for the population at the end of each iteration
# Both vectors are set to zero and have a length equal to the number of iterations
age.wm.vec <- TL.wm.vec <- rep(0, iter)

# STEP #8: stochastic simulations
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Nested Loops: There are two nested loops; the outer loop (e) iterates over the number of simulation runs (iter), and the inner loop (i) steps through each year of the projection period (t).
# Sampling Parameters: Within each year of a simulation run, fertility rates and survival probabilities are sampled using statistical distributions (normal and beta distributions, respectively). This introduces random variations reflecting uncertainty in parameter estimates.
# Catastrophic Events: The model randomly decides whether a catastrophic event occurs each year, which would halve the survival rates, simulating sudden environmental impacts or other disturbances.
# Density Feedback Mechanism: The fertility rate is adjusted based on the current population size relative to a predefined carrying capacity (K). This models the density-dependent effects on population growth, where higher densities might lead to reduced fertility due to limited resources or space.

for (e in 1:iter) {
  # Initialize a vector to store annual growth rates for the current iteration
  r.stoch <- rep(0, t)
  
  # Reset the population matrix to its original configuration at the start of each iteration
  popmat <- popmat.orig
  
  # Loop through each year of the simulation
  for (i in 1:t) {
    # This will pause execution here
    
    # Randomly sample the proportion of mature individuals using a complex transformation involving normal distributions
    p.mat.stoch.a <- rnorm(1, fit.logp.summ$coefficients[1], fit.logp.summ$coefficients[4]) / (1 + (age.vec / rnorm(1, fit.logp.summ$coefficients[2], fit.logp.summ$coefficients[5]))^rnorm(1, fit.logp.summ$coefficients[3], fit.logp.summ$coefficients[6]))
    p.mat.stoch.b <- ifelse(p.mat.stoch.a > 1, 1, p.mat.stoch.a)
    p.mat.stoch <- ifelse(p.mat.stoch.b < 0.001, 0, p.mat.stoch.b)
    p.mat.stoch
    
    # Apply density feedback to the fertility rate based on current population size relative to a starting threshold
    fert.multiplier <- ifelse(sum(n.mat[, i]) >= start.pop, as.numeric(coef(i.K.fit)[1] + coef(i.K.fit)[2] * sum(n.mat[, i])), 1)
    
    # Construct a fertility vector for the population matrix using the sampled fertility and survival rates
    f.fert.stoch <- 0.317 * (p.mat.stoch * litt.pred2) * fert.multiplier  # Account for resting period between calves 
    f.fert.stoch
    
    # Simulate survival rates using a beta distribution parameterized by sampled alpha and beta values
    Sx.sd <- 0.05  # can set to any value
    Sx.alpha <- estBetaParams(surv.vec, Sx.sd^2)$alpha
    Sx.beta <- estBetaParams(surv.vec, Sx.sd^2)$beta
    Sx.stoch <- rep(0, stages)
    Sx.alpha <- ifelse(is.nan(Sx.alpha), 0, Sx.alpha)
    Sx.beta <- ifelse(is.nan(Sx.beta), 0, Sx.beta)
    for (x in 1:stages) {
      Sx.stoch[x] <- rbeta(1, Sx.alpha[x], Sx.beta[x])
    }
    
    # Update the first row of the population matrix 'popmat' with the stochastic fertility vector
    popmat[1, ] <- f.fert.stoch * 0.5 #divide by two to account for females only
    
    # Simulates a catastrophic event occurring with a probability 'cat.pr'. The 'rbinom' function is
    # used here to generate a random number from a binomial distribution, where the number of trials
    # is 1 (i.e., either a catastrophe happens or it doesn't), and the probability of success (catastrophe)
    # is given by 'cat.pr'. This setup allows for a simple yes/no outcome each year.
    catastrophe <- rbinom(1, 1, cat.pr)
    # If a catastrophe occurs (catastrophe equals 1), the survival probabilities for each stage
    # are halved. This reflects a significant detrimental impact on the population, perhaps due to
    # environmental disasters, disease outbreaks, or other sudden changes. The survival probabilities
    # are applied to the main diagonal of the population matrix, excluding the first row, which is
    # dedicated to fertility rates. This operation reduces the survival rate of each stage by 50%,
    # dramatically affecting the population dynamics.
    if (catastrophe == 1) {
      diag(popmat[2:(stages), ]) <- (0.5*Sx.stoch[-stages])}
    # If no catastrophe occurs (catastrophe equals 0), the survival probabilities are applied as
    # originally sampled from the beta distribution, reflecting typical annual survival without
    # external shocks. These probabilities are set to the main diagonal of the population matrix,
    # similarly excluding the first row.
    if (catastrophe == 0) {
      diag(popmat[2:(stages), ]) <- Sx.stoch[-stages]}
    
    # The survival probability for the last stage is set to 0, implying that individuals in this stage
    # do not survive to the next year. This is a modeling choice that can represent the fact that the
    # oldest age class does not contribute to future population growth either due to natural mortality
    # or because they are beyond the reproductive age.
    popmat[stages,stages] <- 0
    
    # Calculate the population vector for the next year by matrix multiplication
    n.mat[,i+1] <- popmat %*% n.mat[,i] 
    
    # Store the logarithm of the ratio of successive total populations to estimate the growth rate
    r.running <- log(sum(n.mat[,i+1], na.rm=T) / sum(n.mat[,i], na.rm=T))
    r.stoch[i] <- ifelse(r.running == -Inf, NA, r.running)
  }
  
  
  # Store the results of this iteration
  r.mat[e, ] <- r.stoch
  n.sums.mat[e, ] <- as.vector(colSums(n.mat))
  
  # Median age & size of final population
  age.wm.vec[e] <- weighted.mean(x = age.vec, w = n.mat[, (t+1)], na.rm = TRUE)
  TL.wm.vec[e] <- 308 / (1 + (3.162162 * exp(-0.146 * age.wm.vec[e])))  # Predict TL from weighted mean age
  
  # Print progress periodically
  if (e %% itdiv == 0) print(e)
}


##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## STEP #9: generation of confidence interval
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# at the end of all iterations, statistical analyses are performed to estimate confidence intervals for the projected population sizes and growth rates.
# Setup plotting parameters to display two plots side by side
par(mfrow=c(1,1))

# Initialize vectors to store the mean, upper, and lower confidence limits of population size
n.up <- n.lo <- n.mn <- rep(0,(t+1))
# Calculate mean and 95% confidence intervals for each year's population size
for (q in 1:(t+1)) {
  n.mn[q] <- mean(n.sums.mat[,q]) # Calculate mean population size for year q
  n.up[q] <- as.vector(quantile(n.sums.mat[,q], probs=0.975)) # 97.5% quantile
  n.lo[q] <- as.vector(quantile(n.sums.mat[,q], probs=0.025)) # 2.5% quantile
}

# Plot mean population size and confidence intervals over time
plot(yrs, n.mn, type="l", xlab="year",ylab="N",xlim=c(yr.now,yr.end),ylim=c(min(n.lo),max(n.up)))
lines(yrs, n.up, lty=2, col="red") # Upper confidence limit
lines(yrs, n.lo, lty=2, col="red") # Lower confidence limit

# plot using classic theme
library(ggplot2)

data_population_projection <- data.frame(yrs = yrs, n.mn = n.mn, n.lo = n.lo, n.up = n.up)

# Create the plot with y-axis starting from 0
population_projection <- ggplot(data_population_projection, aes(x = yrs)) +
  geom_line(aes(y = n.mn), color = "black", size = 0.7) +  # Main line in black
  geom_line(aes(y = n.up), color = "black", size = 0.7, linetype = "dashed") +  # Upper confidence limit with dashed line
  geom_line(aes(y = n.lo), color = "black", size = 0.7, linetype = "dashed") +  # Lower confidence limit with dashed line
  labs(x = "Year", y = "Population Size (N)") +  # Add x and y axis labels
  theme_classic() +  # Apply classic theme
  scale_y_continuous(limits = c(0, max(n.up)),  # Set y limits to include 0
                     breaks = pretty(c(0, max(n.up)), n = 12),  # More y-axis labels
                     expand = expansion(mult = c(0.01, 0.01))) +  # Slightly stretch y-axis
  scale_x_continuous(limits = c(min(yrs), max(yrs)), 
                     breaks = c(seq(2040, max(yrs), by = 20), 2024),  # Keep other tick marks but add 2024
                     expand = expansion(mult = c(0, 0))) +  # Reduce gap between first data point and y-axis
  theme(axis.text = element_text(color = "black", size = 11),  # Adjust axis text size
        axis.title = element_text(color = "black", size = 12),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))  # Increase margin around the plot to prevent labels being cut off

# Display the plot
print(population_projection)


# Calculate mean and 95% confidence intervals for the growth rate 'r' for each year
r.up <- r.lo <- r.mn <- rep(0,(t))
for (q in 1:(t)) {
  r.mn[q] <- mean(r.mat[,q]) # Calculate mean growth rate for year q
  r.up[q] <- as.vector(quantile(r.mat[,q], probs=0.975)) # 97.5% quantile
  r.lo[q] <- as.vector(quantile(r.mat[,q], probs=0.025)) # 2.5% quantile
}
# Plot mean growth rate 'r' and its confidence intervals over time
plot(yrs[-1], r.mn, type="l", xlab="year",ylab="r",xlim=c(yr.now,yr.end),ylim=c(min(r.lo),max(r.up)))
lines(yrs[-1], r.up, lty=2, col="red") # Upper confidence limit
lines(yrs[-1], r.lo, lty=2, col="red") # Lower confidence limit
abline(h=0,lty=3,lwd=2,col="grey") # Reference line at zero growth rate

data_growth_rate_projection <- data.frame(yrs = yrs[-1], r.mn = r.mn, r.lo = r.lo, r.up = r.up)

# Create the plot in ggplot2
growth_rate_projection <- ggplot(data_growth_rate_projection, aes(x = yrs)) +
  geom_line(aes(y = r.mn), color = "black", size = 0.7) +  # Main line for r.mn in black
  geom_line(aes(y = r.up), linetype = "dashed", color = "black", size = 0.7) +  # Upper confidence limit
  geom_line(aes(y = r.lo), linetype = "dashed", color = "black", size = 0.7) +  # Lower confidence limit
  labs(x = "Year", y = "Growth rate (r)") +  # Add x and y axis labels
  theme_classic() +  # Apply classic theme
  scale_y_continuous(limits = c(min(r.lo), max(0.01)),  # Set y limits to fit all data
                     breaks = pretty(c(min(r.lo), max(0.01)), n = 10)) +  # More y-axis labels
  scale_x_continuous(limits = c(yr.now, yr.end), 
                     breaks = c(seq(2040, yr.end, by = 20), 2024), 
                     expand = expansion(mult = c(0, 0))) +  # Reduce space between y-axis and first data point
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.8, color = "red") +  # Reference line at 0
  theme(axis.text = element_text(color = "black", size = 11),  # Adjust axis text size
        axis.title = element_text(color = "black", size = 12), 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))  # Adjust axis title size

# Display the plot
print(growth_rate_projection)


# Restore default plotting parameters
par(mfrow=c(1,1))

