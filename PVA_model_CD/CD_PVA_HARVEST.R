
############################################
## Common dolphin population projection model
############################################
## Harvesting 

source("CD_PVA_MATRIX.R")


######## Harvesting
#calculate probability of being caught based on initial age distribution
init.vec.total.pop <- total_population*stable.stage.dist(popmat)
prob.vec <- init.vec.total.pop / total_population

# Plot the probabilities against the age categories
plot(age.vec, prob.vec, xlab = "age (yrs)", ylab = "prob bycatch",
     type = "l", col = "blue")

# to plot initial age distribution and probabilty on the same graph
plot(age.vec, init.vec.total.pop, xlab = "age (yrs)", ylab = "N", type = "l", col = "blue")

par(new = TRUE)
plot(age.vec, prob.vec, xlab = "", ylab = "", axes = FALSE, type = "l", col = "red", lty = 2)
axis(4, at = pretty(range(prob.vec)), col.axis = "red", col = "red")
mtext("prob bycatch", side = 4, line = 3, col = "red")

sum(prob.vec)



##      Projection setup=
# Stochastic projection simulation
# This loop simulates population changes over time under the specified conditions using the Leslie matrix model.

for (i in 1:t) {
  n.mat[,i+1] <- popmat %*% n.mat[,i]  # Multiply the current population vector by the population matrix to get the next year's population vector.
}

## Create a vector of years for plotting.
yrs <- seq(yr.now,yr.end,1) # Sequence of years from start to end of projection.

# Plot the total population over time
# This plot visualizes how the total population is expected to change from yr.now to yr.end.
plot(yrs,as.vector(colSums(n.mat)),type="l",xlab="year",ylab="N",xlim=c(yr.now,yr.end))
# 'type="l"' creates a line graph, which is helpful for showing trends in population changes over time.


## Time Frame for Projection
yr.now <- 2024 #Define the current year for the start of the projection
yr.end <- 2100 # Define the end year for the projection, setting a long-term analysis

##Time Projection Setup
# Purpose: To establish the timeframe over which the population will be projected.
# Implementation: Calculates the total number of years (t) for the projection based on the difference between the start (yr.now) and end (yr.end) years, and creates a sequence of those years (yrs).
# Define the number of years for the projection based on the start and end years.
t <- (yr.end - yr.now) # Total years for projection
yrs <- seq(yr.now,yr.end,1) # Vector of years from start to end

## Density-Feedback Function for Population Regulation:
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

#Population Matrix + Projection iteration Setup + Storage
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
iter <- 1000 # Number of stochastic iterations
itdiv <- iter/100 # Interval for progress reporting

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

# set harvest loop
# harvest as multiples of age.yr.f, where 0*prop.caught = current harvest, 0.5*prop.caught = 50% more than now, 1*prop.caught = double now
harv.p.vec <- seq(0,1,0.01)
lharv <- length(harv.p.vec)
r.mn.harv <- r.mn.up <- r.mn.lo <- n.mn.harv <- n.mn.up <- n.mn.lo <- n.min.harv <- n.min.up <- n.min.lo <- age.wm.mn <- age.wm.lo <- age.wm.up <- TL.wm.mn <- TL.wm.lo <- TL.wm.up <- rep(0,lharv)


for (h in 1:lharv)  {
  
  # Set storage matrices & vectors
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t + 1))
  r.mat <- matrix(data = 0, nrow = iter, ncol = t)

  
  
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
      popmat[1, ] <- f.fert.stoch * 0.5
      
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
      
      ## Harvest
      n.mat[, i + 1] <- n.mat[, i + 1] - (harv.p.vec[h] * prob.vec * n.mat[, i ]) # Substitute i for 1 in last element for proportional vs. constant harvest, respectively
      n.mat[which(n.mat[, i + 1] < 0), i + 1] <- 0
      
      # Store the logarithm of the ratio of successive total populations to estimate the growth rate
      r.running <- log(sum(n.mat[,i+1], na.rm=T) / sum(n.mat[,i], na.rm=T))
      r.stoch[i] <- ifelse(r.running == -Inf, NA, r.running)
    }
    
    
    # Store the results of this iteration
    r.mat[e, ] <- r.stoch
    n.sums.mat[e, ] <- as.vector(colSums(n.mat))
    
    # Print progress periodically
    if (e %% itdiv == 0) print(e)
  }
  
  
  # 1 generation burn-in
  n.mn <- apply(n.sums.mat[, (gen.l + 1):(t + 1)], MARGIN = 2, mean, na.rm = TRUE)
  n.up <- apply(n.sums.mat[, (gen.l + 1):(t + 1)], MARGIN = 2, quantile, probs = 0.975, na.rm = TRUE)
  n.lo <- apply(n.sums.mat[, (gen.l + 1):(t + 1)], MARGIN = 2, quantile, probs = 0.025, na.rm = TRUE)
  
  plot(yrs[(gen.l + 1):(t + 1)], n.mn, type = "l", xlab = "year", ylab = "N", xlim = c(yrs[(gen.l + 1)], yrs[(t + 1)]), ylim = c(min(n.lo), max(n.up)))
  lines(yrs[(gen.l + 1):(t + 1)], n.up, lty = 2, col = "red")
  lines(yrs[(gen.l + 1):(t + 1)], n.lo, lty = 2, col = "red")

  # r confidence limits with burn-in
  r.mn <- apply(r.mat[, (gen.l + 1):t], MARGIN = 2, mean, na.rm = TRUE)
  r.up <- apply(r.mat[, (gen.l + 1):t], MARGIN = 2, quantile, probs = 0.975, na.rm = TRUE)
  r.lo <- apply(r.mat[, (gen.l + 1):t], MARGIN = 2, quantile, probs = 0.025, na.rm = TRUE)
  
  # Plot with burn-in
  plot(yrs[(gen.l + 1):(t)], r.mn, type = "l", xlab = "year", ylab = "r", xlim = c(yrs[(gen.l)], yrs[(t)]), ylim = c(min(r.lo), max(r.up)))
  lines(yrs[(gen.l + 1):(t)], r.up, lty = 2, col = "red")
  lines(yrs[(gen.l + 1):(t)], r.lo, lty = 2, col = "red")
  abline(h = 0, lty = 3, lwd = 2, col = "grey")
  
  # Store values per h iteration
  r.mn.harv[h] <- median(r.mn, na.rm = TRUE)
  r.mn.up[h] <- quantile(r.mn, probs = 0.975, na.rm = TRUE)
  r.mn.lo[h] <- quantile(r.mn, probs = 0.025, na.rm = TRUE)
  n.mn.harv[h] <- median(n.mn, na.rm = TRUE)
  n.mn.up[h] <- quantile(n.mn, probs = 0.975, na.rm = TRUE)
  n.mn.lo[h] <- quantile(n.mn, probs = 0.025, na.rm = TRUE)
  n.min.harv[h] <- median(n.lo, na.rm = TRUE)
  n.min.up[h] <- quantile(n.lo, probs = 0.975, na.rm = TRUE)
  n.min.lo[h] <- quantile(n.lo, probs = 0.025, na.rm = TRUE)
  TL.wm.mn[h] <- mean(TL.wm.vec, na.rm = TRUE)
  TL.wm.up[h] <- quantile(TL.wm.vec, probs = 0.975, na.rm = TRUE)
  TL.wm.lo[h] <- quantile(TL.wm.vec, probs = 0.025, na.rm = TRUE)
  
  print(paste("harvest increment = ", h, sep = ""))
 }

# Visualizing the results: plotting population trends and changes under different model projections and harvest scenarios.
# Plotting mean population sizes over time.
par(mfrow = c(1, 1))
plot(harv.p.vec, r.mn.harv, type = "l", lwd = 2, ylim = c(-.1, .03), xlab = "bycatch multiplier on current rate", ylab = "mean long-term r")
abline(h = 0, lty = 2, col = "red")
lines(harv.p.vec, r.mn.lo, lty = 2, lwd = 1)
lines(harv.p.vec, r.mn.up, lty = 2, lwd = 1)

plot(harv.p.vec, r.mn.harv, type = "l", lwd = 2, ylim = range(r.mn.harv, r.mn.lo, r.mn.up), xlab = "bycatch multiplier on current rate", ylab = "mean long-term r")
abline(h = 0, lty = 2, col = "red")
lines(harv.p.vec, r.mn.lo, lty = 2, lwd = 1)
lines(harv.p.vec, r.mn.up, lty = 2, lwd = 1)

plot(harv.p.vec, n.mn.harv, type = "l", lwd = 2, ylim = c(min(n.mn.lo), max(n.mn.up)), xlab = "bycatch multiplier on current rate", ylab = "mean long-term N")
lines(harv.p.vec, n.mn.lo, lty = 2, lwd = 1)
lines(harv.p.vec, n.mn.up, lty = 2, lwd = 1)

plot(harv.p.vec, n.min.harv, type = "l", lwd = 2, ylim = c(min(n.min.lo), max(n.min.up)), xlab = "bycatch multiplier on current rate", ylab = "mean long-term N min")
lines(harv.p.vec, n.min.lo, lty = 2, lwd = 1)
lines(harv.p.vec, n.min.up, lty = 2, lwd = 1)
par(mfrow = c(1, 1))



