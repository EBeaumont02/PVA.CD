
######## Harvesting
#calculate probability of being caught based on initial age distribution
total_population <- sum(init.vec)
prob.vec <- init.vec / total_population

# Plot the probabilities against the age categories
plot(age.vec, prob.vec, xlab = "age (yrs)", ylab = "prob bycatch",
     type = "l", col = "blue")

# to plot inital age distiution and probabilty on the same graph
plot(age.vec, init.vec, xlab = "age (yrs)", ylab = "N", type = "l", col = "blue")

par(new = TRUE)
plot(age.vec, prob.vec, xlab = "", ylab = "", axes = FALSE, type = "l", col = "red", lty = 2)
axis(4, at = pretty(range(prob.vec)), col.axis = "red", col = "red")
mtext("prob bycatch", side = 4, line = 3, col = "red")

sum(prob.vec)


# set harvest loop
# harvest as multiples of age.yr.f, where 0*prop.caught = current harvest, 0.5*prop.caught = 50% more than now, 1*prop.caught = double now
harv.p.vec <- seq(0,1,0.01)
lharv <- length(harv.p.vec)
r.mn.harv <- r.mn.up <- r.mn.lo <- n.mn.harv <- n.mn.up <- n.mn.lo <- n.min.harv <- n.min.up <- n.min.lo <- age.wm.mn <- age.wm.lo <- age.wm.up <- TL.wm.mn <- TL.wm.lo <- TL.wm.up <- rep(0,lharv)


for (h in 1:lharv) {
  
  # Set storage matrices & vectors
  n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t + 1))
  r.mat <- matrix(data = 0, nrow = iter, ncol = t)
  age.wm.vec <- TL.wm.vec <- rep(0, iter)
  
  
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
      popmat[1, ] <- f.fert.stoch
      
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
    
    # Median age & size of final population
    age.wm.vec[e] <- weighted.mean(x = age.vec, w = n.mat[, (t+1)], na.rm = TRUE)
    TL.wm.vec[e] <- 308 / (1 + (3.162162 * exp(-0.146 * age.wm.vec[e])))  # Predict TL from weighted mean age
    
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
  age.wm.mn[h] <- mean(age.wm.vec, na.rm = TRUE)
  age.wm.up[h] <- quantile(age.wm.vec, probs = 0.975, na.rm = TRUE)
  age.wm.lo[h] <- quantile(age.wm.vec, probs = 0.025, na.rm = TRUE)
  TL.wm.mn[h] <- mean(TL.wm.vec, na.rm = TRUE)
  TL.wm.up[h] <- quantile(TL.wm.vec, probs = 0.975, na.rm = TRUE)
  TL.wm.lo[h] <- quantile(TL.wm.vec, probs = 0.025, na.rm = TRUE)
  
  print(paste("harvest increment = ", h, sep = ""))
}

# Visualizing the results: plotting population trends and changes under different model projections and harvest scenarios.
# Plotting mean population sizes over time.
par(mfrow = c(3, 1))
plot(harv.p.vec, r.mn.harv, type = "l", lwd = 2, ylim = c(-.1, .03), xlab = "", ylab = "mean long-term r")
abline(h = 0, lty = 2, col = "red")
lines(harv.p.vec, r.mn.lo, lty = 2, lwd = 1)
lines(harv.p.vec, r.mn.up, lty = 2, lwd = 1)

plot(harv.p.vec, n.mn.harv, type = "l", lwd = 2, ylim = c(min(n.mn.lo), max(n.mn.up)), xlab = "", ylab = "mean long-term N")
lines(harv.p.vec, n.mn.lo, lty = 2, lwd = 1)
lines(harv.p.vec, n.mn.up, lty = 2, lwd = 1)

plot(harv.p.vec, n.min.harv, type = "l", lwd = 2, ylim = c(min(n.min.lo), max(n.min.up)), xlab = "harvest multiplier on current rate", ylab = "mean long-term N min")
lines(harv.p.vec, n.min.lo, lty = 2, lwd = 1)
lines(harv.p.vec, n.min.up, lty = 2, lwd = 1)
par(mfrow = c(1, 1))

par(mfrow = c(2, 1))
plot(harv.p.vec, age.wm.mn, type = "l", lwd = 2, ylim = c(min(age.wm.lo), max(age.wm.up)), ylab = "wmn age (yrs) of final pop", xlab = "harvest multiplier")
lines(harv.p.vec, age.wm.lo, lty = 2, lwd = 1)
lines(harv.p.vec, age.wm.up, lty = 2, lwd = 1)

  