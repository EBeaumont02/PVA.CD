############################################
## Common dolphin population projection model
############################################
## Stochastic Simulations - min founding pop

source("CD_PVA_MATRIX.R")
library(foreach)
library(doParallel)
#library(progressr)

##      Projection setup for stochastic simulations
#
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Initial Population Setup
start.pop <- 30000 # half of estimated total pop to represent females only, assuming 1:1 sex ratio # change this to change mvp 
init.vec <- start.pop*stable.stage.dist(popmat.orig) 

# providing a clear visual representation of how the initial population is distributed across age classes.
plot(age.vec,init.vec.total.pop,xlab="age (yrs)", ylab="N", type="l")

## Time Frame for Projection
yr.now <- 2024 #Define the current year for the start of the projection
yr.end <- 2100 # Define the end year for the projection, setting a long-term analysis

##Time Projection Setup
t <- (yr.end - yr.now) # Total years for projection
yrs <- seq(yr.now,yr.end,1) # Vector of years from start to end

## Density-Feedback Function for Population Regulation:
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


###
# Initialize vectors to store the weighted mean age and total length (TL) for the population at the end of each iteration
# Both vectors are set to zero and have a length equal to the number of iterations
age.wm.vec <- TL.wm.vec <- rep(0, iter)

# fertility errors
m.sd.vec <- c(rep(0.35, maxlong+1)) #mean and standard deviations vector, juvenile and adult fertility 

# survival errors
s.sd.vec <- c(rep(0.15, maxlong+1)) #mean and standard deviations vector, juvenile and adult survival

# Define population vector
pop.found <- start.pop 
init.vec <- ssd * pop.found #initial population vector
plot(0:maxlong,ssd,pch=19,type="b")

yr.vec <- seq(yr.now,yr.end) # year vector


# Decrement founding population size and determine MVP & Pr(Qext) -------------------------------------------------

Qthresh <- 500 # quasiextinction threshold. Set at 25 (25f Ne to avoid inbreeding depression), then at whatever cons. managers want 

# sequence vector for founding N
pop.found.vec <- seq(pop.found, 0, -5) # change the increments down to get a smoother line



# Population Matrix + Projection iteration Setup + Storage
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Initialize vectors to store the weighted mean age and total length (TL) for the population at the end of each iteration
# Both vectors are set to zero and have a length equal to the number of iterations
age.wm.vec <- TL.wm.vec <- rep(0, iter)



par(mfrow=c(1,1))

# Set iterations ----------------------------------------------------------
age.max <- maxlong

iter <- 10000 # iterations to run for each projection loop change to 10,000 for final run
itdiv <- iter/1000 #final model rate at iter/1000

nm<-length(pop.found.vec)
#nm<-20

# p loop storage
PrQExt <- minNmd <- minNlo <- minNup <- rep(NA, nm)
s.stor.mat <- matrix(data = NA, nrow = iter, ncol = t)
fert.stor.mat <- matrix(data = NA, nrow = iter, ncol = t)

# Set up parallel backend
num_cores <- detectCores() - 1 # Use all but one core
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Record the start time
start_time <- Sys.time()

# Activate the progress handlers globally
#handlers(global = TRUE)

# Use with_progress to manage the progress bar
#with_progress({
  
  # Create a progressor with the number of steps
#  progress <- progressor(steps = nm)  # nm is the number of iterations
  
outputs <- foreach(p = 1:nm, .combine = 'cbind', .packages = c('base', 'stats')) %dopar% {
  
#for (p in 1:nm) {

  #progress()  # Increment progress bar each iteration
  
  ## initial population vector
  popmat <- popmat.orig
  init.vec <- ssd * pop.found.vec[p] #initial population vector
  
  #print(pop.found.vec[p])
  print(p)
  
  ## stochastic projection with density feedback
  ## set storage matrices & vectors
  n.sums.mat <- qExt.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
  
  # Create a matrix to store the rate of population change ('r') for each year of each iteration
  r.mat <- matrix(data = 0, nrow = iter, ncol = t)
  
  # Initialize a vector to store annual growth rates for the current iteration
  r.stoch <- rep(0, t)
  
  for (e in 1:iter) {
    popmat <- popmat.orig # Reset the population matrix to its original configuration at the start of each iteration
    n.mat <- matrix(0, nrow=age.max+1,ncol=(t+1))
    n.mat[,1] <- init.vec
    
    # Loop through each year of the simulation
    for (i in 1:t) {
      
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
      }  #end x loop
      
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
      
      # print(fert.stoch)
      totN.i <- sum(n.mat[,i], na.rm=T)

      popmat[1,] <- f.fert.stoch * 0.5
      n.mat[,i+1] <- popmat %*% n.mat[,i]
      
      # Store the logarithm of the ratio of successive total populations to estimate the growth rate
      r.running <- log(sum(n.mat[,i+1], na.rm=T) / sum(n.mat[,i], na.rm=T))
      r.stoch[i] <- ifelse(r.running == -Inf, NA, r.running)
      
    } # end i loop
    
    
    # Store the results of this iteration
    r.mat[e, ] <- r.stoch
    n.sums.mat[e,] <- ((as.vector(colSums(n.mat, na.rm=T))))
    qExt.mat[e,] <- ifelse(n.sums.mat[e,] < Qthresh, 1, 0)
    
    if (e %% itdiv==0) print(e)
    
  } # end e loop
  
  n.md <- apply(n.sums.mat, MARGIN=2, median, na.rm=T) # median over all iterations
  n.mn <- apply(n.sums.mat, MARGIN=2, mean, na.rm=T) 
  n.up <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
  n.lo <- apply(n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations
  
  plot(yr.vec,n.md,type="l")
  Qmax <- apply(qExt.mat, MARGIN=1, max, na.rm=T)
  minN <- apply(n.sums.mat, MARGIN=1, min, na.rm=T) # takes the min from n.sums.mat
  
  
  # Return required outputs for this iteration
  list(
    IniPop = pop.found.vec[p],
    PrQExt = sum(Qmax) / iter,
    minNmd = median(minN),
    minNlo = quantile(minN, probs=0.025),
    minNup = quantile(minN, probs=0.975)
  )
  #minNmd[p] <- median(minN)
  #minNlo[p] <- quantile(minN, probs=0.025)
  #minNup[p] <- quantile(minN, probs=0.975)
  #PrQExt[p] <- sum(Qmax)/iter
  
  #plot(yr.vec, n.mn, type="l", xlab="year",ylab="N", xlim=c(yr.now,yr.end),ylim=c(min(n.lo),max(n.up)))
  #lines(yr.vec, n.up, lty=2, col="red") # Upper confidence limit
  #lines(yr.vec, n.lo, lty=2, col="red") # Lower confidence limit
  
  
  #Qmax <- apply(qExt.mat, MARGIN=1, max, na.rm=T)
  #minN <- apply(n.sums.mat, MARGIN=1, min, na.rm=T) # takes the min from n.sums.mat
  #minNmd[p] <- median(minN)
  #minNlo[p] <- quantile(minN, probs=0.025)
  #minNup[p] <- quantile(minN, probs=0.975)
  #PrQExt[p] <- sum(Qmax)/iter
  
  #par(mfrow=c(1,1))
  
  
  # Calculate mean and 95% confidence intervals for each year's population size
#  for (q in 1:(t+1)) {
#    n.mn[q] <- mean(n.sums.mat[,q]) # Calculate mean population size for year q
#    n.up[q] <- as.vector(quantile(n.sums.mat[,q], probs=0.975)) # 97.5% quantile
#    n.lo[q] <- as.vector(quantile(n.sums.mat[,q], probs=0.025)) # 2.5% quantile
    
    ######### Growth rate
    # Calculate mean and 95% confidence intervals for the growth rate 'r' for each year
#    r.mn <- apply(r.mat, MARGIN=2, mean, na.rm=T) 
#    r.up <- apply(r.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
#    r.lo <- apply(r.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
    
#  }
  # Plot mean population size and confidence intervals over time
  #plot(yrs, n.mn, type="l", xlab="year",ylab="N",xlim=c(yr.now,yr.end),ylim=c(min(n.lo),max(n.up)))
  #lines(yrs, n.up, lty=2, col="red") # Upper confidence limit
  #lines(yrs, n.lo, lty=2, col="red") # Lower confidence limit
  
  # Plot mean growth rate 'r' and its confidence intervals over time
  #plot(yrs[-1], r.mn, type="l", xlab="year",ylab="r",xlim=c(yr.now,yr.end),ylim=c(min(r.lo),max(r.up)))
  #lines(yrs[-1], r.up, lty=2, col="red") # Upper confidence limit
  #lines(yrs[-1], r.lo, lty=2, col="red") # Lower confidence limit
  #abline(h=0,lty=3,lwd=2,col="grey") # Reference line at zero growth rate
  
  
  
  #print(" ", quote = FALSE)
  #print("*********************************************", quote = FALSE)
  #print(paste("founding N = ", pop.found.vec[p]), quote = FALSE)
  #print("*********************************************", quote = FALSE)
  #print(" ", quote = FALSE)
} # end p loop

#})

# Record the end time
end_time <- Sys.time()

# Calculate and print the duration
duration <- end_time - start_time
print(paste("Total execution time:", duration, units(duration)))

# Stop the cluster after parallel processing
stopCluster(cl)
# Check the outputs
outputs

#reformating the results into a dataframe
output_matrix<-do.call(rbind,outputs)
nout<-length(output_matrix)
out_rp<-seq(1,nout,5)
nout_rp<-length(out_rp)
output_df<-as.data.frame(matrix(NA,nrow = nout_rp-1,ncol = 5))

output_df[1,]<-output_matrix[out_rp[1]:out_rp[2]-1]
for (i in 2:(nout_rp-1)){
  test<-output_matrix[out_rp[i]:out_rp[i+1]-1]
  output_df[i,]<-test[2:6]
  rm(test)
}

colnames(output_df)<-c("inipop","PrQExt","minNmd",'minNmlo','minMmup')
head(output_df)



plot(output_df$inipop, output_df$PrQExt, type="l", xlab="founding N", ylab="Pr(quasi-extinction)")
#plot(output_df$inipop, output_df$minNmd, type="l", xlab="founding N", ylab="minimum N", ylim=c(min(output_df$minNmlo[is.finite(output_df$minNmlo)],na.rm = T), max(output_df$minMmup[is.finite(output_df$minMmup)],na.rm = T)))
#lines(output_df$inipop, output_df$minNmlo, lty=2, col="red")
#lines(output_df$inipop, output_df$minMmup, lty=2, col="red")
#par(mfrow=c(1,1))

#alldata <- data.frame(pop.found.vec, minNmd, minNlo, minNup, PrQExt)
#saving the results
write.csv(output_df, file = "all_data(FS).csv", row.names = FALSE)

library(ggplot2)

# final plots -------------------------------------------------------------

## plot in ggplot 
output_df<- read.table("all_data(FS).csv",header=T,sep=",")

library(dplyr)
library(ggplot2)

first_below_threshold_x <- output_df %>%
  filter(PrQExt <= 0.1) %>%
  slice_min(order_by = inipop) %>%
  pull(inipop)

Prob.ext.1 <- ggplot(data = output_df, mapping = aes(x=inipop, y=PrQExt)) +
  geom_line(color = "black", linewidth = 1.05) +
  geom_hline(yintercept = 0.1, linetype = 2, color = "red", linewidth = 1.05) +
  geom_vline(xintercept = first_below_threshold_x, linetype = 2, color = "blue",  linewidth = 1.05) +
  scale_x_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 2500), expand = c(0, 0.7)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, by = 0.1), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),   # Remove grid lines
        panel.border = element_blank(), # Remove all borders
        axis.line.x = element_line(linewidth = 0.8),   # Add x-axis line
        axis.line.y = element_line(linewidth = 0.8),   # Add y-axis line
        axis.text = element_text(size = 14, color = "black"),  # Make axis ticks larger and black
        axis.title.x = element_text(size = 16),  # Increase x-axis label text size
        axis.title.y = element_text(size = 16),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm")) + # Increase y-axis label text size
labs(x = "minimum population size (N)", y = "probability of quasi-extinction")

print(Prob.ext.1)


Prob.ext.minN <- ggplot(data = output_df, mapping = aes(x=inipop, y=minNmd)) +
  geom_line(aes(y=minNmd), color = "black") + 
  geom_line(aes(y=minNlo), color = "red", linetype = 2) + 
  geom_line(aes(y=minNup), color = "red", linetype = 2) +
  scale_x_continuous(limits = c(0,520), breaks = seq(0,500, by = 50),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,500), breaks = seq(0,500, by = 50))+
  theme_bw() +
  labs(x = "Founding N", y = "lowest N")
Prob.ext.minN
ggsave("Prob.ext.minN.png")



######### Growth rate
# Calculate mean and 95% confidence intervals for the growth rate 'r' for each year
r.mn <- apply(r.mat, MARGIN=2, mean, na.rm=T) 
r.up <- apply(r.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
r.lo <- apply(r.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)

# Plot mean growth rate 'r' and its confidence intervals over time
plot(yrs[-1], r.mn, type="l", xlab="year",ylab="r",xlim=c(yr.now,yr.end),ylim=c(min(r.lo),max(r.up)))
lines(yrs[-1], r.up, lty=2, col="red") # Upper confidence limit
lines(yrs[-1], r.lo, lty=2, col="red") # Lower confidence limit
abline(h=0,lty=3,lwd=2,col="grey") # Reference line at zero growth rate

# Restore default plotting parameters
par(mfrow=c(1,1))

# Calculate and print confidence limits for the weighted mean age of the population at the end of the projection
# weighted mean age confidence limits
age.wm.mn <- mean(age.wm.vec, na.rm=T)
age.wm.lo <- quantile(age.wm.vec, probs=0.025, na.rm=T)
age.wm.up <- quantile(age.wm.vec, probs=0.975, na.rm=T)
print(c(age.wm.lo, age.wm.mn, age.wm.up))
# Calculate and print confidence limits for the weighted mean total length (TL) of the population
TL.wm.mn <- mean(TL.wm.vec, na.rm=T)
TL.wm.lo <- quantile(TL.wm.vec, probs=0.025, na.rm=T)
TL.wm.up <- quantile(TL.wm.vec, probs=0.975, na.rm=T)
print(c(TL.wm.lo, TL.wm.mn, TL.wm.up))
