############################################
## Common dolphin population projection model
############################################
## Sensitivity Analysis



# Install necessary libraries
install.packages(c("dismo", "gbm", "lhs"))

# Load the libraries
library(dismo)
library(gbm)
library(lhs)


# Initial setup
rm(list = ls())
install.packages("bbmle")
library(bbmle)

#Define functions
#function - estBetaParams: function to estimate the parameters of a beta distribution based on mean and variance.
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2   # Calculate alpha using the rearranged mean and variance formulas.
  beta <- alpha * (1 / mu - 1)   # Calculate beta using the estimated alpha and the mean.
  return(params = list(alpha = alpha, beta = beta))   # Return the estimated parameters as a list.
}

#function - siler model: function to estimate cetacean survival based on age at death distirbution
Si.mod <- function(data, par = c(-0.15, 1.10, 0.15, 0.005, 0.15),
                   rm = 0, method = "Nelder-Mead",
                   control = list(fnscale = -1, maxit = 10000)){
  data$age1 <- data[,1] + 1
  optim(par,
        function(par) {
          a1 <- par[1]
          b1 <- par[2]
          a2 <- par[3]
          a3 <- par[4]
          b3 <- par[5]
          if(rm != 0){
            data <- data[-c(1,rm),]
            data[,1] <- data[,1] - rm
            data[,length(data)] <- data[,length(data)] - rm
          }
          n <- NROW(data)
          S.t <- function(t) {
            return(exp(-a1/b1 * (1 - exp(-b1 * t))) *
                     exp(-a2 * t) * exp(a3/b3 * (1 - exp(b3 * t))))
          }
          dif <- S.t(data[1:n,1]) - S.t(data[1:n,length(data)])
          obs <- data[,2]
          lnlk <- as.numeric(crossprod(obs, log(dif)))
          return(lnlk)
        },
        method = method,
        control = control)
}


# Call data of proportion of each age group from stranding data
age.vec <- seq(0, 29, 1)
surv.data<- read.table("surv.data.csv",header=T,sep=",")
surv.vec <- surv.data$surv
data<-data.frame(age.vec,surv.vec)

#input data into siler model
outfred <- Si.mod(surv.data)

#define parameters
a1 <- outfred$par[1]
b1 <- outfred$par[2]
a2 <- outfred$par[3]
a3 <- outfred$par[4]
b3 <- outfred$par[5]

fred.sil<-exp(-a1/b1 * (1 - exp(-b1 * age.vec))) * exp(-a2 * age.vec) * exp(a3/b3 * (1 - exp(b3 * age.vec)))

# Extracting the Survival Probability Vector

iter <- length(fred.sil)
surv.vec <- numeric(length(fred.sil) - 1)
for (i in 2:iter) {
  surv.vec[i - 1] <- fred.sil[i] / fred.sil[i - 1]
}

# Set the last stage to 0
surv.vec <- c(surv.vec, 0)


#parameters to test:
# carrying capacity?
# starting population?

#### Define parameter ranges 

# density feedback multiplier on fertility between 0.5 and 0.9 (set in the original model arbitrarily to 0.7)
# original model:
# fert.min.mult  <- 0.7
# calculate modifier and add to model:
fert.min.mult.mod  <- seq(0.5, 0.9, 0.1)

#inflection point for age at maturity between 6 and 9 years (set in the original model to 7.345257 years) 
# original model:
# b.new <-  7.345257 
# calculate modifier and add to model:
b.new.mod <- seq(6, 9, 0.1)

# modifier of the catastrophe numerator from 0.1 to 0.2 (set to 0.14 in the original model)
# original model:
# cat.pr <- 0.14/gen.l
# calculate modifier:
cat.num.mod <- seq(0.1, 0.2, 0.01)
# add to model:
# cat.pr <- cat.num.mod/gen.l

# change in the intensity of catastrophic die-offs as a multiplier on survival of between 0.25 and 0.75 (set to 0.5 in the original model)
# original model:
# catastrophe <- rbinom(1, 1, cat.pr)
# if (catastrophe == 1) {
#  diag(popmat[2:(stages), ]) <- (0.5*Sx.stoch[-stages])}
# if (catastrophe == 0) {
#  diag(popmat[2:(stages), ]) <- Sx.stoch[-stages]}
# calculate modifier:
cat.die.off.mod <- seq(0.25, 0.75, 0.01)
# add to model:
# catastrophe <- rbinom(1, 1, cat.pr)
# if (catastrophe == 1) {
#  diag(popmat[2:(stages), ]) <- (cat.die.off.mod*Sx.stoch[-stages])}
# if (catastrophe == 0) {
#  diag(popmat[2:(stages), ]) <- Sx.stoch[-stages]}

# modifier for calving period from 0.25 to 0.5 set to 0.31746 in the original model (breeding once every 3.15 years) 
# original model:
# f.fert.vec <- 0.31746 * (pred.p.mat5*litt.pred2) 
# calculate modifier:
calving.prd.mod <- seq(0.25, 0.5, 0.01)
# add to model:
# f.fert.vec <- calving.prd.mod * (pred.p.mat5*litt.pred2) 

# modifier of the vector of life-table survival estimates:
# original model:
surv.vec
# calculate modifier and add to model:
surv.vec.mod <- c(surv.vec * 0.95, surv.vec * 1.05)
plot (surv.vec.mod)

# LHS SAMPLING

# define vectors 
fert.min.mult.mod  <- seq(0.5, 0.9, 0.1)
b.new.mod <- seq(6, 9, 0.1)
cat.num.mod <- seq(0.1, 0.2, 0.01)
cat.die.off.mod <- seq(0.25, 0.75, 0.01)
calving.prd.mod <- seq(0.25, 0.5, 0.01)
surv.vec.min <- surv.vec * 0.95
surv.vec.max <- surv.vec * 1.05
surv.vec.diff <- surv.vec.max - surv.vec.min

# Define the number of samples
n_samples <- 100
# Combine ranges into a list
bounds <- list(
  fert.min.mult.mod = fert.min.mult.mod,
  b.new.mod = b.new.mod,
  cat.num.mod = cat.num.mod,
  cat.die.off.mod = cat.die.off.mod,
  calving.prd.mod = calving.prd.mod,
  surv.vec.mod = surv.vec.mod
)


# Create the LHS design matrix
lhs_matrix <- lhs::randomLHS(n_samples, length(bounds))

# Function to sample from sequences based on LHS matrix
sample_from_sequences <- function(lhs_matrix, bounds) {
  n_samples <- nrow(lhs_matrix)
  scaled_matrix <- matrix(nrow = n_samples, ncol = length(bounds))
  for (i in seq_along(bounds)) {
    seq_length <- length(bounds[[i]])
    indices <- round(lhs_matrix[, i] * (seq_length - 1)) + 1
    scaled_matrix[, i] <- bounds[[i]][indices]
  }
  return(scaled_matrix)
}

# Apply the sampling function
scaled_matrix <- sample_from_sequences(lhs_matrix, bounds)

# Convert to data frame for easier handling
lhs_df <- as.data.frame(scaled_matrix)
colnames(lhs_df) <- names(bounds)

# Inspect the resulting matrix
print(summary(lhs_df))
print(apply(lhs_df, 2, range))

####MODEL

source("matrixOperators.r") 

# Number of iterations (set this to match the number of rows in your lhs_df)
num_iterations <- 100 # nrow(lhs_df)

# Initialize a list to store the results of each model run
results_pop_size <- results_pop_rate <- vector("list", num_iterations)

# Loop over each set of parameters in the lhs_df matrix
for (r in 1:num_iterations) {

  # Extract the ith set of parameters from lhs_df
  params <- lhs_df[r, ]  # params is a vector of parameters
  
    fert.min.mult <- params[,1]
    b.new <- params[,2]
    cat.num <- params[,3]
    cat.die.off <- params[,4]
    calving.prd <- params[,5]
    surv.temp <- params[,6]*surv.vec.diff
    surv.vec <- surv.vec.min + surv.temp
      
      #defining parameters
      maxlong <- 29 # Set the maximum longevity in the population model at each forecast time step. This value is used to define the upper age limit.
      
      
      # Construct an age vector from 0 to the maximum longevity (maxlong)
      age.vec <- seq(0, maxlong, 1)
      
      # Store the length of the age vector for future use in matrix dimensions and loops.
      lage <- length(age.vec)
      
      # Define the number of life stages based on the length of the age vector.
      stages <- lage
      
      # Initializing the population matrix for modeling.
      # Initially, all elements of the matrix are set to zero, indicating no individuals in any age class.
      popmat <- matrix(0,nrow=stages,ncol=stages)
      colnames(popmat) <- age.vec[1:stages]
      rownames(popmat) <- age.vec[1:stages]
      
      # Maturity data
      # Load maturity data from an external CSV file. 
      prop.mat <- read.table("prop.mat.csv",header=T,sep=",")  # TABLE CREATED USING DATA FROM PALMER # have to change age 8 changed with new data
      
      # Extract the proportion of mature individuals from the dataset. 
      f.prop.mature <- prop.mat$prop.mat
      
      # DEFINE AGE VECTOR (NO ADJUSTEMENT NEEDED)
      f.age.maturity <- prop.mat$Age
      
      # Plot the relationship between age and proportion mature.
      plot(f.age.maturity,f.prop.mature,pch=19,xlab="age (yrs)", ylab="proportion mature")
      
      # Create a data frame from the calculated ages of maturity and the corresponding proportions mature to repare the data for modeling
      mat.dat <- data.frame(f.age.maturity,f.prop.mature)
      
      # Initial parameter values for the logistic power function
      # These are guesses to start the optimization process of fitting the model.
      param.init <- c(1.003150891098860E+00, 1.435062082948057E+01, -3.991451548741554E+01)  
      
      # Fitting a logistic power function to maturity data
      fit.logp <- nls(f.prop.mature ~ a / (1+(f.age.maturity/b)^c), 
                      data = mat.dat,
                      algorithm = "port",# Specifies the optimization algorithm to use.
                      start = c(a = param.init[1], b = param.init[2], c = param.init[3]),# Initial values for parameters.
                      trace = TRUE,# Allows the process to be visible for monitoring convergence.      
                      nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))# Control parameters for the fitting process.
      
      # Summarize the fitting result to assess the quality of the model and parameter estimates.
      fit.logp.summ <- summary(fit.logp)
      
      # Predict the proportion mature over a continuous range of ages based on the fitted model.
      # this is 'bx' in Bradshaw et al 2018  
      age.vec.cont <- seq(1,max(age.vec),0.02) # A continuous age vector for prediction.
      pred.p.mat <- coef(fit.logp)[1] / (1+(age.vec.cont/coef(fit.logp)[2])^coef(fit.logp)[3]) # Predicted values.
      pred.p.matt <- ifelse(pred.p.mat > 1, 1, pred.p.mat)# Ensuring predicted proportions don't exceed 1.
      
      # Add the model predictions to the plot to compare with the observed data.
      lines(age.vec.cont,pred.p.matt,lty=2,lwd=3,col="red")
      
      # Create a data frame of the continuous age vector and the corresponding predicted proportions mature.
      # This can be used for further analysis or exported for documentation.
      mat.fit.out <- data.frame(age.vec.cont, pred.p.matt)
      
      # new b to match the data
      #b.new <- coef(fit.logp)[2] 
      
      # Create a continuous vector of ages for prediction, similar to previous steps.
      age.vec.cont <- seq(1,max(age.vec),0.02)
      
      # Calculate predicted proportions mature using the new 'b' value while keeping 'a' and 'c' from the original fit.
      # This modification allows us to specifically assess the impact of 'b' on the model's predictions.
      pred.p.mat <- coef(fit.logp)[1] / (1+(age.vec.cont/b.new)^coef(fit.logp)[3])
      pred.p.mat2 <- ifelse(pred.p.mat > 1, 1, pred.p.mat) # Ensure that predicted proportions do not exceed 1, as the proportion mature cannot be greater than 100%.
      
      
      # Create a data frame of the continuous age vector and the modified predicted proportions mature.
      out.mat <- data.frame(age.vec.cont,pred.p.mat2)
      
      # Extract a subset of the modified predictions for ages between 6 and 30
      sub.mat <- which(out.mat$age.vec.cont > 6 & out.mat$age.vec.cont < 30)
      out.mat[sub.mat,]
      
      # Finally, create a comprehensive data frame of the age vector and the modified predicted proportions mature
      mat.fit2.out <- data.frame(age.vec.cont, pred.p.mat2)
      
      
      # Litter size data
      
      size.litter<-data.frame(seq(1:maxlong),1)
      colnames(size.litter)<-c("age","litter_size")
      size.litter$litter_size[1:6]<-0
      
      fit.loglit <- nls(litter_size ~ al / (1+(age/bl)^cl), 
                        data = size.litter,
                        algorithm = "port",# Specifies the optimization algorithm to use.
                        start = c(al = param.init[1], bl = param.init[2], cl= param.init[3]),# Initial values for parameters.
                        trace = TRUE,# Allows the process to be visible for monitoring convergence.      
                        nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))# Control parameters for the fitting process.
      
      # Predict the proportion mature over a continuous range of ages based on the fitted model.
      # this is 'bx' in Bradshaw et al 2018  
      fred.age.vec.cont <- seq(1,30,0.02) # A continuous age vector for prediction.
      fred.pred.p.mat <- coef(fit.loglit)[1] / (1+(fred.age.vec.cont/coef(fit.loglit)[2])^coef(fit.loglit)[3]) # Predicted values.
      fred.pred.p.matt <- ifelse(fred.pred.p.mat > 1, 1, fred.pred.p.mat)# Ensuring predicted proportions don't exceed 1.
      
      # Add the model predictions to the plot to compare with the observed data.
      lines(fred.age.vec.cont ,fred.pred.p.matt,lty=2,lwd=3,col="red")
      
      
      ### Building average fertility vector for matrix
      # Calculate predicted proportion mature using the coefficients from the logistic power model (fit.logp)
      pred.p.mat3 <- coef(fit.logp)[1] / (1+(age.vec/b.new)^coef(fit.logp)[3])
      
      # Ensure the predicted proportions do not exceed 1. If the model predicts values over 1,
      # they are adjusted down to 1, reflecting the maximum possible proportion of mature individuals.
      pred.p.mat4 <- ifelse(pred.p.mat3 > 1, 1, pred.p.mat3)
      
      # Ensure that very small predicted values are adjusted to 0 to avoid unrealistic proportions of mature individuals.
      # This step cleans up the data by setting a threshold below which the proportion mature is considered to be effectively zero.
      pred.p.mat5 <- ifelse(pred.p.mat4 < 0.001, 0, pred.p.mat4)
      
      ### FS: Calculate predicted litter sizes based on age using the coefficients from FERTILITY MODEL FITTED L152
      # This represents the expected number of offspring per mature female as a function of age.
      litt.pred2 <- coef(fit.loglit)[1] / (1+(age.vec/coef(fit.loglit)[2])^coef(fit.loglit)[3]) 
      
      # Combine the predicted proportion mature and the predicted litter size to calculate the average fertility vector.
      ##Change this to account for calving period - breeding once every 3.15 years
      f.fert.vec <- calving.prd * (pred.p.mat5*litt.pred2) #  * 0.31746
      
      #### Survival data - removed 
      
      #### Populating the matrix 
      # Updating the First Row of the Population Matrix with Fertility Rates
      popmat[1, ] <- 0.5 * f.fert.vec # * 0.5 for female offspring only
      
      # Updating the Main Diagonal with Survival Probabilities
      diag(popmat[2:(stages), ]) <- surv.vec[-stages]
      
      
      # Setting the Survival Probability in the Last Stage to 0
      popmat[stages,stages] <- 0
      
      popmat.orig <- popmat # Preserve the original matrix setup
      
      max.lambda(popmat) ## 1-yr lambda
      
      max.r(popmat) # rate of population change, 1-yr
      
      ssd <- stable.stage.dist(popmat) ## stable stage distribution
      
      R.val(popmat,stages) # reproductive value
      
      gen.l <- G.val(popmat,stages) # mean generation length
      
      # Calculate the probability of a catastrophic event 
      cat.pr <- cat.num/gen.l
      
      # Initializing the Population Vector
      population_size <- 1440
      init.vec <- 1440*stable.stage.dist(popmat)
      plot(age.vec,init.vec,xlab="age (yrs)", ylab="N", type="l")
      
      
      ##      Projection setup
      ## This section sets up the time frame for the population projection in one-year increments.
      yr.now <- 2022  # Current year for the projection start.
      yr.end <- 2050 # # Target end year for the projection.
      t <- (yr.end - yr.now) # Total number of years over which the projection will run.
      
      # Revert the population matrix to its original state (to make sure we start with a clean slate).
      popmat <- popmat.orig # Use the original population matrix as the baseline for projection.
      
      
      ## Initialize the population storage matrix
      # Create a matrix to store population vectors for each year of the projection.
      n.mat <- matrix(0,nrow=stages,ncol=(t+1)) # The matrix has a column for each year from yr.now to yr.end.
      n.mat[,1] <- init.vec  # The initial population vector is set as the first column of the matrix.
      
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
      
      
      ###########Stochastic simulations
      
      ## Initial Population Setup
      start.pop <- 1440 # Initialize the starting population size at 1000 individuals
      # multiplies this starting population size by the stable stage distribution calculated from the original population matrix (popmat.orig), ensuring that the population structure is balanced across life stages at the outset.
      init.vec <- start.pop*stable.stage.dist(popmat.orig) ## Calculate the initial population vector adjusted by the stable stage distribution of the original population matrix
      
      ## Time Frame for Projection
      yr.now <- 2024 #Define the current year for the start of the projection
      yr.end <- 2100 # Define the end year for the projection, setting a long-term analysis
      
      ##Time Projection Setup
      # Define the number of years for the projection based on the start and end years.
      t <- (yr.end - yr.now) # Total years for projection
      yrs <- seq(yr.now,yr.end,1) # Vector of years from start to end
      
      ## Density-Feedback Function for Population Regulation:
      K <- 3*start.pop # Define a carrying capacity as three times the starting population.
      
      # Create a fertility reduction factor that scales with population size, aiming to simulate density-dependent fertility.
      # sfert.min.mult <- 0.7 # Minimum multiplier for fertility reduction
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
      
      n.mat <- matrix(0,nrow=stages,ncol=(t+1))
      n.mat[,1] <- init.vec  # Set the initial population vector as the first column
      popmat <- popmat.orig # Reset the population matrix to its original state for the simulation
      
      ## Stochastic Projection Iterations
      iter <- 10000 # Number of stochastic iterations
      itdiv <- iter/1000 # Interval for progress reporting
      
      # Set up storage matrices and vectors for simulation results
      n.sums.mat <- matrix(data = 0, nrow = iter, ncol = (t+1))
      
      # Create a matrix to store the rate of population change ('r') for each year of each iteration
      # This matrix has 'iter' rows and 't' columns (one for each year of the projection)
      r.mat <- matrix(data = 0, nrow = iter, ncol = t)
      
      # Initialize vectors to store the weighted mean age and total length (TL) for the population at the end of each iteration
      # Both vectors are set to zero and have a length equal to the number of iterations
      age.wm.vec <- TL.wm.vec <- rep(0, iter)
      
      # STEP #8: stochastic simulations
      
      for (e in 1:iter) {
        # Initialize a vector to store annual growth rates for the current iteration
        r.stoch <- rep(0, t)
        
        # Reset the population matrix to its original configuration at the start of each iteration
        popmat <- popmat.orig
        
        # Loop through each year of the simulation
        for (i in 1:t){
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
          #######problems here
          for (x in 1:stages) {
            Sx.stoch[x] <- rbeta(1, Sx.alpha[x], Sx.beta[x])
          }
        
          # Update the first row of the population matrix 'popmat' with the stochastic fertility vector
          popmat[1, ] <- f.fert.stoch * 0.5
          
          # Simulates a catastrophic event occurring with a probability 'cat.pr'.
          catastrophe <- rbinom(1, 1, cat.pr)
          if (catastrophe == 1) {
            diag(popmat[2:(stages), ]) <- (cat.die.off*Sx.stoch[-stages])}
          if (catastrophe == 0) {
            diag(popmat[2:(stages), ]) <- Sx.stoch[-stages]}
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
      
      # Initialize vectors to store the mean, upper, and lower confidence limits of population size
      n.up <- n.lo <- n.mn <- rep(0,(t+1))
      # Calculate mean and 95% confidence intervals for each year's population size
      for (q in 1:(t+1)) {
        n.mn[q] <- mean(n.sums.mat[,q]) # Calculate mean population size for year q
        n.up[q] <- as.vector(quantile(n.sums.mat[,q], probs=0.975)) # 97.5% quantile
        n.lo[q] <- as.vector(quantile(n.sums.mat[,q], probs=0.025)) # 2.5% quantile
      }
      
      # Calculate mean and 95% confidence intervals for the growth rate 'r' for each year
      r.up <- r.lo <- r.mn <- rep(0,(t))
      for (q in 1:(t)) {
        r.mn[q] <- mean(r.mat[,q]) # Calculate mean growth rate for year q
        r.up[q] <- as.vector(quantile(r.mat[,q], probs=0.975)) # 97.5% quantile
        r.lo[q] <- as.vector(quantile(r.mat[,q], probs=0.025)) # 2.5% quantile
      }
  
  # Store the results of the model output
  mean.long.term.r <- mean(r.mn)
  mean.long.term.n <- mean(n.mn)
  
  results_pop_size[[r]] <- mean.long.term.n
  results_pop_rate[[r]] <- mean.long.term.r
} 



SA.pop.size <- cbind(lhs_df, pop_size = unlist(results_pop_size))
SA.pop.rate <- cbind(lhs_df, pop_rate = unlist(results_pop_rate))              

## BRT
library(dismo)
library(gbm)
brt.fit.n <- gbm.step(SA.pop.size, gbm.x = attr(SA.pop.size, "names")[1:6], gbm.y = attr(SA.pop.size, "names")[7], family="gaussian", tolerance = 0.0001, learning.rate = 0.01, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit.n)
gbm.plot(brt.fit.n)
gbm.plot.fits(brt.fit.n)


brt.fit.r <- gbm.step(SA.pop.rate, gbm.x = attr(SA.pop.rate, "names")[1:6], gbm.y = attr(SA.pop.rate, "names")[7], family="gaussian", tolerance = 0.0001, learning.rate = 0.01, bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit.r)
gbm.plot(brt.fit.r)
gbm.plot.fits(brt.fit.r)



dev.off()

brt.summary <- summary(brt.fit.n, plotit = FALSE)

# Load ggplot2 for custom plotting
library(ggplot2)

# Create a custom bar plot with ggplot2, removing grid lines and adding black axis lines
ggplot(brt.summary, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flips the bars horizontally
  labs(x = "Variables", y = "Relative influence (%)") +
  theme_minimal(base_size = 14) +  # Minimal theme with larger font size
  theme(panel.grid = element_blank(),  # Remove grid lines
        axis.line = element_line(color = "black"),  # Add black axis lines
        axis.ticks = element_line(color = "black"),  # Add black tick marks
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))


# Fit the model for population rate data
brt.fit.r <- gbm.step(SA.pop.rate, gbm.x = attr(SA.pop.rate, "names")[1:6], gbm.y = attr(SA.pop.rate, "names")[7], 
                      family = "gaussian", tolerance = 0.0001, 
                      learning.rate = 0.01, bag.fraction = 0.75, 
                      tree.complexity = 2)

# Extract the summary data from the model without plotting
brt.summary.r <- summary(brt.fit.r, plotit = FALSE)

# Load ggplot2 for custom plotting
library(ggplot2)

# Create a custom bar plot with ggplot2, removing grid lines and adding black axis lines
ggplot(brt.summary.r, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flips the bars horizontally
  labs(x = "Variables", y = "Relative influence (%)") +
  theme_minimal(base_size = 14) +  # Minimal theme with larger font size
  theme(panel.grid = element_blank(),  # Remove grid lines
        axis.line = element_line(color = "black"),  # Add black axis lines
        axis.ticks = element_line(color = "black"),  # Add black tick marks
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

