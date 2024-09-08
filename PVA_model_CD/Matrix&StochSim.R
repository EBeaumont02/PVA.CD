############################################
## Common dolphin population projection model
############################################
## Matrix

# Initial setup
rm(list = ls())
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




# Data Setup for Analysis - create matrix
#
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set the working directory and source necessary R scripts

# Load custom functions for matrix operations from an external script.
source("matrixOperators.r") 

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
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
# The logistic power function form is: y = a / (1 + (x/b)^c), where 'y' is the proportion mature, 
# 'x' is the age, and 'a', 'b', 'c' are parameters to be estimated.
# This model attempts to capture how the proportion of mature individuals varies with age.
fit.logp <- nls(f.prop.mature ~ a / (1+(f.age.maturity/b)^c), 
                data = mat.dat,
                algorithm = "port",# Specifies the optimization algorithm to use.
                start = c(a = param.init[1], b = param.init[2], c = param.init[3]),# Initial values for parameters.
                trace = TRUE,# Allows the process to be visible for monitoring convergence.      
                nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))# Control parameters for the fitting process.

# Summarize the fitting result to assess the quality of the model and parameter estimates.
fit.logp.summ <- summary(fit.logp)

# Plot the observed data to visualize the original relationship between age and proportion mature.
plot(f.age.maturity,f.prop.mature,pch=19,xlab="age (yrs)", ylab="proportion mature")

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
b.new <- coef(fit.logp)[2] 

# Re-plot the original data to visualize the relationship between age and proportion mature.
# This serves as a baseline comparison for the modified model predictions.
plot(f.age.maturity,f.prop.mature,pch=19, xlab="age (yrs)", ylab="proportion mature")

# Create a continuous vector of ages for prediction, similar to previous steps.
age.vec.cont <- seq(1,max(age.vec),0.02)

# Calculate predicted proportions mature using the new 'b' value while keeping 'a' and 'c' from the original fit.
# This modification allows us to specifically assess the impact of 'b' on the model's predictions.
pred.p.mat <- coef(fit.logp)[1] / (1+(age.vec.cont/b.new)^coef(fit.logp)[3])
pred.p.mat2 <- ifelse(pred.p.mat > 1, 1, pred.p.mat) # Ensure that predicted proportions do not exceed 1, as the proportion mature cannot be greater than 100%.

# Add the modified model predictions to the plot. The use of a different color or line type helps
# distinguish these predictions from the original data and model fit.
lines(age.vec.cont,pred.p.mat,lty=2,lwd=3,col="red")
lines(age.vec.cont,pred.p.mat2,lty=2,lwd=3,col="green")

# Create a data frame of the continuous age vector and the modified predicted proportions mature.
out.mat <- data.frame(age.vec.cont,pred.p.mat2)

# Extract a subset of the modified predictions for ages between 6 and 30
sub.mat <- which(out.mat$age.vec.cont > 6 & out.mat$age.vec.cont < 30)
out.mat[sub.mat,]

# Finally, create a comprehensive data frame of the age vector and the modified predicted proportions mature
# The first column contains a continuous sequence of ages, starting from 1 year to the maximum age present in the age vector (age.vec), incremented by 0.02 years. This fine granularity allows for a smooth curve when plotting the model's predictions.
# The second column contains the predicted proportions of mature individuals corresponding to each age in age.vec.cont. These predictions are calculated using the modified logistic power function with the new value of b. The predictions are adjusted to ensure they fall within the range of 0 to 1, representing realistic proportions (from 0% to 100% maturity).
mat.fit2.out <- data.frame(age.vec.cont, pred.p.mat2)


# Litter size data
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

size.litter<-data.frame(seq(1:maxlong),1)
colnames(size.litter)<-c("age","litter_size")
size.litter$litter_size[1:6]<-0

fit.loglit <- nls(litter_size ~ al / (1+(age/bl)^cl), 
                  data = size.litter,
                  algorithm = "port",# Specifies the optimization algorithm to use.
                  start = c(al = param.init[1], bl = param.init[2], cl= param.init[3]),# Initial values for parameters.
                  trace = TRUE,# Allows the process to be visible for monitoring convergence.      
                  nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))# Control parameters for the fitting process.

# Summarize the fitting result to assess the quality of the model and parameter estimates.
fit.loglit.summ <- summary(fit.loglit)

# Plot the observed data to visualize the original relationship between age and proportion mature.
plot(litter_size ~ age, data = size.litter, pch=19,xlab="age (yrs)", ylab="litter size")

# Predict the proportion mature over a continuous range of ages based on the fitted model.
# this is 'bx' in Bradshaw et al 2018  
litt.age.vec.cont <- seq(1,30,0.02) # A continuous age vector for prediction.
litt.pred.p.mat <- coef(fit.loglit)[1] / (1+(litt.age.vec.cont/coef(fit.loglit)[2])^coef(fit.loglit)[3]) # Predicted values.
litt.pred.p.matt <- ifelse(litt.pred.p.mat > 1, 1, litt.pred.p.mat)# Ensuring predicted proportions don't exceed 1.

# Add the model predictions to the plot to compare with the observed data.
lines(litt.age.vec.cont ,litt.pred.p.matt,lty=2,lwd=3,col="red")


## Building average fertility vector for matrix
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Calculate predicted proportion mature using the coefficients from the logistic power model (fit.logp)
# adjusted for a new parameter value b.new. This model predicts the proportion of mature individuals
# based on age, adjusted to the new understanding of parameter 'b'
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
#####Change this to account for calving period - breeding once every 3.15 years
f.fert.vec <- 0.31746 * (pred.p.mat5*litt.pred2) #  * 0.31746

plot(f.fert.vec, 
     pch = 16,             # Solid black points
     col = "black",        # Black color for points
     xlab = "Age",         # Label for x-axis
     ylab = "Fertility",   # Label for y-axis
     las = 1)              # Orientation of axis labels


# Adding the original Fert vec line to the plot
lines(f.fert.vec, col = "red", lwd = 2, lty = 2)  # Red dashed line

# Survival data
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

plot(age.vec,surv.data$surv,pch=19,xlab="age (yrs)", ylab="survival")
lines(age.vec,fred.sil,lty=2,lwd=3,col="red")

# Plot the initial data with the primary y-axis
plot(age.vec, surv.data$surv, pch=19, xlab="age (yrs)", ylab="survival", type="p",ylim=c(min(data$surv), max(data$surv)), col="blue")
par(new=TRUE)# Overlay a new plot for the secondary data
plot(age.vec, fred.sil, type="n", axes=FALSE, xlab="", ylab="")  # Type "n" to not plot the data points
axis(side=4, at=pretty(range(fred.sil)))  # Add a secondary y-axis on the right
mtext("survival probability", side=4, line=3)  # Label the secondary y-axis
lines(age.vec, fred.sil, lty=2, lwd=3, col="red") # Add lines to the secondary plot
axis(side=1) # Ensure the x-axis is correctly labeled


# Extracting the Survival Probability Vector
# Assuming the CSV file has a column named 'Sx' representing survival probabilities (the proportion of individuals
# surviving to the next age class or year), 

iter <- length(fred.sil)
surv.vec <- numeric(length(fred.sil) - 1)
for (i in 2:iter) {
  surv.vec[i - 1] <- fred.sil[i] / fred.sil[i - 1]
}

# Set the last stage to 0
surv.vec <- c(surv.vec, 0)

# Plot to visualise the age vector
plot(age.vec, surv.vec, pch=19, xlab = "Age", ylab = "surv.vec", type="p", col = "blue")
surv.vec

plot(age.vec, surv.vec, pch = 19, col = "black", 
     xlab = "Age", ylab = "Survival Probability", 
     type = "p")
lines(age.vec[-length(age.vec)], surv.vec[-length(surv.vec)], col = "black")

## Populating the matrix
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Updating the First Row of the Population Matrix with Fertility Rates
popmat[1, ] <- 0.5 * f.fert.vec # * 0.5 for female offspring only

# Updating the Main Diagonal with Survival Probabilities
# The last stage's survival rate is omitted because it typically transitions out of the modeled population.
diag(popmat[2:(stages), ]) <- surv.vec[-stages]

# Setting the Survival Probability in the Last Stage to 0
popmat[stages,stages] <- 0

# Saving a Copy of the Original Population Matrix
popmat.orig <- popmat # Preserve the original matrix setup

# Displaying Portions of the Population Matrix for Inspection
popmat[1,]  # Show the first row, representing fertility rates adjusted for female offspring.
popmat[1:10,1:10] # Show the top left corner of the matrix, useful for smaller matrices or a quick check.
popmat[19:29,19:29] # Show the bottom right corner, often where terminal stages and survival probabilities are of interest.

# Calculate the dominant eigenvalue (lambda) of the population matrix.
# The dominant eigenvalue, often called 'lambda', represents the growth rate of the population over one time unit (e.g., one year).
# A lambda > 1 indicates growth, lambda = 1 indicates stability, and lambda < 1 indicates decline.
max.lambda(popmat) ## 1-yr lambda

# Calculate the intrinsic rate of increase (r), which is the natural logarithm of lambda.
# The intrinsic rate of increase provides another measure of population growth rate per time unit.
max.r(popmat) # rate of population change, 1-yr

# Calculate the stable stage distribution, which shows the proportional distribution of individuals across different stages
# in a long-term stable population. This distribution is independent of the initial population structure.
ssd <- stable.stage.dist(popmat) ## stable stage distribution

# Calculate the reproductive value (R0) for the population.
# Reproductive value gives the relative contribution of individuals in different stages to future generations of the population.
R.val(popmat,stages) # reproductive value

# Calculate the mean generation length (G), which is the average time between when mothers give birth and when their offspring do.
gen.l <- G.val(popmat,stages) # mean generation length

# Calculate the probability of a catastrophic event affecting the population, based on the mean generation length.
# The example uses a specific formula from Reed et al. (2003), where the catastrophe probability is inversely related
# to the mean generation length, illustrating a simplistic model where longer generations imply less frequent catastrophes.
cat.pr <- 0.14/gen.l # probability of catastrophe (Reed et al. 2003)

# Initializing the Population Vector
# The 'stable.stage.dist(popmat)' function calculates the stable stage distribution for the population,
# which represents the proportion of individuals expected in each age class when the population reaches equilibrium.
# This distribution is multiplied by 1000 to scale up the proportions to actual numbers, creating an initial population vector.
# This vector, 'init.vec', thus represents the initial number of individuals in each age class for the modeled population.
effective_population <- 1400
population_size <- 1400
init.vec <- 1400*stable.stage.dist(popmat)

# Visualizing the Initial Population Distribution
# The 'plot' function is used to create a line graph displaying the initial population distribution across age classes.
# 'age.vec' on the x-axis represents the different age classes, while 'init.vec' on the y-axis shows the number of individuals
# in each age class. The labels for the x and y axes are set to "age (yrs)" and "N" respectively, indicating the age in years
# and the number of individuals. The 'type="l"' argument specifies that the data should be plotted as a line rather than points,
# providing a clear visual representation of how the initial population is distributed across age classes.
plot(age.vec,init.vec,xlab="age (yrs)", ylab="N", type="l")


##      Calculating total population
#
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


effective_population <- 1500

# plot initial age distribution vector
plot(age.vec,init.vec,xlab="age (yrs)", ylab="N", type="l")

# calculate the proportion mature in the inital vector
proportion_init_mature <- sum(init.vec[age.vec > 6]) / sum(init.vec)

# calculate population using proportion_init_mature
total_population <- effective_population / proportion_init_mature


#############################################################################################################################
## Population projection & growth rate 

##      Projection setup
#
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Revert the population matrix to its original state (to make sure we start with a clean slate).
popmat <- popmat.orig # Use the original population matrix as the baseline for projection.


## Initialize the population storage matrix
# Create a matrix to store population vectors for each year of the projection.
#n.mat <- matrix(0,nrow=stages,ncol=(t+1)) # The matrix has a column for each year from yr.now to yr.end.
#n.mat[,1] <- init.vec  # The initial population vector is set as the first column of the matrix.

# Stochastic projection simulation
# This loop simulates population changes over time under the specified conditions using the Leslie matrix model.

#for (i in 1:t) {
#  n.mat[,i+1] <- popmat %*% n.mat[,i]  # Multiply the current population vector by the population matrix to get the next year's population vector.
#}

## Create a vector of years for plotting.
#yrs <- seq(yr.now,yr.end,1) # Sequence of years from start to end of projection.

# Plot the total population over time
# This plot visualizes how the total population is expected to change from yr.now to yr.end.
#plot(yrs,as.vector(colSums(n.mat)),type="l",xlab="year",ylab="N",xlim=c(yr.now,yr.end))
# 'type="l"' creates a line graph, which is helpful for showing trends in population changes over time.


####################################################Stochastic simulations

## Initial Population Setup
start.pop <- 2800 # Initialize the starting population size at 1000 individuals
# multiplies this starting population size by the stable stage distribution calculated from the original population matrix (popmat.orig), ensuring that the population structure is balanced across life stages at the outset.
init.vec <- start.pop*stable.stage.dist(popmat.orig) ## Calculate the initial population vector adjusted by the stable stage distribution of the original population matrix

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
n.mat[,1] <- init.vec  # Set the initial population vector as the first column
popmat <- popmat.orig # Reset the population matrix to its original state for the simulation

## Stochastic Projection Iterations
# Purpose: To perform multiple simulations of population projections under varying conditions to capture uncertainty and variability in population dynamics.
# Implementation: Specifies the number of iterations (iter) for the stochastic projections and sets up for detailed simulation processes in the subsequent lines 
# Define the number of iterations for the stochastic projection to simulate variability.
iter <- 10000 # Number of stochastic iterations
itdiv <- iter/1000 # Interval for progress reporting

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
    popmat[1, ] <- f.fert.stoch ######### here i think it needs to be halved to account for just females
    
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

