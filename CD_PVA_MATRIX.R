############################################
## Common dolphin population projection model
############################################
## Matrix

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
fred.age.vec.cont <- seq(1,30,0.02) # A continuous age vector for prediction.
fred.pred.p.mat <- coef(fit.loglit)[1] / (1+(fred.age.vec.cont/coef(fit.loglit)[2])^coef(fit.loglit)[3]) # Predicted values.
fred.pred.p.matt <- ifelse(fred.pred.p.mat > 1, 1, fred.pred.p.mat)# Ensuring predicted proportions don't exceed 1.

# Add the model predictions to the plot to compare with the observed data.
lines(fred.age.vec.cont ,fred.pred.p.matt,lty=2,lwd=3,col="red")



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


f.fert.vec

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

