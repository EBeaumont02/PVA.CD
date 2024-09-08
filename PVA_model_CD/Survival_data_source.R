
# Survival data - Bycatch, Stranding, Mass stranding
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Siler model:

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


# Strandings

# Create data frame of proportion of each age group from stranding data
age.vec <- seq(0, 29, 1)
surv.s<-c(10,6,8,3,2,1,3,2,6,2,0,6,3,2,2,3,4,0,4,0,3,0,1,0,1,0,0,0,0,1)
data.s<-data.frame(age.vec,surv.s)

outfred.s <- Si.mod(data.s)

a1 <- outfred.s$par[1]
b1 <- outfred.s$par[2]
a2 <- outfred.s$par[3]
a3 <- outfred.s$par[4]
b3 <- outfred.s$par[5]

fred.sil.s<-exp(-a1/b1 * (1 - exp(-b1 * age.vec))) * exp(-a2 * age.vec) * exp(a3/b3 * (1 - exp(b3 * age.vec)))

plot(age.vec,data.s$surv.s,pch=19,xlab="age (yrs)", ylab="survival")
lines(age.vec,fred.sil.s,lty=2,lwd=3,col="red")

# Plot the initial data with the primary y-axis
plot(age.vec, data.s$surv.s, pch=19, xlab="age (yrs)", ylab="survival", type="p",ylim=c(min(data.s$surv.s), max(data.s$surv.s)), col="blue")
par(new=TRUE)# Overlay a new plot for the secondary data
plot(age.vec, fred.sil.s, type="n", axes=FALSE, xlab="", ylab="")  # Type "n" to not plot the data points
axis(side=4, at=pretty(range(fred.sil.s)))  # Add a secondary y-axis on the right
mtext("survival probability", side=4, line=3)  # Label the secondary y-axis
lines(age.vec, fred.sil.s, lty=2, lwd=3, col="red") # Add lines to the secondary plot
axis(side=1) # Ensure the x-axis is correctly labeled

# Extracting the Survival Probability Vector
# Assuming the CSV file has a column named 'Sx' representing survival probabilities (the proportion of individuals
# surviving to the next age class or year), 

iter <- length(fred.sil.s)
surv.vec.s <- numeric(length(fred.sil.s) - 1)
for (i in 2:iter) {
  surv.vec.s[i - 1] <- fred.sil.s[i] / fred.sil.s[i - 1]
}

surv.vec.s <- c(surv.vec.s, 0)

plot(age.vec, surv.vec.s, pch=19, xlab = "Age", ylab = "surv.vec.s", type="p", col = "blue")



# Mass Strandings

# Create data frame of proportion of each age group from stranding data
age.vec <- seq(0, 29, 1)
surv.m<-c(2,0,0,2,2,1,0,2,2,1,1,1,4,1,2,2,1,0,2,0,0,0,1,0,0,1,0,1,0,0)
data.m<-data.frame(age.vec,surv.m)

outfred.m <- Si.mod(data.m)

a1 <- outfred.m$par[1]
b1 <- outfred.m$par[2]
a2 <- outfred.m$par[3]
a3 <- outfred.m$par[4]
b3 <- outfred.m$par[5]

fred.sil.m<-exp(-a1/b1 * (1 - exp(-b1 * age.vec))) * exp(-a2 * age.vec) * exp(a3/b3 * (1 - exp(b3 * age.vec)))

plot(age.vec,data.m$surv.m,pch=19,xlab="age (yrs)", ylab="survival")
lines(age.vec,fred.sil.m,lty=2,lwd=3,col="red")

# Plot the initial data with the primary y-axis
plot(age.vec, data.m$surv.m, pch=19, xlab="age (yrs)", ylab="survival", type="p",ylim=c(min(data.m$surv.m), max(data.m$surv.m)), col="blue")
par(new=TRUE)# Overlay a new plot for the secondary data
plot(age.vec, fred.sil.m, type="n", axes=FALSE, xlab="", ylab="")  # Type "n" to not plot the data points
axis(side=4, at=pretty(range(fred.sil.m)))  # Add a secondary y-axis on the right
mtext("survival probability", side=4, line=3)  # Label the secondary y-axis
lines(age.vec, fred.sil.m, lty=2, lwd=3, col="red") # Add lines to the secondary plot
axis(side=1) # Ensure the x-axis is correctly labeled

# Extracting the Survival Probability Vector
# Assuming the CSV file has a column named 'Sx' representing survival probabilities (the proportion of individuals
# surviving to the next age class or year), 

iter <- length(fred.sil.m)
surv.vec.m <- numeric(length(fred.sil.m) - 1)
for (i in 2:iter) {
  surv.vec.m[i - 1] <- fred.sil.m[i] / fred.sil.m[i - 1]
}

surv.vec.m <- c(surv.vec.m, 0)

plot(age.vec, surv.vec.m, pch=19, xlab = "Age", ylab = "surv.vec.m", type="p", col = "blue")



# Bycatch

# Create data frame of proportion of each age group from stranding data
age.vec <- seq(0, 29, 1)
surv.b<-c(0,1,0,1,0,0,2,0,2,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0)
data.b<-data.frame(age.vec,surv.b)

outfred.b <- Si.mod(data.b)

a1 <- outfred.b$par[1]
b1 <- outfred.b$par[2]
a2 <- outfred.b$par[3]
a3 <- outfred.b$par[4]
b3 <- outfred.b$par[5]

fred.sil.b<-exp(-a1/b1 * (1 - exp(-b1 * age.vec))) * exp(-a2 * age.vec) * exp(a3/b3 * (1 - exp(b3 * age.vec)))

plot(age.vec,data.b$surv.b,pch=19,xlab="age (yrs)", ylab="survival")
lines(age.vec,fred.sil.b,lty=2,lwd=3,col="red")

# Plot the initial data with the primary y-axis
plot(age.vec, data.b$surv.b, pch=19, xlab="age (yrs)", ylab="survival", type="p",ylim=c(min(data.b$surv.b), max(data.b$surv.b)), col="blue")
par(new=TRUE)# Overlay a new plot for the secondary data
plot(age.vec, fred.sil.b, type="n", axes=FALSE, xlab="", ylab="")  # Type "n" to not plot the data points
axis(side=4, at=pretty(range(fred.sil.b)))  # Add a secondary y-axis on the right
mtext("survival probability", side=4, line=3)  # Label the secondary y-axis
lines(age.vec, fred.sil.b, lty=2, lwd=3, col="red") # Add lines to the secondary plot
axis(side=1) # Ensure the x-axis is correctly labeled

# Extracting the Survival Probability Vector
# Assuming the CSV file has a column named 'Sx' representing survival probabilities (the proportion of individuals
# surviving to the next age class or year), 

iter <- length(fred.sil.b)
surv.vec.b <- numeric(length(fred.sil.b) - 1)
for (i in 2:iter) {
  surv.vec.b[i - 1] <- fred.sil.b[i] / fred.sil.b[i - 1]
}

surv.vec.b <- c(surv.vec.b, 0)

plot(age.vec, surv.vec.b, pch=19, xlab = "Age", ylab = "surv.vec.b", type="p", col = "blue")





