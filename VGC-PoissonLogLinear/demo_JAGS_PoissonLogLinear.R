# Poisson Log Linear Regression example for the paper 
# "Variational Gaussian Copula Inference", 
# Shaobo Han, Xuejun Liao, David. B. Dunson, and Lawrence Carin, 
# The 19th International Conference on Artificial Intelligence and Statistics (AISTATS 2016)
# -----------------------------
# Download and install JAGS as per operating system requriements.
# http://mcmc-jags.sourceforge.net/ 
# -----------------------------
# JAGS (Implemented via RJAGS) 
# For more information about this dataset, see
# https://github.com/petrkeil/Statistics/blob/master/Lecture%203%20-%20poisson_regression/poisson_regression.Rmd
# -----------------------------
# Shaobo Han, Duke University
# shaobohan@gmail.com
# 10/01/2015

library(R.matlab)
library(rjags)

path <- getwd()
pathname <- file.path(path, "TreeData50.mat")
data <- readMat(pathname)
attach(data)

##  JAGS
jags.data <- list(N.cells = length(y), n50 = y, elev50 = x)
cat("\n      model\n      {\n        # priors\n
    beta0 ~ dnorm(0,sigma)\n 
    beta1 ~ dnorm(0,sigma)\n
    beta2 ~ dnorm(0,sigma)\n 
    sigma <- pow(tau, -1)\n
    tau ~ dgamma(1, 1)\n
    \n        # likelihood\n
    for(i in 1:N.cells)\n        
    {\n          n50[i] ~ dpois(lambda[i])\n          
    log(lambda[i]) <- beta0 + beta1*elev50[i] + beta2*pow(elev50[i],2)\n         
    # this part is here in order to make nice prediction curves:\n          
    prediction[i] ~ dpois(lambda[i])\n        } \n      }\n  ", 
    file = "PoissonLogLinearModel.txt")
params <- c("beta0", "beta1", "beta2", "tau", "prediction")
jm <- jags.model("PoissonLogLinearModel.txt", data = jags.data, n.chains = 10, n.adapt = 1000)
update(jm, n.iter = 10000)
jm.sample <- jags.samples(jm, variable.names = params, n.iter = 10000, thin = 1)

summary(as.mcmc.list(jm.sample$beta0))
summary(as.mcmc.list(jm.sample$beta1))
summary(as.mcmc.list(jm.sample$beta2))
summary(as.mcmc.list(jm.sample$tau))

# Save MCMC samples to mat file
filename <- paste("PoissonLogMCMC", ".mat", sep="")
writeMat(filename, beta0=jm.sample$beta0, beta1=jm.sample$beta1, beta2=jm.sample$beta2,
         tau=jm.sample$tau)