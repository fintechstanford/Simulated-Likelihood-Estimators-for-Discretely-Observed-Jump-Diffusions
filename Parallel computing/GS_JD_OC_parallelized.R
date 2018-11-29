##############################################################################################################
### This code computes the optimal configuration for the density estimator of the paper "Simulated Likelihood Estimators for Discretely Observed Jump-Diffusions" by Kay Giesecke and Gustavo Schwenkler. More specifically, the code solves the optimization problem (22) of the paper.
### The codes were implemented by the authors under special assistance of Francois Guay at Boston University. The copyright to these codes remains with the authors.
##############################################################################################################

source("GS_JD_algorithms_parallelized.R")
source("GS_JD_randomelement.R")
source("GS_JD_modelinput.R")

# Enter the parameter at which to evaluate the optimal configuration
theta = c(0.8542, 0.0330, 0.0173, 0.2162*250, 0, 0.0004, 0.0058)

# Enter the value of Delta at which to evaluate the density estimator from Theorem 4.1
Delta = 1/12

# Enter the number of cores you want to use for parallel computing
numworkers = 4

# Enter the number of times you want to run the maximization with different starting points
N=5

# Enter an initial guess for l, L and rho
ell_L_rho0 = c(56.25, 3.9, 0)

# Enter some upper and lower bounds on l, L and rho
lbounds = c(0.1,1,-10)
ubounds = c(100,10,10)

### DO NOT CHANGE ANYTHING BELOW ###

# L, ell and rho are defined globally
ell = ell_L_rho0[1]
L = ell_L_rho0[2]
rho = ell_L_rho0[3]

# Run a first estimate for the value function
minsofar = max_var_rt(ell_L_rho0, Delta, theta, lbounds, ubounds, numworkers)
ellLrhosofar = ell_L_rho0

# Try to estimate the OC parameters with 20 different starting points
for (n in 1:N) {

  # initial randomized guess
	ellLrho = c(lbounds[1],lbounds[2],lbounds[3]) + runif(3)*c(ubounds[1]-lbounds[1], ubounds[2]-lbounds[2], ubounds[3]-lbounds[3])

  # Estimate the optimal configuration
	Try = opt_ell_L_rho(ellLrho, Delta, theta, lbounds, ubounds, numworkers)
  # Keep the optima if the value function is lower
	if (Try[1] < minsofar) {
		minsofar = Try[1]
		ellLrhosofar = Try[2:4]
	}
}

best = c(minsofar, ellLrhosofar)

print("The optimal configuration parameters are : ")
print(best[2:4])
print("The optimal value is : ")
print(best[1])

