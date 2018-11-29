##############################################################################################################
### This code computes the density estimator of the paper "Simulated Likelihood Estimators for Discretely Observed Jump-Diffusions" by Kay Giesecke and Gustavo Schwenkler. More specifically, the code compute hat_p_K as in equation (23), which is a Monte-Carlo average over K samples of the density estimator computed in accordance with Theorem 4.1.
### The codes were implemented by the authors under special assistance of Francois Guay at Boston University. The copyright to these codes remains with the authors.
##############################################################################################################

source("GS_JD_algorithms.R")
source("GS_JD_randomelement.R")
source("GS_JD_modelinput.R")

# Enter the parameter at which to evaluate the density at
theta = c(0.8542, 0.0330, 0.0173, 0.2162*250, 0, 0.0004, 0.0058)

# Enter the value of "v" at which to evaluate the density estimator of Theorem 4.1.
v = theta[2]

# Enter the range for values of w at which to evaluate the density estimator of Theorem 4.1. "lb" is the smallest value of w to consider, ub is the largest value of w to consider, and le is the length of the vector of possible values of w
lb = 0
ub = 0.1
le = 100
w = seq(from = lb, to = ub, length = le)

# Enter the value of Delta at which to evaluate the density estimator from Theorem 4.1
Delta = 1/12

# Enter the values of L, ell, and rho to use
L = 3.65311279
ell = 58.17718506
rho = 0.02044678

# Enter the number of Monte-Carlo samples to use
K = 5000

### DO NOT CHANGE ANYTHING BELOW ###

# Generate the random element R
PT = gen_P_T(Delta, K, ell)
all_E = gen_E(PT, L)
# this functions checks if maxk is high enough and stops the code if not
checkmaxk(PT, all_E)
all_U = gen_U(PT)
all_V = gen_V(PT)
all_N1 = gen_N(PT)
all_N2 = gen_N(PT)
all_N3 = gen_N(PT)
all_D = gen_D(PT)

# Compute the density estimator. 
temp = log_hat_p_K_vec(rep(v, times = le), Delta, w, theta)
UD = exp(temp)

# Store the density estimator in a txt file
write(c(theta, UD), "UD.txt", ncolumns = length(UD) + length(theta), sep = "\t", append = TRUE)
