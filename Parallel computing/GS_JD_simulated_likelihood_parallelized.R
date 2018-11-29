##############################################################################################################
### This code solves the maximum likelihood inference problem for a one-dimensional jump diffusion based on the density estimator of the paper "Simulated Likelihood Estimators for Discretely Observed Jump-Diffusions" by Kay Giesecke and Gustavo Schwenkler. The codes were implemented by the authors under special assistance of Francois Guay at Boston University. The copyright to these codes remains with the authors.
##############################################################################################################

source("GS_JD_algorithms_parallelized.R")
source("GS_JD_modelinput.R")
source("GS_JD_randomelement.R")

# Enter the values of L, ell, and rho to use
L = 3.65311279
ell = 58.17718506
rho = 0.02044678

# Enter the data
data = read.table("exact_paths_das_600.txt", header = FALSE)
data = as.matrix(data)
data = data[1,]
datalength = length(data)

# Enter an initial parameter guess
theta0 = c(0.8542, 0.0330, 0.0173, 54.05, 0, 0.0004, 0.0058)

# Enter lower and upper bounds for the parameter space
lower_b = c(1e-4, -1, 1e-4, 1e-4, -10, -0.1, 1e-4)
upper_b = c(3   ,  1, 1   , 100 ,  10,  0.1, 0.1)

# Enter the time period between consecutive observations of the data
Delta = 1/12

# Enter the number of cores you want to use for parallel computing
numworkers = 4

# Enter the maximum number of iterations for the Nelder-Mead algorithm
maxIter = 5000

# Choose the number of Monte-Carlo simulations to use to compute the simulated density estimator hat_p_K in equation (23). This has an impact on the asymptotic distribution of the SMLE as indicated in Theorem 6.2
K = 10*datalength

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
	
# Solve the maximum simulated likelihood problem.
result = max_llf(theta0, data, lower_b, upper_b, maxIter, numworkers)
	
# Compute the standard errors
sterr = st_errors(result[4:10], data, Delta, numworkers)

# Return the maximum value of the log-likelihood, as well as the SMLE that maximizes the simulated likelihood.
cat("The value of the likelihood function is: ", result[1])
cat("The number of iterations in the search algorithm is: ", result[2])
cat("The optimal parameter is: ", result[4:10])
cat("Its standard errors are: ", sterr)

