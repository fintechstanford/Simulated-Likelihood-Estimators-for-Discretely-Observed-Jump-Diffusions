##############################################################################################################
### These functions implement the density estimator of the paper "Simulated Likelihood Estimators for Discretely Observed Jump-Diffusions" by Kay Giesecke and Gustavo Schwenkler. The codes were implemented by the authors under special assistance of Francois Guay at Boston University. The copyright to these codes remains with the authors.
##############################################################################################################

# maxk and maxp should be changed if the algorithm stops and tells you to do so.
maxk = 1
maxp = 400
pi = 3.14159265

# You need this library to compute the gradient of the function numerically
library("stats")

# this package is used to evaluate numerical derivatives and integrals when some functions are not known in closed-form
library("numDeriv")

# The prefix "cmpfun" in front of every function is used to compile the functions only once at the beginning of program. This speeds up the running time of the functions. In order to use the compiling options, we need to install and call the "compiler" library
#install.packages("compiler", repos = "http://cran.r-project.org")
library("compiler")

# you need this library to use paralel computing
#install.packages("parallel", repos = "http://cran.r-project.org")
library("parallel")

##############################################################################################################
#   Do not change the functions below; these implement the density estimator from Giesecke & Schwenkler (2014)
##############################################################################################################


# mu_Y(y;theta) : drift function of Y
mu_Y = cmpfun({function(y, theta, x0) {
	
  mu_X(F_inv(y,theta,x0),theta)/sigma(F_inv(y,theta,x0),theta) - 0.5*sigma_prime(F_inv(y,theta,x0),theta)
  
}
})


# mu_Y'(y;theta) : First order Derivative of mu_Y(y;theta) with respect to y
mu_Y_prime = cmpfun({function(y, theta, x0) {

  F_inv_prime(y,theta,x0) *  (  (mu_X_prime(F_inv(y,theta,x0),theta)*sigma(F_inv(y,theta,x0),theta)-mu_X(F_inv(y,theta,x0),theta)*sigma_prime(F_inv(y,theta,x0),theta)   )  /  sigma(F_inv(y,theta,x0),theta)^2  - 0.5* sigma_prime_prime(F_inv(y,theta,x0),theta))
  
}
})


# Gamma_Y(y,d;theta) : Jump Size Function of the process Y
Gamma_Y = cmpfun({function(y,d, theta,x0) {
	
  F(F_inv(y,theta,x0) + Gamma_X(F_inv(y,theta,x0),d,theta),theta,x0)-y
	
}
})


# Jump Size Function of the process Y but vectorized
Gamma_Y_vec = cmpfun({function(y, d, theta,x0) {
	
	rep(Gamma_Y(y,d,theta,x0), times = length(y))
	
}
})


# a(y;theta) = "integral of (mu_y(u;theta)-rho) for values of u between 0 and y"
a_alg = cmpfun({function(y, theta, x0) {
  
  integral_mu_Y(rep(0,length(y)),y,x0,theta) - rho*y
  
}
})

# b(y;theta) = Lambda(F^(-1)(y;theta);theta) - l + 0.5*(mu_y(y;theta)^2-rho^2+mu_Y'(y;theta))
b_alg = cmpfun({ function(y, theta, x0) {
	
	temp1 = Lambda_X(F_inv(y,theta,x0), theta) - ell
	temp2 = 0.5*((mu_Y(y, theta, x0)^2) - (rho^2))	
	temp3 = 0.5*mu_Y_prime(y,theta,x0)

	temp1 + temp2 + temp3
}
})

# log(c(y,d;theta)) as in page 9 
log_c_alg = cmpfun({function(y, d, theta, x0) {
  
  D = Gamma_Y(y,d, theta,x0)
  
  temp1 = Lambda_X(F_inv(y,theta,x0), theta)/ell
  # -Integral (mu_y(u;theta)-rho)du between y and y+Gamma_Y
  temp2 = -integral_mu_Y(y,y+D,x0,theta)
  temp3 = rho*D
  
  log(temp1) + temp2 + temp3
}
})

#############################################################################################################
# Optimization

# This function performs the optimization via the Neldermead algorithm of the function optim
# theta: initial parameter
# data: time series of data
# lbounds: lower bound of paramter set
# ubounds: upper bound of parameter set
# numworkers: the number of cores to use if you can do parallel computing
max_llf = cmpfun({function(theta, data, lbounds, ubounds, maxIter, numworkers) {
	
	mll = optim(par = theta, fn = log_llf, data = data, lb = lbounds, ub = ubounds, nw = numworkers, method = "Nelder-Mead", control = list(fnscale = -1, maxit = maxIter))
	
  # returns the value of the likelihood function, the number of iterations, the optimal theta 
	c(mll$value, mll$counts, mll$par)

}
})


# return the parallelized log likelihood function of the data (scalar)
#theta : parameter at which the llf is evaluated
#data : time series of the data
#lb : lower bound on the parameter
#ub : upper bound on the parameter
#nw : the number of cores used to compute the llf with paralel computing 
log_llf = cmpfun({function(theta, data, lb, ub, nw) {
	
	m = length(data);
  
	hit_ub = min(ub - theta)
	hit_lb = min(theta - lb)
	
  # If you you hit the lb or ub, then the llf is set to be equal to -Inf
  if ((hit_lb >= 0) && (hit_ub >= 0)) {
		logllf = log_hat_p_K_vec(data[1:(m-1)], Delta, data[2:m], theta, nw)
		logllf = sum(logllf)
	} else {
		logllf = -Inf

}
	
  logllf
}
})

# return the log likelihood function of the data between each t->t+1 (vector)
#theta : parameter at which the llf is evaluated
#data : time series of the data
#lb : lower bound on the parameter
#ub : upper bound on the parameter
#nw : the number of cores used to compute the llf with paralel computing 
log_llf_vec = cmpfun({function(theta, data, lb, ub, nw) {
  
  m = length(data);
  
  hit_ub = min(ub - theta)
  hit_lb = min(theta - lb)
  
  # If you you hit the lb or ub, then the llf is set to be equal to -Inf
  if ((hit_lb >= 0) && (hit_ub >= 0)) {
    logllf = log_hat_p_K_vec(data[1:(m-1)], Delta, data[2:m], theta, nw)
  } else {
    logllf = -Inf
  }
  
  logllf
}
})

#########################################################################################################

# Monte-Carlo density estimator in equation (23)
log_hat_p_K_vec = cmpfun({function(v, t, w, theta, nw, PT1 = PT, all_E1 = all_E, all_U1 = all_U, all_V1 = all_V, all_N11 = all_N1, all_N21 = all_N2, all_N31 = all_N3, all_D1 = all_D) {
  
# v: observations at t
# w: observations at t+1
# t: frequency of the observations
# theta: parameter at which the log likelihood function is computed  
# nw: number of cores used to compute the log likelihod function  
  
  # Transform v and w using the Lamperti Transform  
  x_vec = F(v, theta, v[1])
  y_vec = F(w, theta, v[1])
  
  # Maximum number of jumps generated in the ramdom numbers
  maxP = max(PT1[,1])
  # Number of MC simulations
  K = nrow(PT1)
  
  # this vector is used for the parallel processing function
  values <- 1:K
  ## Number of workers to use:
  numWorkers <- nw
  ## Parallel calculation (mclapply):
  # return a list of J vectors, each returned by log_alg_5_1_vec
  # log_alg_5_1_vec compute log_psi according to formula (15) in the paper
    paral <- mclapply(values, function(z) log_alg_5_1_vec(t, x_vec, y_vec, theta, v[1], PT1[z,], all_E1[z,], all_U1[z,], all_V1[z,], all_N11[z,], all_N21[z,], all_N31[z,], all_D1[z,], maxP), mc.cores = numWorkers)

  # make the list a matrix, by taking each element of the list as a row of the matrix
  log_psi = do.call(rbind, paral)
  
  # test wether there are any errors returned 
  Test = rep(0,K)
  for(k in 1:K){
    # if makk is too low  
    if(is.character(unlist(paral[k]))==TRUE)
    {
      stop("You need to increase maxp - the number of uniforms generated in allV, all_N1, all_N2 and all_N3 is too low")
    }
  }
  
  # log_td is the log of the transition density according to formula (12) in the paper
  temp = a_alg(y_vec, theta, v[1]) - a_alg(x_vec, theta, v[1]) - log(sigma(w, theta)) - 0.5*log(2*pi)
  log_td = log_psi + matrix(rep(temp,each=K), nrow=K)

  # compute the expectation of the likelihood function
  res = colSums(exp(log_td))/K
  
  # take the log of that expectation to get the log likelihood function
  log(res)

}
})

# Steps 2 and 3 of algorithm 5.1 for one sample of the random element R when v and w are vectors. It returns the logarithm of the end result
log_alg_5_1_vec = cmpfun({function(t, x, y, theta, x0, P_T, allE, allU, allV, allN1, allN2, allN3, allD, maxP) {
  
  # t: frequency of the observations  
  # x: Lamperti Transform of the observations at t (vector form)
  # y: Lamperti Transform of the observations at t+1 (vector form)
  # theta: parameter at which the log likelihood function is computed  
  # x0: initial point (v(1)) used to compute the Lamperti Transform 
  
  # Number of jumps between [0,t]
  P = P_T[1]
  # Jump times, including 0 and t ([0,T1,T2,...,t])
  T = c(0, P_T[2:(P+2)])  
  
  # loading the random variables
  # Ecal : Eta's from Chen and Huang (2013) - Exit times
  Ecal = matrix(allE[1:((maxP+1)*maxk)], nrow = maxP+1, ncol = maxp)
  # U : uniforms [0,1] to compute the Brownian Motion at exit times eta 
  U = matrix(allU[1:((maxP+1)*maxp)], nrow = maxP+1, ncol = maxp)
  # V : uniforms [0,1] to compute the poisson process of intensity M 
  V = matrix(allV[1:((maxP+1)*(maxk+1)*maxp)], nrow = maxP+1, ncol = maxp*(maxk+1))
  # Random normal variables to compute the brownian bridges - see Beskos et al. (2009)
  N1 = matrix(allN1[1:((maxP+1)*(maxk+1)*maxp)], nrow = maxP+1, ncol = maxp*(maxk+1))
  N2 = matrix(allN2[1:((maxP+1)*(maxk+1)*maxp)], nrow = maxP+1, ncol = maxp*(maxk+1))
  N3 = matrix(allN3[1:((maxP+1)*(maxk+1)*maxp)], nrow = maxP+1, ncol = maxp*(maxk+1))
  
  K = rep(0, times = P+1)
  w = matrix(nrow = P+1, ncol = 1)
  # For each jump time
  for (i in 2:(P+2)) {
    # Compute the number of exit times within the interval [T[i-1],T[i]]
    K[i-1] = findInterval(T[i]-T[i-1], cumsum(Ecal[i-1,]))
    
    # If you have more exit times than w allow, make w bigger
    if (K[i-1] + 1 > ncol(w)) {
      w = cbind(w, matrix(nrow = P+1, ncol = K[i-1] + 1 - ncol(w)))	
    };
    
    # Each time we reach an exit time, determine wether the BM hit -L (U < 0.5) or L (U > 0.5)
    for (j in 1:(K[i-1]+1)) {
      w[i-1,j] = if(U[i-1,j] > 0.5) L else -L;
    }
  }
  
  le_y = length(y)
  # Y_minus and Y_plus are equal to x (vector form) - the starting point of the process Y at t=0
  Y_minus = Y_plus = matrix(x, nrow = le_y, ncol = 1)
  
  # log_Psi will contain the logarithm of Psi_hat according to formula (16) in the paper
  log_Psi = rep(0, times = le_y)
  
  # for each jump time (or on [0,t] is n=1) we can compute log Psi_hat 
  for (n in 1:(P+1)) {
    
    # load the randon variables for n
    temp_V = matrix(V[n,], nrow = maxk+1, ncol = maxp, byrow = FALSE)
    temp_N1 = matrix(N1[n,], nrow = maxk+1, ncol = maxp, byrow = FALSE)
    temp_N2 = matrix(N2[n,], nrow = maxk+1, ncol = maxp, byrow = FALSE)
    temp_N3 = matrix(N3[n,], nrow = maxk+1, ncol = maxp, byrow = FALSE)
    
    # log_E = log f_hat (see expression (20)
    log_E = rep(0, times = le_y)
    
    # y_alg is either x or YTi (Y at the jump time Ti, right after the jump)
    y_alg = matrix(Y_plus[,n], nrow = le_y, ncol = 1)
    
    # for each exit time - K[n] = I in the paper
    for (k in 1:(K[n]+1)) {
      
      # btime is minimum of first jump time and first localization time
      # etime is first localization time
      if (K[n] == 0) {
        # If there is no exit time, then we need to evaluate values of Y on [T[n],T[n+1]]
        # However we still need etime to simulate Y until the BM leave [-L,L]
        btime = T[n+1]-T[n]
        etime = Ecal[n,1]
      } else {
        # Else, then we need to evaluate values of Y until the first exit time
        btime = if (k < K[n]+1) Ecal[n,k] else T[n+1]-sum(Ecal[n,1:(k-1)]) - T[n]
        etime = Ecal[n,k]
      }
      
      # Evaluate y_alg at the exit time - for it to be our starting point at the next exit time
      y_alg = cbind(y_alg, y_alg[,k] + rep(w[n,k], times = le_y) + rho*btime)
      
      # Find m_i(theta) - algorithm 5-2 step (i)
      m = m_nk(y_alg[,k], theta, L + abs(rho*btime), x0)
      
      # Find the max of b(y;theta) -  algorithm 5-2 step (i)
      M = M_nk(y_alg[,k], theta, L + abs(rho*btime), x0)
      
      # Generate the taus - algorithm 5-2 step (ii)
      # They are jump times of a dominating poisson process on "[0, eta]" with rate M_i(theta) (M-m in the code)
      tau = matrix(0, nrow = le_y, ncol = 1)
      p = 1
      # This is the method used to draw Poisson jump times
      iat = -log(temp_V[k,p]) / (M - m)
      cont = rep(1, times = le_y)
      
      # Generate Poisson jump times until we reach the exit time
      while (max(tau[,p] + iat <= btime) == 1) {
        p = p + 1	
        tau = cbind(tau, tau[,p-1] + iat)
        if (p<=maxp) {
          iat = -log(temp_V[k,p]) / (M - m) 
        } else {
          # Not enough random variables are generated
          # maxp is too small, it needs to be increased
          stop("error with maxp - need more uniforms")
        } 
      }	
      
      # Combine the jump times tau, btime and the first localization if this happens after btime
      tau2 = cbind(tau, rep(btime, times = le_y), rep(etime, times = le_y))
      # Sort so that the time is ascending
      tau2 = t(apply(tau2, 1, sort, decreasing = FALSE))
      
      # Detect the times that we will keep when computing log_E
      Lambda_mat = (tau2 < btime)*1
      
      # Compute the Brownian Bridges at the different taus. algorithm 5-2 step (iii)
      # See Beskos et al. (2009) for the method.
      le = ncol(tau2)
      while (le-1>ncol(temp_N1)) {
        temp_N1 = cbind(temp_N1, rnorm(nrow(temp_N1)))
        temp_N2 = cbind(temp_N2, rnorm(nrow(temp_N2)))
        temp_N3 = cbind(temp_N3, rnorm(nrow(temp_N3)))
      }
      
      b1 = matrix(0, nrow = le_y, ncol = le)
      b2 = matrix(0, nrow = le_y, ncol = le)
      b3 = matrix(0, nrow = le_y, ncol = le)
      
      b1[,2] = temp_N1[k,1] * sqrt( ((etime - tau2[,2])^2)*(tau2[,2] - tau2[,1]) / ( etime*(etime - tau2[,2])*(etime - tau2[,1]) ) )
      b2[,2] = temp_N2[k,1] * sqrt( ((etime - tau2[,2])^2)*(tau2[,2] - tau2[,1]) / ( etime*(etime - tau2[,2])*(etime - tau2[,1]) ) )
      b3[,2] = temp_N3[k,1] * sqrt( ((etime - tau2[,2])^2)*(tau2[,2] - tau2[,1]) / ( etime*(etime - tau2[,2])*(etime - tau2[,1]) ) )
      
      if (le > 2) {
        for (i in 3:le) {
          temp = matrix(temp_N1[k,1:(i-1)], ncol = i-1, nrow = le_y, byrow = TRUE)
          b1[,i] = rowSums(temp * sqrt( ((etime - tau2[,i])^2)*(tau2[,2:i] - tau2[,1:(i-1)]) / ( etime*(etime - tau2[,2:i])*(etime - tau2[,1:(i-1)]) ) ), na.rm = TRUE)
          
          temp = matrix(temp_N2[k,1:(i-1)], ncol = i-1, nrow = le_y, byrow = TRUE)
          b2[,i] = rowSums(temp * sqrt( ((etime - tau2[,i])^2)*(tau2[,2:i] - tau2[,1:(i-1)]) / ( etime*(etime - tau2[,2:i])*(etime - tau2[,1:(i-1)]) ) ), na.rm = TRUE)
          
          temp = matrix(temp_N3[k,1:(i-1)], ncol = i-1, nrow = le_y, byrow = TRUE)
          b3[,i] = rowSums(temp * sqrt( ((etime - tau2[,i])^2)*(tau2[,2:i] - tau2[,1:(i-1)]) / ( etime*(etime - tau2[,2:i])*(etime - tau2[,1:(i-1)]) ) ), na.rm = TRUE)
        }
      }
      
      
      b1[,1] = b1[,le] = b2[,1] = b2[,le] = b3[,1] = b3[,le] = rep(0, times = le_y);
      
      temp = abs(w[n,k])*(etime - tau2)/etime
      B = ((temp + sqrt(etime)*b1)^2) + etime*(b2^2 + b3^2)
      
      # Take into account that the BM is equal to L or -L at exit time
      W = if (w[n,k] < 0) w[n,k] + sqrt(B) else w[n,k] - sqrt(B);
      W[,1] = rep(0, times = nrow(W))
      
      # Step (iv) of Algorithm (5.2) 
      log_e = -m*btime
      
      # Computing the remaining of the expression (20) (sum over j of log(1 - (b-m)/M)) 
      log_prod = rep(0, times = le_y)
      # sum over j - the jump times
      for (j in 2:ncol(W)) {
        test = 1 - ( (b_alg(y_alg[,k] + W[,j] + rho*tau2[,j], theta, x0) - m)/(M - m) )
        test[test<=0]=1
        log_prod = log_prod + log(test)*Lambda_mat[,j]
      }
      
      # add log_e and log_prod to the log_E computed for the precedent exit time
      log_E = log_E + log_e + log_prod
    }
    
    # Find the position of the jump time T[n]
    pos_T = Lambda_mat[,1:(ncol(Lambda_mat)-1)] -  Lambda_mat[,2:ncol(Lambda_mat)]
    pos_T = cbind(0, pos_T)
    # compute the Y_T[n]-, the value of the process Y right before the jump at T[n]
    Y_minus = cbind(Y_minus, y_alg[,k] + rowSums(pos_T*W, na.rm = TRUE) + rowSums(pos_T*tau2, na.rm = TRUE)*rho)
    
    # Compute the value of Y if we evaluate it at a jump time 
    # Y_T[n] = Y_T[n]- + Gamma_Y(Y_T[n]-,D_n,theta)
    Y_plus = cbind(Y_plus, Y_minus[,n+1] + 1*(n <= P)*Gamma_Y(Y_minus[,n+1], allD[n], theta,x0))
    
    # Update log_Psi - do not add log_c_alg(Y_minus[,n+1], allD[n], theta, x0) if we are not at a jump time
    log_Psi = if (n <= P) log_Psi + log_E + log_c_alg(Y_minus[,n+1], allD[n], theta, x0) else log_Psi + log_E; 
  }
  
  # Add the second term in expression (21) - the Gaussian Term - and -0.5*log(Delta-T_N_delta)
  # T[P+1] is the time of the last jump
  # Y_plus[,P+1] is the value of Y at the last jump time
  log_Psi = log_Psi - (((y - rho*(t - T[P+1]) - Y_plus[,P+1])^2)/(2*(t - T[P+1]))) - 0.5*log(t - T[P+1])
  
  # return log_Psi as a double
  as.double(log_Psi)	
}
})


######################################################################################################
# Standard errors for SMLE in accordance with Theorem 6.2 

st_errors = cmpfun({function(theta, data, t, nw) {
		
	m = length(data);
	
	v = data[1:(m-1)]
	w = data[2:m]
	
  # Compute the gradient of the  log likelihood function, but in vector form
	myenv <- new.env()
	assign("theta", theta, envir = myenv)
	assign("w", w, envir = myenv)
	assign("v", v, envir = myenv)
	assign("t", t, envir = myenv)
	assign("nw", nw, envir = myenv)
	grad_log_td = numericDeriv(quote(log_hat_p_K_vec(v, t, w, theta, nw)), c("theta"), myenv)
	grad_log_td = attr(grad_log_td, "gradient")
	
  # Compute the Hessian matrix using the information equality
  mat = 1/(m-1)*t(grad_log_td)%*%grad_log_td
  # Compute the inverse of the Hessian
  avcm = solve(mat)

  # The standard errors are just the square root  of the diagonal of the information matrix, divided by the number of
  # observations
	se = diag(avcm)
	sqrt(se/m)	
}
})

#############################################################################################################
# these functions are used to compute the density from t to t+1 (not the whole vector)

# Logarithm of the density estimator of Theorem 4.1 when v and w are one-dimensional
log_hat_p_K = cmpfun({function(v, t, w, theta, x0, nw, PT2 = PT, all_E2 = all_E, all_U2 = all_U, all_V2 = all_V, all_N12 = all_N1, all_N22 = all_N2, all_N32 = all_N3, all_D2 = all_D) {
	
	x = F(v,theta, x0)
	y = F(w,theta, x0)
	
	K = nrow(PT2)
	maxP = max(PT2[,1])
	
	values <- 1:K
	## Number of workers (R processes) to use:
	numWorkers <- nw
	## Parallel calculation (mclapply):
	# return a list of J vectors, each returned by log_alg_5_1
	paral <- mclapply(values, function(z) log_alg_5_1(t, x, y, theta, x0,  PT2[z,], all_E2[z,], all_U2[z,], all_V2[z,], all_N12[z,], all_N22[z,], all_N32[z,], all_D2[z,], maxP), mc.cores = numWorkers)
	
	# make the list a matrix, by taking each element of the list as a row of the matrix
	log_psi = do.call(rbind, paral)
	temp = a_alg(y, theta, x0) - a_alg(x, theta, x0) - log(sigma(w, theta)) - 0.5*log(2*pi)
	log_td = log_psi + temp
	
	res = sum(exp(log_td))/K
	
	log(res)

}
})

# Steps 2 and 3 of Algorithm 5.1 for one sample of the random element R when v and w are one-dimensional. It returns the logarithm of the end result
log_alg_5_1 = cmpfun({function(t, x, y, theta, x0, P_T, allE, allU, allV, allN1, allN2, allN3, allD, maxP) {
  
  # Number of jumps between [0,t]
  P = P_T[1]
  # Jump times, including 0 and t ([0,T1,T2,...,t])
  T = c(0, P_T[2:(P+2)])  
  
  # loading the random variables
  # Ecal : Eta's from Chen and Huang (2013) - Exit times	
  Ecal = matrix(allE[1:((maxP+1)*maxk)], nrow = maxP+1, ncol = maxp)
  # U : uniforms [0,1] to compute the Brownian Motion at exit times eta
  U = matrix(allU[1:((maxP+1)*maxp)], nrow = maxP+1, ncol = maxp)
  # V : uniforms [0,1] to compute the poisson process of intensity M
  V = matrix(allV[1:((maxP+1)*(maxk+1)*maxp)], nrow = maxP+1, ncol = maxp*(maxk+1))
  # Random normal variables to compute the brownian bridges - see Beskos et al. (2009)
  N1 = matrix(allN1[1:((maxP+1)*(maxk+1)*maxp)], nrow = maxP+1, ncol = maxp*(maxk+1))
  N2 = matrix(allN2[1:((maxP+1)*(maxk+1)*maxp)], nrow = maxP+1, ncol = maxp*(maxk+1))
  N3 = matrix(allN3[1:((maxP+1)*(maxk+1)*maxp)], nrow = maxP+1, ncol = maxp*(maxk+1))
  
  if (P > 0) {
    # Scale the D's adequately
    D = allD[1:P]	
  }
  
  K = rep(0,times = P + 1)
  w = matrix(nrow = P+1, ncol = 1)
  # For each jump time
  for (i in 2:(P+2)) {
    # Compute the number of exit times within the interval [T[i-1],T[i]]
    K[i-1] = findInterval(T[i]-T[i-1], cumsum(Ecal[i-1,]))
    
    # If you have more exit times than w allow, make w bigger
    if (K[i-1] + 1 > ncol(w)) {
      w = cbind(w, matrix(nrow = P+1, ncol = K[i-1] + 1 - ncol(w)))	
    };
    
    # Each time we reach an exit time, determine wether the BM hit -L (U < 0.5) or L (U > 0.5)
    for (j in 1:(K[i-1]+1)) {
      w[i-1,j] = if(U[i-1,j] > 0.5) L else -L;
    }
  }
  
  # Y_minus and Y_plus are equal to x - the starting point of the process Y at t=0
  Y_minus = Y_plus = c(x)
  
  # log_Psi will contain the logarithm of Psi_hat according to formula (15) in the paper
  log_Psi = 0
  
  # for each jump time (or on [0,t] is n=1) we can compute log Psi_hat 
  for (n in 1:(P+1)) {
    
    # load the randon variables for n
    temp_V = matrix(V[n,], nrow = maxk+1, ncol = maxp, byrow = FALSE)
    temp_N1 = matrix(N1[n,], nrow = maxk+1, ncol = maxp, byrow = FALSE)
    temp_N2 = matrix(N2[n,], nrow = maxk+1, ncol = maxp, byrow = FALSE)
    temp_N3 = matrix(N3[n,], nrow = maxk+1, ncol = maxp, byrow = FALSE)
    
    # log_E = log f_hat (see expression (19)
    log_E = 0
    
    # y_alg is either x or YTi (Y at the jump time Ti, right after the jump)
    y_alg = c(Y_plus[n])
    
    # for each exit time - K[n] = I in the paper
    for (k in 1:(K[n]+1)) {
      
      # Evaluate y_alg at the exit time - for it to be our starting point at the next exit time
      y_alg[k+1] = y_alg[k] + w[n,k]
      
      # Find m_i(theta) - algorithm 5-2 step (i)
      m = m_nk(y_alg[k], theta, L, x0)
      
      # Find the max of b(y;theta) - algorithm 5-2 step (i)
      M = M_nk(y_alg[k], theta, L, x0)
      
      # btime is minimum of first jump time and first localization time
      # etime is first localization time
      if (K[n] == 0) {
        # If there is no exit time, then we need to evaluate values of Y on [T[n],T[n+1]]
        # However we still need etime to simulate Y until the BM leave [-L,L]
        btime = T[n+1]-T[n]
        etime = Ecal[n,1]
      } else {
        # Else, then we need to evaluate values of Y until the first exit time
        btime = if (k < K[n]+1) Ecal[n,k] else T[n+1]-sum(Ecal[n,1:(k-1)]) - T[n]
        etime = Ecal[n,k]
      }
      
      # Generate the taus - algorithm 5-2 step (ii)
      # They are jump times of a dominating poisson process on "[0, eta]" with rate M_i(theta) (M-m in the code)
      p = 1
      # This is the method used to draw Poisson jump times
      iat = -log(temp_V[k,p]) / (M - m)
      tau = c(0)
      # Generate Poisson jump times until we reach the exit time
      while (tau[p] + iat <= btime) {
        p = p + 1	
        
        tau[p] = tau[p-1] + iat
        if(p<=maxp){
          iat = -log(temp_V[k,p]) / (M - m)
        } else{
          # In case we did not generate enough random variables, we then need to generate more N1, N2 and N3 to generate the Brownian Bridge at the tau's
          # maxp should be increased then, so the algorithm has to stop
          stop("You need to increase maxp - the number of uniforms generated in allV, all_N1, all_N2 and all_N3 is too low")
        } 
      }	
      
      # Combine the jump times tau, btime and the first localization if this happens after btime
      # Sort so that the time is ascending
      tau2 = if (btime!=etime) sort(c(tau, btime, etime)) else sort(c(tau, btime))
      
      # Detect the times that we will keep when computing log_E
      Lambda = (tau2 < btime)*1
      
      # Compute the Brownian Bridges at the different taus - algorithm 5-2 step (iii)
      # See Beskos et al. (2009) for the method.
      le = length(tau2)
      
      b1 = c(0)
      b2 = c(0)
      b3 = c(0)
      
      for (i in 2:le) {		
        b1[i] = sum(temp_N1[k,1:(i-1)] * sqrt( ((etime - tau2[i])^2)*(tau2[2:i] - tau2[1:(i-1)]) / ( etime*(etime - tau2[2:i])*(etime - tau2[1:(i-1)]) ) ))
        b2[i] = sum(temp_N2[k,1:(i-1)] * sqrt( ((etime - tau2[i])^2)*(tau2[2:i] - tau2[1:(i-1)]) / ( etime*(etime - tau2[2:i])*(etime - tau2[1:(i-1)]) ) ))
        b3[i] = sum(temp_N3[k,1:(i-1)] * sqrt( ((etime - tau2[i])^2)*(tau2[2:i] - tau2[1:(i-1)]) / ( etime*(etime - tau2[2:i])*(etime - tau2[1:(i-1)]) ) ))
      }
      
      b1[1] = b1[le] = b2[1] = b2[le] = b3[1] = b3[le] = 0;
      
      temp = abs(w[n,k])*(etime - tau2)/etime
      B = ((temp + sqrt(etime)*b1)^2) + etime*(b2^2 + b3^2)
      
      # Take into account that the BM is equal to L or -L at exit time
      W = if (w[n,k] < 0) w[n,k] + sqrt(B) else w[n,k] - sqrt(B);
      W[1] = 0
      
      # Step (iv) of Algorithm (5.2) 
      log_e = -m*btime
      
      # Computing the remaining of the expression (19) (sum over j of log(1 - (b-m)/M)) 
      log_E = log_E + log_e
      temp = log(1 - (( b_alg(y_alg[k] + W, theta, x0) - m)/(M - m)))*Lambda
      temp[is.na(temp)] = 0 
      log_E = log_E + sum(temp)
      
    }
    
    # compute the Y_T[n]-, the value of the process Y right before the jump at T[n]
    k = K[n] + 1
    p = max(which(Lambda > 0)) + 1
    Y_minus[n+1] = y_alg[k] + W[p] + rho*tau2[p]
    
    # Compute the value of Y if we evaluate it at a jump time 
    # Y_T[n] = Y_T[n]- + Gamma_Y(Y_T[n]-,D_n,theta)
    Y_plus[n+1] = if (n <= P) Y_minus[n+1] + Gamma_Y(Y_minus[n+1], D[n], theta, x0) else Y_minus[n+1]
    
    # Update log_Psi - do not add log_c_alg(Y_minus[,n+1], allD[n], theta, x0) if we are not at a jump time
    log_Psi = if (n <= P) log_Psi + log_E + log_c_alg(Y_minus[n+1], D[n], theta, x0) else log_Psi + log_E; 
  }
  
  # Add the second term in expression (20) - the Gaussian Term - and -0.5*log(Delta-T_N_delta)
  # T[P+1] is the time of the last jump
  # Y_plus[,P+1] is the value of Y at the last jump time
  log_Psi = log_Psi - (((y - Y_plus[P+1])^2)/(2*(t - T[P+1]))) - 0.5*log(t - T[P+1])
  
  log_Psi
  
}
})

#############################################################################################################
# These functions compute the optimal allocation of the algorithm as in Section 5

# Solves the optimization problem in equation (22)
opt_ell_L_rho = cmpfun({function(ell_L_rho, Delta, theta0, lbounds, ubounds, numworkers) {
  
  # Compute maximum variance times running time
  res = optim(ell_L_rho, max_var_rt, t2 = Delta, theta = theta0, lb=lbounds, ub = ubounds, nw2 = numworkers, method = "Nelder-Mead")
  
  c(res$value, res$par)
  
}
})


# Objective function of the optimization problem in equation (22) based on K = 1000 Monte Carlo samples of the random element R
max_var_rt = cmpfun({function(ellLrho, t2, theta, lb, ub, nw2) {
  
  ell1 = ellLrho[1]
  L1 = ellLrho[2]
  rho1 = ellLrho[3]
  
  hit_ub = min(ub[1] - ell1, ub[2] - L1, ub[3] - rho1)
  hit_lb = min(ell1 - lb[1], L1 - lb[2], rho1 - lb[3])
  
  if ((hit_lb >= 0) && (hit_ub >= 0)) {  
    
    # Make global the values of L, l and rho to compute the second moment of the transition density
    assign("ell", ell1, envir = .GlobalEnv)
    assign("L", L1, envir = .GlobalEnv)
    assign("rho", rho1, envir = .GlobalEnv)
    
    # Generate random elements
    K = 1000
    PT2 = gen_P_T(t2, K, ell)
    all_E2 = gen_E(PT2, L)
    # this functions checks if maxk is high enough and stops the code if not
    checkmaxk(PT2,all_E2)
    all_U2 = gen_U(PT2)
    all_V2 = gen_V(PT2)
    all_N12 = gen_N(PT2)
    all_N22 = gen_N(PT2)
    all_N32 = gen_N(PT2)
    all_D2 = gen_D(PT2)
    
    par = c(theta[2], theta)
    
    # Compute maximum variance times running time
    mvt = max_var_rt_loc(par, t2, nw2, PT2, all_E2, all_U2, all_V2, all_N12, all_N22, all_N32, all_D2) 
    
  } else {
    mvt = Inf
  }
  
  mvt
  
}
})

# This function computes the logarithm of p_hat^2(v,w;theta)*R(v,w;theta) in equation (22)
max_var_rt_loc = cmpfun({function(pars, t, nw, PT, all_E, all_U, all_V, all_N1, all_N2, all_N3, all_D) {
  
  # generate random variables used as the discrete process estimated
  v = pars[1]
  theta = pars[2:length(pars)]
  w = seq(from = v-theta[3]*qnorm(0.995), to = v+theta[3]*qnorm(0.995), length = 100)
  v = rep(v, times = length(w))
  
  # if you hit lower bound on the v's generated
  hit_ub = min(0.1-v)
  hit_lb = min(v)
  
  if ((hit_lb >= 0) && (hit_ub >= 0)) {
    
    # Compute run time to generate the logarithm of the simulated density squared and its maximum
    time = system.time({mlv = max_log_var(w, t, v, theta, nw, PT, all_E, all_U, all_V, all_N1, all_N2, all_N3, all_D)})
    
    # add the log of the second moment and the log of the time it takes to compute it
    res = mlv + log(as.double(time[3]))
  } else {
    res = -Inf
  }
  
  res
}
})

# Compute the maximum log second moment over all considered values of w in equation (22)
max_log_var = cmpfun({function(wvec, Delta, vvec, theta, nw, PT2 = PT, all_E2 = all_E, all_U2 = all_U, all_V2 = all_V, all_N12 = all_N1, all_N22 = all_N2, all_N32 = all_N3, all_D2 = all_D) {
  
  var_vec = log_hat_V_vec2(vvec, Delta, wvec, theta, nw, PT2, all_E2, all_U2, all_V2, all_N12, all_N22, all_N32, all_D2)
  
  max(var_vec)
  
}
})

# This and the functions are a slight modification of Algorithm 5.1 that computes the second moment of the density estimator needed in equation (22)
log_hat_V_vec2 = cmpfun({function(v, t, w, theta, nw, PT1 = PT, all_E1 = all_E, all_U1 = all_U, all_V1 = all_V, all_N11 = all_N1, all_N21 = all_N2, all_N31 = all_N3, all_D1 = all_D) {
  
  # Transform v and w using the Lamperti Transform
  x_vec = F(v, theta, v[1])
  y_vec = F(w, theta, v[1])
  # Maximum number of jumps generated in the ramdom numbers
  maxP = max(PT1[,1])
  # Number of MC simulations
  K = nrow(PT1)
  
  # this vector is used for the parallel processing function
  values <- 1:K
  ## Number of workers to use:
  numWorkers <- nw
  ## Parallel calculation (mclapply):
  # return a list of K vectors, each returned by log_alg_5_1_vec
  # log_alg_5_1_vec compute log_psi according to formula (15) in the paper
  paral <- mclapply(values, function(z) 2*log_alg_5_1_vec(t, x_vec, y_vec, theta, v[1], PT1[z,], all_E1[z,], all_U1[z,], all_V1[z,], all_N11[z,], all_N21[z,], all_N31[z,], all_D1[z,], maxP), mc.cores = numWorkers)
  
  # make the list a matrix, by taking each element of the list as a row of the matrix
  log_psi_square = do.call(rbind, paral)
  
  # test wether there are any errors returned 
  Test = rep(0,K)
  for(k in 1:K){
    # if makk is too low  
    if(is.character(unlist(paral[k]))==TRUE)
    {
      stop("You need to increase maxp - the number of uniforms generated in allV, all_N1, all_N2 and all_N3 is too low")
    }
  }
  
  # log_td is the log of the transition density according to formula (12) in the paper
  temp = 2*a_alg(y_vec, theta, v[1]) - 2*a_alg(x_vec, theta, v[1]) - 2*log(sigma(w, theta)) - log(2*pi)
  log_td_square = log_psi_square + matrix(rep(temp,each=K), nrow=K)
  
  # compute the expectation of the likelihood function
  res_square = colSums(exp(log_td_square))/K
  
  # take the log of that expectation to get the log likelihood function of the density squared
  log(res_square)
  
}
})
