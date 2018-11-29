##############################################################################################################
### These functions implement the random variables needed to compute the density estimator of the paper "Simulated Likelihood Estimators for Discretely Observed Jump-Diffusions" by Kay Giesecke and Gustavo Schwenkler. The codes were implemented by the authors under special assistance of Francois Guay at Boston University. The copyright to these codes remains with the authors.
##############################################################################################################

gen_P_T = function(tau, K, ell) {
  
  # Generate K poisson random variable with intensity l*tau
  P_here = rpois(K, ell*tau)
  
  # Generate the jump times between [0,tau]
  res = matrix(tau, nrow = K, ncol = max(P_here)+2)
  res[,1] = P_here
  for (k in 1:K) {
    # p is the number of jumps, drawn from Poisson(l*tau)
    p = res[k,1]
    if (p > 0) {
      # generate p numbers from uniform distribution on [0,tau], and sort them
      res[k,2:(p+1)] = sort(runif(p,0,tau))
    }
  }
  
  res
}

# Compute the exit times - See Chen & Huang (2013)
gen_E = function(PT, L, max_k = maxk) {
  
  # maximum number of jumps within the K MC samples
  sample_P = max(PT[,1]) 
  # number of MC samples
  K = nrow(PT)
  
  # Returns a vector of exit times, of length (sample_P+1)*(max_k+1)
  sample = matrix((L^2)*gen_zeta1(K*(sample_P+1)*(max_k+1)), nrow = K, ncol = (sample_P+1)*(max_k+1), byrow = TRUE)
  sample
  
}

# Generate the exit times following section 4.1 of Chen & Huang (2013)
gen_zeta1 = function(hm) {
  zetas = rep(0, times = hm)
  for (i in 1:hm) {
    accepted = 0;
    while (!accepted) {
      X = rgamma(1, shape = 1.088870, scale = 0.810570);
      Y = runif(1)*1.243707*pgamma(X, shape = 1.088870, scale = 0.810570);
      sqrt2piX3 = sqrt(2*pi*X^3);
      N = max(ceiling(0.275*X), 3)
      K = 1+2*(-N:N);
      fN0 = sum(((-1)^(-N:N))*K*exp(-K^2/(2*X)))/sqrt2piX3;
      N = N + 1;
      fN1 = fN0 + ((-1)^N)*((1-2*N)*exp(-(1-2*N)^2/(2*X))+(1+2*N)*exp(-(1+2*N)^2/(2*X)))/sqrt2piX3;
      while (sign((Y - fN0)*(Y - fN1)) == -1) {
        fN0 = fN1;
        N = N + 1;
        fN1 = fN0 + (-1)^N*((1-2*N)*exp(-(1-2*N)^2/(2*X)) + (1+2*N)*exp(-(1+2*N)^2/(2*X)))/sqrt2piX3;
      }
      if (Y <= fN1) {
        accepted = 1;	
        zetas[i] = X
      }
    }
  }	
  zetas
}

# Check if max_k is big enough
checkmaxk = function(PT, allE, max_k = maxk ){
  
  K = nrow(PT)
  maxP = max(PT[,1])
  count=0
  
  for(k in 1:K){
    
    TJ = c(0,PT[k,2:(maxP+2)])
    ET = matrix(allE[1:((maxP+1)*(maxk+1))], nrow = maxP+1, ncol = maxk+1)
    
    for(i in 1:(maxP+1)){
      
      b = cumsum(ET[i,])
      # check if the last exit time generated happens before the next jump
      # If yes, then we need to generate more exit times
      if(b[maxk+1]<TJ[i+1]-TJ[i]){
        count=count+1
      }
    }
  } 
  # count is the number of times the violation happens
  if(count>0){
    print(count)
    stop('Not enough exit times generated - please increase maxk')
  }
}

# Generate uniforms to determine the value of the Brownian motion when it exits [-L,L] - See section 4.2 of Chen & Huang (2013)
gen_U = function(PT, max_k = maxk) {
  
  sample_P = max(PT[,1])
  K = nrow(PT)
  
  sample = matrix(runif(K*(sample_P+1)*(max_k+1)), nrow = K, ncol = (sample_P+1)*(max_k+1), byrow = TRUE)
  sample
  
}

# Generate uniforms random variables to compute the jump times of intensity M_i (for algorithm 4.4 part (ii))
gen_V = function(PT, max_k = maxk, max_p = maxp) {
  
  sample_P = max(PT[,1])
  K = nrow(PT)
  
  sample = matrix(runif(K*(sample_P+1)*(max_k+1)*max_p), nrow = K, ncol = (sample_P+1)*(max_k+1)*max_p, byrow = TRUE)
}

# Generation the standard normal random variables to compute the brownian bridges (for algorithm 4.4 part (iii))
gen_N = function(PT, max_k = maxk, max_p = maxp) {
  
  sample_P = max(PT[,1])
  K = nrow(PT)
  
  sample = matrix(rnorm(K*(sample_P+1)*(max_k+1)*max_p), nrow = K, ncol = (sample_P+1)*(max_k+1)*max_p, byrow = TRUE)
  sample
}

# Generate the mark variables from the density pi (standard normal here)
gen_D = function(PT) {
  
  sample_P = max(PT[,1])
  K = nrow(PT)
  sample = matrix(rnorm(K*(sample_P+1), 0, 1), nrow = K, ncol = sample_P+1, byrow = TRUE)
  sample
  
}

