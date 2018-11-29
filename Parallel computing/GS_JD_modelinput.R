#####################################################################################################
#   Model specification to be given as an input by the end-user for the density estimator of the paper "Simulated Likelihood Estimators for Discretely Observed Jump-Diffusions" by Kay Giesecke and Gustavo Schwenkler. The original version of this file implemented by the authors gives the model specification for the numerical case study in Section 7 of the paper. 
#####################################################################################################


# Enter the drift function mu_X(x;theta) of X 
mu_X = cmpfun({function(x, theta) {
  ka = theta[1]
  Xb = theta[2]
  si = theta[3]
  la0 = theta[4]
  la1 = theta[5]
  ga1 = theta[6]
  ga2 = theta[7]
  
  ka*(Xb - x)
  
}
})


# Enter the first derivative of the drift function mu_X(x;theta) with respect to x if available in closed form
mu_X_prime = cmpfun({function(x, theta) {
  ka = theta[1]
  Xb = theta[2]
  si = theta[3]
  la0 = theta[4]
  la1 = theta[5]
  ga1 = theta[6]
  ga2 = theta[7]
  
  -ka 
}
})


# Enter the volatility function sigma(x;theta) of X
sigma = cmpfun({function(x, theta) {
  ka = theta[1]
  Xb = theta[2]
  si = theta[3]
  la0 = theta[4]
  la1 = theta[5]
  ga1 = theta[6]
  ga2 = theta[7]
  
  rep(si, times = length(x))
  
}
})


# Enter the first derivative of the volatility function sigma(x;theta) with respect to x if available in closed form. This is needed to compute mu_Y(y,theta)
sigma_prime = cmpfun({function(x, theta) {
  ka = theta[1]
  Xb = theta[2]
  si = theta[3]
  la0 = theta[4]
  la1 = theta[5]
  ga1 = theta[6]
  ga2 = theta[7]
  
  rep(0, times = length(x))
  
}
})


# Enter the second derivative of the volatility function sigma(x;theta) with respect to x if available in closed form. This is needed to compute b(y;theta)
sigma_prime_prime = cmpfun({function(x, theta) {
  ka = theta[1]
  Xb = theta[2]
  si = theta[3]
  la0 = theta[4]
  la1 = theta[5]
  ga1 = theta[6]
  ga2 = theta[7]
  
  rep(0, times = length(x))
  
}
})


# Enter the intensity function Lambda(x,theta) driving the jump times of the process X
Lambda_X = cmpfun({function(x, theta) {
  ka = theta[1]
  Xb = theta[2]
  si = theta[3]
  la0 = theta[4]
  la1 = theta[5]
  ga1 = theta[6]
  ga2 = theta[7]
  
  res = rep(0, times = length(x))
  # res return the jump intensity, which cannot be negative or null
  for (j in 1:length(x)) {
    res[j] = max(la0 + la1*x[j],0)
  }
  
  res	
  
}
})


# Enter the function Gamma_X(x,d;theta) specifying the jump size of X when the value of X right before the jump is x and the mark variable is d
Gamma_X = cmpfun({function(x,d, theta) {
  ka = theta[1]
  Xb = theta[2]
  si = theta[3]
  la0 = theta[4]
  la1 = theta[5]
  ga1 = theta[6]
  ga2 = theta[7]
  
  
  (ga1 + ga2*d)
  
}
})


# Enter a closed-form expression for the Lamperti transform x = F(w;theta,x0) in equation (6) if available
# Otherwise you can replace the expression by the following numerical integral:
# sig <- function(x) sigma(x,theta)
# integrate(sig,lower=x0,upper=w)
F = cmpfun({function(w, theta, x0) {
  ka = theta[1]
  Xb = theta[2]
  si = theta[3]
  la0 = theta[4]
  la1 = theta[5]
  ga1 = theta[6]
  ga2 = theta[7]
  
  # if available in closed-form
  (w - x0)/si
  
}
})

# Enter a closed-form expresion for the derivative of the Lamperti transform F(w;theta,x0) with respect to w if available
# Otherwise, you can replace the expression by the numerical derivative:
# Ffunc <- function(v) F(v,theta,x0)
# grad(Ffunc,w)
F_prime = cmpfun({function(w,theta,x0) {
  ka = theta[1]
  Xb = theta[2]
  si = theta[3]
  la0 = theta[4]
  la1 = theta[5]
  ga1 = theta[6]
  ga2 = theta[7]
  
  # if available in closed-form
  1/si
  
}
})


# Enter a closed-form expression for the inverse w = F^(-1)(y;theta) of the Lamperti transform in equation (6) if available
# Otherwise, you can find the expression by doing the following:
# solve for x in y-F(x,theta,x0)=0
# f <- function(x) y-F(x,theta,x0)
# unitroot(f, c(bound1, bound2), tol = 0.0001)
# (bound1 and bound2 need to be specified by the end user)
F_inv = cmpfun({function(y, theta, x0) {
  ka = theta[1]
  Xb = theta[2]
  si = theta[3]
  la0 = theta[4]
  la1 = theta[5]
  ga1 = theta[6]
  ga2 = theta[7]
  
  # if available in closed-form
  si*y + x0
  
  
}
})


# Enter the first partial derivative F^(-1)'(y;theta) of the inverse of the Lamperti transform in equation (6) with respect to y  if available in closed form. This is needed to compute b(y;theta)
# Otherwise, you can replace this expression by:
# 1/F_prime(F_inv(y,theta,x0),theta,x0)
F_inv_prime = cmpfun({function(y, theta, x0) {
  ka = theta[1]
  Xb = theta[2]
  si = theta[3]
  la0 = theta[4]
  la1 = theta[5]
  ga1 = theta[6]
  ga2 = theta[7]
  
  # if available in closed-form
  si
  
}
})

# IF available in closed-form, enter the integral of muY(y;theta) - It is used to compute a(y;theta) and log c(y,d;theta)
# Otherwise, you can do the following:
# if not available in closed-form
# muY <- function(y) mu_Y(y,theta,x0)
# integ = rep(0,length(bound2))
# for(i in 1:length(bound2)) {
#  integ[i] = unlist(integrate(muY,lower=bound1[i],upper=bound2[i])$value)
# }
integral_mu_Y = cmpfun({ function(bound1,bound2, x0, theta) {
  ka = theta[1]
  Xb = theta[2]
  si = theta[3]
  la0 = theta[4]
  la1 = theta[5]
  ga1 = theta[6]
  ga2 = theta[7]
  
  # if available in closed-form
  (ka*(Xb-x0)/si*bound2 - ka*0.5*bound2^2)-(ka*(Xb-x0)/si*bound1 - ka*0.5*bound1^2)
  
}
})


# m_nk computes the value of m_i in step (i) of Algorithm 5.2
m_nk = cmpfun({function(y, theta, L_here, x0) {
  ka = theta[1]
  Xb = theta[2]
  si = theta[3]
  la0 = theta[4]
  la1 = theta[5]
  ga1 = theta[6]
  ga2 = theta[7]
  
  # L_here is the upper bound for |x| 
  # wstar is a stationary point of b
  wstar = ((Xb - x0)/si) - (la1*si/(ka^2)) - y 
  
  # b(y;theta) is computed on the left(y-L_here), on the right(y+L_here) and at its stationary point wstar
  left = b_alg(y - L_here, theta, x0)
  right = b_alg(y + L_here, theta, x0)
  opt = b_alg(y + wstar, theta, x0)
  
  res = rep(0, times = length(y))
  # We want to find the minimum of b: it is reached at its boundaries or at its stationnary point
  for (j in 1:length(y)) {
    if (abs(wstar[j]) <= L_here) {
      res[j] = min(left[j], right[j], opt[j])
    } else {
      res[j] = min(left[j], right[j])	
    }
  }
  
  res
  
}
})


# M_nk computes the value of M_i in step (i) of Algorithm 5.2
M_nk = cmpfun({function(y, theta, L_here, x0) {
  ka = theta[1]
  Xb = theta[2]
  si = theta[3]
  la0 = theta[4]
  la1 = theta[5]
  ga1 = theta[6]
  ga2 = theta[7]
  
  # L_here is the upper bound for |x|
  # wstar is a stationary point of b
  wstar = ((Xb - x0)/si) - (100*la1*si/(ka^2)) - y 
  
  # b(y;theta) is computed on the left(y-L_here), on the right(y+L_here) and at its stationary point wstar
  left = b_alg(y - L_here, theta, x0)
  right = b_alg(y + L_here, theta, x0)
  opt = b_alg(y + wstar, theta, x0)
  
  res = rep(0, times = length(y))
  # We want to find the maximum of b: it is reached at its boundaries or at its stationnary point
  for (j in 1:length(y)) {
    if (abs(wstar[j]) <= L_here) {
      res[j] = max(left[j], right[j], opt[j])
    } else {
      res[j] = max(left[j], right[j])	
    }
  } 
  
  res
}
})
