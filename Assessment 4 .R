## Things to do
## 1. Wolfe conditions
## 2. Reduce complexity for approximating Hessian
## 3. Errors or warnings
## 4. Gradient attribute - meaning and if not suppled
## 5. eps value
## 6. Comments






## Work Group 1 - Sanika Baxi - s2159255		Ayush Oza - s2184992		Gowtham Palepu - s2113890

## BFGS Optimizer

## Overview

##*****************************************************************************************************

## The BFGS Optimizer is an iterative method for solving non-linear optimization problems. The BFGS determines the 
## direction of descent by preconditioning the gradient with curvature information. This is achieved by gradually improving 
## an approximation to the Hessian matrix of the loss function obtained by gradient evaluation. This algorithm is named after
## its inventors Broyden–Fletcher–Goldfarb–Shanno (BFGS). 

## Details of the BFGS Optimizer below: 

## The function bfgs is used to optimize the given objective function by taking in various arguments as described below: 
## The bfgs function - bfgs(theta, f, ..., tol = 1e-5, fscale = 1, maxit = 100)
## theta is the vector of initial values for our optimization problem
## f is the objective function argument. The first argument of function is a vector of optimization parameters and second
## argument is a logical indicating if gradient is to be calculated
## ... is for any arguments of f that are to be passed to it after the above two
## fscale indicates a rough estimate of magnitude of f at the optimum
## maxit is the maximum number of BFGS iterations taken by the function before stopping further optimization

## The BFGS function above returns the below values upon calling it with valid arguments: 

## f being the scalar value of the objective function at the minimum
## theta is the vector of values of parameters at the minimum
## iter denoting the number of iterations taken by the BFGS to reach minimum
## g denoting the gradient vector at the minimum
## H denoting the approximated Hessian matrix at the minimum

## *********************************************************************************************************

## Defining and calculating the finite differencing function to calculate the gradient/differentiation of the objective function 
## It takes in the objective function and theta as its arguments and returns the value after calculating the finite difference/first order derivative

fd1d <- function(f=f , theta) {          ## Declaring and defining the function fd1d with necessary arguments
  eps <- 1e-7                            ## Using epsilon value as 1e-7 as delta for differencing
  
  ## Initiating the differencing with the starting point as obtained from theta
  
  f0 <- f(theta)                         ## Initiating f0 starting point of the function            
  fd1 <- theta                           ## Initiating fd1 with theta vector of starting parameters

  ## The below loop is to calculate the finite difference/first order derivative of the objective function f
  ## First, we create a vector of 0s and then perturbate with the epsilon value while performing the differencing
  ## The finite differencing formula can be stated as fd = (f1 - f0)/eps where f1 is next value of the curve obtained 
  ## by using the theta and eps values and f0 is the current point on the function
  
  for (i in 1:length(theta)) {          ## Looping over the length of theta to calculate finite differences
    pet <- rep(0, length(theta))        ## Creating a vector of 0s to perturbate with eps values 
    pet[i] <- eps                       ## Performing perturbation
    x1 <- theta + pet                   ## Determining/approximating to the next point on the function
    f1 <- f(x1)                         ## Calculating function value at the next point
    fd1[i] <- (f1 - f0)/eps             ## Calculating finite difference between new point on the function with the old point
  } ## End of loop 
  return(fd1)                           ## Finally returning the obtained value of the finite difference
}


## By the end of above function, we will have successfuly calculated the first order derivative for a given objective function
## with given theta parameters as starting points. We now use the above function in BFGS algorithm below for various approximations 

## Start of BFGS algorithm taking in arguments as descibed in the overview

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100) {     ## Declaring and defining the BFGS function 
  
  ## Since calculating Hessian matrices is complex and costly, we use the first order gradients obtained above from fd1d to approximate
  ## our Hessian matrix. This is achieved by combining gradient information over several steps to build up an approximate to Hessian. 
  
  ## Here eps is again used to perform the differencing. Then we use the theta details to perform further algorithm iterations. 
  ## c1 and c2 are the step length conditions chosen astisfying the Wolfe conditions. 

  eps <- 1e-7                                       ## Using eps as 1e-07 for differencing           
  thetaold <- theta                                 ## Using thetaold variable to store the previous value of the function. Initially given theta
  n <- length(theta)                                ## Initializing n to be length of theta to create a diagonal matrix
  B <- I <- diag(n)                                 ## Initializing B and I diagonal matrices of n where B is inverse Hessian matrix 
  c1 <- 0.5                                         ## Initializing c1 and c2 step lengths chosen satisfying Wolfe conditions
  c2 <- 0.9
  counter <- 0                                      ## Initializing a counter to return the number of iterations for optimization 

  ## Now that we have all variables required in place, we start to loop over the maxit value and start the optimization process
  
  for (i in 1:maxit) {                              ## Start of the loop over maxit = 100 condition to stop optimizing if maxit is reached before
    
	## gradold initialized to store gradient of the function at previous point on the function
    
	gradold <- fd1d(f, thetaold)                    ## gradient is obatained from fd1d function above              
    
	## The Quasi-Newton step where step = -InvHessian * gradient 
	
    step <- -1 *(B %*% gradold)                     ## Determining the step using gradient value obtained from finite difference
    
	## Determining the next point on the function using the step obtained above 
	
    thetanew <- thetaold + drop(step)               ## Storing the new position in thetanew variable 
    
    # if (f(thetanew) > f(thetaold) + c1 * (t(gradold) %*% step) 
    #  t(fd1d(thetaold + step)) %*% step <  c2 * (t(fd1d(f, thetaold)) %*% step) ) {
      
    #}
	
	## From here all the steps below are calculated to approimate the Hessian using the gradient values. We calculate all necessary parameters 
	## use them in performing the iterations of BGS
    
	## Calculating difference between the new and old theta values
	
    s <- thetanew - thetaold                       
    
	## Calculating the new gradient with new theta as obtained from step above using fd1d function again
    	
    gradnew <- fd1d(f, thetanew) 

    ## Calculating difference between the new and old gradients
     	
    y <- gradnew - gradold
	
	## Calculaitng rho for BFGS update condition where rho = transpost(s) * y
 
    rho <- drop(1/(t(s) %*% y))                 ## Calculating rho value from s and y as per definition 
	
	## Approximating the Hessian B where B is given below as per BFGS algorithm: 
	
	## B = (I - rho * s * transpose(y)) * B * (I - rho*(y * transpose(s))) + (rho * s * transpose(s))
	## Determining B as defined above using matrix multiplications as applicable

    B <- ((I - rho*(s %*% t(y))) %*% B %*% (I - rho*(y %*% t(s)))) + (rho * s %*% t(s))    
  
    ## Updated B 
	
    counter = counter + 1                     ## Updating counter for each iteration
    
	## Checking for the first Wolfe condition
	
    if (max(abs(gradnew)) < (abs(f(thetanew))+fscale)*tol) {
      break
    }
    
	## Updating thetaold with thetanew for next iteration 

    thetaold <- thetanew
    
  }
  

  ## Calculating approximate Hessian matrix 

  Hfd <- matrix(0,n,n)                                ## Initiating all 0s finite difference Hessain matrix 
  for (i in 1:n) {                                    ## Looping over the parameter length value to iterate over the matrix and calculate entries
    th1 <- thetanew; th1[i] <- th1[i] + eps           ## Updating th1 to thetanew obtained from above BFGS update and the eps value  
    g1 <- fd1d(f, th1)                                ## Calculating the finite difference of the function with th1 value 
    Hfd[i,] <-(g1 - gradnew)/eps                      ## Calculating/approximating second order derivative Hessian matrix  
  }
  
  Hfd <-0.5 * (t(Hfd) + Hfd)                          ## Obtaining the final finite difference approximated Hessian matrix
  
  ## Returning all the calculated values as required by the BFG function as a list - function value at thetanew, thetanew, 
  ## number of iterations taken for approximation, newgradient, Hessian 
  
  return(list(f = f(thetanew), theta = thetanew, iter = counter, g = gradnew, H = Hfd))        ## Return list
  
}


## Test Code for BFGS as suggested with Rosenbrick objective function

bfgs(theta = c(-1,2), f = rb, maxit = 100, tol = 1e-5)     ## Calling bfgs function with required arguments
fd1d(f = rb, theta = c(-1,2))                              ## Calculating finite difference 

## Defining Rosenbrock function


rb <-function(theta, getg=FALSE,k=10) {
  ## Rosenbrock objective function, suitable for use by ’bfgs’
  z <-theta[1]; x <-theta[2]
  f <-k*(z-x^2)^2 + (1-x)^2 + 1
  if (getg) {
    attr(f,"gradient") <-c(2*k*(z-x^2),
    -4*k*x*(z-x^2) -2*(1-x)) 
    } 
    f
} ## rb





## Rubbish messing around/workings


new <-function(theta,getg=FALSE,k=10) {
  ## Rosenbrock objective function, suitable for use by ’bfgs’
  z <-theta[1]
  f <- exp(z) - 2 * z^2 + 3 * z - 1
  } 
  f
} ## rb


rb(c(-1,2))

bfgs(theta = c(-1,2), rb, maxit=20)

rb(theta = c(-1,2))

fd1d(rb, c(-1,2))

theta <- c(-1,2)
thetaold <- theta 
B <- I <- diag(length(theta))
c1 <- 0.5
c2 <- 0.9

for (i in 1:20) {
  
  step <- -B %*% fd1d(rb, thetaold)
  
  #if (f(thetaold + step) > f(thetaold) + c1 * (t(fd1d(f, thetaold)) %*% step) 
  #  | t(fd1d(thetaold + step)) %*% step <  c2 * (t(fd1d(f, thetaold)) %*% step) ) {
  
  #}
  cat("step ",step)
  thetanew <- thetaold + step
  cat("thetanew ", thetanew)
  s <- thetanew - thetaold
  cat("s ", s)
  y <- fd1d(f, thetanew) - fd1d(rb, thetaold)
  cat("y ", y)
  rho <- drop(1/(t(s) %*% s))
  cat("rho ", rho)
  B <- ((I - rho*(s %*% t(y))) %*% B %*% (I - rho*(s %*% t(y)))) - drop((rho * t(s) %*% s)) 
  cat("B ", B)
  thetaold <- thetanew
}
print(thetanew)


f <- function(z, x){
  f <-10*(z-x^2)^2 + (1-x)^2 + 1
}

f <- function(x) {
  return(exp(x) - 2 * x^2 + 3 * x - 1)
}


add <- function(a=a,b, ...){
  print(e)
  return(a(...)+b)
}

d <- function(e,f){
  return(e*f)
}


add(a=d, e = 2, f = 3, b = 4)


d()

for (i in 1:2) {
  print(i)
}
