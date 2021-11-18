## Work Group 1 - Sanika Baxi - s2159255		Ayush Oza - s2184992		Gowtham Palepu - s2113890
## GitHub repo link: https://github.com/ayushoza1/SP-Assessment-4.git

## BFGS Optimizer

## Overview
## ********

## The BFGS Optimizer is an iterative method for solving non-linear optimization problems. The BFGS determines the 
## direction of descent by preconditioning the gradient with curvature information. This is achieved by gradually improving 
## an approximation to the Hessian matrix of the loss function obtained by gradient evaluation. We also Define and calculate
## the finite differencing function to calculate the gradient/differentiation of the objective function if gradient is not 
## provided. Wolfe conditions are also checked and the step length updated to ensure or objective function is decaresing and 
## that Hessian is positive definite.. 

get_grad <- function(f ,theta, ...) {
  
  ## Get_grad function 
  ## *****************
  ## The loop below is used to finite difference/first order derivative of the objective function f
  ## First, we create a vector of 0s and then perturb with the epsilon value while performing the differencing
  ## The finite differencing formula can be stated as fd = (f1 - f0)/eps where f1 is next value of the curve obtained 
  ## by using the theta and eps values and f0 is the current point on the function. f = function of gradient to be found
  ## theta = parameter of values where gradient is to be found.
  
  ## INPUT: objective function=f and theta=parameter values
  ## OUTPUT: returns gradient of the objective function calculated by finite difference or using gradient attribute = grad
  
  if (formals(f)[2] == FALSE | is.null(attr(f(theta, ...), "gradient"))) { ## Check to see if gradient calculation is needed by finite differencing
    
    eps <- sqrt(.Machine$double.eps) ## Declare peterbation value
    
    f0 <- f(theta, ...)  ## Initiating the vector of paramters where gradient is calculated            
    grad <- theta   ## Initiating a storage vector for the gradient to be passed into

    for (i in 1:length(theta)) {  ## Loop through parameters
      pet <- rep(0, length(theta))  ## Creating a vector of 0s to perturb with eps values 
      pet[i] <- eps ## Peturb the value of parameter i 
      x1 <- theta + pet ## new paramter value   
      f1 <- f(x1, ...) ## Calculate function at the new parameter value
      grad[i] <- (f1 - f0)/eps  
    } 

  } else {  ## If finite differncing is not required gradient function is provided
    
    grad <- attr(f(theta, ...), "gradient") 

  }

  return(grad)                        
}


wolf_check <- function(f ,thetaold, step, c2, ...) {
  
  ## Wolfe_check function 
  ## *******************
  ## Function checks that the objective function has been reduced and that the Hessian matrix is positive
  ## definite for the next run of the bfgs iteration. A nested while loop is used to ensuer both conditions
  ## are met and the wolfe conditions rechecked once one condition has been met. If the first condition is 
  ## not met the step length is halved and if the second condition is not met the step length is multiplied
  ## by 1.5. The function 'breaks' if no step length is found satifying both conditions after 50 iterations.
  
  ## INPUT: objective function = f, theta at old position on the function = thetaold, step length= step, wolfe constants = c2
  ## OUTPUT: step satisfying Wolfe conditions = step
  
  wolfe_1 <- f(thetaold + drop(step), ...) <= f(thetaold, ...)  ## Logical to see if objective function has been reduced
  wolfe_2 <- t(get_grad(f,thetaold + drop(step), ...)) %*% step >=  c2 * (t(get_grad(f,thetaold, ...))) %*% step ## Logical to see if updated Hessian is positive deffinite 
  wolfe_counter <- 0 ## Initialize counter to count updates of step length if wolfe conditions not met
    
  while (wolfe_1 == FALSE | wolfe_2 == FALSE) { ## Continuous loop if wolfe conditions do not hold 
    
    while (wolfe_1 == FALSE) { ## Function not being reduced
      step <- step * 0.25 ## Reduce step length by 1/4
      wolfe_counter <- wolfe_counter + 1 ## Update counter 
      wolfe_1 <- f(thetaold + drop(step), ...) <= f(thetaold, ...) ## Update logical if first wolfe condition satisfied
      if (wolfe_counter > 5) { ## break out of function if finding step length is taking too long
        stop("Stuck in loop with Wolfe conditions 1")
      } 
    }
    
    while (wolfe_2 == FALSE) { ## Hessian not positive definite
      step <- step * 1.25 ## Increase step length by 1.25
      wolfe_counter <- wolfe_counter + 1 ## Update counter 
      wolfe_2 <- t(get_grad(f,thetaold + drop(step), ...)) %*% step >=  c2 * (t(get_grad(f,thetaold, ...))) %*% step ## Update logical if second wolfe condition satisfied
      if (wolfe_counter > 5) { ## break out of function if finding step length is taking too long
        stop("Stuck in loop with Wolfe conditions 2")
        
        wolfe_1 <- f(thetaold + drop(step), ...) <= f(thetaold, ...) ## Update logicals to test both conditions have been met
        wolfe_2 <- t(get_grad(f,thetaold + drop(step), ...)) %*% step >=  c2 * (t(get_grad(f,thetaold, ...))) %*% step ## Udpate logicals to test both conditions have been met      
        
      } 
    }
    
  }
  
  return(step)
  
}




bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100) {   
  
  ## bfgs function 
  ## *************
  ## The BFGS Optimizer is an iterative method for solving non-linear optimization problems. The BFGS determines the 
  ## direction of descent by preconditioning the gradient with curvature information. This is achieved by gradually improving 
  ## an approximation to the Hessian matrix of the loss function obtained by gradient evaluation.  
  
  ## INPUT:
  ## theta is the vector of initial values for our optimization problem
  ## f is the objective function argument. The first argument of function is a vector of optimization parameters and second
  ## argument is a logical indicating if gradient is to be calculated
  ## ... is for any arguments of f that are to be passed to it after the above two
  ## tol is the tolerance value to calculate step difference 
  ## fscale indicates a rough estimate of magnitude of f at the optimum
  ## maxit is the maximum number of BFGS iterations taken by the function before stopping further optimization
  
  ## OUTPUT: 
  ## The BFGS function above returns the below values upon calling it with valid arguments: 
  ## f being the scalar value of the objective function at the minimum
  ## theta is the vector of values of parameters at the minimum
  ## iter denoting the number of iterations taken by the BFGS to reach minimum
  ## g denoting the gradient vector at the minimum
  ## H denoting the approximated Hessian matrix at the minimum


  eps <- sqrt(.Machine$double.eps) ## Declare perturbation value for use in finite differencing          
  thetaold <- theta  ## Declaring theta^k-1 
  n <- length(theta)  ## Number of paramenters 
  B <- I <- diag(n) ## Initializing by the the identity matrix                                     
  c2 <- 0.9 ## Initialize c2 for use in second Wolfe condition check
  counter <- 0  ## Initializing a counter to store number of iterations of bfgs 
  
  if (any(is.finite(f(theta, ...))) == FALSE | any(is.finite(get_grad(f, theta, ...))) == FALSE) { ## Check to see if the gradient or function value is non finite
    stop("Gradient or function is non-finite, please try different parameter values") ## Stop function if so 
  } 
  
  
  for (i in 1:maxit) {  ## Loop over the max itterations (loop is broken if convergence is found before)
  
    
	  gradold <- get_grad(f, theta = thetaold, ...) ## Gradient of old or k-1 paramter vector stored    
	
    step <- -1 *(B %*% gradold) ## Determine the Quasi-Newton step = -InvHessian * gradient or the descent direction
    
    step <- wolf_check(f,thetaold, step, c2, ...) ## Check Wolfe conditions and return new step length
    
    thetanew <- thetaold + step   ## Calculate the new parameter vector, or theta k
	
    s <- thetanew - thetaold  ## Calculate step length vector                   
    	
    gradnew <- get_grad(f,thetanew, ...) ## Calculate gradient of the new parameter vector
     	
    y <- gradnew - gradold ## Calculate difference between old and new gradients
 
    rho <- drop(1/(t(s) %*% y)) ## rho = transpose(s) * y         

    ## B evaluated to reduce computational time by  matrix/vector multiplication rather than matrix/matrix
    
    B <- B - (rho*B%*%y%*%t(s)) - rho*s%*%(t(y)%*%B) + (rho^2)*s%*%(t(y)%*%B)%*%y%*%t(s) + (rho * s %*% t(s)) 	## Approximating the Hessian B where B is given below as per BFGS algorithm:
	
    counter = counter + 1 ## Update counter for each iteration
	
    if (max(abs(gradnew)) < (abs(f(thetanew, ...))+fscale)*tol) { ## Check for convergence
      break ## break function if convergence is found
    } else if (counter == maxit) { ## If we reach max number of iterations and convergence not found then issue warning
      warning("Max number of itterations reached and convergance has not occured")
    } else if (f(thetanew, ...) - f(thetaold, ...) > 0){ ## If objective function has not been reduced issue warning 
      warning(paste("Objective function not being reduced for itteration", i))
    } 

    thetaold <- thetanew ## Update theta so theta k is now theta k-1
    
  }
  

  ## Calculating approximate Hessian matrix 

  Hfd <- matrix(0,n,n)  ## Initiate the Hessian matrix to be calculated by finite differencing
  for (i in 1:n) {  ## Looping over the parameter length value to iterate over the matrix and calculate entries
    th1 <- thetanew; th1[i] <- th1[i] + eps ## Perturb gradient vector 
    g1 <- get_grad(f,th1, ...)  ## Calculating the gradient vector 
    Hfd[i,] <-(g1 - gradnew)/eps  
  }
  
  Hfd <-0.5 * (t(Hfd) + Hfd)  ## Ensure the Hessian is symmetric
  
  ## Returning all the calculated values as required by the BFG function as a list - function value at thetanew, thetanew, 
  ## number of iterations taken for approximation, newgradient, Hessian 
  
  return(list(f = f(thetanew, ...), theta = drop(thetanew), iter = counter, g = gradnew, H = Hfd)) ## Return list
  
}