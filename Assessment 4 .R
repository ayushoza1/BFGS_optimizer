fd1d <- function(f=f , theta) {
  eps <- 1e-7
  f0 <- f(theta)
  fd1 <- theta
  for (i in 1:length(theta)) {
    pet <- rep(0, length(theta))
    pet[i] <- eps
    x1 <- theta + pet
    f1 <- f(x1)
    fd1[i] <- (f1 - f0)/eps
  }
  return(fd1)
}

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100) {

  eps <- 1e-7
  thetaold <- theta 
  n <- length(theta)
  B <- I <- diag(n)
  c1 <- 0.5
  c2 <- 0.9
  counter <- 0

  for (i in 1:maxit) {
    
    gradold <- fd1d(f, thetaold) 
    
    step <- -1 *(B %*% gradold)
    
    thetanew <- thetaold + drop(step)
    
    if (f(thetanew) > f(thetaold) + c1 * (t(gradold) %*% step) 
      | t(fd1d(thetaold + step)) %*% step <  c2 * (t(fd1d(f, thetaold)) %*% step) ) {
      
    }
    
    s <- thetanew - thetaold
    
    gradnew <- fd1d(f, thetanew) 

    y <- gradnew - gradold

    rho <- drop(1/(t(s) %*% y))

    B <- ((I - rho*(s %*% t(y))) %*% B %*% (I - rho*(y %*% t(s)))) + (rho * s %*% t(s)) 
  
    counter = counter + 1
    
    if (max(abs(gradnew)) < (abs(f(thetanew))+fscale)*tol) {
      break
    }
    
    
    thetaold <- thetanew
    
  }
  

  Hfd <- matrix(0,n,n) ## finite diference Hessian 
  for (i in 1:n) { ## loop over parameters 
    th1 <- thetanew; th1[i] <- th1[i] + eps ## increase th0[i] by eps 
    g1 <- fd1d(f, th1) ## compute resulting nll 
    Hfd[i,] <-(g1 - gradnew)/eps ## approximate second derivs 
  }
  
  Hfd <-0.5 * (t(Hfd) + Hfd)
  
  
  return(list(f = f(thetanew), theta = thetanew, iter = counter, g = gradnew, H = Hfd))
  
}



bfgs(theta = c(-1,2), f = rb, maxit = 100, tol = 1e-5)



theta <- c(-1,2)
fd1d(f = rb, theta = c(-1,2))








?optim


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
