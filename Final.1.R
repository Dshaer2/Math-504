## Part A
library(Matrix)
nls.data <- read.table("C:/Users/562644/Documents/~ Education/Georgetown/2015 Spring/504 - Numerical Methods/HW6/nonlinear_least_squares.txt", header=TRUE, quote="\"")

x <- nls.data$x
y <- nls.data$y

## Normalizing function
norm.func <- function(x){
  norm <- sqrt(sum(x^2))
  return(norm)
}

## Loss Function
fx = function(C){
  d <- C[1]
  r <- C[2]  
  fx <- y - d*exp(-r*x)
  return(sum(fx^2)) ## CHANGE *** not sqrt too 
}

## Gradient Function
grad.fx = function(C){
  d <- C[1]
  r <- C[2]
  grad.d <- sum(-2*(y*exp(-r*x) - d*exp(-r*x)*exp(-r*x))) ## CHANGE
  grad.r <- sum(2*(y*d*x*exp(-r*x) - d^2 * x * exp(-r*x)*exp(-r*x))) ## CHANGE
  grad <- c(grad.d, grad.r)
  return(grad)
}

## Hessian Function
hess.fx <- function(C){
  d <- C[1]
  r <- C[2]
  drr <- sum(-2*y*d*(x^2)*exp(-r*x) + 4*(d^2)*(x^2)*exp(-2*r*x))
  ddr <- sum(2*y*x*exp(-r*x) - 4*d*x*exp(-2*r*x))
  ddd <- sum(2*exp(-2*r*x))
  hess.fx <- matrix(c(ddd,ddr,ddr,drr), nrow=2, ncol=2) ## CHANGED Order
  return(hess.fx)
}

## Modified Newton Method

newton.mod = function(C, eps, adjustment){
  path <- matrix(c(0,0), nrow=1)
  iterations <- 0
  while (norm.func(grad.fx(C)) > eps) {
    step <- 1
    M <- hess.fx(C)
    evs <- eigen(M)$values
    lambda <- min(evs)
    
    ## If eigenvalues are negative, shift by adjustment
    if (lambda < 0) M <- (abs(lambda) + adjustment)*diag(rep(1,nrow(M))) + M
   
    d <- -solve(M)%*%grad.fx(C)
    
    ## Backtracking
    while(fx(C + step*d) > fx(C)){
      step = step/2
    }
    
    path <- rbind(t(C), path) 
    C <- as.vector(C + step*d)
    iterations <- iterations + 1
  }
  return(path[1,])
}


C <- c(1,1)
eps <- .001
adjustment <- .0001

results <- newton.mod(C, eps, adjustment)      ## Didn't even get this far in the homework
d <- results[1]
r <- results[2]

fx.1 = function(x,y){
  y <- d*exp(-r*x)
  return(y)
}

summary(nls.data)
curve(fx.1, from = 0, to = 15, col = "red", ylab = "y", main = "Newton Method Fit", ylim=c(0,1700))
points(x, y, pch=1, cex = .5, col = "dark gray")


debug(newton.mod)
undebug(newton.mod)

## b

nls.fit <- nls(y ~ d*exp(-r*x), start=list(d=2000,r=.2))

plot(x, predict(nls.fit, newx=seq(0,15,0.01)), ylim=c(0,1700), main = "nls Fit", ylab = "y")
points(x, y, pch=1, cex = .5, col = "dark gray")

nls.fit
results

## Note the number of iterations to convergence 4 vs 16