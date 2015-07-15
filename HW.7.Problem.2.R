## Part A -- 
## Let L(t) be the probability the 40 year old lives past the age 40 + t where t is any positive real number. 
## Estimate L(t) by first considering t = 0, 1, 2, .... These values of L(t) can be computed using the spreadsheet data.
## (For example, for 40-41 and 41-42). For other t values, interpolate using a cubic spline. In R you can use the spline
## and splinefun commands to construct cubic splines, see the help documentation. Graph teh interpolating cubic spline


## Import data and organize data
life <- read.csv("C:/Users/562644/Documents/~ Education/Georgetown/2015 Spring/504 - Numerical Methods/HW7/US Life Expectancy 2003.csv", 
                           header=TRUE, quote="\"")
life <- life[1:101, 1:2]

## L(0) = 1
## L(1) = 1 - P(40 <= X <= 41) = .9979962
## L(2) = L(1) - P(41 <= X <= 42) * P(40 <= X <= 41) ...

L.integer <- function(t){
  p.death <- c()
  for (i in 1:length(t)){
    if(t[i]==0){
      p.death[i]=0 ## as opposed to the p.live, which is 1
    }
    else{
      p.death[i]=life$qx[41]
      if(t[i]==1){
        p.death[i]=p.death[i]
      }
      if (t[i] > 60){
        p.death[i]=1
      }
      else{
        for (j in 1:(t[i]-1)){
          p.death[i]=p.death[i]+prod(1-life$qx[41:(40+j)])*life$qx[(41+j)]
        }
      }
    }
  }
  p.live <- 1 - p.death
  return(p.live)
}

L <- function(t){
  l=c()
  for (i in 1:length(t)){
    ## If it's an integer
    if (t[i]-ceiling(t[i])==0){
      l[i] = L.integer(t[i])
    }
    else{
      x <- c(floor(t[i]), ceiling(t[i]))
      y <- c(L.integer(x[1]), L.integer(x[2]))
      f <- splinefun(x,y)
      l[i] <- f(t[i])
    }
  }
  return(l)
}


L(1)
L(2)

t <- seq(0, 100, 0.2)
plot(t, L(t))


life.matrix <- cbind(c(0:100), life$qx)
x <- c(41:101)
y <- cumprod(1 - life.matrix[x,2])

i <- c(1:720)
L <- splinefun(x,y)(40+i/12)

sum(200*L*exp(-0.05*i/12))
