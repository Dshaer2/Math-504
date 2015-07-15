## Objective #1: Segment the user data such that users are divided into 2 classes – 
## those that like watching Oprah and those that don’t like watching Oprah (I call them “NOprahs”). 
## Then run PCA on the data, plot the data on the first two principal components then create an LDA classifier. 
## Note the effectiveness of the classifier on the training data, then try it on a test observation.

## Import data
users <- read.table("C:/Users/562644/Documents/~ Education/Georgetown/2015 Spring/504 - Numerical Methods/HW12/user-shows.txt", header=FALSE, quote="\"")
shows <- read.table("C:/Users/562644/Documents/~ Education/Georgetown/2015 Spring/504 - Numerical Methods/HW12/shows.txt", header=FALSE, quote="\"")

## Merge into one data frame
shows <- as.vector(shows[,1])
colnames(users) <- shows

## Create new class
users$class <- "red" ## Red for Oprah
users$class[users[,266]==0] <- "blue" ## Blue for NOprah

sum(users$class=="red")/9985 ## 30% likes Oprah

## We need a square matrix to run PCA, so we sample 562 users from the full 9985
s <- sample(1:9985, 562, replace = FALSE)
user.matrix <- cbind(users[s,1:265], users[s,267:563])
user.class <- users[s,564]

## Center Data and build covariance matrix
user.matrix.cen <- sapply(1:562, function(x) user.matrix[,x] - mean(user.matrix[,x]))
X <- cov(user.matrix.cen)

## PCA using the power method:
## normalizing function
norm <- function(x){
  norm <- sqrt(sum(x^2))
  return(norm)
}
## power method
power_method <- function(M, v1_old, v2_old, eps) {
  v1_new <- M %*% v1_old
  v1_new <- v1_new / norm(v1_new)
  
  v2_new <- M %*% v2_old
  v2_new <- v2_new - as.numeric( (t(v1_new) %*% v2_new) / (t(v1_new) %*% v1_new) ) * v1_new
  v2_new <- v2_new/norm(v2_new)
  steps <- 1
  
  while(norm(v1_new - v1_old) > eps) {
    v1_old <- v1_new
    v2_old <- v2_new
    
    v1_new <- M %*% v1_old
    v1_new <- v1_new / norm(v1_new)
    
    v2_new <- M %*% v2_old
    v2_new <- v2_new <- v2_new - as.numeric( (t(v1_new) %*% v2_new) / (t(v1_new) %*% v1_new) ) * v1_new
    v2_new <- v2_new/norm(v2_new)
    
    steps <- steps + 1
  }
  
  lambda1 = t(v1_new) %*% M %*% v1_new
  lambda2 = t(v2_new) %*% M %*% v2_new
  
  list(v1_new = v1_new, lambda1 = lambda1, topnode1 = which.max(v1_new),
       v2_new = v2_new, lambda2 = lambda2, topnode2 = which.max(v2_new),
       iterations = steps)
}

## Now we run the power method to find the eigens. We will initialize with 2 vectors of length 562
eps <- 1e-10
v1 <- rep(1,562)
v2 <- rep(2,562)

results <- power_method(X, v1, v2, eps)

q1 <- results$v1_new
q2 <- results$v2_new
e1 <- results$lambda1
e2 <- results$lambda2

## Eigens using R
eVals <- eigen(X)$values
eVecs <- eigen(X)$vectors

## As you can see, R's Eigens confirm my results for first component
head(q1)
head(eVecs[,1])

e1
eVals[1]

## We project the centered data onto the first and second principal components
C <- t(sapply(1:562, function(x) t(user.matrix.cen[x,]) %*% eVecs[,1:2]))

plot(C[,1], C[,2], pch=18, ylab = "Second Principal Component", xlab = "First Principle Component", col = user.class)
plot(C[,1], pch=18, xlab = "First Principle Component", col = user.class)
plot(C[,2], pch=18, xlab = "Second Principle Component", col = user.class)

## Isolate the principal components for each class (Oprahs and NOprahs)
Oprah.1pc  <- C[user.class == "red", 1]
NOprah.1pc <- C[user.class == "blue", 1]

Oprah.2pc  <- C[user.class == "red", 2]
NOprah.2pc <- C[user.class == "blue", 2]

## Let's make sure the Variances are close because inconsistent variances among classes will screw up LDA. 
var(Oprah.1pc)
var(NOprah.1pc)

var(Oprah.2pc)
var(NOprah.2pc)

## These are close enough. Let's define the parameters, mean, sigma, pi_i
mean.Oprah <- c(mean(Oprah.1pc), mean(Oprah.2pc))
mean.NOprah <- c(mean(NOprah.1pc), mean(NOprah.2pc))
SIGMA <- cov(C)
pi_Oprah <- 3/10 ## recall 30% of our sample liked Oprah

## Write LDA function
LDA.func <- function(x, SIGMA, pi_Oprah){
  score.Oprah  <- pi_Oprah*exp(-t(x - mean.Oprah) %*% solve(SIGMA) %*% (x - mean.Oprah) / 2)
  score.NOprah <- (1-pi_Oprah)*exp(-t(x - mean.NOprah) %*% solve(SIGMA) %*% (x - mean.NOprah) / 2)
  matrix(c(score.Oprah, score.NOprah), nrow=1)
}

## Run LDA
result <- t(sapply(1:562, function(x) LDA.func(C[x,], SIGMA, pi_Oprah)))

## The resulting matrix contains the "scores" from the discriminant functions.
## Users in the "Oprah" class will have a larger column 1 score and "NOprahs" will have a larger column 2 score
head(result)

## Test on Training Data
classification.result <- ifelse(result[,1] > result[,2], "red", "blue")

## I built a function to write a confusion matrix so that we can better understand the results of the classification
confusion.matrix <- function(classification.result, user.class){
  truepositive <- sum((classification.result == "red") & (user.class == "red"))
  falsepositive <- sum((classification.result == "red") & (user.class == "blue"))
  truenegative <- sum((classification.result == "blue") & (user.class == "blue"))
  falsenegative <- sum((classification.result == "blue") & (user.class == "red"))
  
  confusion <- rbind(c(truepositive, falsenegative),c(falsepositive,truenegative))
  overall <- (truepositive+truenegative) / (truepositive+falsepositive+truenegative+falsenegative)
  specificity <- truenegative/(truenegative + falsepositive)
  sensitivity <- truepositive/(truepositive+falsenegative)
  list(confusion, overall, specificity, sensitivity)
}

Oprah.confusion <- confusion.matrix(classification.result, user.class)
Oprah.confusion


## Plot the line that divides the classes
## Step 1: Identify centroids of each cluster
mean.Oprah
mean.NOprah

## Step 2: compute line that connects centroid (simple algebra derived from y = mx + b)
s1 <- (mean.Oprah[2]-mean.NOprah[2])/(-mean.NOprah[1]+mean.Oprah[1])
b1 <- mean.Oprah[1]*s1 - mean.Oprah[2]

plot(C[,1], C[,2], pch=".", ylab = "Second Principal Component", xlab = "First Principle Component", col = user.class, cex = 3)
points(x=mean.Oprah[1], y=mean.Oprah[2], pch=3, cex=3)
points(x=mean.NOprah[1], y=mean.NOprah[2], pch=3, cex=3)
lines (x=c(mean.Oprah[1], mean.NOprah[1]), y=c(mean.Oprah[2],mean.NOprah[2]), lty=3, col = "black")

## Step 3: compute orthogonal line that crosses the midpoint
## negative reciprocal of slope
s2 <- -1/s1

## midpoint
mid <- c((mean.Oprah[1]+mean.NOprah[1])/2, (mean.Oprah[2]+mean.NOprah[2])/2)
b2 <- mid[2] - mid[1]*s2
plot(C[,1], C[,2], pch=".", ylab = "Second Principal Component", xlab = "First Principle Component", col = user.class, cex = 3)
abline(a=b2, b=s2)

## Oversimplified linear boundary, I know...
## but hey, it classified over 76% of the training set correctly



## Objective #2: Run SVD on the same data
## Plot the television space and identify programs with projections similar to The Oprah Winfrey Show
## Plot the user space and make recommendations for a new user.

## Even though SVD can handle rectangular matrices, we are going to stick with the same square user matrix
M <- user.matrix

## Run SVD
svd.M <- svd(M)
S <- diag(svd.M$d)
U <- svd.M$u
V <- svd.M$v

## Plot the TV Space by projecting onto U
C. <- t(sapply(1:562, function(x) matrix(as.numeric(M[,x]),nrow=1) %*% U[,1:2]))
plot(C.[,1], C.[,2], ylab = "b_2", xlab = "b_1", main="TV Space", pch=18)

## Let's plot Oprah on the TV Space
Oprah <- users[s,266] ## remember we omitted "The Oprah Winfrey Show" (266) from the user.matrix

O <- t(matrix(as.numeric(Oprah),nrow=1) %*% U[,1:2])
points(x=O[1], y=O[2], col = "dark orange", pch = 17, cex = 1.5)
## She is on the boundary of the data

## Let's examine all of Oprah's nearest neighbor projections
b_1.peers <- which( (C.[,1] > O[1]-1) & (C.[,1] < O[1]+1) ) ## These programs are all within 1 b_1 unit of Oprah
Oprah.peers <- b_1.peers[(C.[b_1.peers,2] > O[2]-1) & (C.[b_1.peers,2] < O[2]+1)] ## Of those programs, these are all within 1 b_2 unit of Oprah
Oprah.peers.df <- user.matrix[, Oprah.peers]
colnames(Oprah.peers.df)

points(x=C.[Oprah.peers,1], y=C.[Oprah.peers,2], pch = 17, col = "gray") ## Visually confirm that these are the nearest neighbors

## Now let's plot User Space by projecting onto V
C <- t(sapply(1:562, function(x) matrix(as.numeric(M[x,]),nrow=1) %*% V[,1:2]))
plot(C[,1], C[,2], ylab = "a_2", xlab = "a_1", main="User Space", pch=18, col=user.class)

## You can see the Oprah watchers towards the left boundary

## I had my girlfriend, Kristianne select programs that she enjoyed watching from the list of 562 (omitting Oprah)
kristianne <- c(10, 11, 13, 28, 29, 85, 87, 98, 121, 139, 143, 150, 151, 175, 228, 234, 273, 339, 355, 358, 362, 397, 405, 429, 435, 438, 441, 442, 447, 456, 470, 497)
matrix(colnames(user.matrix), ncol=1)[kristianne] ## check to make sure I recorded these right

ka <- rep(0, 563)
ka[kristianne] <- 1
ka <- c(ka[1:265], ka[267:563])

## Plot Kristianne's projection
ka.userspace <- matrix(as.numeric(ka), nrow=1) %*% V[,1:2]
plot(C[,1], C[,2], ylab = "a_2", xlab = "a_1", main="User Space", pch=18)
points(x=ka.userspace[1], y=ka.userspace[2], col = "dark orange", pch = 17, cex = 1.5)

## Identify Kristianne's nearest neighbors
a_1.peers <- which( (C[,1] > ka.userspace[1]-.5) & (C[,1] < ka.userspace[1]+.5) ) ## These programs are all within 0.5 a_1 units of Kristianne
ka.peers <- a_1.peers[(C[a_1.peers,2] > ka.userspace[2]-.5) & (C[a_1.peers,2] < ka.userspace[2]+.5)] ## Of those programs, these are all within 0.5 1_2 units of Kristianne
ka.peers.df <- user.matrix[ka.peers,]
rownames(ka.peers.df) ## These are users with similar taste to Kristianne

points(x=C[ka.peers,1], y=C[ka.peers,2], pch = 17, col = "gray") ## Visually confirm that these are the nearest neighbors
points(x=ka.userspace[1], y=ka.userspace[2], col = "dark orange", pch = 17, cex = 1.5)

ka.peers.df <- t(ka.peers.df)
ka.recs <- sort(rowSums(ka.peers.df), decreasing = TRUE)

matrix(names(ka.recs[1:15]), ncol=1)
matrix(ka.recs[1:15],ncol=1)

## Of the Top 15 recommendations, Kristianne says that she is only interested in watching 8. This is barely better than 50%
## I do not know how the sports and news stuff got mixed in there?

## Finally, let's try to predict if Kristianne is an Oprah or a NOprah
## we just need to project her tastes onto the principal components and plot
ka.pca <- t(ka) %*% eVecs[,1:2]
C <- t(sapply(1:562, function(x) t(user.matrix.cen[x,]) %*% eVecs[,1:2]))

plot(C[,1], C[,2], pch=".", ylab = "Second Principal Component", xlab = "First Principle Component", col = user.class, cex = 3)
abline(a=b2, b=s2)
points(x=ka.pca[1], y=ka.pca[2], pch = 17, cex = 1.5, col = "dark orange")

## Thus, Kristianne is classified as an Oprah fan (red) according to my LDA classifier.
## This is an accurate assessment. Kristianne likes watching Oprah








kristianne <- c(2,6,13,16,21,26,27,28,32,44,65,73,98,111,112,119,140,143,149,150,
  151,157,167,169,175,256,271,275,278,293,298,308,339,340,388,417,478,500,524)                          

