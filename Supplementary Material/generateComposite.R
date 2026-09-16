
#set parameters: 
p = n = 5 #number of variables (aka p)
cor <- .6 #desired correlation between the two variables in the sum score

#replace next 2 lines with however you generate your original p*p matrix
#including convert to pcor, add sparsity, convert back, etc. 
S <- matrix(.5, p, p) 
diag(S) <- 1


#to turn the latent into a composite: 

#step 1: add another variable that is orthogonal to everything else: 
Splus <- diag(p+1)
Splus[1:p,1:p] <- S #latent.sparse

#step 2: generate item-level data using lambda and theta matrices and a funky measurement model
lambda <- diag(x = .9, (n+1), (n+1))
lambda[n, n] <- lambda[n,(n+1)] <- sqrt((cor + 1)/2)
lambda[(n+1), n] <- sqrt(((cor + 1)/2) - cor)
lambda[(n+1),(n+1)] <- -(sqrt(((cor + 1)/2) - cor))

theta <- matrix(0, n+1, n+1)
diag(theta)[1:(n-1)] <- .19

item.mat <- t(lambda) %*% Splus %*% lambda + theta

## TEST TO SEE THAT IT WORKS (don't put this part in your code) ##
# generate data from this item.mat and create sum score out of items 5 and 6 
# to make sure it works: 
d <- data.frame(MASS::mvrnorm(Sigma = item.mat, n = 1000, mu = rep(0, n+1), empirical = TRUE))
d$sum <- d$X5 + d$X6

cov(d)
cor(d) 
########