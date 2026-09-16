cov2cor = function(mat){
  inv.sd = diag(1/sqrt(diag(mat)))
  matrix = inv.sd %*% mat %*% inv.sd
  return(matrix)
}
cor2pcor <- function(mat){
  precision.mat <- solve(mat)
  pcor.mat <- qgraph::wi2net(precision.mat)
  pcor.mat <- as.matrix(pcor.mat)
  return(pcor.mat)
}

cor.gen = function(nvar = NULL,
                   bg.cor = 0.30,
                   mn.cor = 0.50,
                   prop.cor = 1,
                   iter.lim = 100,
                   sd = 0.1,
                   tol.eig = 1e-1,
                   mean.tol = 0.005,
                   clip = 0.99){
  
  
  
  flag = 0
  iter = 0
  
  
  n_lower = (nvar^2 - nvar) / 2 # Unique off-diagonal elements 
  n.big = floor(n_lower * prop.cor) # floor will never be greater than itself, 10.5 = 10
  # take # of unique off-diag elements and multiply proportion of correlations for n.big so this is how many will be larger correlations 
  while(flag < 1){
    iter = iter + 1
    
    temp = matrix(rnorm(nvar^2, mean = bg.cor, sd = sd), nvar, nvar) # Initial correlation matrix 
    diag(temp) = 1 # Set diagonals to 1 to be correlations 
    
    idx.lower = which(lower.tri(temp)) # R column-major indexing 
    
    if (n.big > 0) {
      idx.big = sample(idx.lower, n.big)
      temp[idx.big] = rnorm(n.big, mean = mn.cor, sd = sd) # Replace whatever proportion of elements should have mn.cor  size correlations
    } else {
      idx.big = integer(0)
    }
    
    temp[lower.tri(temp)] = pmax(pmin(temp[lower.tri(temp)], clip), -clip)
    # any values less than or greater than clip will be replaced by clip, making the correlations bounded between -1 and 1
    temp[upper.tri(temp)] = t(temp)[upper.tri(temp)]
    # Replace lower elements with upper
    diag(temp) = 1
    
    # Eigen check
    E = eigen(temp, symmetric = TRUE) # Provides eigenvalues and vectors, eigen decomposition 
    # t(E$vectors) %*% E$vectors will give the identity matrix, all eigenvalues must be >= 0
    # E$vectors %*% diag(E$values, nrow = 5, ncol = 5) %*% solve(E$vectors) = temp (although symmetric)
    vals = E$values # Eignevalues
    vals[vals < tol.eig] = tol.eig # Replace any with a really small positive value that are less than 0
    temp.psd = E$vectors %*% (diag(vals, nvar, nvar)) %*% t(E$vectors) # eigenvectors are orthonormal so that its inverse is equal to its transpose
    # Rebuild the matrix that is psd  
    d = sqrt(diag(temp.psd)) # Ensure it is correlation matrix
    temp_cor = temp.psd / outer(d, d)
    diag(temp_cor) = 1
    
    if (length(idx.big) > 0) {
      big.mean = mean(temp_cor[idx.big]) # Average correlation size of larger correlations 
    } else {
      big.mean = NA_real_
    }
    setdiff(idx.lower, idx.big) # Returns elements present in first object but not in 2nd
    idx.bg = setdiff(idx.lower, idx.big)
    if (length(idx.bg) > 0) {
      bg.mean = mean(temp_cor[idx.bg])
    } else {
      bg.mean = NA_real_
    }
    # Write what iteration and the average correlations of either baseline or larger
    cat("\r", paste0("Iter ", iter,
                     " | big mean=", ifelse(is.na(big.mean), "NA", round(big.mean, 3)),
                     " | bg mean=", ifelse(is.na(bg.mean), "NA", round(bg.mean, 3))), sep = "")
    flush.console()
    # Tolerance, if less than whatever mean.tol is then move on(how off we are willing to be from set correlation)
    ok_big = if (length(idx.big) > 0) abs(big.mean - mn.cor) <= mean.tol else TRUE
    ok_bg  = if (length(idx.bg) > 0)  abs(bg.mean - bg.cor) <= mean.tol else TRUE
    ok = ok_big && ok_bg
    # If TRUE send loop and return PSD matrix 
    if (ok) {
      flag = 1
      cat("\n")
      return(temp_cor)
    }
    
    if (iter >= iter.lim) {
      flag = 1
      cat("\nReached Iteration Limit. Returning last PSD correlation matrix.\n")
      return(temp_cor)
    }
  }
  
  cat("\n")
  return(NULL)
}


sparse.prec.mat = function(mat, prop = .3, tol = 1e-1, max.iter = 10000){
  
  
  p = nrow(mat)
  fail.count = 0
  prec.mat = solve(mat)
  idx.lower = which(lower.tri(prec.mat)) # Get indices of lower triangle
  p.zero = floor(prop * length(idx.lower)) # Number of 0s 
  
   for (i in 1:max.iter) {
  idx.zero = sample(idx.lower, p.zero) # Sample from indices that should be 0
  prec.mat[idx.zero] = 0 # Change the sampled indices to 0 
    
  prec.mat[upper.tri(prec.mat)] = t(prec.mat)[upper.tri(prec.mat)] # Make upper triangle the same
    
  # Check for PSDness
  E = eigen(x = prec.mat, symmetric = TRUE)
  vals = E$values
  vals[vals < 0] = tol
  psd = E$vectors %*% (diag(vals, p, p)) %*% t(E$vectors)
  
  psd[idx.zero] = 0 # Change the sampled indices to 0 
  psd[upper.tri(psd)] = t(psd)[upper.tri(psd)] # Make upper triangle the same
  
  #psd[abs(psd) < tol] = 0 # Change noisy 0s to straight up 0s
  
  #zeros = sum(psd[lower.tri(psd)] == 0)
  #proportion = round(zeros/(p*(p-1)/2), 1) # Tried making sure proportion of 0s match what sparsity I want .3, .5, etc 
  
  if(min(eigen(psd)$values) > tol){
    message("Found correct PSD matrix after ", fail.count, " failed attempt(s)")
    return(psd)
   }
  
  fail.count = fail.count + 1
  
   }
  message("No valid matrix found after ", max.iter, " iterations")
  return(NULL)
   
}

# Produces really strange matrices (inflated values)
# Fails a lot, I tried different sparsity levels (.3, .5., .7) at 10 and 20 network sizes and saw consistently more failures
# and strange matrices
lat.mat = cor.gen(nvar = 10) # default avg corr is .5
mat = sparse.prec.mat(mat = lat.mat, prop = .3)
mat


# # Ok what if I make a copy straight from latent.sparse?
# S = sparse.latent.mat
# copy = sparse.latent.mat[1:(p-1), p]
# Splus <- diag(p+1)
# Splus[1:p,1:p] <- S 
# Splus[1:(p-1), p+1] = copy
# Splus[p+1, 1:(p-1)] = copy
# Splus[p, p+1] = .81
# Splus[p + 1, p] = .81
# 
# # Get item correlations now                         
# lambda = diag(x = loadings, nrow = nrow(Splus) , ncol = ncol(Splus))
# item.mat = lambda %*% Splus %*% t(lambda)
# theta = diag(x = 1 - (loadings)^2,
#              nrow = nrow(Splus),
#              ncol = ncol(Splus))
# 
# sigma.nothing = item.mat + theta # Sigma (p+1) X (p+1)
# 
# vals = eigen(sigma.nothing, symmetric = TRUE)$values
# if(!all(vals >= 0)) {
#   return(NULL)
# }
# return(sigma.nothing)
# 
# 
# cor2pcor(sigma.nothing)
# 
# p = nvar + 1
# mat = cor.gen(nvar = p) # Start with (p+1) X (p+1) latent correlation matrix 
# prec.mat = sparse.prec.mat(mat = mat, 
#                            prop = sparsity) # Eventually change for the correct function to get sparsity 
# latent.sparse = cov2cor(solve(prec.mat))
# 
# 
# # Get item correlations now                         
# lambda = diag(x = loadings, nrow = nrow(latent.sparse) , ncol = ncol(latent.sparse))
# item.mat = lambda %*% latent.sparse %*% t(lambda)
# theta = diag(x = 1 - (loadings)^2,
#              nrow = nrow(latent.sparse),
#              ncol = ncol(latent.sparse))
# 
# sigma.nothing = item.mat + theta # Sigma (p+1) X (p+1)
# # Copy 2nd to last col and make it the same as last col
# copy = sigma.nothing[1:(nvar - 1), nvar]
# sigma.nothing[1:(nvar - 1), p] = copy
# sigma.nothing[p, 1:(nvar - 1)] = copy
# 
# # Match correlation between target and pseudo-clone so that it looks like redundancy
# fake.redun = desired.cor * desired.cor
# sigma.nothing[p, nvar] = fake.redun
# sigma.nothing[nvar, p] = fake.redun
# 
# # rownames(sigma.nothing) = colnames(sigma.nothing) = c(paste0("p", 1:(nvar-1)), "Target", "Pseudo-Clone")
# vals = eigen(sigma.nothing, symmetric = TRUE)$values
# if(!all(vals >= 0)) {
#   return(NULL)
# }
# return(sigma.nothing)


# Add desired correlation between two items
mat = cor.gen(nvar = 5)
sparsity = .3
p = 5
peripheral.loadings = .9
clone.loading = .9
mat = latent.gen(mat)

mat = mat$latent.sparse
loadings = .9
# 
# lambda[p+1, p+1] = .7 # Just need to manipulate the redundancy here

donothing.gen = function(mat, loadings = .9, clone.loading = .9){
  S = mat
  copy = mat[1:(p-1), p]
  Splus <- diag(p+1)
  Splus[1:p,1:p] <- S 
  Splus[1:(p-1), p+1] = copy
  Splus[p+1, 1:(p-1)] = copy
  Splus[p, p+1] = .999 # So the matrix is still invertible
  Splus[p + 1, p] = .999

  # Get item correlations now                         
  lambda = diag(x = loadings, nrow = nrow(Splus) , ncol = ncol(Splus))
  lambda[nrow(lambda), ncol(lambda)] = clone.loading # Change level of redundancy between nodes
  item.mat = lambda %*% Splus %*% t(lambda)
  theta = diag(x = 1 - (loadings)^2,
               nrow = nrow(Splus),
               ncol = ncol(Splus))
  theta[nrow(theta), ncol(theta)] = 1 - (clone.loading)^2
  
  sigma.nothing = item.mat + theta # Sigma (p+1) X (p+1)
  
  vals = eigen(sigma.nothing, symmetric = TRUE)$values
  if(!all(vals >= 0)) {
    return(NULL)
  }
  return(sigma.nothing)
}







# Generating Issue --------------------------------------------------------

#what we expect to happen
A <- matrix(.5, 20, 20)
diag(A) <- 1
Ainv <- solve(A)
Anet <- -1 * cov2cor(Ainv)
diag(Anet) = 0



#what actually happens
X <- cor.gen(nvar = 20, mn.cor =.5)
Xnet = cor2pcor(X)
# Xinv <- solve(X)
# Xnet <- -1 * cov2cor(Xinv)
# diag(Xnet) = 0

par(mfrow = c(1, 2))
hist(X[lower.tri(X)], xlim = c(-1, 1), breaks = seq(-1.1, 1.1, by = 0.1))
# ?hist
hist(Xnet[lower.tri(Xnet)], xlim = c(-1, 1), breaks = seq(-1.1, 1.1, by = .1))
mean(A[lower.tri(A)])
sd(A[lower.tri(A)])
qgraph::qgraph(Xnet, theme = 'colorblind', layout = "circle",
               edge.labels = TRUE, minimum = 0.5)
# Original study
# Needed values
nv = 20
mn.cor = .5
loading = .9
clone.loading = .9
ind.corr <- function(matrix, loadings, clone.loading = .9, redundancy = TRUE){
  # Lambda
  dimensions <- dim(matrix)
  if (redundancy == FALSE) {
    nv = nv + 1
    matrix = cor.gen(nvar = nv, mn.cor = mn.cor)
    dimensions = dim(matrix)
    lambda_matrix = diag(x = loadings, nrow = dimensions[1], ncol = dimensions[2])
    last_element = nrow(lambda_matrix) * ncol(lambda_matrix)
    lambda_matrix[last_element] = clone.loading
    
    # Lambda * Psi * t(Lambda)
    random_mat <- lambda_matrix %*% matrix %*% t(lambda_matrix)
    
    # Theta
    dimensions_theta <- dim(random_mat)
    theta <- diag((x = 1 - (loadings)^2), nrow = dimensions_theta[1], ncol = dimensions_theta[2])
    theta[dimensions_theta[1],dimensions_theta[2]] <- (1-(clone.loading)^2)
    random_mat <- random_mat + theta
    return(random_mat)
  }
  lambda_matrix <- diag(x = loadings, nrow = dimensions[1] + 1 , ncol = dimensions[2])
  last_element <- nrow(lambda_matrix) * ncol(lambda_matrix)
  lambda_matrix[last_element] <- clone.loading
  
  # Lambda * Psi * t(Lambda)
  redun_mat <- lambda_matrix %*% matrix %*% t(lambda_matrix)
  
  # Theta
  dimensions_theta <- dim(redun_mat)
  theta <- diag((x = 1 - (loadings)^2), nrow = dimensions_theta[1], ncol = dimensions_theta[2])
  theta[dimensions_theta[1],dimensions_theta[2]] <- (1-(clone.loading)^2)
  redun_mat <- redun_mat + theta
  cat("\n Mean =", mean(redun_mat[lower.tri(redun_mat)]), "\n")
  return(redun_mat)
}

lat.corr = cor.gen(nvar = nv, mn.cor = mn.cor)
redun.corr = ind.corr(matrix = lat.corr,
                      loadings = loading,
                      clone.loading = clone.loading,
                      redundancy = TRUE)

true.corr = redun.corr[1:nv, 1:nv]

true.pcor = cor2pcor(true.corr)
redun.pcor = cor2pcor(redun.corr)


par(mfrow = c(1, 2))
hist(true.pcor[lower.tri(true.pcor)], xlim = c(-1, 1))
hist(redun.pcor[lower.tri(redun.pcor)], xlim = c(-1, 1), breaks = seq(-1.1, 1.1, by = .1))









# NaN vals and the mystery of 0s  -----------------------------------------

# Need result matrices and specific condition/iteration



# NaN vals, I found out that sometimes at a sparsity of .5 all of the edges to
# the target are wiped out causing the target to have no closeness because the 
# nature of the equation

# Get matrices of iteration where I saw the vals of NaN 
# Condition is   p = 10, lat.cor = .5, sparsity = .5, lvl.redun = .9

true.net = matrices[[1]][[10]]$omega.pseudo
true.net = matrices[[1]][[10]]$omega.pseudo # Collapse so we are comparing to the average of the last two

# Get the sigma used in that iteration and apply correct solution
est.net = apply.composite(matrices[[1]][[10]]$pseudo.sigma, p = 10)



# qgraph(true.net, edge.labels = TRUE)
# qgraph(est.net, edge.labels = TRUE)
true.net2 = pcor.to.igraph(true.net)
est.net2 = pcor.to.igraph(est.net)

# Closeness becomes NaN cuz it has no connections
closeness(graph = true.net2, normalized = FALSE,
          weights = 1/abs(E(true.net2)$weight))




# Dealing with inflated vals for apl.bias 
true.net = matrices[[1]][[10]]$omega.latent

est.net = apply.removal(matrices[[1]][[10]]$redundant.sigma, p = 10)

true.net2 = pcor.to.igraph(true.net)
est.net2 = pcor.to.igraph(est.net)



E(est.net2)$weight[abs(E(est.net2)$weight) < 1e-4] = 0 # worked for other ocasions
E(true.net2)$weight[E(true.net2)$weight < 1e-4] = 0


# Test fix using APL function
mean_distance(graph = est.net2,
              weights = 1/abs(E(est.net2)$weight))


mean_distance(graph = true.net2,
              weights = 1/abs(E(true.net2)$weight))

calc.apl.bias(est.net = est.net2, true.net = true.net2)




# Pseudo problem

true.net = matrices[[1]][[2]]$omega.pseudo
true.net = collapse.pseudo(true.net, p =10)
est.net = apply.lnm(matrices[[1]][[2]]$pseudo.sigma, p = 10)

true.net2 = pcor.to.igraph(true.net)
E(true.net2)$weight[abs(E(true.net2)$weight) < 1e-4] = 0
est.net2 = pcor.to.igraph(est.net)
E(est.net2)$weight[abs(E(est.net2)$weight) < 1e-4] = 0 


mean_distance(graph = est.net2,
              weights = 1/abs(E(est.net2)$weight))

mean_distance(graph = true.net2,
              weights = 1/abs(E(true.net2)$weight))


# Weird Bias at 0 Sparsity ------------------------------------------------

# Example of the problem Condition 9 iteration 42
# avg lat corr = 0.9, sparsity = 0, lvl.redun 0.7
# Goal network: Composite
# Solution: lnm




prob.sigma = mats[[69]][[70]]$sparse.latent.mat
# prob.omega = apply.lnm(mats[[9]][[1]]$redundant.sigma,p = 10)
prob.omega = mats[[69]][[70]]$omega.pseudo

hist(prob.omega)




# mean(prob.omega[(lower.tri(prob.omega))])
# mean(prob.sigma[(lower.tri(prob.sigma))])
# hist(prob.omega[(lower.tri(prob.omega))])
# hist(prob.sigma[(lower.tri(prob.sigma))])
# avg partial correlation = 0.11

# qgraph::qgraph(prob.omega, edge.labels = TRUE)

# library(psychonetrics)
omega = apply.removal(sigma = prob.sigma, p = 10)

prob.omega = collapse.pseudo(prob.omega, p = 10)


true.ig = delete_edges(true, E(true)[abs(E(true)$weight) < 1e-2])
est.ig = delete_edges(est, E(est)[abs(E(est)$weight) < 1e-2])




# apply.composite(prob.sigma, p = 10)

# mean(abs(omega_lnm[(lower.tri(omega_lnm))]))
# hist(omega_lnm[(lower.tri(omega_lnm))])

# library(igraph)
est = pcor.to.igraph(omega)
true = pcor.to.igraph(prob.omega)



strength(est)
strength(true)
calc.strgth.bias(est.net = est, true.net = true)
calc.RMSE(est.net = est, true.net = true)
calc.apl.bias(est.net =  est, true.net = true)


sigma.apl = mean_distance(graph = est.ig,
                         weights = 1/abs(E(est.ig)$weight))
true.apl = mean_distance(graph = true.ig,
                         weights = 1/abs(E(true.ig)$weight))








apply.lnm = function(sigma, p, nobs = 10000) {
  lambda = diag(x = 1, nrow = p + 1, ncol = p)
  lambda[p + 1, p] = 0.9
  mod = lnm(cors = prob.sigma, lambda = lambda, nobs = 10000,
            omega_zeta = "full", identification = "variance") |>
    runmodel()
  getmatrix(mod, "omega_zeta")
}








# Weird LNM? Goal Network: Composite --------------------------------------

# Find condition, looking for something lvl.redun = 0.7, p = 20, solution = lnm,
# sparsity = 0

# Need to see that exact condition
# Isolate where the issue is arising and find any bias greater or less than +/- 1
perf.df %>% 
  filter(p == 20 & sparsity == 0 & lvl.redun == 0.7 & sigma.type == "composite" &
           truth == "composite" & solution == "lnm") %>% 
  filter(abs(strgth.bias) >= 1) # 37 cases when sparsity = 0, 9 cases when sparsity = 0.3, 2 cases when sparsity = 0.5, and 0 when sparsity = 0.7
  # Of those 37 cases 5 are greater than 10

# Let's check the matrices we are getting when they are greater than 1 

sigma.6.4 = mats[[6]][[4]]$composite.sigma
omega.6.4 = mats[[6]][[4]]$omega.composite

# Avg partial correlation in omega.6.4 
mean(omega.6.4[(lower.tri(omega.6.4))]) # 0.05067
hist(omega.6.4[(lower.tri(omega.6.4))]) # Nothing out of the ordinary

# Condition 6 iteraion 4 apply solution lnm
solution.6.4 = apply.lnm(sigma = sigma.6.4, p =  20)

# Avg partial correlation in solution.6.4 
mean(solution.6.4[(lower.tri(solution.6.4))]) # 0.0479
hist(solution.6.4[(lower.tri(solution.6.4))])
# Warning from model indicated this run may have a Heywood case, optimization problem

# Compute bias

# convert to igraph objects
est.net = pcor.to.igraph(solution.6.4)
true.net = pcor.to.igraph(omega.6.4)

true.strength.all = strength(graph = true.net, weights = abs(E(true.net)$weight))
true.strength.all[p]
strength(graph = est.net, weights = abs(E(est.net)$weight))[p]

# Let's check the matrices we are getting when they are greater than 10 

sigma.6.39 = mats[[49]][[1]]$redundant.sigma
omega.6.39 = mats[[49]][[1]]$omega.latent
library(qgraph)
qgraph(omega.6.39, edge.labels = TRUE)

# Avg partial correlation in omega.6.13 
mean(omega.6.39[(lower.tri(omega.6.39))]) # 0.05059
hist(omega.6.39[(lower.tri(omega.6.39))]) # Nothing out of the ordinary

# Condition 6 iteraion 13 apply solution lnm
solution.6.39 = apply.composite(sigma = sigma.6.39, p =  10)

solution.6.39 - omega.6.39
# 2: In runmodel(lnm(cors = sigma, lambda = lambda, nobs = nobs, omega_zeta = "full",  :
# One or more parameters were estimated to be near its bounds. This may be indicative of, for example, a Heywood case, but also of an optimization problem. Interpret results and fit with great care. For unconstrained estimation, set bounded = FALSE.
# 3: In addSEs_cpp(x, verbose = verbose, approximate_SEs = approximate_SEs) :
#  Exact standard errors could not be obtained because the Fischer information matrix could not be inverted. Falling back to approximate standard errors. This can occur with zero cells in crosstables (common in multi-group Ising models) or near-boundary estimates. Interpret with care.
# 4: In runmodel(lnm(cors = sigma, lambda = lambda, nobs = nobs, omega_zeta = "full",  :
# Model might not have converged properly: mean(abs(gradient)) > 1.


# Avg partial correlation in solution.6.4 
mean(solution.6.39[(lower.tri(solution.6.39))]) # 0.2612
hist(solution.6.39[(lower.tri(solution.6.39))])
# Warning from model indicated this run may have a Heywood case, optimization problem

# Compute bias

# convert to igraph objects
est.net = pcor.to.igraph(solution.6.39)
true.net = pcor.to.igraph(omega.6.39)

true.strength.all = strength(graph = true.net, weights = abs(E(true.net)$weight))
true.strength.all[p] # 1.7753
strength(graph = est.net, weights = abs(E(est.net)$weight))[p] # 11.8058











perf.df %>% 
  filter(p == 20 & sparsity == 0 & lvl.redun == 0.7 & sigma.type == "composite" &
           truth == "composite" & solution == "lnm") %>% 
  summarize(
    n_total   = n(),
    n_dropped = sum(abs(strgth.bias) > 2, na.rm = TRUE),
    mean_raw  = mean(strgth.bias, na.rm = TRUE),
    mean_trim = mean(strgth.bias[abs(strgth.bias) <= 2], na.rm = TRUE)
  )



true.net = matrices[[1]][[10]]$omega.pseudo
true.net = matrices[[1]][[10]]$omega.pseudo 



# Sparsity Systematic Bias? -----------------------------------------------

# Turn correlations into binaries

# Get matrices
mats = readRDS("/Users/aginigue/Desktop/simulation_matricesP2_100v3.rds")

p          <- 10
target     <- mats_0.7_10
zero.tol   <- 1e-10

idx         <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE) # T/F matrix 45/190
edge.labels <- paste0(idx[, "row"], "-", idx[, "col"])

zero.list <- lapply(names(target), function(cond) {
  flags <- sapply(target[[cond]], function(it) {
    pc <- cor2pcor(it$sparse.latent.mat)
    abs(pc[upper.tri(pc)]) < zero.tol
  })
  data.frame(condition = cond, edge = edge.labels,
             count = rowSums(flags), n.iter = ncol(flags), row.names = NULL)
})

zero.df <- do.call(rbind, zero.list)
zero.df$condition <- factor(zero.df$condition, levels = names(target))







jpeg("sparsity_0.7_10.jpeg",
     width = 12,
     height = 6, 
     res = 800, 
     units = "in")
ggplot(zero.df, aes(x = edge, y = count)) +
  geom_col() +
  facet_wrap(~ condition, ncol = 3) +
  geom_hline(yintercept = 0.7 * 125, linetype = "dashed", colour = "red") +
  labs(x = "Edge Pairs", y = "Count") +
  theme_light() +
  theme(axis.text.x  = element_blank(),
              axis.ticks.x = element_blank())
dev.off()



# Singular Matrices -------------------------------------------------------


apply.lnm = function(sigma, p, nobs = 10000) {
  lambda = diag(x = 1, nrow = p + 1, ncol = p)
  lambda[p + 1, p] = 1
  mod = lnm(cors = sigma, lambda = lambda, nobs = nobs,
            omega_zeta = "full", identification = "variance") |>
    runmodel()
  getmatrix(mod, "omega_zeta")
  parameters(mod)
}


mat = mats[[36]][[1]]$composite.sigma


sparse.latent.mat = mats[[36]][[1]]$sparse.latent.mat


apply.lnm(sigma = mat, 20, nobs = 10000)



composite.gen = function(mat, clone.loading = .85, peripheral.loadings = .9) {
  p = nrow(mat) #number of variables
  cor = 0.9 * clone.loading #desired correlation between the two variables that form a sum score
  
  S = mat # Latent.sparse
  
  #step 1: add another variable that is orthogonal to everything else: 
  Splus <- diag(p+1)
  Splus[1:p,1:p] <- S 
  
  #step 2: generate item-level data using lambda and theta matrices and a funky measurement model
  lambda <- diag(x = peripheral.loadings, (p+1), (p+1))
  lambda[p, p] <- lambda[p,(p+1)] <- sqrt((cor + 1)/2)
  lambda[(p+1), p] <- sqrt(((cor + 1)/2) - cor)
  lambda[(p+1),(p+1)] <- -(sqrt(((cor + 1)/2) - cor))
  
  theta <- matrix(0, p+1, p+1)
  diag(theta)[1:(p-1)] <- .19
  
  sigma.composite <- t(lambda) %*% Splus %*% lambda + theta
  
  colnames(sigma.composite) = rownames(sigma.composite) = c(paste0("p", 1:(p -1)),
                                                            "target",
                                                            "composite_clone")
  
  return(sigma.composite)
  
}

new.sigma = composite.gen(mat = sparse.latent.mat)

apply.lnm(new.sigma, p = 20, nobs = 100)


# Incorrect one 

mat2 = mats[[30]][[90]]$pseudo.sigma
sparse.latent.mat2 = mats[[30]][[90]]$sparse.latent.mat

apply.lnm(sigma = mat2, 20, nobs = 10000)


donothing.gen = function(mat, loadings = .9, clone.loading = .9){
  S = mat
  copy = mat[1:(p-1), p]
  Splus <- diag(p+1)
  Splus[1:p,1:p] <- S 
  Splus[1:(p-1), p+1] = copy
  Splus[p+1, 1:(p-1)] = copy
  Splus[p, p+1] = .999 # So the matrix is still invertible
  Splus[p + 1, p] = .999
  
  # Get item correlations now                         
  lambda = diag(x = loadings, nrow = nrow(Splus) , ncol = ncol(Splus))
  lambda[nrow(lambda), ncol(lambda)] = clone.loading # Change level of redundancy between nodes
  item.mat = lambda %*% Splus %*% t(lambda)
  theta = diag(x = 1 - (loadings)^2,
               nrow = nrow(Splus),
               ncol = ncol(Splus))
  theta[nrow(theta), ncol(theta)] = 1 - (clone.loading)^2
  
  sigma.nothing = item.mat + theta # Sigma (p+1) X (p+1)
  
  vals = eigen(sigma.nothing, symmetric = TRUE)$values
  if(!all(vals >= 0)) {
    return(NULL)
  }
  return(sigma.nothing)
}
p = 20
new.sigma2 = donothing.gen(mat = sparse.latent.mat2)
apply.lnm(sigma = new.sigma2, 20)





# Check 
mat3 = mats[[1]][[1]]$pseudo.sigma
apply.lnm(sigma = mat3, 10, nobs = 10000)






# Sparsity Change ---------------------------------------------------------

p = 5
lat.cor = 0.5
prop = 0.3
lat.mat = cor.gen(nvar = 5, mn.cor = lat.cor)

latent.mats = latent.gen(mat = lat.mat)

cor2pcor(latent.mats$latent.sparse)
latent.mats$prec.mat



redundant.sigma = sigma.gen(mat = latent.mats$latent.sparse, peripheral.loadings = .9, clone.loading = 0.9)

cor2pcor(redundant.sigma)




sparse.prec.mat = function(mat, prop = .3, tol = 1e-1, max.iter = 10000){
  
  
  p = nrow(mat)
  fail.count = 0
  prec.mat.full = solve(mat)
  idx.lower = which(lower.tri(prec.mat.full)) # Get indices of lower triangle
  p.zero = floor(prop * length(idx.lower)) # Number of 0s 
  
  
  
  
  
  
  for (i in 1:max.iter) {
    prec.mat = prec.mat.full
    idx.zero = sample(idx.lower, p.zero) # Sample from indices that should be 0
    prec.mat[idx.zero] = 0 # Change the sampled indices to 0 
    
    prec.mat[upper.tri(prec.mat)] = t(prec.mat)[upper.tri(prec.mat)] # Make upper triangle the same
    
    # Check for PSDness
    E = eigen(x = prec.mat, symmetric = TRUE)
    vals = E$values
    vals[vals < 0] = tol
    psd = E$vectors %*% (diag(vals, p, p)) %*% t(E$vectors)
    
    psd[idx.zero] = 0 # Change the sampled indices to 0 
    psd[upper.tri(psd)] = t(psd)[upper.tri(psd)] # Make upper triangle the same
    
    #psd[abs(psd) < tol] = 0 # Change noisy 0s to straight up 0s
    
    #zeros = sum(psd[lower.tri(psd)] == 0)
    #proportion = round(zeros/(p*(p-1)/2), 1) # Tried making sure proportion of 0s match what sparsity I want .3, .5, etc 
    
    if(min(eigen(psd)$values) > 1e-8){
      message("Found correct PSD matrix after ", fail.count, " failed attempt(s)")
      return(psd)
    }
    
    fail.count = fail.count + 1
    
  }
  message("No valid matrix found after ", max.iter, " iterations")
  return(NULL)
  
}




# cov2cor(solve(sparse.prec.mat(mat = redundant.sigma)))






# Get "sparse" latent correlation matrix, true sparse latent network, and precision matrix of the latent variables
latent.gen = function(mat){
  prec.mat = sparse.prec.mat(mat = mat,
                             prop = sparsity)
  
  
  
  # Take "prec.mat" and turn back into a sparse latent correlation matrix
  latent.sparse = cov2cor(solve(prec.mat))
  colnames(latent.sparse) = colnames(prec.mat) = c(paste0("p", 1:(p -1)),
                                                   "target")
  rownames(latent.sparse) = rownames(prec.mat) = c(paste0("p", 1:(p -1)),
                                                   "target")
  return(list(latent.sparse = latent.sparse,
              prec.mat = prec.mat))
}

