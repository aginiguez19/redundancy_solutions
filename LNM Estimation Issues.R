
library(psychonetrics)
library(qgraph)
library(dplyr)

# Results 
perf.list = readRDS("/Users/aginigue/Desktop/simulation_resultsP2_100v3.rds")
perf.df = do.call(rbind, perf.list)
row.names(perf.df) = NULL

# Matrices 
mats = readRDS("/Users/aginigue/Desktop/simulation_matricesP2_100v3.rds")

# Conditions

conditions = expand.grid(
  p = c(10, 20),
  lat.cor = c( .3, .5, .7), # Removed .7 and .9 
  sparsity = c(0, .3, .5, .7), # Removed .7
  lvl.redun = c(.7, .8, .9) # Removed .85
)


# Functions Needed --------------------------------------------------------

apply.lnm = function(sigma, p, nobs = 10000) {
  lambda = diag(x = 1, nrow = p + 1, ncol = p)
  lambda[p + 1, p] = 1
  mod = lnm(cors = sigma, lambda = lambda, nobs = nobs,
            omega_zeta = "full", identification = "variance") |>
    runmodel()
  getmatrix(mod, "omega_zeta")
  parameters(mod)
}


# Composite ---------------------------------------------------------------



mat = mats[[36]][[1]]$composite.sigma # Condition and iteration that returns bad matrix

# P = 20, lat.cor = 0.7, sparsity = 0.3, redundancy level = 0.8

apply.lnm(sigma = mat, 20, nobs = 10000) # Check lambdas

# Change loading of peripheral nodes 0.9 --> 0.85
composite.gen = function(mat, clone.loading = .9, peripheral.loadings = .85) {
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

# Get generating latent correlation matrix from same condition and iteration
sparse.latent.mat = mats[[36]][[1]]$sparse.latent.mat

# Get new composite matrix with changed peripheral loadings 
new.sigma = composite.gen(mat = sparse.latent.mat)

# Apply lnm solution again and see results differ
apply.lnm(new.sigma, p = 20, nobs = 100)



# Unique Variable ---------------------------------------------------------

# Same process as in composite 
# Incorrect one 
mat2 = mats[[30]][[90]]$pseudo.sigma # Check size of matrix in case p needs to change 
# p = 20, lat.cor = 0.7, sparsity = 0, redundancy = 0.8

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


# Modify to have each item load on one latent varible
apply.lnm.new = function(sigma, p, nobs = 10000) {
  lambda = diag(x = 1, nrow = p + 1, ncol = p + 1)
  mod = lnm(cors = sigma, lambda = lambda, nobs = nobs,
            omega_zeta = "full", identification = "variance") |>
    runmodel()
  getmatrix(mod, "omega_zeta")
  parameters(mod)
} # Rerun to see if result is improved 

apply.lnm.new(sigma = mat2, 20, nobs = 10000) # Result does improve, could the issue stem from being forced to think two items measure one latent var when they don't?




# No issue check 
mat3 = mats[[1]][[1]]$pseudo.sigma
apply.lnm(sigma = mat3, 10, nobs = 10000)

# p = 10, lat.cor = 0.3, sparsity = 0, redundancy = 0.7



