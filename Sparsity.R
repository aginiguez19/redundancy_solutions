
# Functions ---------------------------------------------------------------


# Get latent correlation matrix
cor.gen = function(nvar = NULL,
                   bg.cor = 0.30,
                   mn.cor = 0.50,
                   prop.cor = 1,
                   iter.lim = 100,
                   sd = 0.05,
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

# Get sigma correlation matrix 
sigma.gen = function(mat, peripheral.loadings = .9, clone.loading = .9){
  dim.mat = dim(mat)
  lambda = diag(x = peripheral.loadings, nrow = dim.mat[1] + 1, ncol = dim.mat[2])
  lambda[nrow(lambda), ncol(lambda)] = clone.loading
  item.mat = lambda %*% mat %*% t(lambda)
  dim.theta = dim(item.mat)
  theta = diag(x = 1 - (peripheral.loadings)^2,
               nrow = dim.theta[1],
               ncol = dim.theta[2])
  theta[nrow(theta), ncol(theta)] = 1 - (clone.loading)^2
  
  sigma = item.mat + theta # Sigma (p+1 X p+1) with redundancy
  
  
  # Add sparsity to sigma
  # Add same pattern of sparsity as sparse precision matrix of latent variables 
  # sigma.prec = solve(sigma) # Get precision matrix of sigma
  # block = sigma.prec[1:p, 1:p] # Get (p X p)  block
  # block[pattern] = 0 # Add matching pattern of sparsity
  # sigma.prec[1:p, 1:p] = block # Add block back into full "sigma.prec"
  # 
  # target.zeros = which(sigma.prec[p, ] == 0) # Find which of the target's edges are 0
  # sigma.prec[p+1, target.zeros] = 0 # Match those 0 edges with the clone's edges
  # sigma.prec[target.zeros, p+1] = 0
  # 
  # vals = eigen(sigma.prec, symmetric = TRUE, only.values = TRUE)$values
  # if(!all(vals >= 1e-8)) warning("Precision matrix is not positive semi-definite!")
  # 
  # sigma.sparse = cov2cor(solve(sigma.prec))
  colnames(sigma) = rownames(sigma) = c(paste0("p", 1:(p -1)),
                                        "target",
                                        "clone")
  # No longer need the code above. Tried something that was not actually related to our design 
  
  
  # # Get sparse sigma true (p X p) and turn into network
  # sigma.true = sigma[1:p, 1:p]
  # true.prec = solve(sigma.true)
  # true.prec[pattern] = 0
  # 
  # sigma.true.net = -1*cov2cor(true.prec)
  # diag(sigma.true.net) = 0 
  # colnames(sigma.true.net) = rownames(sigma.true.net) = c(paste0("p", 1:(p -1)),
  #                                                        "target")
  # No longer need the code above. Tried something that was not actually related to our design 
  
  return(sigma)
  # true.removal.net = sigma.true.net))
}

# Correlation to partial correlation matrix
cor2pcor <- function(mat){
  precision.mat <- solve(mat)
  pcor.mat <- qgraph::wi2net(precision.mat)
  pcor.mat <- as.matrix(pcor.mat)
  return(pcor.mat)
}

# Covariance Matrix to correlation 

cov2cor = function(mat){
  inv.sd = diag(1/sqrt(diag(mat)))
  matrix = inv.sd %*% mat %*% inv.sd
  return(matrix)
}



# Example -----------------------------------------------------------------

lat.cor = 0.5
p = 10

# 1. Get latent correlation matrix 
lat.mat = cor.gen(nvar = 10, mn.cor = lat.cor)

# 2. Get sigma correlation matrix with clone node

redundant.sigma = sigma.gen(mat = lat.mat, 
                            peripheral.loadings = 0.9,
                            clone.loading = 0.9)


# 3 Rebuild sparsity function 


sparse.prec.mat = function(mat, prop = .3, tol = 1e-1, max.iter = 10000){
  
  p = nrow(mat)
  
  fail.count = 0
  prec.mat.full = solve(mat)
  
  
  idx.lower = which(lower.tri(prec.mat.full)) # Get indices of lower triangle
  idx.pair = (p - 2) * p + p # the target and clone index
  idx.target = (1:(p-2)-1)*p + (p-1) # target's index
  idx.clone = (1:(p-2)-1) * p + p # clone index
  
  idx.pool = setdiff(idx.lower, c(idx.pair,idx.clone)) # elements that are not the pair and clone's index
  
  p.zero = floor(prop * length(idx.lower)) # Number of 0s
  p.zero = min(p.zero, length(idx.pool)) # safeguard if sparsity is set extremely high

 
  
  for (i in 1:max.iter) {
    prec.mat = prec.mat.full
    idx.zero = sample(idx.pool, p.zero) # Sample from indices that should be 0
    prec.mat[idx.zero] = 0 # Change the sampled indices to 0
  
    prec.mat[idx.clone] = prec.mat[idx.target] # Match clone and target
    
    prec.mat[upper.tri(prec.mat)] = t(prec.mat)[upper.tri(prec.mat)] # Make upper triangle the same
    
    # p = nrow(prec.mat)
    # prec.mat[p, 1:(p-1)-1] = prec.mat[p-1,1:(p-1)-1]
    # prec.mat[1:(p-1)-1, p] = prec.mat[1:(p-1)-1, p-1]
    
    # Preserve pattern of 0s 
    
    idx.zero.restore = which(prec.mat == 0 & lower.tri(prec.mat))
  
    # Check for PSDness
    E = eigen(x = prec.mat, symmetric = TRUE)
    vals = E$values
    vals[vals < 0] = tol
    psd = E$vectors %*% (diag(vals, p, p)) %*% t(E$vectors)
    
    psd[idx.zero.restore] = 0 # Change the sampled indices to 0 
    psd[upper.tri(psd)] = t(psd)[upper.tri(psd)] # Make upper triangle the same
    
    
    if(min(eigen(psd)$values) > 1e-8 && psd[p, p-1] != 0){
      message("Found correct PSD matrix after ", fail.count, " failed attempt(s)")
      return(psd)
    }
    
    fail.count = fail.count + 1
    
  }
  message("No valid matrix found after ", max.iter, " iterations")
  return(NULL)
  
}

# 4 Test 

# Parameters 
p = 10 # Will need to change if nvar is changed
lat.cor = 0.5
prop = 0.3
lat.mat = cor.gen(nvar = 10, mn.cor = lat.cor)
redundant.sigma = sigma.gen(mat = lat.mat, 
                            peripheral.loadings = 0.9,
                            clone.loading = 0.9)


mat = sparse.prec.mat(mat = redundant.sigma, prop = 0.3)

partial.mat = -1*cov2cor(mat)
diag(partial.mat) = 0
qgraph::qgraph(partial.mat, edge.labels = TRUE,
                theme = "colorblind")


sim = 1:3
conds[sim]
small.sim = matrix(NA, 500, 3)
conds = c(0.27, 0.43, 0.61)
conds.t = c(0.30, 0.50, 0.70)
for(sims in 1:3){
  for(i in 1:500){
    p = 10
    lat.cor = 0.7
    prop = conds[sims]
    lat.mat = cor.gen(nvar = 10, mn.cor = lat.cor)
    
    redundant.sigma = sigma.gen(mat = lat.mat,
                                peripheral.loadings = 0.9,
                                clone.loading = 0.9)
    
    
    mat = sparse.prec.mat(mat = redundant.sigma, prop = prop)
    
    partial.mat = -1*cov2cor(mat)
    diag(partial.mat) = 0
    # qgraph::qgraph(partial.mat, edge.labels = TRUE,
    #                theme = "colorblind")
    small.sim[i,sims] = sum(partial.mat[lower.tri(partial.mat)] != 0 )/length(partial.mat[lower.tri(partial.mat)])
  }
}

(1 - colMeans(small.sim)) - conds.t
apply(small.sim, 2, sd)




