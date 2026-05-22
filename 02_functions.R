

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


# Generating Functions ----------------------------------------------------

## For testing only ##
# mat = lat.mat 
# prop = .3
# max.iter = 100



# Get latent correlation matrix
cor.gen = function(nvar = NULL,
                   bg.cor = 0.30,
                   mn.cor = 0.50,
                   prop.cor = 1,
                   iter.lim = 100,
                   sd = 0.1,
                   tol.eig = 1e-8,
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
    temp.psd = E$vectors %*% (diag(vals, nvar, nvar)) %*% t(E$vectors) # eigenvectors are orthonormal so that its inverse is equal to its tranpose
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


# Get sparse precision matrix
sparse.prec.mat = function(mat, prop = .3, max.iter = 10000){
  mat = solve(mat) # Turn into a precision matrix first
  idx.lower = which(lower.tri(mat)) # Get indices of lower triangle
  p.zero = floor(prop * length(idx.lower)) # Number of 0s 
  
  for (i in 1:max.iter){
    mat.new = mat  
    idx.zero = sample(idx.lower, p.zero) # Sample from indices that should be 0
    mat.new[idx.zero] = 0 # Change the sampled indices to 0 
    
    mat.new[upper.tri(mat.new)] = t(mat.new)[upper.tri(mat.new)] # Make upper triangle the same
    
    vals = eigen(mat.new, symmetric = TRUE, only.values = TRUE)$values # Check eigenvalues
    
    if (all(vals >= 1e-8)){
      return(mat.new)
    }
  }
  return(NULL) # Add `if (is.null(result)) next` to simulation to skip failed iterations
}



# Get "sparse" latent correlation matrix, true sparse latent network, and precision matrix of the latent variables
latent.gen = function(mat){
  prec.mat = sparse.prec.mat(mat = mat,
                             prop = sparsity)
  
  true.net = -1 * cov2cor(prec.mat) 
  diag(true.net) = 0 # Obtain a true sparse latent network
  
  # Take "prec.mat" and turn back into a sparse latent correlation matrix
  latent.sparse = cov2cor(solve(prec.mat))
  colnames(latent.sparse) = colnames(true.net) = colnames(prec.mat) = c(paste0("p", 1:(p -1)),
                                                                        "target")
  rownames(latent.sparse) = rownames(true.net) = rownames(prec.mat) = c(paste0("p", 1:(p -1)),
                                                                        "target")
  return(list(latent.sparse = latent.sparse,
              true.latent.net = true.net,
              prec.mat = prec.mat))
}



# Get model-implied correlation matrix (sigma) with redundancy and generating mechanism for removal solution (true (item)removal network) 
sigma.gen = function(mat, peripheral.loadings = .9, clone.loading = .9, pattern){
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
  sigma.prec = solve(sigma) # Get precision matrix of sigma
  block = sigma.prec[1:p, 1:p] # Get (p X p)  block
  block[pattern] = 0 # Add matching pattern of sparsity
  sigma.prec[1:p, 1:p] = block # Add block back into full "sigma.prec"
  
  target.zeros = which(sigma.prec[p, ] == 0) # Find which of the target's edges are 0
  sigma.prec[p+1, target.zeros] = 0 # Match those 0 edges with the clone's edges
  sigma.prec[target.zeros, p+1] = 0
  
  vals = eigen(sigma.prec, symmetric = TRUE, only.values = TRUE)$values
  if(!all(vals >= 1e-8)) warning("Precision matrix is not positive semi-definite!")

  sigma.sparse = cov2cor(solve(sigma.prec))
  colnames(sigma.sparse) = rownames(sigma.sparse) = c(paste0("p", 1:(p -1)),
                                        "target",
                                        "clone")
  
  
  # Get sparse sigma true (p X p) and turn into network
  sigma.true = sigma[1:p, 1:p]
  true.prec = solve(sigma.true)
  true.prec[pattern] = 0
  
  sigma.true.net = -1*cov2cor(true.prec)
  diag(sigma.true.net) = 0 
  colnames(sigma.true.net) = rownames(sigma.true.net) = c(paste0("p", 1:(p -1)),
                                                         "target")

  return(list(sigma.sparse = sigma.sparse,
              true.removal.net = sigma.true.net))
}


## For testing only ## 
# mat = sparse.latent.mats$latent.sparse


# Get "sparse" composite item correlation matrix and true sparse composite network
composite.gen = function(mat, clone.loading = .9, peripheral.loadings = .9, pattern) {
  p = nrow(mat) #number of variables
  cor = clone.loading * clone.loading #desired correlation between the two variables that form a sum score
 
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
   
   
   # Add sparsity to sigma.composite
   
   sigma.composite.prec = solve(sigma.composite) # Precision matrix
   block.comp = sigma.composite.prec[1:p, 1:p] # Get (p X p)  composite block
   block.comp[pattern] = 0 # Add matching pattern of sparsity
   sigma.composite.prec[1:p, 1:p] = block.comp # Add block back into full "sigma.composite.prec"
   
   target.zeros = which(sigma.composite.prec[p, ] == 0) # Find which of the target's edges are 0
   sigma.composite.prec[p+1, target.zeros] = 0 # Match those 0 edges with the clone's edges
   sigma.composite.prec[target.zeros, p+1] = 0
   
   
   vals = eigen(sigma.composite.prec, symmetric = TRUE, only.values = TRUE)$values
   if(!all(vals >= 1e-8)) warning("Composite precision matrix is not positive semi-definite!")
   
   
   sigma.composite.sparse = cov2cor(solve(sigma.composite.prec))
   colnames(sigma.composite.sparse) = rownames(sigma.composite.sparse) = c(paste0("p", 1:(p -1)),
                                                       "target",
                                                       "clone")
  
   
   # Get sparse composite sigma true (p X p) and turn into network
   I = diag(x = 1, nrow = p, ncol = p + 1)
   I[p, p + 1] = 1
   
   sum.mat = I %*% sigma.composite.sparse %*% t(I)
   sum.mat = cov2cor(sum.mat)
   
   true.composite.net = -1*cor2pcor(sum.mat)
   diag(true.composite.net) = 0
   colnames(true.composite.net) = rownames(true.composite.net) = c(paste0("p", 1:(p-1)), "target")
   
   return(list(sigma.composite.sparse = sigma.composite.sparse,
               true.composite.net = true.composite.net))
}




# Performance Measures ----------------------------------------------------

 
calc.correlation = function(est.net, true.net) {
  true.edges = true.net[upper.tri(true.net)]
  est.edges = est.net[upper.tri(est.net)]
  
  correlation = cor(true.edges, est.edges)
  return(correlation)
}


# Must be an igraph object!

weighted.density = function(igrph.object, abs = TRUE){
  if (abs == FALSE){
    density = sum(E(igrph.object)$weight)/(length(igrph.object)*(length(igrph.object-1)/2))
  }
  else {
    density = sum(abs(E(igrph.object)$weight))/(length(igrph.object)*(length(igrph.object-1)/2))
  }
  
  return(density)
}


calc.wd.rbias = function(est.net, true.net){
  if (is_igraph(est.net) & is_igraph(true.net)) {
    # True weighted density
    true.density = weighted.density(true.net)
    # Est weighted density
    est.density = weighted.density(est.net)
    # Relative bias
    return((est.density - true.density)/true.density)
  }
}

calc.apl.rbias = function(est.net, true.net){
  if (is_igraph(est.net) & is_igraph(true.net)) {
    # True APL
    true.apl = mean_distance(graph = true.net,
                             weights = 1/abs(E(true.net)$weight))
    # Est APL
    est.apl = mean_distance(graph = est.net,
                            weights = 1/abs(E(est.net)$weight))
    # Relative bias
    return((est.apl - true.apl)/true.apl)
  }
}

calc.strgth.rbias = function(est.net, true.net){
  if (is_igraph(est.net) & is_igraph(true.net)) {
    # True strength of target  
    true.strength = strength(graph = true.net,
             weights = abs(E(true.net)$weight))[p]
    # Est strength of target
    est.strength = strength(graph = est.net,
                            weights = abs(E(est.net)$weight))[p]
    # Relative bias
    return((est.strength - true.strength)/true.strength)

  }
}

calc.expctinflu.rbias = function(est.net, true.net) {
  if (is_igraph(est.net) & is_igraph(true.net)) {
    # True expected influence of target
    true.expected = strength(graph = true.net, 
                             weights = E(true.net)$weight)[p]
    # Est expected influence of target 
    est.expected = strength(graph = est.net,
                            weights = E(est.net)$weight)[p]
    # Relative bias
    return((est.expected - true.expected)/true.expected)
    
  }
}


calc.closeness.rbias = function(est.net, true.net){
  if (is_igraph(est.net) & is_igraph(true.net)) {
    # True closeness of target

    true.closeness = closeness(graph = true.net,
                              normalized = TRUE,
                             weights = 1/abs(E(true.net)$weight))[p]
    # Est closeness of target
    est.closeness = closeness(graph = est.net,
                             normalized = TRUE,
                             weights = 1/abs(E(est.net)$weight))[p]
    # Relative bias
    return((est.closeness - true.closeness)/true.closeness)
  }
}



calc.betweenness.rbias = function(est.net, true.net){
  if (is_igraph(est.net) & is_igraph(true.net)) {
    
    # True betweenness of target
    true.betweenness = betweenness(graph = true.net,
                               normalized = TRUE,
                               weights = 1/abs(E(true.net)$weight))[p]
    # Est betweenness of target
    est.betweenness = betweenness(graph = est.net,
                              normalized = TRUE,
                              weights = 1/abs(E(est.net)$weight))[p]
    # Relative bias
    return((est.betweenness - true.betweenness)/true.betweenness)
  }
}


calc.RMSE = function(est.net, true.net) {
  if (is_igraph(est.net) & is_igraph(true.net)) {
    est.net = as_adjacency_matrix(est.net,
                                  attr = "weight",
                                  sparse = FALSE)
    true.net = as_adjacency_matrix(true.net,
                                   attr = "weight",
                                   sparse = FALSE)
    }
  idx = 1:(p - 1)
  return(sqrt(mean((est.net[idx, p] - true.net[idx, p])^2)))
  }


