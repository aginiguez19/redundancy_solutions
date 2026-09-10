

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


# Get sparse precision matrix
sparse.prec.mat.sim2 = function(mat, prop = .3, tol = 1e-1, max.iter = 10000){
  
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



sparse.prec.mat.sim1 = function(mat, prop = .3, tol = 1e-1, max.iter = 10000){
  
  
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
    
    if(min(eigen(psd)$values) > 1e-8){
      message("Found correct PSD matrix after ", fail.count, " failed attempt(s)")
      return(psd)
    }
    
    fail.count = fail.count + 1
    
  }
  message("No valid matrix found after ", max.iter, " iterations")
  return(NULL)
  
}


# Generate true latent network 

write.cfa <- function(p, cloneloading) {
  lineend <- "\n "
  cfa.mod <- ""
  v <- 0
  for (j in 1:(p-1))   {
    v <- v + 1
    cfa.mod <- paste0(cfa.mod, "F", j, " =~ 1*p", v, sep = "")
    cfa.mod <- paste0(cfa.mod, lineend, sep = " ")
  }
  cfa.mod <- paste0(cfa.mod, "F", p, " =~ .9*target + ", cloneloading, "*clone")
  return(cfa.mod)
}

# cfamod <- write.cfa(p = p, cloneloading = lvl.redun)
# fitcfa <- cfa(model = cfamod, sample.cov = redundant.sigma, sample.nobs = 100000, std.lv = TRUE)
# true.R <- lavInspect(fitcfa, "cov.lv")
# true.net <- cor2pcor(true.R)

# cfamod = write.cfa(p = 10, cloneloading = 0.9)
# 
# fitcfa <- cfa(model = cfamod, sample.cov = redundant.sigma, sample.nobs = 100000, std.lv = TRUE)
# true.R <- lavInspect(fitcfa, "cov.lv")
# true.net = cor2pcor(true.R)
# 
# apply.lnm(sigma = redundant.sigma, p = 10, nobs = 10000)

# Get sparse precision matrix then convert back into correlation matrix
sigma.sparse.sim2 = function(mat){
  prec.mat = sparse.prec.mat.sim2(mat = mat,
                             prop = sparsity)
  
  # Take "prec.mat" and turn back into a sparse latent correlation matrix
  sigma.corr = cov2cor(solve(prec.mat))
  colnames(sigma.corr) = colnames(prec.mat) = c(paste0("p",
                                                       1:(p -1)),
                                                "target",
                                                "clone")
  rownames(sigma.corr) = rownames(prec.mat) = c(paste0("p", 1:(p -1)),
                                                "target",
                                                "clone")
  return(list(sigma.corr = sigma.corr,
              prec.mat = prec.mat))
}



latent.gen.sim1 = function(mat){
  prec.mat = sparse.prec.mat.sim1(mat = mat,
                             prop = sparsity)
  
  # Take "prec.mat" and turn back into a sparse latent correlation matrix
  latent.sparse = round(cov2cor(solve(prec.mat)), 3)
  colnames(latent.sparse) = colnames(prec.mat) = c(paste0("p",
                                                       1:(p -1)),
                                                "target")
  rownames(latent.sparse) = rownames(prec.mat) = c(paste0("p", 1:(p -1)),
                                                "target")
  return(list(latent.sparse = latent.sparse,
              prec.mat = prec.mat))
}


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


## For testing only ## 
# mat = sparse.latent.mats$latent.sparse


# Get "sparse" composite item correlation matrix and true sparse composite network
composite.gen = function(mat, clone.loading = .9, peripheral.loadings = .9) {
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
   
   
   # Add sparsity to sigma.composite
   
   # sigma.composite.prec = solve(sigma.composite) # Precision matrix
   # block.comp = sigma.composite.prec[1:p, 1:p] # Get (p X p)  composite block
   # block.comp[pattern] = 0 # Add matching pattern of sparsity
   # sigma.composite.prec[1:p, 1:p] = block.comp # Add block back into full "sigma.composite.prec"
   # 
   # target.zeros = which(sigma.composite.prec[p, ] == 0) # Find which of the target's edges are 0
   # sigma.composite.prec[p+1, target.zeros] = 0 # Match those 0 edges with the clone's edges
   # sigma.composite.prec[target.zeros, p+1] = 0
   # 
   # 
   # vals = eigen(sigma.composite.prec, symmetric = TRUE, only.values = TRUE)$values
   # if(!all(vals >= 1e-8)) warning("Composite precision matrix is not positive semi-definite!")
   # 
   # 
   # sigma.composite.sparse = cov2cor(solve(sigma.composite.prec))
   colnames(sigma.composite) = rownames(sigma.composite) = c(paste0("p", 1:(p -1)),
                                                       "target",
                                                       "composite_clone")
   # 
   # 
   # # Get sparse composite sigma true (p X p) and turn into network
   # I = diag(x = 1, nrow = p, ncol = p + 1)
   # I[p, p + 1] = 1
   # 
   # sum.mat = I %*% sigma.composite.sparse %*% t(I)
   # sum.mat = cov2cor(sum.mat)
   # 
   # true.composite.net = -1*cor2pcor(sum.mat)
   # diag(true.composite.net) = 0
   # colnames(true.composite.net) = rownames(true.composite.net) = c(paste0("p", 1:(p-1)), "target")
   # 
   return(sigma.composite)
               # true.composite.net = true.composite.net))
}
# No longer need the code above. Tried something that was not actually related to our design 






# Get pseudo redundancy sigma matrix


# Start with single item measurement model
# Similar to random network except one is an exact copy just loading on a different factor
# Mat is then a (p+1) x (p+1) matrix generated from cor.gen
# Gonna need latent sparse mat (latent.gen()) when sparse.prec.mat is finalized
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



# Performance Measures ----------------------------------------------------


# Must be an igraph object!
calc.correlation = function(est.net, true.net) {
  true.edges = true.net[upper.tri(true.net)]
  est.edges = est.net[upper.tri(est.net)]
  
  correlation = cor(true.edges, est.edges)
  return(correlation)
}



weighted.density = function(igrph.object, abs = TRUE){
  if (abs == FALSE){
    density = sum(E(igrph.object)$weight)/(length(igrph.object)*(length(igrph.object-1)/2))
  }
  else {
    density = sum(abs(E(igrph.object)$weight))/(length(igrph.object)*(length(igrph.object-1)/2))
  }
  
  return(density)
}




calc.wd.bias = function(est.net, true.net){
  if (is_igraph(est.net) & is_igraph(true.net)) {
    # True weighted density
    true.density = weighted.density(true.net)
    # Est weighted density
    est.density = weighted.density(est.net)
    # Relative bias
    return((est.density - true.density))
  }
}

calc.apl.bias = function(est.net, true.net){
  if (is_igraph(est.net) & is_igraph(true.net)) {
    # True APL
    true.apl = mean_distance(graph = true.net,
                             weights = 1/abs(E(true.net)$weight))
    # Est APL
    est.apl = mean_distance(graph = est.net,
                            weights = 1/abs(E(est.net)$weight))
    # Relative bias
    return((est.apl - true.apl))
  }
}

calc.strgth.bias = function(est.net, true.net, avg.last.two = FALSE){
  if (is_igraph(est.net) & is_igraph(true.net)) {
    # True strength of target  
    true.strength.all = strength(graph = true.net, weights = abs(E(true.net)$weight))
    true.strength = true.strength.all[p]
    # Est strength of target
    est.strength = strength(graph = est.net,
                            weights = abs(E(est.net)$weight))[p]
    # Relative bias
    return((est.strength - true.strength))

  }
}

calc.expctinflu.bias = function(est.net, true.net, avg.last.two = FALSE) {
  if (is_igraph(est.net) & is_igraph(true.net)) {
    # True expected influence of target
    true.expected.all = strength(graph = true.net, weights = E(true.net)$weight)
    true.expected = true.expected.all[p]
    # Est expected influence of target 
    est.expected = strength(graph = est.net,
                            weights = E(est.net)$weight)[p]
    # Relative bias
    return((est.expected - true.expected))
    
  }
}



calc.closeness.bias = function(est.net, true.net, avg.last.two = FALSE){
  if (is_igraph(est.net) & is_igraph(true.net)) {
    # True closeness of target
    true.closeness.all = closeness(graph = true.net, normalized = TRUE,
                                   weights = 1/abs(E(true.net)$weight))
    true.closeness = true.closeness.all[p]
    # Est closeness of target
    est.closeness = closeness(graph = est.net,
                             normalized = TRUE,
                             weights = 1/abs(E(est.net)$weight))[p]
    # Relative bias
    return((est.closeness - true.closeness))
  }
}



calc.betweenness.bias = function(est.net, true.net, avg.last.two = FALSE){
  if (is_igraph(est.net) & is_igraph(true.net)) {
    
    # True betweenness of target
    true.betweenness.all = betweenness(graph = true.net,
                                       normalized = TRUE,
                                       weights = 1/abs(E(true.net)$weight))
    true.betweenness = true.betweenness.all[p]
    # Est betweenness of target
    est.betweenness = betweenness(graph = est.net,
                              normalized = TRUE,
                              weights = 1/abs(E(est.net)$weight))[p]
    # Bias
    return((est.betweenness - true.betweenness))
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


# Helper Functions --------------------------------------------------------


# Partial correlation matrix → weighted igraph (zero diagonal first)
pcor.to.igraph = function(mat) {
  graph_from_adjacency_matrix(mat, mode = "undirected", weighted = TRUE, diag = FALSE)
}


# Three solutions: each takes (p+1)×(p+1) sigma, returns p×p network

apply.removal = function(sigma, p) {
  cor2pcor(sigma[1:p, 1:p]) # Subset by p×p 
}

apply.composite = function(sigma, p) {
  I = diag(x = 1, nrow = p, ncol = p + 1)
  I[p, p + 1] = 1
  cor2pcor(cov2cor(I %*% sigma %*% t(I)))
}

apply.lnm = function(sigma, p, nobs = 10000) {
  lambda = diag(x = 1, nrow = p + 1, ncol = p)
  lambda[p + 1, p] = 1
  mod = lnm(cors = sigma, lambda = lambda, nobs = nobs,
             omega_zeta = "full", identification = "variance") |>
    runmodel()
  getmatrix(mod, "omega_zeta")
}


# Est.mat and true.mat are both p×p partial correlation matrices

calc.all.metrics = function(est.mat, true.mat, avg.last.two = FALSE) {
  est.ig  <- pcor.to.igraph(est.mat)
  # E(est.ig)$weight[abs(E(est.ig)$weight) < 1e-4] = 0 # Added to deal with inflated vals
  # est.ig = delete_edges(est.ig, E(est.ig)[abs(E(est.ig)$weight) < 1e-3])
  true.ig <- pcor.to.igraph(true.mat)
  # E(true.ig)$weight[abs(E(true.ig)$weight) < 1e-4] = 0
  # true.ig = delete_edges(true.ig, E(true.ig)[abs(E(true.ig)$weight) < 1e-3])
  
  if (avg.last.two){
    true.mat.final = collapse.pseudo(true.mat, p)
    true.ig.final = pcor.to.igraph(true.mat.final)
    # E(true.ig.final)$weight[abs(E(true.ig.final)$weight) < 1e-4] = 0
    # true.ig.final = delete_edges(true.ig.final, E(true.ig.final)[abs(E(true.ig.final)$weight) < 1e-3])
  } else {
    true.mat.final = true.mat
    true.ig.final = true.ig
  }
  
  list(
    correlation       = calc.correlation(est.mat, true.mat.final),
    wd.bias          = calc.wd.bias(est.ig, true.ig.final),
    apl.bias         = calc.apl.bias(est.ig, true.ig.final),
    strgth.bias      = calc.strgth.bias(est.ig, true.ig.final, avg.last.two),
    expctinflu.bias  = calc.expctinflu.bias(est.ig, true.ig.final, avg.last.two),
    closeness.bias   = calc.closeness.bias(est.ig, true.ig.final, avg.last.two),
    betweenness.bias = calc.betweenness.bias(est.ig, true.ig.final, avg.last.two),
    RMSE              = calc.RMSE(est.mat, true.mat.final)
  )
}

# Collapse Pseudo so it can be compared when using global measures and correlations

collapse.pseudo = function(mat, p) {
  out = mat[1:p, 1:p]                                  
  out[1:(p-1), p] = (mat[1:(p-1), p] + mat[1:(p-1), p + 1]) / 2    # avg peripheral -> {target, clone}
  out[p, 1:(p-1)] = out[1:(p-1), p]                               # keep symmetric                                        
  return(out)
}


# Test for results why is everything overestimating 
# Why is removal doing bad too? 

# p = 10
# sparsity = 0
# mat  = cor.gen(nvar = 10)
# mat2 = latent.gen(mat)
# sigma = donothing.gen(mat2$latent.sparse)
# apply.removal(sigma, p =4) # Removal solution
# true = cor2pcor(sigma) # True pseudo
# collapse.pseudo(mat = true, p =4) # Works correctly
# 
# sigma = sigma.gen(mat2$latent.sparse, peripheral.loadings = 0.9,
#           clone.loading =  0.7)





# Test why I saw LNM perform poorly even when it is the correct solution
# p = 10
# sparsity = 0
# mat  = cor.gen(nvar = 10)
# mat2 = latent.gen(mat)
# 
# sigma = sigma.gen(mat2$latent.sparse, peripheral.loadings = 0.9,
#                   clone.loading =  0.7)
# 
# solution = apply.lnm(sigma, p = 10)
# lambda.true = diag(x = 1, nrow = p + 1, ncol = p)
# lambda.true[p + 1, p] = 1
# mod.lnm = lnm(cors = sigma,
#     nobs = 10000,
#     lambda = lambda.true,
#     omega_zeta = "full",
#     identification = "loadings") %>% 
#   runmodel()
# omega.latent = getmatrix(mod.lnm, "omega_zeta")
# 
# true = pcor.to.igraph(omega.latent)
# est = pcor.to.igraph(solution)
# 
# calc.strgth.bias(est, true)
# strength(true)
# strength(est)

