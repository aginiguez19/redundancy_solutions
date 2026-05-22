# Required packages
library(igraph)
library(psychonetrics)
library(qgraph)


# Additional functions
source(file = "02_functions.R")


####################################
# Simulation network recovery 
#
###################################

# Conditions:

# Test:
# p = 5
# sparsity = .3
# lat.cor = .5
# lvl.redun = .9

p = c(10, 20)
lat.cor = c(.1, .3, .5, .7, .9)
sparsity = c(0, .3, .5, .7)
lvl.redun = c(.7, .85, 9)

method = c("removal", "sum", "RNM")

# sparse.deltas.c = matrix(NA, nrow = 100, ncol = p)
for (i in 1:10) {
  # Generate latent correlation matrix 
  lat.mat = cor.gen(nvar = p,
                    mn.cor = lat.cor,
                    sd = .05)

  

  if (method == "removal") {
    remove.mat = item.mat[1:p, 1:p]
    
    net = cor2pcor(remove.mat)
 
    
  } else if (method == "sum") {
    I = diag(x = 1,
             nrow = p,
             ncol = p + 1)
    I[p, p + 1] = 1
    sum.mat = I %*% item.mat %*% t(I)
    sum.mat = cov2cor(sum.mat)
    colnames(sum.mat) = c(paste0("p", 1:(p -1)),
                                "target")
    net = cor2pcor(sum.mat)
    
  } else if (method == "RNM") {
    
    nobs = 1000
    vars = colnames(item.mat)
    lambda = matrix(0, 6, 1, dimnames = list(vars, "LV"))
    lambda[p:(p + 1), "LV"] = 1
    mod.rnm = rnm(cors = item.mat,
                  lambda = lambda, 
                  vars = vars,
                  nobs = nobs,
                  omega_epsilon = "full",
                  identification = "variance") %>%
      runmodel()
    
    net = getmatrix(mod.rnm, "omega_epsilon")
  }
}






for (i in 1:5) {
  lat.mat = cor.gen(nvar = p,
                    mn.cor = lat.cor)
  
  sparse.latent.mats = latent.gen(mat = lat.mat)
  zero.pattern = which(sparse.latent.mats$prec.mat == 0)
  
  sparse.sigma.mats = sigma.gen(mat = sparse.latent.mats$latent.sparse,
                      pattern = zero.pattern)
  
  sparse.composite.mats = composite.gen(mat = sparse.latent.mats$latent.sparse,
                                        pattern = zero.pattern)
}

# Latents
sparse.latent.mats$latent.sparse #Matrix that is input into all other generating mechanisms
sparse.latent.mats$true.latent.net

# Sigmas
sparse.sigma.mats$sigma.sparse # Matrix in which all solutions are applied to 
sparse.sigma.mats$true.removal.net

# Composites
sparse.composite.mats$sigma.composite.sparse # Other matrix in which all solutions are applied to 
sparse.composite.mats$true.composite.net




## Testing ##

# results = vector("logical", 1000)
# results = NA
# failed = vector("logical", 1000)
# failed.mats = list()
# 
# for (i in 1:1000){ 
#   lat.mat = cor.gen(nvar = p, mn.cor = lat.cor)
#   
#   tryCatch({
#     sparse.latent.mats = sparse.lat.mats(mat = lat.mat)
#     zero.pattern = which(sparse.latent.mats$prec.mat == 0)
#     item.prec = solve(item.mat$sigma.true)
#     item.prec[zero.pattern] = 0
#     results[i] = any(eigen(item.prec)$values < 0)
#   }, error = function(e) {
#     failed[i] <<- TRUE
#     failed.mats[[as.character(i)]] <<- lat.mat  # store the matrix that caused the failure
#   })
# }
# 
# which(failed)
# sum(failed)
# View(as.data.frame(failed))
# View(as.data.frame(results))

