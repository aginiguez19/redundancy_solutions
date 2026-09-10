# Required packages
library(igraph)
library(lavaan)
library(psychonetrics)
library(qgraph)


# Additional functions
source(file = "02_functions.R")


####################################
# Simulation Network Recovery
###################################

# Conditions:

# Test:
# p = 5
# sparsity = .3
# lat.cor = .5
# lvl.redun = .9





conditions = expand.grid(
  p = c(10, 20),
  lat.cor = c(.3, .5, .7),  
  sparsity = c(0, .27, .43, .61),
  lvl.redun = c(.7, .8, .9) 
)




ncond = nrow(conditions)
n.iter = 5

solutions = list(
  removal   = apply.removal,
  composite = apply.composite,
  lnm       = apply.lnm
)

perf.list = list()

matrices = vector(mode = "list", length = ncond)



for (c in 1:ncond) {
  p = conditions$p[c]
  lat.cor = conditions$lat.cor[c]
  sparsity = conditions$sparsity[c]
  lvl.redun = conditions$lvl.redun[c]
  
  matrices[[c]] = vector(mode = "list", length = n.iter)
  
  cat(sprintf("\nCondition %d / %d | p=%d  lat.cor=%.1f  sparsity=%.1f  lvl.redun=%.2f\n", # Needs to be edited
              c, ncond, p, lat.cor, sparsity, lvl.redun))
  
  for (i in 1:n.iter){
    
    tryCatch({
      
    # Generate
    lat.mat = cor.gen(nvar = p, mn.cor = lat.cor) # Start by generating a latent pXp correlation mat
    
    # Get (p+1) X (p+1) redundant sigma
    redundant.sigma = sigma.gen(mat = lat.mat, peripheral.loadings = .9, clone.loading = lvl.redun)
    redundant.sigma = sigma.sparse.sim2(mat = redundant.sigma)$sigma.corr
    redundant.sigma = round(redundant.sigma, 3)
  
 
    # True
    
    # lambda.true = diag(x = 1, nrow = p + 1, ncol = p)
    # lambda.true[p + 1, p] = 2
    # mod.lnm = lnm(cors = redundant.sigma,
    #               nobs = 10000,
    #               lambda = lambda.true,
    #               omega_zeta = "full",
    #               identification = "variance") |>
    #   runmodel()
    # omega.latent = getmatrix(mod.lnm, "omega_zeta")
    
    
    cfamod <- write.cfa(p = p, cloneloading = lvl.redun)
    fitcfa <- cfa(model = cfamod, sample.cov = redundant.sigma, sample.nobs = 100000, std.lv = TRUE)
    true.R <- lavInspect(fitcfa, "cov.lv")
    omega.latent <- round(cor2pcor(true.R), 3)
    
    

   # sum(abs(omega.latent[lower.tri(omega.latent)]) < 1e-4)

    composite.sigma <- composite.gen(mat = lat.mat,
                                     clone.loading = lvl.redun,
                                     peripheral.loadings = .9)
    composite.sigma = sigma.sparse.sim2(mat = composite.sigma)$sigma.corr
    composite.sigma = round(composite.sigma, 3)
    
    I.mat <- diag(x = 1, nrow = p, ncol = p + 1)
    I.mat[p, p + 1] <- 1
    omega.composite <- round(cor2pcor(cov2cor(I.mat %*% composite.sigma %*% t(I.mat))), 3)
    
    # Get (p+1) X (p+1) pseudo-redundancy sigma 
    pseudo.sigma = donothing.gen(mat = lat.mat, clone.loading = lvl.redun)
    pseudo.sigma = sigma.sparse.sim2(pseudo.sigma)$sigma.corr
    pseudo.sigma = round(pseudo.sigma, 3)
    
    # True
    
    omega.pseudo = round(cor2pcor(pseudo.sigma), 3)
    
    mats = list(
      latent.corr.mat = lat.mat,
      redundant.sigma   = redundant.sigma,
      composite.sigma   = composite.sigma,
      pseudo.sigma      = pseudo.sigma,
      omega.latent      = omega.latent,
      omega.composite   = omega.composite,
      omega.pseudo      = omega.pseudo
    )
    matrices[[c]][[i]] = mats
      
      # Each sigma paired with its own truth 
    sigmas = list(
      redundant = mats$redundant.sigma,
      composite = mats$composite.sigma,
      pseudo = mats$pseudo.sigma
      )
      
    truths = list(
      redundant = mats$omega.latent,                  # p×p
      composite = mats$omega.composite,               # p×p
      pseudo = mats$omega.pseudo                   # p+1×p+1
      )
      
      for (sigma.name in names(sigmas)) { # Pulls redundant, composite, then pseudo
        for (sol.name in names(solutions)) { # Pulls each solution and applies for each sigma above
          # Solutions gets the correct solution function and the corresponding sigma
          est.mat  <- solutions[[sol.name]](sigmas[[sigma.name]], p) # Need p so functions can run 
          true.mat <- truths[[sigma.name]]              # matched truth
          metrics  <- calc.all.metrics(est.mat, true.mat,
                                       avg.last.two = (sigma.name == "pseudo")) # Will always be the correct truth that all solutions are tested against
          # Should match the sigma name so for example the composite sigma has all solutions applied to it and then tested against the true composite
          
          perf.list[[length(perf.list) + 1]] = data.frame(
            condition = c,
            p = p, 
            lat.cor = lat.cor, 
            sparsity = sparsity,
            lvl.redun = lvl.redun,
            iteration  = i,
            sigma.type = sigma.name,
            truth = sigma.name, # Decide whether to add additional column or change sigma.type to truth 
            solution   = sol.name,
            as.data.frame(metrics),
            row.names = NULL
          )
          
        }
      }
    }, error = function(e) {cat(sprintf("  Iter %d FAILED: %s\n", i, conditionMessage(e)))
      
    })
  }
  
  saveRDS(object = perf.list, 
          file = "/Users/aginigue/Desktop/simulation_resultsP2_100v5.rds")
  saveRDS(object = matrices,
          file = "/Users/aginigue/Desktop/simulation_matricesP2_100v5.rds")
}

perf.df <- do.call(rbind, perf.list) # call rbind and connect perf.list together 
# Need to save results as well and the matrices 
row.names(perf.df) = NULL
View(perf.df)

# mat = matrices[[46]][[1]]$pseudo.sigma
# mat.sol = apply.removal(mat, p = 10)
# est.ig  <- pcor.to.igraph(mat.sol)
# est.ig = delete_edges(est.ig, E(est.ig)[abs(E(est.ig)$weight) < 1e-3])

# true.mat.final = collapse.pseudo(matrices[[46]][[1]]$omega.pseudo, 10)
# true.ig.final = pcor.to.igraph(true.mat.final)
# true.ig.final = delete_edges(true.ig.final, E(true.ig.final)[abs(E(true.ig.final)$weight) < 1e-4])
# 
# calc.apl.bias(est.net = est.ig, true.net = true.ig.final)
# 
# mean_distance(graph = true.net,
#               weights = 1/abs(E(true.net)$weight))
# 20.45
# E(est.net)$weight
# 
# mean_distance(graph = est.net,
#               weights = 1/abs(E(est.net)$weight))
# 
# calc.apl.bias = function(est.net, true.net){
#   if (is_igraph(est.net) & is_igraph(true.net)) {
#     # True APL
#     true.apl = mean_distance(graph = true.net,
#                              weights = 1/abs(E(true.net)$weight))
#     # Est APL
#     est.apl = mean_distance(graph = est.net,
#                             weights = 1/abs(E(est.net)$weight))
#     # Relative bias
#     return((est.apl - true.apl))
#   }
# }
# 
# 
# 
# perf.list = readRDS("/Users/aginigue/Desktop/simulation_resultsP2_100v4.rds")
# perf.df = do.call(rbind, perf.list)
# row.names(perf.df) = NULL
# 
# 
# mats = readRDS("/Users/aginigue/Desktop/simulation_matricesP2_100v4.rds")
# View(mats)
# View(mats)
# true = mats[[118]][[56]]$omega.composite
# 
# est = apply.lnm(mats[[118]][[56]]$composite.sigma, p = 20)
# 
# est = pcor.to.igraph(est)
# est = delete_edges(est, E(est)[abs(E(est)$weight) < 1e-4])
# 
# true = pcor.to.igraph((true))
# true = delete_edges(true, E(true)[abs(E(true)$weight) < 1e-4]) # Important
# calc.betweenness.bias(est, true)
# 
# 
# betweenness(graph = est,
#             normalized = TRUE,
#             weights = 1/abs(E(est)$weight))[20]
# 
# betweenness(graph = true,
#             normalized = TRUE,
#             weights = 1/abs(E(true)$weight))[20]
# 
# 
# idx = 1:(p-1)
# 
# sqrt(mean((est[idx, p] - true[idx, p])^2))


# true.net = matrices[[1]][[10]]$omega.pseudo
# true.net = collapse.pseudo(true.net, 10)
# est.net = apply.composite(matrices[[1]][[10]]$pseudo.sigma, p = 10)
# 
# 
# 
# # qgraph(true.net, edge.labels = TRUE)
# # qgraph(est.net, edge.labels = TRUE)
# true.net2 = pcor.to.igraph(true.net)
# est.net2 = pcor.to.igraph(est.net)
# 
# # Closeness
# closeness(graph = true.net2, normalized = FALSE,
#           weights = 1/abs(E(true.net2)$weight))
# 
# 
# E(est.net2)$weight[abs(E(est.net2)$weight) < 1e-4] <- 0
# E(true.net2)$weight[E(true.net2)$weight < 1e-4] = 0
# 
# mean_distance(graph = est.net2,
#               weights = 1/abs(E(est.net2)$weight))
# 
# strength(graph = est.net2, weights = abs(E(est.net2)$weight))
# 
# mean_distance(graph = true.net2,
#               weights = 1/abs(E(true.net2)$weight))
# calc.apl.rbias(est.net = est.net, true.net = true.net)
# for (i in 1:5){
#  lat.mat = cor.gen(nvar = p, mn.cor = lat.cor) # Start by generating a latent pXp correlation mat
#  sparse.latent.mat = latent.gen(mat = lat.mat)$latent.sparse # Then make it sparse, now we can use 
#  
#  # Get (p+1) X (p+1) redundant sigma
#  redundant.sigma = sigma.gen(mat = sparse.latent.mat, peripheral.loadings = .9, clone.loading = .9)
#  
#  # True
# 
#  lambda.true = diag(x = 1, nrow = p + 1, ncol = p)
#  lambda.true[p + 1, p] = 1
#  mod.lnm = lnm(cors = redundant.sigma,
#                nobs = 10000,
#                lambda = lambda.true,
#                omega_zeta = "full",
#                identification = "loadings") %>% 
#    runmodel()
#  omega.latent = getmatrix(mod.lnm, "omega_zeta") # Hard to tell where da 0s at but matches relatively well with
#  
#  # Get (p+1) X (p+1) composite sigma
#  
#  composite.sigma = composite.gen(mat = sparse.latent.mat, clone.loading = lvl.redun, peripheral.loadings = .9)
#  
#  # True
#  
#  I = diag(x = 1,
#           nrow = p,
#           ncol = p + 1)
#  I[p, p + 1] = 1
#  sum.mat = I %*% composite.sigma %*% t(I)
#  sum.mat = cov2cor(sum.mat)
#  omega.composite = cor2pcor(sum.mat) # Same as omega_latent. Because they have the same spare latent correlation matrix?
#  
#  
#  
#  # Get (p+1) X (p+1) pseudo-redundancy sigma 
#  pseudo.sigma = donothing.gen(mat = sparse.latent.mat)
# 
#  
#  # True
# 
#  omega.pseudo = cor2pcor(pseudo.sigma)
#  
#  # Save all matrices in each iteration
#  matrices[[i]] = list(
#    latent.corr.mat = lat.mat, 
#    sparse.latent.mat = sparse.latent.mat,
#    redundant.sigma   = redundant.sigma,
#    composite.sigma   = composite.sigma,
#    pseudo.sigma      = pseudo.sigma,
#    omega.latent      = omega.latent,
#    omega.composite   = omega.composite,
#    omega.pseudo      = omega.pseudo
#  )
# }
#  
# 
# 
#  
# # Comparison loop
# 
# solutions = list(
#   removal   = apply.removal,
#   composite = apply.composite,
#   lnm       = apply.lnm
# )
# 
# perf.list = list()
# 
# for (i in seq_along(matrices)) {
#   
#   mats = matrices[[i]]
#   
#   # Each sigma paired with its own truth 
#   sigmas = list(
#     redundant = mats$redundant.sigma,
#     composite = mats$composite.sigma,
#     pseudo    = mats$pseudo.sigma
#   )
#   
#   truths = list(
#     redundant = mats$omega.latent,                  # p×p
#     composite = mats$omega.composite,               # p×p
#     pseudo    = mats$omega.pseudo                   # p+1×p+1
#   )
#   
#   for (sigma.name in names(sigmas)) { # Pulls redundant, composite, then pseudo
#     for (sol.name in names(solutions)) { # Pulls each solution and applies for each sigma above
#       # Solutions gets the correct solution function and the corresponding sigma
#       est.mat  <- solutions[[sol.name]](sigmas[[sigma.name]], p) # Need p so functions can run 
#       true.mat <- truths[[sigma.name]]              # matched truth
#       metrics  <- calc.all.metrics(est.mat, true.mat) # Will always be the correct truth that all solutions are tested against
#       # Should match the sigma name so for example the composite sigma has all solutions applied to it and then tested against the true composite
#       
#       perf.list[[length(perf.list) + 1]] = data.frame(
#         iteration  = i,
#         sigma.type = sigma.name,
#         solution   = sol.name,
#         as.data.frame(metrics),
#         row.names = NULL
#       )
#     }
#   }
# }
# 
# perf.df <- do.call(rbind, perf.list) # call rbind and connect perf.list together 
# # Need to save results as well and the matrices 
# row.names(perf.df) = NULL # Remove weird row numbering 
# View(perf.df)

 

# # Betweenness issue
# 
# est.net = pcor.to.igraph(apply.composite(matrices[[2]]$pseudo.sigma, p))
# E(est.net)$weight
# true.net = pcor.to.igraph(matrices[[2]]$omega.pseudo) # Betweenness for target and pseudoclone in omega.pseudo is 0 in every iteration, this makes the relative bias for betweenness INF
# E(true.net)$weight    # I think writing the methods will help clarify what's going on and how I can solve this problem 
# calc.betweenness.rbias = function(est.net, true.net){
#   if (is_igraph(est.net) & is_igraph(true.net)) {
#     
#     # True betweenness of target
#     true.betweenness = betweenness(graph = true.net,
#                                    normalized = TRUE,
#                                    weights = 1/abs(E(true.net)$weight))[p]
#     # Est betweenness of target
#     est.betweenness = betweenness(graph = est.net,
#                                   normalized = TRUE,
#                                   weights = 1/abs(E(est.net)$weight))[p]
#     # Relative bias
#     return((est.betweenness - true.betweenness)/true.betweenness)
#   }
# }

# I want to take each sigma and apply a solution to it and compare it to the truth
#  
#  removal.mat = redundant.sigma[1:p, 1:p]
#  net1 = cor2pcor(removal.mat)
#  
#  
#  I = diag(x = 1,
#           nrow = p,
#           ncol = p + 1)
#  I[p, p + 1] = 1
#  sum.mat = I %*% redundant.sigma %*% t(I)
#  sum.mat = cov2cor(sum.mat)
#  colnames(sum.mat) = c(paste0("p", 1:(p -1)),
#                        "target")
#  net2 = cor2pcor(sum.mat)
#  
#  
# 
#  lambda = diag(x = 1, nrow = p + 1, ncol = p)
#  lambda[p + 1, p] = 1
#  mod.rnm = lnm(cors = redundant.sigma,
#                lambda = lambda,
#                nobs = 10000,
#                omega_zeta = "full",
#                identification = "variance") %>%
#    runmodel()
#  
# net3 = getmatrix(mod.lnm, "omega_zeta") 
 
# Then compare each matrix to the truth or the omega_latent 


# Need to take care of the average strength, etc comparison...
# mat = matrices[[1]]$omega.pseudo
# mat2 = matrices[[1]]$omega.pseudo[1:p, 1:p]
# igraph.ob = pcor.to.igraph(mat)
# igraph.ob2 = pcor.to.igraph(mat2)
# strength(igraph.ob)
# strength(igraph.ob2)
# betweenness(igraph.ob, normalized = TRUE,
#             weights = 1/abs(E(igraph.ob)$weight))





# # Latents
# sparse.latent.mats$latent.sparse #Matrix that is input into all other generating mechanisms
# sparse.latent.mats$true.latent.net
# 
# # Sigmas
# sparse.sigma.mats$sigma.sparse # Matrix in which all solutions are applied to 
# sparse.sigma.mats$true.removal.net
# 
# # Composites
# sparse.composite.mats$sigma.composite.sparse # Other matrix in which all solutions are applied to 
# sparse.composite.mats$true.composite.net

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






# Add 07/15, maybe


# apply.lnm <- function(sigma, p, nobs = 10000) {
#   lambda <- diag(x = 1, nrow = p + 1, ncol = p)
#   lambda[p + 1, p] <- 1
#   mod <- lnm(cors = sigma, lambda = lambda, nobs = nobs,
#              omega_zeta = "full", identification = "variance") %>% runmodel()
#   om <- getmatrix(mod, "omega_zeta")
#   attr(om, "improper") <- any(abs(om[upper.tri(om)]) >= 1) || !all(is.finite(om))
#   om
# }


