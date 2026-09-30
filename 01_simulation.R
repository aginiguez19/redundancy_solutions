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
  lat.cor = c(0.3, 0.5, 0.7),
  sparsity = c(0, 0.3, 0.5), 
  lvl.redun = c(0.7, 0.8, 0.9) 
)



ncond = nrow(conditions)
n.iter = 50


solutions = list(
  removal = apply.removal,
  composite = apply.composite,
  lnm = apply.lnm)

perf.list = list()
edge.list = list()

matrices = vector(mode = "list", length = ncond)


start.time = Sys.time()
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
    
    latent.mat = cor.gen(nvar = p, mn.cor = lat.cor) # Start by generating a latent pXp correlation matrix
    sparse.latent.mat = sparse.lat(mat = latent.mat)$lat.mat # Then make the precision matrix sparse, now we can use 


  
    # Get (p+1) X (p+1) redundant sigma
    redundant.sigma = sigma.gen(mat = sparse.latent.mat, peripheral.loadings = .9, clone.loading = lvl.redun)
    
    # True
    
    cfamod = write.cfa(p = p, cloneloading = lvl.redun)
    fitcfa = cfa(model = cfamod, sample.cov = redundant.sigma, sample.nobs = 100000, std.lv = TRUE)
    true.R = lavInspect(fitcfa, "cov.lv")
    omega.latent = round(cor2pcor(true.R), 2)
    omega.latent[omega.latent < 0.05] = 0
    
    composite.sigma = composite.gen(mat = sparse.latent.mat,
                                     clone.loading       = lvl.redun,
                                     peripheral.loadings = .9)
    
    I.mat = diag(x = 1, nrow = p, ncol = p + 1)
    I.mat[p, p + 1] = 1
    omega.composite = round(cor2pcor(cov2cor(I.mat %*% composite.sigma %*% t(I.mat))), 2)
    omega.composite[omega.composite < 0.05] = 0
    
    # Get (p+1) X (p+1) pseudo-redundancy sigma 
    pseudo.sigma = donothing.gen(mat = sparse.latent.mat)
    
    # True
    
    omega.pseudo = round(cor2pcor(pseudo.sigma), 2)
    omega.pseudo[omega.pseudo < 0.05] = 0
    
    
    mats = list(
      latent.corr.mat = latent.mat, 
      sparse.latent.mat = sparse.latent.mat,
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
          est.mat  = solutions[[sol.name]](sigmas[[sigma.name]], p) # Need p so functions can run 
          true.mat = truths[[sigma.name]]              # matched truth
          metrics  = calc.all.metrics(est.mat, true.mat,
                                       avg.last.two = (sigma.name == "pseudo"))
          edge.metrics = calc.edge.metrics(est.mat, true.mat,
                                           avg.last.two = (sigma.name == "pseudo"))
          # Will always be the correct truth that all solutions are tested against
          # Should match the sigma name so for example the composite sigma has all solutions applied to it and then tested against the true composite
          
          perf.list[[length(perf.list) + 1]] = data.frame(
            condition = c,
            p = p, 
            lat.cor = lat.cor, 
            sparsity = sparsity,
            lvl.redun = lvl.redun,
            iteration  = i,
            sigma.type = sigma.name,
            truth = sigma.name, 
            solution   = sol.name,
            as.data.frame(metrics),
            row.names = NULL
          )
          
          
          edge.list[[length(edge.list) + 1]] = data.frame(
            condition = c,
            p = p, 
            lat.cor = lat.cor, 
            sparsity = sparsity,
            lvl.redun = lvl.redun,
            iteration  = i,
            sigma.type = sigma.name,
            truth = sigma.name, 
            solution   = sol.name,
            as.data.frame(edge.metrics),
            row.names = NULL)
          
        }
      }
    }, error = function(e) {cat(sprintf("  Iter %d FAILED: %s\n", i, conditionMessage(e)))
      
    })
  }
  
  saveRDS(object = perf.list, 
          file = "/Users/aginigue/Desktop/simulation_resultsR2.rds")
  saveRDS(object = matrices,
          file = "/Users/aginigue/Desktop/simulation_matricesR2.rds")
}

perf.df = do.call(rbind, perf.list) # call rbind and connect perf.list together 
edge.df = do.call(rbind, edge.list)
# Need to save results as well and the matrices 
row.names(perf.df) = NULL
row.names(edge.df) = NULL

end.time = Sys.time()
time.taken = end.time - start.time
print(time.taken)
xfun::session_info()





# 
# sparsity = 0.3
# p = 5
# lvl.redun = 0.9
# nsim = 500
# ratios = matrix(data = NA, nrow = 100, ncol = 1, byrow = TRUE)
# for (i in 1:nsim){
# latent.mat = cor.gen(nvar = p, mn.cor = 0.5)
# sparse.latent.mat = sparse.lat(mat = latent.mat)$lat.mat
# redundant.sigma = sigma.gen(mat = sparse.latent.mat, peripheral.loadings = .9, clone.loading = lvl.redun)
# cfamod = write.cfa(p = p, cloneloading = lvl.redun)
# fitcfa = cfa(model = cfamod, sample.cov = redundant.sigma, sample.nobs = 100000, std.lv = TRUE)
# true.R = lavInspect(fitcfa, "cov.lv")
# omega.latent = round(cor2pcor(true.R), 2)
# omega.latent[omega.latent < 0.05] = 0
# 
# sparse.ratio = sum(omega.latent[lower.tri(omega.latent)] == 0)/((p * (p-1))/2)
# ratios[i] = sparse.ratio
# }
# 
# 
# mean(ratios)
# ratios
