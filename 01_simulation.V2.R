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

