pxp_mat = cor.gen(nvar = 5) # start with pxp matrix
uncor_mat = cbind(pxp_mat, c(pxp_mat[1:4, 5], rnorm(1, mean = .3, sd = .01)))
copy = uncor_mat[,6]
uncor_mat = rbind(uncor_mat, c(copy, 1)) # Get different form of redundancy

eigen(uncor_mat)$values # every so often negative


cor2pcor(uncor_mat)

