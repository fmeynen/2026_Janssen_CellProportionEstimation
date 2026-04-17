
# Idea ------------------------------------------------------------------------------------------------------------

#Start small: one cell type in one donor (later: increase cell types and/or increase number of donors)
#Use binomial distribution to model cell type (later: beta binomial or more difficult frameworks)
#Go for precision instead of hypothesis testing (no hypothesis testing possible in one cell type on donor)
#Compare based on:
  #* Scopit-like where c = 1 (probability of at least observing 1 count)
  #* ARE: absolute relative error
  #* Coefficient of variation (simple since easy formula for binomial distribution)

# Setup -----------------------------------------------------------------------------------------------------------
rm(list = ls())

#params
N <- c(10, 50, 100, 500, 1000) #number of cells
p <- c(0.01,0.05,0.1) #proportions

param.grid <- expand.grid(N = N,p = p)

#repetitions
B <- 1000

#simulation
#switch the loops around.
lapply(1:B,function(b){
  set.seed(2609 + b)
  res.b <- lapply(1:nrow(param.grid), function(i){
    N <- param.grid$N[i]
    p <- param.grid$p[i]
    
    counts <- rbinom(1, N, p)
    p_hat <- counts/N
    
  })
})