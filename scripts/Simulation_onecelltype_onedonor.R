
# Idea ------------------------------------------------------------------------------------------------------------

#Start small: one cell type in one donor (later: increase cell types and/or increase number of donors)
#Use binomial distribution to model cell type (later: beta binomial or more difficult frameworks)
#Go for precision instead of hypothesis testing (no hypothesis testing possible in one cell type on donor)
#Compare based on:
  #* Scopit-like where c = 1 (probability of at least observing 1 count)
  #* ARE: absolute relative error
  #* Coefficient of variation (simple since easy formula for binomial distribution)

#Idea of sampling:
  #* Choose sampling parameters
  #* Construct parameter grid
  #* Choose number of samples B
  #* For each combination of sampling parameters (=row in parameter grid) sample B times
  #* Calculate the comparison statistics

#Presentation
  #* Present summary stats
  #* Graphs of summary stats
  #* Distribution of ARE and ALOR

# Setup -----------------------------------------------------------------------------------------------------------
rm(list = ls())

#params
N <- c(10, 50, 100, 1000, 10000) #number of cells
p <- c(0.01, 0.05, 0.1) #proportions
c <- c(1) #Scopit limit

B <- c(1000) #repetitions

param.grid <- expand.grid(N = N, p = p, c = c, B = B)


# Simulation ------------------------------------------------------------------------------------------------------

res.sim <- lapply(1:nrow(param.grid), function(i){
  set.seed(2609*i)
  #recover parameters
  with(param.grid[i,], {
  #get samples
  counts <- rbinom(B, N, p)
  p_hat <- counts / N
  
  #return stats
  list("params" = param.grid[i,],
       # "counts" = counts, "p_hat" = p_hat,
       "scopit" = mean(counts >= c),
       "ARE" = abs(p_hat - p) / p,
       "ALOR" = abs(log(p_hat / p)),
       "CV" = sd(p_hat) / mean(p_hat))
  })
})
names(res.sim) <- paste0("sim_", 1:nrow(param.grid))

# Summary stats ---------------------------------------------------------------------------------------------------

matrix(unlist(lapply(res.sim, function(l){
  c(l$params, l$scopit, mean(l$ARE), mean(l$ALOR), l$CV)
})), nrow = nrow(param.grid), byrow = TRUE)

# Visualization ---------------------------------------------------------------------------------------------------


