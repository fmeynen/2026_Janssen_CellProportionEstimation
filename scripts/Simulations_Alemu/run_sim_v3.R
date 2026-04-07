library(BiocParallel)
library(tidyverse)  
library(DirichletReg) 

options(dplyr.summarise.inform = FALSE)


# simulation setting
K.set = c(5, 10, 20, 30) # number of cell types
n.set = c(1, 5, 10, 30)  # the number of biological samples  
N.set <- c(500, 1000,2500, 5000, 10000, 25000, 50000, 100000, 500000, 1000000, 5000000) # the number of cells to sequence 
gma.set <- c(0.25, 1, 4) # level of biological variability

R = 200 # number of replicates per sample
B = 200 # number of repitition for sampling donors

param.mat <- expand_grid(K=K.set, N=N.set, n=n.set, gma=gma.set) %>%
  subset(!(n==1 & gma!=1))

# simulate counts
bplapply(1:B,function(b){
  cat(b,"- ")
  set.seed(9551+b)
  res.b <- lapply(1:nrow(param.mat), function(i){
    #cat(i, "_")
    K <- param.mat$K[i]
    N <- param.mat$N[i]
    n <- param.mat$n[i]
    gma <- param.mat$gma[i]
    
    ## dirichlet parameters 
    B0  <- rnorm(n=K, mean=0, sd=0.5)  
    theta0   = gma*exp(B0)
    
    ## multinomial probabilities
    pi0   = rdirichlet(n, theta0) 
    
    ## sim counts
    if(n==1){
      counts.donor <- lapply(1:n, function(j){
        count = t(rmultinom(R, size = N, prob = pi0[j,]))
        sample_id = paste0("sample_", 1:R)
        rownames(count) <- sample_id 
        colnames(count) =  paste0("ppl_", 1:K)
        
        count0.longdf <- count %>% as.data.frame() %>%
          rownames_to_column("sample_id") %>%
          pivot_longer(cols = -sample_id,  names_to = "cellType", 
                       values_to = "count") 
        
        ## add true frequencies
        true.freq <- data.frame(freq.true= as.numeric(pi0[j,]), 
                                cellType=paste0("ppl_", 1:K))
        
        count0.longdf <-  count0.longdf %>%
          left_join(true.freq, by="cellType") %>%
          mutate(donor_id = paste0("donor_", j)) %>%
          dplyr::select(donor_id, sample_id, cellType, count, freq.true)
        
        count0.longdf
      }) %>%bind_rows() %>% 
        mutate(freq.est = count/N,
               error.est = freq.est - freq.true,
               ARE= abs(error.est)/freq.true,
               ALOR= abs(log(freq.est/freq.true))) %>%
        group_by(donor_id, sample_id) %>%
        summarise(all_pass1p_ARE=all(ARE <= 0.01),
                  all_pass5p_ARE=all(ARE <= 0.05),
                  all_pass10p_ARE=all(ARE <= 0.1),
                  all_pass20p_ARE=all(ARE <= 0.2), 
                  
                  all_pass1p_ALOR=all(ALOR <= 0.01),
                  all_pass5p_ALOR=all(ALOR <= 0.05),
                  all_pass10p_ALOR=all(ALOR <= 0.1),
                  all_pass20p_ALOR=all(ALOR <= 0.2), 
                  nchk=n()) %>%
        ungroup() %>%
        #group_by(donor_id) %>%
        summarise(prob1p_ARE=mean(all_pass1p_ARE),
                  prob5p_ARE=mean(all_pass5p_ARE),
                  prob10p_ARE=mean(all_pass10p_ARE),
                  prob20p_ARE=mean(all_pass20p_ARE),
                  
                  prob1p_ALOR=mean(all_pass1p_ALOR),
                  prob5p_ALOR=mean(all_pass5p_ALOR),
                  prob10p_ALOR=mean(all_pass10p_ALOR),
                  prob20p_ALOR=mean(all_pass20p_ALOR),
                  nchk=n()) %>%
        #ungroup() %>%
        mutate(cond=i, N = N, K=K, n=n, gma=gma)
    }else{
      counts.donor <- lapply(1:n, function(j){
        count = t(rmultinom(R, size = N, prob = pi0[j,]))
        sample_id = paste0("sample_", 1:R)
        rownames(count) <- sample_id 
        colnames(count) =  paste0("ppl_", 1:K)
        
        count0.longdf <- count %>% as.data.frame() %>%
          rownames_to_column("sample_id") %>%
          pivot_longer(cols = -sample_id,  names_to = "cellType", 
                       values_to = "count") 
        
        ## add true frequencies
        true.freq <- data.frame(freq.true= as.numeric(pi0[j,]), 
                                cellType=paste0("ppl_", 1:K))
        
        count0.longdf <-  count0.longdf %>%
          left_join(true.freq, by="cellType") %>%
          mutate(donor_id = paste0("donor_", j)) %>%
          dplyr::select(donor_id, sample_id, cellType, count, freq.true)
        
        count0.longdf
      }) %>%bind_rows() %>% 
        mutate(freq.est = count/N) %>%
        group_by(sample_id, cellType) %>%
        summarise(mfreq.est = mean(freq.est), mfreq.true=mean(freq.true)) %>%
        ungroup() %>%
        mutate(error.est = mfreq.est - mfreq.true,
               mARE= abs(error.est)/mfreq.true,
               mALOR= abs(log(mfreq.est/mfreq.true))) %>%
        # group_by(sample_id, cellType) %>%
        # summarise(mARE=mean(ARE)) %>%
        # ungroup() %>%
        group_by(sample_id) %>%
        summarise(all_pass1p_ARE=all(mARE <= 0.01),
                  all_pass5p_ARE=all(mARE <= 0.05),
                  all_pass10p_ARE=all(mARE <= 0.1),
                  all_pass20p_ARE=all(mARE <= 0.2), 
                  
                  all_pass1p_ALOR=all(mALOR <= 0.01),
                  all_pass5p_ALOR=all(mALOR <= 0.05),
                  all_pass10p_ALOR=all(mALOR <= 0.1),
                  all_pass20p_ALOR=all(mALOR <= 0.2), 
                  nchk=n()) %>%
        ungroup() %>%
        #group_by(sample_id) %>%
        summarise(prob1p_ARE=mean(all_pass1p_ARE),
                  prob5p_ARE=mean(all_pass5p_ARE),
                  prob10p_ARE=mean(all_pass10p_ARE),
                  prob20p_ARE=mean(all_pass20p_ARE),
                  
                  prob1p_ALOR=mean(all_pass1p_ALOR),
                  prob5p_ALOR=mean(all_pass5p_ALOR),
                  prob10p_ALOR=mean(all_pass10p_ALOR),
                  prob20p_ALOR=mean(all_pass20p_ALOR),
                  nchk=n()) %>%
        #ungroup() %>%
        mutate(cond=i, N = N, K=K, n=n, gma=gma)
    }
     
    counts.donor 
  }) %>%
    bind_rows() %>%
    mutate(sim.rep=b)
  
  saveRDS(res.b, paste0("Results/simResults/fixedMarginError/sim.count_itr_", b, ".rds"))
},BPPARAM = BiocParallel::MulticoreParam(workers = 30)) 

