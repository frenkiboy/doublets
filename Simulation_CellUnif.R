```{r}
sl = list()
for(i in c(0,.01,.03,.05,0.1,0.15,0.2)){
  
  message(i)
  if(i > 0){  
    g1 = sample(which(sim$Group == 'Group1'), ncol(sim)*i)
    g2 = sample(which(sim$Group == 'Group2'), ncol(sim)*i)
    
    sim_doub = sim[,g1]
    assays(sim_doub)$TrueCounts = assays(sim_doub)$TrueCounts + assays(sim[,g2])$TrueCounts
    assays(sim_doub)$counts     = assays(sim_doub)$counts     + assays(sim[,g2])$counts
    cname = paste0(colData(sim_doub)$Cell, 'drop')
    colData(sim_doub)$Cell = cname
    rownames(colData(sim_doub)) = cname
    assays(sim_doub) = lapply(assays(sim_doub), function(x){colnames(x)=cname;x})
    sim_doub$type = 'doublet'
    
    sim_sing = sim[,-sample(c(g1,g2), length(g1))]
    sim_sing$type = 'singlet'
    
    sim_all = cbind(sim_doub, sim_sing)
    sim_all = normalise(sim_all)
  }else{
    sim_all = sim
    sim_all$type = 'singlet'
  }
  
  param = newSplatParams(dropout.type = 'experiment', 
                         nGenes       = nrow(sim_all), 
                         batchCells   = ncol(sim_all))
  # param = setParam(param, 'dropout.mid', 1)
  # param = setParam(param, 'dropout.shape', -1)
  param = setParam(param, 'dropout.mid', rpois(ncol(sim_all),2))
  param = setParam(param, 'dropout.shape', rep(-1, ncol(sim_all)))
  param = setParam(param, 'dropout.type', 'cell')
  
  sim_drop = splatter:::splatSimDropout(
    sim_all, 
    params=param)
  sim_drop$type = sim_all$type[match(sim_drop$Cell, sim_all$Cell)]
  
  
  name = paste0('d',i)
  sl[[name]] = sim_drop
  
  #   for(dropout in c('unif','cell','doub')){
  #     param = setParam(param, 'dropout.mid', rpois(ncol(sim_all),10))
  #   param = setParam(param, 'dropout.shape', rep(-1, ncol(sim_all)))
  #   param = setParam(param, 'dropout.type', 'cell')
  #   message(dropout)
  # }                     
  
  
  
}
sl = sl[order(names(sl))]
```