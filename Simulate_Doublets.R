# ---------------------------------------------------------------------------- #
Simulate_Expression = function(
  nGenes     = 500,
  batchCells = 1000,
  group.prob = c(0.5, 0.5),
  de.prob    = 0.1
){
  suppressPackageStartupMessages(library(splatter))
  sim = splatSimulate(
    nGenes     = nGenes,
    batchCells = batchCells,
    group.prob = group.prob,
    method     = "groups",
    de.prob    = de.prob
  )
  sim$type = 'singlet'
  return(sim)
}


# ---------------------------------------------------------------------------- #
Create_Doublets = function(
  sim,
  doublet_percentage = c(.01,.03,.05,0.1,0.15,0.2)
){

  if(is.null(sim))
      stop('sim not provided')

  suppressPackageStartupMessages(library(scater))

  sl = list()
  for(i in doublet_percentage){
    # determines the number of cells to sample
    nsamp = round(ncol(sim)/(1-i) - ncol(sim))
    g1 = sample(which(sim$Group == 'Group1'), nsamp)
    g2 = sample(which(sim$Group == 'Group2'), nsamp)

    sim_doub = sim[,g1]
    assays(sim_doub)$TrueCounts   = assays(sim_doub)$TrueCounts + assays(sim[,g2])$TrueCounts
    assays(sim_doub)$counts       = assays(sim_doub)$counts     + assays(sim[,g2])$counts
    assays(sim_doub)$CellMeans    = assays(sim_doub)$CellMeans  + assays(sim[,g2])$CellMeans
    cname = paste0(colData(sim_doub)$Cell, 'drop')
    colData(sim_doub)$Cell      = cname
    rownames(colData(sim_doub)) = cname
    assays(sim_doub) = lapply(assays(sim_doub), function(x){colnames(x)=cname;x})
    sim_doub$type = 'doublet'

    sim_all = cbind(sim_doub, sim)
    sim_all = normalise(sim_all)
    name = paste0('d',i)
    sl[[name]] = sim_all
  }
  sl[['d0']] = sim
  sl = sl[order(names(sl))]
  return(sl)
}


# ---------------------------------------------------------------------------- #
Simulate_Dropout = function(
  sl         = NULL,
  mids       = c(1,2,3),
  shape_vars = c(0, 0.15,0.25,0.35)
){
  if(is.null(sl))
    stop('sim not provided')

  ldrop = list()
  for(name in names(sl)){
    sim = sl[[name]]
    param = newSplatParams(dropout.type = 'experiment',
                           nGenes       = nrow(sim),
                           batchCells   = ncol(sim))

    for(mid in mids){
      for(shape_var in shape_vars){

        message(paste(mid, shape_var))
        param = setParam(param, 'dropout.mid',   rep(mid, ncol(sim)))
        param = setParam(param, 'dropout.shape', rnorm(ncol(sim), -1, shape_var))
        param = setParam(param, 'dropout.type', 'cell')
        sim_drop = splatter:::splatSimDropout(
          sim,
          params=param)
        sim_drop$type = sim$type[match(sim_drop$Cell, sim$Cell)]

        dname = paste(name, 'mid',mid, 'shape',shape_var, sep='.')
        ldrop[[dname]] = sim_drop
      }
    }
  }
  ldrop = ldrop[order(names(ldrop))]
  return(ldrop)

}

# ---------------------------------------------------------------------------- #
# ss = Simulate_Single_Cell()
Simulate_Single_Cell = function(){

  message('Expression ...')
    sim   = Simulate_Expression()

  message('Doublets ...')
    sl    = Create_Doublets(sim)

  message('Dropout ...')
    ldrop = Simulate_Dropout(sl)

  return(ldrop)
}
