suppressWarnings(suppressMessages(require(GSVA)))
suppressWarnings(suppressMessages(require(snow)))
suppressWarnings(suppressMessages(require(rlecuyer)))
suppressWarnings(suppressMessages(require(reshape2)))
suppressWarnings(suppressMessages(require(parallel)))
suppressWarnings(suppressMessages(require(viper)))
no_cores = 2#detectCores()
# if(no_cores>=4)
#   no_cores =  no_cores - 1
op = options(digits.secs = 6)




# ------------------------------------------------------------------------------------------------------------------------------------
# SLEA: Sample Level Enrichment Analysis
# ------------------------------------------------------------------------------------------------------------------------------------
SLEA = function(E, genesets, method, M = NULL, permutations = 1000, filter_E = F, min_geneset_size = 3, max_shared = 10){
  # E = Expression matrix (rows=genes; columns=samples) scaled and recentered (use gene_expression_statistic).
  # genesets = list with two elemens: 1) NAME = vector with the names of the genes set; 2) GENES = named list of gene sets
  # method = name of the method used for the SLEA
  # M = Methylation binary matrix (optional). If null will be ignored
  # permutations. Number of permutations needed for the method=MEAN
  # filter_E. Logical indicating if genes in E not in the gene set should be removed/ Default FALSE.  
  E = as.matrix(E)
  E = E[ ! duplicated(rownames(E)) , ]
  E = E[ order(rownames(E)), ]
  genesets = SLEA.clean_genesets(genesets, E, min_geneset_size = min_geneset_size, max_shared = max_shared)
  message('Getting Enrichment Scores')
  # cat('Please be aware that this test should be applied on expression values that are on the same scale.\nIt assumes scaled and recentered expression data\n')
  # cat('If not scaled and recentered, please use the \'gene_expression_statistic()\' function\n')
  if (filter_E){
    E = E[ rownames(E) %in% unlist(lapply(genesets$GENES, names)), ]
    M = M[ rownames(M) %in% unlist(lapply(genesets$GENES, names)), ]
  }
  if ( ! is.null(M) ){
    clME = SLEA.clean_ME(M, E)
    M = clME$M
    E = clME$E
  }
  if( method == 'GSVA' ){
    cat('Calculating SLEA scores using GSVA\n')
    RESULTS = SLEA.GSVA(E, genesets, M, method = 'GSVA')
  }else if( method == 'ssGSEA' ){
    cat('Calculating SLEA scores using ssGSEA\n')
    RESULTS = SLEA.GSVA(E, genesets, M, method = 'ssgsea')
  }else if( method == 'MEAN' ){
    cat('Calculating SLEA scores using MEAN with', sum(permutations), 'permutations\n')
    RESULTS = SLEA.MEAN(E, genesets, permutations, M)
  }else if( method == 'GSEAlm' ){
    cat('Calculating SLEA scores using GSEAlm\n')
    RESULTS = SLEA.GSEAlm(E, genesets, M)
  }else if( method == 'ManW' ){
    cat('Calculating SLEA scores using ManW\n')
    RESULTS = SLEA.ManW(E, genesets, M)
  }else if( method == 'VIPER' ){
    cat('Calculating SLEA scores using VIPER\n')
    RESULTS = SLEA.VIPER(E, genesets, M)
  }else{
    cat('Invalid method. The available methods are: \'ssGSEA\', \'GSVA\', \'GSEAlm\', \'VIPER\', \'ManW\' and \'MEAN\' \n' )
    RESULTS = NULL
  }
  RESULTS$regulons = sapply(genesets$GENES, length)
  message('\nDone!\n')
  return(RESULTS)
}





SLEA.GSEAlm = function(E, genesets, M = NULL){ # NOTE: method modified from Shahzad Asif 2010
  gsm = geneset.list2binaryMat(gsl = genesets$GENES, E = E)
  gsm = gsm[ , ! duplicated(t(gsm))]
  ES = matrix(-1, ncol = ncol(E), nrow = ncol(gsm) )
  colnames(ES) = colnames(E)
  rownames(ES) = colnames(gsm)
  PVAL = ES
  FDR = ES
  pb = txtProgressBar(min = 0, max = ncol(ES), style = 3)
  for( s in colnames(ES)){
      x = gsealm(E, sa = s, gsm = gsm, M = M)
      ES[,s] = x$score
      PVAL[,s] = x$pval
      setTxtProgressBar(pb, which(colnames(ES) == s))
  }
  for( gs in rownames(PVAL) )
    FDR[gs,] = p.adjust(PVAL[gs,], method = 'fdr')
  RESULTS = list(ES=ES, PVAL=PVAL, FDR=FDR, method='lm')

  return(RESULTS)
}


gsealm = function(E, sa, gsm, M = NULL){ 
  expression_pattern = E[, sa]
  if ( ! is.null(M) ){ # If methylation, remove methylated genes
    methylation_pattern = M[ , sa]
    hidden_genes = names(methylation_pattern)[ methylation_pattern == 1 ]
    expression_pattern = expression_pattern[ ! names(expression_pattern) %in% hidden_genes ]
  }
    X = gsm[ names(expression_pattern) , ]
    fit = summary(lm(expression_pattern~X))
    if ( ! all(colnames(gsm) %in% gsub('^X', '', rownames(fit$coefficients))) ){
      rownames(fit$coefficients) = unlist(lapply( rownames(fit$coefficients) , function(x) strsplit(x, ']')[[1]][2] ))
    }else{
      rownames(fit$coefficients) = gsub('^X', '', rownames(fit$coefficients))
    }
    return(list(score = fit$coefficients[ -1, 'Estimate'], pval = fit$coefficients[ -1, 'Pr(>|t|)' ]))
}



SLEA.ManW = function(E, genesets, M = NULL){ # NOTE: method proposed by Claudia I. Hernández-Armenta in Beltrao's Group
  ES = matrix(-1, ncol = ncol(E), nrow = length(genesets$NAME) )
  colnames(ES) = colnames(E)
  rownames(ES) = genesets$NAME
  PVAL = ES
  FDR = ES
  pb = txtProgressBar(min = 0, max = length(genesets$NAME), style = 3)
  if( min(unlist(genesets$GENES)) < 0)
    E = apply(E, 2, function(x) (x-min(x)) / (max(x)-min(x)) ) # rescale gene expression 
  for( gs in genesets$NAME ){
    Eset = E[ rownames(E) %in% names(genesets$GENES[[gs]]), ]
    Eset[ rownames(Eset) %in% names(genesets$GENES[[gs]])[ genesets$GENES[[gs]] < 0 ] , ]  =  1 - Eset[ rownames(Eset) %in% names(genesets$GENES[[gs]])[ genesets$GENES[[gs]] < 0 ] , ] # invert values of repressed genes
    Enonset = E[ ! rownames(E) %in% rownames(Eset), ]
    for( s in colnames(ES)){
      x = ManW(s, Eset[,s], Enonset[,s], M)
      ES[gs,s] = x$score
      PVAL[gs,s] = x$pval
    }
    FDR[gs,] = p.adjust(PVAL[gs,], method = 'fdr')
    setTxtProgressBar(pb, which(genesets$NAME == gs))
  }
  RESULTS = list(ES=ES, PVAL=PVAL, FDR=FDR, method='ManW')
  return(RESULTS)
}


ManW = function(sa, setvalues, nonsetvalues, M = NULL){ 
  if( ! is.null(M) ){
    methylation_pattern = M[,sa]
    hidden_genes = names(methylation_pattern)[ methylation_pattern == 1 ]
    test = wilcox.test(setvalues[ ! names(setvalues) %in% hidden_genes ], nonsetvalues[ ! names(nonsetvalues) %in% hidden_genes ])
  }else{
    test = wilcox.test(setvalues, nonsetvalues)
  }
  pval = test$p.value
  score = ifelse(mean(setvalues) < mean(nonsetvalues), log10(pval), -log10(pval))
  return(list(score = score, pval = pval))
}




SLEA.VIPER = function(E, genesets, M = NULL){
  vpres = new_aREA(eset = E, genesets = genesets, M = M)
  ES = vpres$es
  NES = vpres$nes
  PVAL = vpres$pval
  FDR = t(apply(PVAL, 1, p.adjust, method = 'fdr'))
  RESULTS = list(ES=ES, NES=NES, PVAL=PVAL, FDR=FDR, method='VIPER', permutations=NULL)
  return(RESULTS)
}


new_aREA = function(eset, genesets, M = NULL) {
  if ( ! is.null(M) ){
    hidden_genes = apply(M, 2, function(x) names(x)[x==1] )
    temp = lapply(1:length(genesets$GENES), 
                  function(i) {
                    gsw = genesets$GENES[[i]]
                    gsn = names(gsw)
                    es = lapply(colnames(eset), function(s){
                      genes = setdiff(rownames(eset), hidden_genes[[s]])
                      t2 = rank(eset[ genes ,s], na.last="keep") / (length(genes)+1)
                      t2 = qnorm(t2)
                      gsgenes = intersect(gsn, genes)
                      xst2 = t2[ gsgenes ]
                      wts = gsw[ gsgenes ]/sum(gsw[ gsgenes ]) #new
                      sum1 = matrix(wts, nrow = 1, ncol = length(gsgenes)) %*% xst2 #new
                      # sum1 = matrix(gsw[ gsgenes ], nrow = 1, ncol = length(gsgenes)) %*% xst2 
                      nes = sum1 * sqrt(sum(gsw[ gsgenes ]^2))
                      pval = pnorm(nes, mean = 0)#, sd = sqrt(sum(wts^2)))
                      pval = ifelse(pval<0.5, pval, 1-pval)*2
                      return(list(es=nes, nes=nes, pval = pval))
                    } )
                    return( es ) 
                  })
    es = t(sapply(temp, function(x) sapply(x, function(xx) xx$es)))
    nes = t(sapply(temp, function(x) sapply(x, function(xx) xx$nes)))
    pval = t(sapply(temp, function(x) sapply(x, function(xx) xx$pval)))
    es = apply(es, 2, as.numeric)
    nes = apply(nes, 2, as.numeric)
    pval = apply(pval, 2, as.numeric)
    colnames(es) = colnames(eset)
    rownames(es) = genesets$NAME
    colnames(nes) = colnames(eset)
    rownames(nes) = genesets$NAME
    colnames(pval) = colnames(eset)
    rownames(pval) = genesets$NAME
  }else{
    t2 = t(t(apply(eset, 2, rank, na.last="keep"))/(nrow(eset)+1))
    t2 = qnorm(t2)
    temp = lapply(1:length(genesets$GENES), 
                    function(i) {
                      x = genesets$GENES[[i]]
                      pos = match(names(x), rownames(t2))
                      # sum1 = matrix(x, nrow = 1, ncol = length(x)) %*% filterRowMatrix(t2, pos)
                      sum1 = matrix(x/length(x), nrow = 1, ncol = length(x)) %*% filterRowMatrix(t2, pos) #new
                      return( sum1 )
                    })
    es = t(sapply(temp, function(x) x))
    es = apply(es, 2, as.numeric)
    colnames(es) = colnames(t2)
    rownames(es) = genesets$NAME   
    w = sapply(genesets$GENES, function(x) sqrt(sum(x^2)) )
    nes = es*w# normalize enrihment score
    pval = t(sapply(genesets$NAME, function(g)
      pnorm(nes[g,], mean = 0)#, sd = sqrt(sum((genesets$GENES[[g]]/sum(genesets$GENES[[g]]))^2)))
      ))
    pval = apply(pval, 2, function(x) ifelse(x<0.5, x, 1-x)*2 )
    colnames(pval) = colnames(nes)
  }
  res =  list(es=es, nes=nes, pval = pval)
  return( res  )
}



SLEA.GSVA = function(E, genesets, M = NULL, method='GSVA'){ # NOTE: GSVA  method developed by S Hänzelmann - 2013
                                                            # NOTE: ssGSEA method developed by S Hänzelmann 2013 as proposed by D A BArbie 2009
  if( is.null(M) ){
    genesets$GENES = lapply(genesets$GENES, names)
    if (method=='GSVA'){
      ES = gsva(expr = E, genesets$GENES, method = 'gsva', verbose = F, kernel = F, parallel.sz=no_cores, parallel.type='FORK')$es.obs
      NES = ES
      PVAL = t(sapply(genesets$NAME, function(g)  pnorm(NES[g, ], mean = 0, sd = sqrt(  1/length(genesets$GENES[[g]]) ))  )) # use empirical approximation
    }
    if (method %in% c('ssgsea', 'ssGSEA')){
      ES = gsva(expr = E, genesets$GENES, method = 'ssgsea', ssgsea.norm=T, verbose = F, kernel = T, parallel.sz=no_cores, parallel.type='FORK')#, tau=0)
      NES = ES
      PVAL = t(sapply(genesets$NAME, function(g)  pnorm(ES[g, ], mean = mean(ES), sd = sqrt(  1/length(genesets$GENES[[g]]) ))  )) # use empirical approximation
    }
    PVAL = ifelse(PVAL > .5, 1-PVAL, PVAL)*2
  }else{
    ES = matrix(-1, ncol = ncol(E), nrow = length(genesets$NAME) )
    colnames(ES) = colnames(E)
    rownames(ES) = genesets$NAME
    NES = ES
    PVAL = ES
    FDR = ES
    pb = txtProgressBar(min = 0, max = ncol(ES), style = 3)
    res = lapply(colnames(E), function(s){
      setTxtProgressBar(pb, which(colnames(ES) == s))
      SLEA.myssGSVA(E, s, genesets, M, method)
      })
    names(res) = colnames(E)
    for (s in colnames(ES) ){
      ES[,s] = res[[s]]$es
      NES[,s] = res[[s]]$nes
      PVAL[,s] = res[[s]]$pval
    }
  }
  FDR = PVAL
  for (gs in genesets$NAME )
    FDR[gs,] = p.adjust(PVAL[gs,], method = 'fdr')
  RESULTS = list(ES=ES, NES=NES, PVAL=PVAL, FDR=FDR, method=method, permutations=NULL)
  return(RESULTS)
}


SLEA.myssGSVA = function(E, s, genesets, M = NULL, method = 'GSVA'){ # runs gsva per sample for all genesets. Only when M != NULL
  sE = E[,s]
  if( ! is.null(M) ){
    sM = M[,s]
    hidden_genes = names(sM)[ sM == 1 ]
    sE = sE[ ! names(sE) %in% hidden_genes ]
  }
  if( all(unlist(genesets$GENES) > 0)){
    genesets$GENES = lapply(genesets$GENES, names)
    if( ! is.null(M) ){
      genesets$GENES = lapply(genesets$GENES, setdiff, y = hidden_genes)
    }
    if( method == 'GSVA'){
      sE = cbind(sE, sE)# trick to avoid errors in gsva. Expr must have more than one column
      ssES = gsva(expr = sE, genesets$GENES, method = 'gsva', verbose = F, kernel = F, parallel.sz=no_cores, parallel.type='FORK')$es.obs
      es = ssES[,1]
      nes = es
      pval = sapply(genesets$NAME, function(g) pnorm(es[g], mean = 0, sd = sqrt(  1/length(genesets$GENES[[g]]) )) ) # use analytical approximation
    }
    if( method %in% c('ssgsea', 'ssGSEA')){
      sE = cbind(sE, sapply(1:2, function(i) sample(sE) )) # randomize the data for empiricaly computed null distribution
      ssES = gsva(expr = sE, genesets$GENES, method = 'ssgsea', ssgsea.norm=T, verbose = F, kernel = T, parallel.sz=no_cores, parallel.type='FORK', tau=0)
      es = ssES[,1]
      nes = es
      pval = sapply(genesets$NAME, function(g)  pnorm(nes[g], mean = 0, sd = sqrt(  1/length(genesets$GENES[[g]]) ))  ) # use empirical approximation
    }
  }else{
    sE = (sE - min(sE)) / (max(sE) - min(sE))
    ssES = lapply(genesets$NAME, function(gsn){
      gs = genesets$GENES[[ gsn ]]
      if( ! is.null(M) ){
        gs = gs[setdiff(names(gs), hidden_genes)]
      }
      sE[ names(sE) %in% names(gs)[gs<0] ] = 1 - sE[ names(sE) %in% names(gs)[gs<0] ]
      gs = list(names(gs))
      names(gs) = gsn
      if( method == 'GSVA'){
        sE = cbind(sE, sapply(1:2, function(i) sample(sE) )) # trick to avoid errors in gsva function
        ssES = gsva(expr = sE, gs, method = 'gsva', verbose = F, kernel = F, parallel.sz=no_cores, parallel.type='FORK')$es.obs
        es = ssES[,1]
        nes = es
        pval = pnorm(es, mean = 0, sd = sqrt(1/length(genesets$GENES[[g]]))  )# use analytical approximation
      }
      if( method == 'ssgsea'){
        sE = cbind(sE, sapply(1:100, function(i) sample(sE) )) # randomize the data for empiricaly computed null distribution
        ssES = gsva(expr = sE, gs, method = 'ssgsea', ssgsea.norm=T, verbose = F, kernel = F, parallel.sz=no_cores, parallel.type='FORK')
        es = ssES[,1]
        nes = es
        pes = ssES[,-1]
        pval = pnorm(es, mean = mean(pes), sd = sd(pes)) # use empirical approximation
      }
      return(list(es=es, nes=nes, pval = pval))
    })
    names(ssES) = genesets$NAME
    es = sapply(ssES, function(x) x$es )
    nes = sapply(ssES, function(x) x$nes  )
    pval = sapply(ssES, function(x) x$pval  )
  }
  pval = ifelse(pval > .5, 1-pval, pval)*2
  names(pval) = names(nes)
  return(list(es=es, nes=nes, pval=pval))
}





SLEA.MEAN = function(E, genesets, permutations = 1000, M = NULL){ # NOTE: method as proposed by G Gundem - 2012
  ES = matrix(-1, ncol = ncol(E), nrow = length(genesets$NAME) )
  colnames(ES) = colnames(E)
  rownames(ES) = genesets$NAME
  PVAL = ES
  NES = ES
  FDR = ES
  cl = makeCluster(no_cores)
  clusterExport(cl, "SLEA.MEAN.sample_geneset_level")
  pb = txtProgressBar(min = 0, max = length(genesets$NAME), style = 3)
  sample_geneset_level = vector('list', ncol(ES))
  for( gs in genesets$NAME ){
    genes = genesets$GENES[[gs]]
    sample_geneset_level = parLapply(cl, colnames(ES), function(s, E, genes, permutations, M){ SLEA.MEAN.sample_geneset_level(s = s, E = E, genes = genes, permutations = permutations, M = M)}, E, genes, permutations, M)
    ES[gs,] = unlist(lapply(sample_geneset_level, function(x) x$es))
    NES[gs,] = unlist(lapply(sample_geneset_level, function(x) x$nes))
    PVAL[gs,] = unlist(lapply(sample_geneset_level, function(x) x$pval))
    FDR[gs,] = p.adjust(PVAL[gs,], method = 'fdr')
    setTxtProgressBar(pb, which(genesets$NAME == gs))
  }
  RESULTS = list(ES=ES, NES=NES, PVAL=PVAL, FDR=FDR, method='mean', permutations=permutations)
  stopCluster(cl)

  return(RESULTS)
}




SLEA.MEAN.sample_geneset_level = function(s, E, genes, permutations = 10000, M = NULL){
  desc = 'median'
  sE = E[,s]
  if( any(genes < 0) ){
    sE = (sE - min(sE)) / (max(sE) - min(sE))
    sE[ names(sE) %in% names(genes)[ genes < 0 ] ] = 1 - sE[ names(sE) %in% names(genes)[ genes < 0 ] ]
  }
  genes = names(genes)
  if( ! is.null(M)){
    sM = M[,s]
    hidden_genes = names(sM)[ sM == 1 ]
    sE = sE[ ! names(sE) %in% hidden_genes ]
    genes = setdiff(genes, hidden_genes)
    if( length(genes) == 0 )
      return(list(es=NA, nes=NA, pval=NA))  
  }
  # obtain the observed Enrichment Scores as the mean/median expression of the genes in the gene set
  observed = sE[ names(sE) %in% genes ]
  o = do.call(desc, list(observed))
  # obtain the expected Enrichment Score on permuted genesets
  expected = lapply(1:permutations, function(i) sample(sE, length(genes), replace = T)) # get random sets of the same size
  e = sapply(expected, function(e) do.call(desc, list(e)) )
  # compute z-score
  e_mean = mean(e) # expected mean
  e_sd = sd(e)*sqrt((length(e)-1)/(length(e))) # expected sd
  ezscore = (e - e_mean) / e_sd # expected z-scores
  ozscore = (o - e_mean) / e_sd  # observed z-score
  # get pvalue
  f = ecdf(ezscore) # empirical cummulative distribution of e
  fq = f(ozscore)  # pvalue ~ observed percentile in ref distribution
  pval = ifelse(fq < 0.5, fq, 1 - fq)*2 # Correct fot two-tail test
  pval
  return(list(es=ozscore, nes=fq, pval=pval))  
}
# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------




# ------------------------------------------------------------------------------------------------------------------------------------
# GSVA package: modifying functions to improve performance
# ------------------------------------------------------------------------------------------------------------------------------------
rndWalk <- function(gSetIdx, geneRanking, j, R, alpha) {
  # indicatorFunInsideGeneSet <- match(geneRanking, gSetIdx)
  # indicatorFunInsideGeneSet[!is.na(indicatorFunInsideGeneSet)] <- 1
  # indicatorFunInsideGeneSet[is.na(indicatorFunInsideGeneSet)] <- 0
  indicatorFunInsideGeneSet <- (geneRanking %in% gSetIdx) + 0
  jR <- R[geneRanking, j]
  stepCDFinGeneSet <- cumsum((jR * indicatorFunInsideGeneSet)^alpha) /
    sum((jR * indicatorFunInsideGeneSet)^alpha)
  stepCDFoutGeneSet <- cumsum(!indicatorFunInsideGeneSet) /
    sum(!indicatorFunInsideGeneSet)
  walkStat <- stepCDFinGeneSet - stepCDFoutGeneSet
  sum(walkStat) 
}
assignInNamespace("rndWalk",rndWalk, ns="GSVA")
# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------





# ------------------------------------------------------------------------------------------------------------------------------------
# Getting gene expression statistics over samples. 
# This brings distinct expression profiles to a common scale to avoid gene specific effects.
# ------------------------------------------------------------------------------------------------------------------------------------
gene_expression_statistic = function(expression, method, rnaseq=FALSE){
  if( method == 'FiNorm' ){
    message('Scaling and recentering gene expression using Fi method.\nPlease note that this step may take several hours')
    recentered_m = t(apply(expression, 1, GES.Fi.KDEpdf))
  }else if( method == 'GSVAnorm' ){
    message('Scaling and recentering gene expression using GSVA method')
    recentered_m = GES.gsva.GaussianKernel(as.matrix(expression), rnaseq=rnaseq)
  }else if( method == 'scale' ){
    message('Scaling and recentering gene expression using R scale method')
    recentered_m = GES.scale(as.matrix(expression))
  }else{
    message('Invalid method. The available methods are: \'GSVAnorm\' and FiNorm' )
  }
  return(recentered_m)
}



GES.gsva.GaussianKernel  =  function(E, sample.idxs=NULL, rnaseq=FALSE){
  # As GSVA (Hanzelmann et al 2013)
  # library(GSVA)
  # a kernel non-parametric estimation of the empirical cumulative distribution function 
  samples = colnames(E)
  genes = rownames(E)
  n.test.samples  =  ncol(E)
  n.genes  =  nrow(E)
  if( is.null(sample.idxs) ){
    sample.idxs = colnames(E)
  }
  sample.idxs = intersect(colnames(E), sample.idxs)
  n.density.samples  =  length(sample.idxs)
  A = .C("matrix_density_R",
         as.double(t(E[ ,sample.idxs, drop=FALSE])),
         as.double(t(E)),
         R = double(n.test.samples * n.genes),
         n.density.samples,
         n.test.samples,
         n.genes,
         as.integer(rnaseq))$R
  
  gene.density  =  t(matrix(A, n.test.samples, n.genes))
  colnames(gene.density) = samples
  rownames(gene.density) = genes
  return (gene.density)  
}


GES.gsva.ecdf = function(E, sample.idxs=NULL, rnaseq=FALSE){
  # As GSVA (Hanzelmann et al 2013)
  # library(GSVA)
  # empirical cumulative distribution function directly estimated from the observed data
  samples = colnames(E)
  genes = rownames(E)
  n.test.samples = ncol(E)
  n.genes = nrow(E)
  if( is.null(sample.idxs) ){
    sample.idxs = colnames(E)
  }
  sample.idxs = intersect(colnames(E), sample.idxs)
  n.density.samples  =  length(sample.idxs)
  gene.density <- t(apply(E, 1, function(x, sample.idxs) {
    f <- ecdf(x[sample.idxs])
    f(x)
  }, sample.idxs))
  gene.density <- log(gene.density / (1-gene.density))
  return(gene.density)
}


GES.Fi.KDEpdf  = function(signal){
  signal[ is.na(signal) ] = median(signal, na.rm = T)# impute missing values with the median
  if( min(signal, na.rm = T) < 0 ){ #avoid negative values
    signal = signal + abs(min(signal, na.rm = T)) + 1 
  }
  PDF = density(signal, width=sd(signal)/4)
  PDF = approxfun(PDF$x, PDF$y, yleft=0, yright=0) # the area under the curve of a density function represents the probability of getting an x value between a range of x-a/x+b values
  PLL = lapply(signal, function(value) quad(PDF, xa = 0, xb = value)) # computed integral: Area before the point
  CLL = lapply(signal, function(value) quad(PDF, xa = value, xb = max(signal)+1)) # computed integral: Area after the point
  PLL = unlist(PLL)
  CLL = unlist(CLL)
  PLL[ which( PLL > 1 )]  = 1 
  CLL[ which( CLL == 0 )]  = .Machine$double.eps
  norm_signal = log(PLL/CLL)
  return(norm_signal)
}


GES.scale = function(E){
  rnames = rownames(E)
  cnames = colnames(E)
  E = t(scale(t(E), center = T, scale = T))
  rownames(E) = rnames
  colnames(E) = cnames
  return(E)
}
# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------






# ------------------------------------------------------------------------------------------------------------------------------------
# SLEA Data managment
# ------------------------------------------------------------------------------------------------------------------------------------
geneset.list2binaryMat = function(gsl, E){
  mgs = melt(lapply(gsl, names))
  names(mgs) = c('tar', 'TF')
  mgs$value = melt(gsl)$value
  bmgs = dcast(mgs, tar~TF, fill = 0)
  rownames(bmgs) = bmgs$tar
  bmgs = as.matrix(bmgs[,-1])
  genes2add = setdiff(rownames(E), rownames(bmgs))
  bmgs = rbind(bmgs, 
               matrix(0, ncol = ncol(bmgs), nrow = length(genes2add), dimnames = list(genes2add, colnames(bmgs)) ) )
  return(bmgs)
}

SLEA.clean_genesets = function(genesets, E, min_geneset_size = 3, max_shared = 10){
  if( ! is.null(max_shared) ){
    cat('Removing targets under more than', max_shared, ' TF\n')
    rem = names(which(table(unlist(sapply(genesets$GENES, names) ))>max_shared))
    g = names(which(table(unlist(sapply(genesets$GENES, names) ))<=max_shared))
    genesets$GENES = lapply(genesets$GENES, function(x) x[ intersect(names(x), g) ] )
    cat('\t', length(g), ' targets keept\n')
    cat('\t', length(rem), ' targets removed\n')
  }
  cat('Filtering genesets: removing targets not in the expression matrix\n')
  genesets$GENES = lapply(genesets$GENES, function(x) x[ intersect(names(x), rownames(E)) ] )
  cat('Removing gene sets with less than ', min_geneset_size, ' genes\n')
  n = sum(unlist(lapply(genesets$GENES, length)) < min_geneset_size )
  genesets$NAME = genesets$NAME[ sapply(genesets$GENES, length) >= min_geneset_size ]
  genesets$GENES = genesets$GENES[ sapply(genesets$GENES, length) >= min_geneset_size ]
  cat('\t', n, ' gene sets removed\n')
  cat('\t', length(genesets$GENES), ' gene sets used covering ', sum(rownames(E) %in% unlist(lapply(genesets$GENES, names)) ), ' genes in the expression matrix\n')
  genesets$GENES = genesets$GENES[ order(genesets$NAME) ]
  genesets$NAME = genesets$NAME[ order(genesets$NAME) ]
  return(genesets)
}

SLEA.clean_ME = function(M, E){
  cat('Methylatyion matrix provided\n')
  cat('\t', ncol(E), 'samples in E\n')
  cat('\t', ncol(M), 'samples in M\n')
  cat('\t', nrow(M), 'genes in M\n')
  M = M[ rownames(M) %in% rownames(E), ]
  samples = intersect(colnames(M), colnames(E))
  E = E[ , samples ]
  M = M[ , samples ]
  M = M[ order(rownames(M)), ]
  cat('\t', ncol(E), 'samples in E after cleaning\n')
  cat('\t', ncol(M), 'samples in M after cleaning\n')
  cat('\t', nrow(M), 'genes in M after cleaning\n')
  return(list(M=M, E=E))
}

# SLEA.VIPER.formatgenset = function(genesets){
#   reg = lapply(genesets$NAME, function(g) {
#     x = genesets$GENES[[g]]
#     likelihood = rep(1, length(x))
#     tfmode = likelihood
#     names(tfmode) = x
#     list("tfmode" = tfmode,
#          "likelihood"= likelihood)
#   })
#   names(reg) = genesets$NAME
#   return(reg)
# }
# 
# SLEA.VIPER.expressionSet = function(E){
#   pdata =  data.frame(randomclass = rep('B', ncol(E)) ,  samples = colnames(E), stringsAsFactors = F)
#   pdata$randomclass[ sample(ncol(E), 1) ] = 'A'
#   pdata$randomclass = factor(pdata$randomclass)
#   rownames(pdata) = colnames(E)
#   phenoData = new("AnnotatedDataFrame", pdata)
#   E = ExpressionSet(E, phenoData = phenoData)
#   return(E)
# }
# 
# ------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------
