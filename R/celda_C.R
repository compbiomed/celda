# -----------------------------------
# Variable description
# -----------------------------------
# C = Cell
# S or s = Sample
# G = Gene
# CP = Cell population
# n = counts of transcripts
# m = counts of cells
# K = Total number of cell populations
# nM = Number of cells
# nG = Number of genes
# nS = Number of samples

# -----------------------------------
# Count matrices descriptions
# -----------------------------------

# All n.* variables contain counts of transcripts
# n.CP.by.TS = Number of counts in each Cellular Population per Transcriptional State
# n.TS.by.C = Number of counts in each Transcriptional State per Cell 
# n.CP.by.G = Number of counts in each Cellular Population per Gene
# n.by.G = Number of counts per gene (i.e. rowSums)
# n.by.TS = Number of counts per Transcriptional State

## All m.* variables contain counts of cells
# m.CP.by.S = Number of cells in each Cellular Population per Sample

# nG.by.TS = Number of genes in each Transcriptional State


#' celda Cell Clustering Model
#' 
#' @param counts A numeric count matrix
#' @param sample.label A vector indicating the sample for each cell (column) in the count matrix
#' @param K An integer or range of integers indicating the desired number of cell clusters (for celda_C / celda_CG models)
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
#' @param algorithm Use 'EM' or 'Gibbs' sampling for clustering of cells into subpopulations. EM is much faster for larger datasets. Default 'EM'.
#' @param stop.iter Number of iterations without improvement in the log likelihood to stop the inference algorithm. Default 10.
#' @param max.iter Maximum iterations of inference algorithm to perform regardless of convergence. Default 200.
#' @param split.on.iter On every 'split.on.iter' iteration, a heuristic will be applied to determine if a gene/cell cluster should be reassigned and another gene/cell cluster should be split into two clusters. Default 10.
#' @param split.on.last After the the chain has converged according to 'stop.iter', a heuristic will be applied to determine if a gene/cell cluster should be reassigned and another gene/cell cluster should be split into two clusters. If a split occurs, then 'stop.iter' will be reset. Default TRUE.
#' @param random.state.order Whether to sample cells in a random order when performing Gibbs sampling. Defaults to TRUE.
#' @param count.checksum An MD5 checksum for the provided counts matrix
#' @param seed Parameter to set.seed() for random number generation
#' @param z.init Initial values of z. If NULL, z will be randomly sampled. Default NULL.
#' @param logfile If NULL, messages will be displayed as normal. If set to a file name, messages will be redirected messages to the file. Default NULL.
#' @return An object of class celda_C with clustering results and Gibbs sampling statistics
#' @export
celda_C = function(counts, sample.label=NULL, K, alpha=1, beta=1,
					         algorithm = c("EM", "Gibbs"), stop.iter = 10, max.iter=200, 
					         split.on.iter=10, split.on.last=TRUE,
                 	 random.state.order=TRUE, count.checksum=NULL, seed=12345,
                 	 z.init = NULL, logfile=NULL) {
  
  ## Error checking and variable processing
  counts = processCounts(counts)  
    
  s = processSampleLabels(sample.label, ncol(counts))
  if (is.null(sample.label)) {
    sample.label = s
  }
  
  if(is.null(count.checksum)) {
    count.checksum = digest::digest(counts, algo="md5")
  }

  algorithm <- match.arg(algorithm)
  algorithm.fun <- ifelse(algorithm == "Gibbs", "cC.calcGibbsProbZ", "cC.calcEMProbZ")
  
  ## Randomly select z and y or set z/y to supplied initial values
  z = initialize.cluster(K, ncol(counts), initial = z.init, fixed = NULL, seed=seed)
  z.best = z
  
  # Global variables for decomposeCounts
  setGlobalVariables.celda_C(counts, K, s, z)

  
  ## Calculate counts one time up front
  p = cC.decomposeCounts(counts, s, z, K)
  nS = p$nS
  nG = p$nG
  nM = p$nM
  m.CP.by.S = p$m.CP.by.S
  n.G.by.CP = p$n.G.by.CP
  n.CP = p$n.CP
  n.by.C = p$n.by.C
  
  ll = cC.calcLL(m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, s=s, K=K, nS=nS, nG=nG, alpha=alpha, beta=beta)
  
  set.seed(seed)
  logMessages(date(), "... Starting celda_C to cluster cells.", logfile=logfile, append=FALSE)
  
  iter = 1L
  num.iter.without.improvement = 0L
  do.cell.split = TRUE
  while(iter <= max.iter & num.iter.without.improvement <= stop.iter) {
    
    next.z = do.call(algorithm.fun, list(counts=counts, m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, n.by.C=n.by.C, n.CP=n.CP, z=z, s=s, K=K, nG=nG, nM=nM, alpha=alpha, beta=beta))

    m.CP.by.S = next.z$m.CP.by.S
    n.G.by.CP = next.z$n.G.by.CP
    n.CP = next.z$n.CP
    z = next.z$z

    ## Perform split on i-th iteration of no improvement in log likelihood
    if(K > 2 & (((iter == max.iter | num.iter.without.improvement == stop.iter) & isTRUE(split.on.last)) | (split.on.iter > 0 & iter %% split.on.iter == 0 & isTRUE(do.cell.split)))) {

      logMessages(date(), " ... Determining if any cell clusters should be split.", logfile=logfile, append=TRUE, sep="")
	  res = cC.splitZ(counts, m.CP.by.S, n.G.by.CP, s, z, K, nS, nG, alpha, beta, z.prob=t(next.z$probs), max.clusters.to.try=10, min.cell=3)
      logMessages(res$message, logfile=logfile, append=TRUE)

	  # Reset convergence counter if a split occured
	  if(!isTRUE(all.equal(z, res$z))) {
		num.iter.without.improvement = 0L
		do.cell.split = TRUE
	  } else {
		do.cell.split = FALSE
	  }
            
      ## Re-calculate variables
      z = res$z
      #m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
      #n.G.by.CP = colSumByGroup(counts, group=z, K=K)
      #n.CP = as.integer(colSums(n.G.by.CP))
      m.CP.by.S = res$m.CP.by.S
      n.G.by.CP = res$n.G.by.CP
      n.CP = res$n.CP
    }
    
    
    ## Calculate complete likelihood
    temp.ll = cC.calcLL(m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, s=s, K=K, nS=nS, nG=nG, alpha=alpha, beta=beta)
    if((all(temp.ll > ll)) | iter == 1) {
      z.best = z
      ll.best = temp.ll
      num.iter.without.improvement = 1L
    } else {  
      num.iter.without.improvement = num.iter.without.improvement + 1L   
    }
    ll = c(ll, temp.ll)
    
    #logMessages(date(), "... Completed iteration:", iter, "| logLik:", temp.ll, logfile=logfile, append=TRUE)
    iter = iter + 1    
  }
  
  names = list(row=rownames(counts), column=colnames(counts), sample=levels(sample.label))
  
  result = list(z=z.best, completeLogLik=ll,  
                finalLogLik=ll.best, seed=seed, K=K, 
                sample.label=sample.label, alpha=alpha, 
                beta=beta, count.checksum=count.checksum, 
                names=names)
  
  class(result) = "celda_C"
  result = reorder.celda_C(counts = counts, res = result)
  
  return(result)
}


# Gibbs sampling for the celda_C Model
cC.calcGibbsProbZ = function(counts, m.CP.by.S, n.G.by.CP, n.by.C, n.CP, z, s, K, nG, nM, alpha, beta, do.sample=TRUE, random.state.order=TRUE) {
  
  ## Set variables up front outside of loop  
  probs = matrix(NA, ncol=nM, nrow=K)
  #temp.n.G.by.CP = n.G.by.CP
  #temp.n.CP = n.CP
  
  if(isTRUE(random.state.order)) {
    ix = sample(1:nM)
  } else {
    ix = rev(1:nM)
  } 
  
  #start_time = Sys.time()
  for(i in ix) {
    ## Subtract current cell counts from matrices
    m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] - 1L
    n.G.by.CP[,z[i]] = n.G.by.CP[,z[i]] - counts[,i]
    n.CP[z[i]] = n.CP[z[i]] - n.by.C[i]
    
    g_col_sums_vector = vector(length=K)
    g2_elements_vector = vector(length=K)
    for(j in 1:K) {
      g_col_sums_vector[j] = sum(lgamma(n.G.by.CP[,j] + beta)) # Calculate lgamma for every column of n.G.by.CP
      g2_elements_vector[j] = lgamma(n.CP[j] + (nG * beta)) # Calculate lgamma for every j element in n.CP
    }
    g_sum = sum(g_col_sums_vector) # sum of n.G.by.CP gammas
    g2_sum = sum(g2_elements_vector) # sum of n.CP gammas
    
    for(j in 1:K) {
      new_col_g_sum = sum(lgamma((n.G.by.CP[,j] + counts[,i]) + beta)) # calculate new sum of COLUMN j
      new_g_sum = g_sum - g_col_sums_vector[j] + new_col_g_sum # old sum - sum of old column + sum of new column

      new_g2_elememt = lgamma((n.CP[j] + n.by.C[i]) + (nG * beta)) # calculate new sum of ELEMENT j
      new_g2_sum = g2_sum - g2_elements_vector[j] + new_g2_elememt # old sum - old element + new element

      probs[j,i] = log(m.CP.by.S[j,s[i]] + alpha) + new_g_sum - new_g2_sum
    }
    
    ## Sample next state and add back counts
    if(isTRUE(do.sample)) z[i] = sample.ll(probs[,i])
    
    m.CP.by.S[z[i],s[i]] = m.CP.by.S[z[i],s[i]] + 1L
    n.G.by.CP[,z[i]] = n.G.by.CP[,z[i]] + counts[,i]
    n.CP[z[i]] = n.CP[z[i]] + n.by.C[i]
  }
  return(list(m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, n.CP=n.CP, z=z, probs=probs))
}


cC.calcEMProbZ = function(counts, m.CP.by.S, n.G.by.CP, n.by.C, n.CP, z, s, K, nG, nM, alpha, beta, do.sample=TRUE) {

  ## Expectation given current cell population labels
  theta = log(normalizeCounts(m.CP.by.S + alpha, scale.factor=1))
  phi = log(normalizeCounts(n.G.by.CP + beta, scale.factor=1))
  
  ## Maximization to find best label for each cell
  probs = eigenMatMultInt(phi, counts) + theta[, s]  
  #probs = (t(phi) %*% counts) + theta[, s]  
  
  z.previous = z
  z = apply(probs, 2, which.max)

  ## Recalculate counts based on new label
  #p = cC.decomposeCounts(counts, s, z, K)
  p = cC.reDecomposeCounts(counts, s, z, z.previous, n.G.by.CP, K)
  m.CP.by.S = p$m.CP.by.S
  n.G.by.CP = p$n.G.by.CP
  n.CP = p$n.CP

  return(list(m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, n.CP=n.CP, z=z, probs=probs))
}

#' Simulate cells from the cell clustering generative model
#' 
#' @param model Celda model to use for simulation. One of 'available_models'. 
#' @param S Total number of samples
#' @param C.Range Vector of length 2 given the range (min,max) of number of cells for each sample to be randomly generated from the uniform distribution
#' @param N.Range Vector of length 2 given the range (min,max) of number of counts for each cell to be randomly generated from the uniform distribution
#' @param G Total number of Genes to be simulated
#' @param K An integer or range of integers indicating the desired number of cell clusters (for celda_C / celda_CG models)
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
#' @param seed starting point used for generating simulated data
#' @param ... Other arguments
#' @export
simulateCells.celda_C = function(model, S=10, C.Range=c(10, 100), N.Range=c(100,5000), 
                                 G=500, K=5, alpha=1, beta=1, seed=12345, ...) {
  
  cC.global_simCellsFlag <<- TRUE
  set.seed(seed) 
  
  phi <- rdirichlet(K, rep(beta, G))
  theta <- rdirichlet(S, rep(alpha, K))
  
  ## Select the number of cells per sample
  nC <- sample(C.Range[1]:C.Range[2], size=S, replace=TRUE)  
  cell.sample.label <- rep(1:S, nC)
  
  ## Select state of the cells  
  z <- unlist(lapply(1:S, function(i) sample(1:K, size=nC[i], prob=theta[i,], replace=TRUE)))
  
  ## Select number of transcripts per cell
  nN <- sample(N.Range[1]:N.Range[2], size=length(cell.sample.label), replace=TRUE)
  
  ## Select transcript distribution for each cell
  cell.counts <- sapply(1:length(cell.sample.label), function(i) stats::rmultinom(1, size=nN[i], prob=phi[z[i],]))
  
  rownames(cell.counts) = paste0("Gene_", 1:nrow(cell.counts))
  colnames(cell.counts) = paste0("Cell_", 1:ncol(cell.counts)) 
  cell.sample.label = paste0("Sample_", 1:S)[cell.sample.label]
  
  ## Peform reordering on final Z and Y assigments:
  names = list(row=rownames(cell.counts), column=colnames(cell.counts), 
               sample=unique(cell.sample.label))
  result = list(z=z, completeLogLik=NULL, 
                finalLogLik=NULL, K=K, 
                alpha=alpha, beta=beta, seed=seed, 
                sample.label=cell.sample.label, names=names,
                count.checksum=NULL)
  class(result) = "celda_C" 
  result = reorder.celda_C(counts = cell.counts, res = result)
  
  cC.global_simCellsFlag <<-FALSE
  return(list(z=result$z, counts=processCounts(cell.counts), sample.label=cell.sample.label, K=K, alpha=alpha, beta=beta, C.Range=C.Range, N.Range=N.Range, S=S))
}

#' Generate factorized matrices showing each feature's influence on the celda_C model clustering 
#' 
#' @param counts A numeric count matrix
#' @param celda.mod Object return from celda_C function
#' @param type A character vector containing one or more of "counts", "proportions", or "posterior". "counts" returns the raw number of counts for each entry in each matrix. "proportions" returns the counts matrix where each vector is normalized to a probability distribution. "posterior" returns the posterior estimates which include the addition of the Dirichlet concentration parameter (essentially as a pseudocount).
#' @export
factorizeMatrix.celda_C = function(counts, celda.mod, type=c("counts", "proportion", "posterior")) {

  K = celda.mod$K
  z = celda.mod$z
  alpha = celda.mod$alpha
  beta = celda.mod$beta
  sample.label = celda.mod$sample.label
  s = processSampleLabels(sample.label, ncol(counts))

  p = cC.decomposeCounts(counts, s, z, K)
  m.CP.by.S = p$m.CP.by.S
  n.G.by.CP = p$n.G.by.CP
  
  K.names = paste0("K", 1:K)
  rownames(n.G.by.CP) = celda.mod$names$row
  colnames(n.G.by.CP) = K.names
  rownames(m.CP.by.S) = K.names
  colnames(m.CP.by.S) = celda.mod$names$sample
  
  counts.list = c()
  prop.list = c()
  post.list = c()
  res = list()
  
  if(any("counts" %in% type)) {
    counts.list = list(sample.states=m.CP.by.S, gene.states=n.G.by.CP)
    res = c(res, list(counts=counts.list))
  }
  if(any("proportion" %in% type)) {
    ## Need to avoid normalizing cell/gene states with zero cells/genes
    unique.z = sort(unique(z))
    temp.n.G.by.CP = n.G.by.CP
    temp.n.G.by.CP[,unique.z] = normalizeCounts(temp.n.G.by.CP[,unique.z], scale.factor=1)
    
    prop.list = list(sample.states = normalizeCounts(m.CP.by.S, scale.factor=1),
                     gene.states = temp.n.G.by.CP)
    res = c(res, list(proportions=prop.list))
  }
  if(any("posterior" %in% type)) {
    post.list = list(sample.states = normalizeCounts(m.CP.by.S + alpha, scale.factor=1),
                     gene.states = normalizeCounts(n.G.by.CP + beta, scale.factor=1))
    res = c(res, posterior = list(post.list))                           
  }
  
  return(res)
}


# Calculate log-likelihood for celda_C model
cC.calcLL = function(m.CP.by.S, n.G.by.CP, s, z, K, nS, nG, alpha, beta) {
  
  ## Calculate for "Theta" component
  a = nS * lgamma(K * alpha)
  b = sum(lgamma(m.CP.by.S + alpha))
  c = -nS * K * lgamma(alpha)
  d = -sum(lgamma(colSums(m.CP.by.S + alpha)))
  
  theta.ll = a + b + c + d
  
  ## Calculate for "Phi" component
  a = K * lgamma(nG * beta)
  b = sum(lgamma(n.G.by.CP + beta))
  c = -K * nG * lgamma(beta)
  d = -sum(lgamma(colSums(n.G.by.CP + beta)))
  
  phi.ll = a + b + c + d
  
  final = theta.ll + phi.ll
  return(final)
}


#' Calculate the celda_C log likelihood for user-provided cluster assignments
#' 
#' @param counts A numeric count matrix
#' @param sample.label A vector indicating the sample label for each cell (column) in the count matrix
#' @param z A numeric vector of cluster assignments
#' @param K The total number of clusters in z
#' @param alpha Non-zero concentration parameter for sample Dirichlet distribution
#' @param beta Non-zero concentration parameter for gene Dirichlet distribution
#' @param ... Additional parameters
#' @export
calculateLoglikFromVariables.celda_C = function(counts, sample.label, z, K, alpha, beta) {
  s = processSampleLabels(sample.label, ncol(counts))
  p = cC.decomposeCounts(counts, s, z, K)  
  final = cC.calcLL(m.CP.by.S=p$m.CP.by.S, n.G.by.CP=p$n.G.by.CP, s=s, z=z, K=K, nS=p$nS, nG=p$nG, alpha=alpha, beta=beta)
  return(final)
}

setGlobalVariables.celda_C = function(counts, K, s, z){
  cC.global_previousZ <<- integer(length(z)) # vector of 0s
  cC.global_previousS <<- 0
  cC.global_zChanged <<- TRUE
  cC.global_sChanged <<- TRUE
  cC.global_nS <<- 0
  cC.global_nG <<- 0
  cC.global_nM <<- 0
  cC.global_m.CP.by.S <<- matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=length(unique(s)))
  cC.global_n.G.by.CP <<- colSumByGroup(counts, group=z, K=K)
  cC.global_n.CP <<- 0
  cC.global_n.by.C <<- 0
  cC.globalFlag <<- FALSE
  cC.global_variables_set <<- TRUE
}

#' Takes raw counts matrix and converts it to a series of matrices needed for log likelihood calculation
#' @param counts A numeric count matrix
#' @param s An integer vector indicating the sample label for each cell (column) in the count matrix
#' @param z A numeric vector of cluster assignments
#' @param K The total number of clusters in z
cC.decomposeCounts = function(counts, s, z, K) {
  if (!exists('cC.global_variables_set')){
    setGlobalVariables.celda_C(counts, K, s, z)
  }
  
  cC.global_zChanged <<- if(identical(cC.global_previousZ, z)) FALSE else TRUE
  cC.global_previousZ <<- z
  cC.global_sChanged <<- if(identical(cC.global_previousS, s)) FALSE else TRUE
  cC.global_previousS <<- s
  
  if(!cC.globalFlag){
    cC.global_n.by.C <<- as.integer(colSums(counts))
    cC.global_nG <<- nrow(counts)
    cC.global_nM <<- ncol(counts)
    cC.globalFlag = TRUE
  }
  n.by.C = cC.global_n.by.C
  nG = cC.global_nG 
  nM = cC.global_nM 
  
  if(cC.global_sChanged){
    cC.global_nS <<- length(unique(s))
  }
  nS = cC.global_nS
  
  if(cC.global_zChanged || cC.global_sChanged){
    cC.global_m.CP.by.S <<- matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
  }
  m.CP.by.S = cC.global_m.CP.by.S
  
  if(cC.global_zChanged){
    cC.global_n.G.by.CP = colSumByGroup(counts, group=z, K=K)
    cC.global_n.CP = as.integer(colSums(cC.global_n.G.by.CP))
  }
  n.G.by.CP = cC.global_n.G.by.CP
  n.CP = cC.global_n.CP
  
  return(list(m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, n.CP=n.CP, n.by.C=n.by.C, nS=nS, nG=nG, nM=nM))
}

cC.reDecomposeCounts = function(counts, s, z, previous.z, n.G.by.CP, K) {

  ## Recalculate counts based on new label
  n.G.by.CP = colSumByGroupChange(counts, n.G.by.CP, z, previous.z, K)
  nS = length(unique(s))
  m.CP.by.S = matrix(as.integer(table(factor(z, levels=1:K), s)), ncol=nS)
  n.CP = as.integer(colSums(n.G.by.CP))

  return(list(m.CP.by.S=m.CP.by.S, n.G.by.CP=n.G.by.CP, n.CP=n.CP))  
}


#' Calculates the conditional probability of each cell belong to each cluster given all other cluster assignments
#'
#' @param celda.mod A model returned from the 'celda_C' function
#' @param counts The original count matrix used in the model
#' @param log If FALSE, then the normalized conditional probabilities will be returned. If TRUE, then the unnormalized log probabilities will be returned.  
#' @param ... Other arguments
#' @return A list containging a matrix for the conditional cell cluster probabilities. 
#' @export
clusterProbability.celda_C = function(celda.mod, counts, log=FALSE, ...) {
  
  z = celda.mod$z
  sample.label = celda.mod$sample.label
  s = processSampleLabels(sample.label, ncol(counts))
  
  K = celda.mod$K
  alpha = celda.mod$alpha
  beta = celda.mod$beta
  
  p = cC.decomposeCounts(counts, s, z, K)  
  
  next.z = cC.calcGibbsProbZ(counts=counts, m.CP.by.S=p$m.CP.by.S, n.G.by.CP=p$n.G.by.CP, n.by.C=p$n.by.C, n.CP=p$n.CP, z=z, s=s, K=K, nG=p$nG, nM=p$nM, alpha=alpha, beta=beta, do.sample=FALSE)  
  z.prob = t(next.z$probs)
  
  if(!isTRUE(log)) {
    z.prob = normalizeLogProbs(z.prob)
  }
  
  return(list(z.probability=z.prob))
}


#' @export
calculatePerplexity.celda_C = function(counts, celda.mod) {
  
  factorized = factorizeMatrix(counts = counts, celda.mod = celda.mod, "posterior")
  theta = log(factorized$posterior$sample.states)
  phi = log(factorized$posterior$gene.states)
  sl = celda.mod$sample.label
  
  inner.log.prob = (t(phi) %*% counts) + theta[, sl]  
  log.px = sum(apply(inner.log.prob, 2, matrixStats::logSumExp))
  
  perplexity = exp(-(log.px/sum(counts)))
  return(perplexity)
}  


reorder.celda_C = function(counts, res){
  if(res$K > 2 & isTRUE(length(unique(res$z)) > 1)) {
    res$z = as.integer(as.factor(res$z))
    fm <- factorizeMatrix(counts = counts, celda.mod = res)
    unique.z = sort(unique(res$z))
    d <- cosineDist(fm$posterior$gene.states[,unique.z])
    h <- hclust(d, method = "complete")
    res <- recodeClusterZ(res, from = h$order, to = 1:length(h$order))
  }  
  return(res)
}


#' finalClusterAssignment for celda Cell clustering funciton 
#' @param celda.mod A celda model object of class "celda_C"
#' @export
finalClusterAssignment.celda_C = function(celda.mod) {
  return(celda.mod$z)
}


#' getK for celda Cell clustering function 
#' @param celda.mod A celda model object of class "celda_C"
#' @export
getK.celda_C = function(celda.mod) {
  return(celda.mod$K)
}


#' getL for celda Cell clustering function 
#' @param celda.mod A celda model object of class "celda_C"
#' @export
getL.celda_C = function(celda.mod) { return(NA) }


#' celdaHeatmap for celda Cell clustering function 
#' @param celda.mod A celda model object of class "celda_C"
#' @param counts A numeric count matrix
#' @param ... extra parameters passed onto the renderCeldaHeatmap
#' @export
celdaHeatmap.celda_C = function(celda.mod, counts, ...) {
  renderCeldaHeatmap(counts, z=celda.mod$z, ...)
}
