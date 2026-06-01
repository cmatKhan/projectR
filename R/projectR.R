#' @importFrom stats hclust kmeans prcomp
setOldClass("kmeans")
setOldClass("hclust")
setOldClass("prcomp")

#######################################################################################################################################
#' @import limma
#' @importFrom stats model.matrix
#' @param NP vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @param model Optional arguements to choose method for projection
#' @param bootstrapPval logical to indicate whether to generate p-values using bootstrap, not available for prcomp and rotatoR objects
#' @param bootIter number of bootstrap iterations, default = 1000
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",loadings="matrix"),function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames = NULL, # a vector with names of data rows
  loadingsNames = NULL, # a vector with names of loadings rows
  NP=NA, # vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  bootstrapPval=FALSE, # logical to indicate whether to generate p-values using bootstrap
  bootIter=1e3 # No of bootstrap iterations
  ){

  ifelse(!is.na(NP),loadings<-loadings[,NP],loadings<-loadings)

  #match genes in data sets
  if(is.null(dataNames)){
    dataNames <- rownames(data)
  }
  if(is.null(loadingsNames)){
    loadingsNames <- rownames(loadings)
  }
  dataM<-geneMatchR(data1=data, data2=loadings, data1Names=dataNames, data2Names=loadingsNames, merge=FALSE)
  print(paste(as.character(dim(dataM[[2]])[1]),'row names matched between data and loadings'))
  print(paste('Updated dimension of data:',as.character(paste(dim(dataM[[2]]), collapse = ' '))))
  # do projection
  Design <- model.matrix(~0 + dataM[[1]])
  colnames(Design) <- colnames(dataM[[1]])
  projection <- lmFit(as.matrix(t(dataM[[2]])),Design)
  projectionPatterns <- t(projection$coefficients)
  projection.ts<-t(projection$coefficients/projection$stdev.unscaled/projection$sigma)

  #For limma
  pval.matrix<-2*pnorm(-abs(projection.ts))

  if(bootstrapPval){
  boots <- lapply(1:bootIter,function(x){
  rows <- sample(nrow(Design),nrow(Design),replace = T)
  projection <- lmFit(as.matrix(t(dataM[[2]][rows,])),Design[rows,])
  return(projection$coefficients)
    })
  bootPval <- compareBoots(projection$coefficients,boots)
  }

  if(full & bootstrapPval){
      projectionFit <- list('projection'=projectionPatterns, 'pval'=pval.matrix, 'bootstrapPval' = bootPval)
      return(projectionFit)
  } else if(full){
      projectionFit <- list('projection'=projectionPatterns, 'pval'=pval.matrix)
      return(projectionFit)
  }
  else{return(projectionPatterns)}
})

#######################################################################################################################################
#' @import MatrixModels
#' @importFrom stats model.matrix
#' @param NP vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
#' @param full logical indicating whether to return the full model solution. By default only the new pattern object is returned.
#' @param model Optional arguements to choose method for projection
#' @param chopBy number of columns to chop the data into (chopping helps runnning large datasets)
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="dgCMatrix",loadings="matrix"),function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames = NULL, # a vector with names of data rows
  loadingsNames = NULL, # a vector with names of loadings rows
  NP=NULL, # vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  chopBy=1000 # number of columns to chop the data into
  ){

  if(!is.null(NP)) {
    loadings<-loadings[,NP]
  }

  print("dgCMatrix detected, projecting in chunks.")
  #columns of dgcMatrix are LHS for stats::lm, and columns of loadings are the
  #dense RHS (predictors). sometimes dgcMatrix is too big to fit RAM, so we
  #just fit chunks of lm models as supported by stats::lm/limma::lmFit

  chop <- function(sparsematrix) {
    coln <- ncol(sparsematrix)
    bins <- seq(1, coln, by = chopBy)
    lapply(seq_along(bins), function(i) {
      start <- bins[i]
      end <- ifelse(i < length(bins), bins[i + 1] - 1, coln)
      return(start:end)
    })
  }
  #discard print statements projectR generates each time a chunk is called
  w <- invisible(capture.output(
    projectionList <- lapply(chop(data), function(i) {
      projectR(as.matrix(data[,i]), loadings, full=full)
    })
  ))
  #since chopping by columns, it's enough to print matching rows only once
  print(w[1])

  if(full==TRUE) {
      if(length(projectionList)==1) {#if only one chunk all OK
        res <- projectionList[[1]]
      } else {#if multiple chunks - gather pvalues and projections
        pvalues <- do.call(cbind,
          lapply(projectionList, function(x) x[["pval"]]))
        projections <- do.call(cbind,
          lapply(projectionList, function(x) x[["projection"]]))
        res <- list(projection=projections, pval=pvalues)
      }
  } else {
    res <- do.call(cbind, projectionList)
  }
  return(res)
})


#######################################################################################################################################
#' @import limma
#' @import SingleCellExperiment
#' @importFrom NMF fcnnls
#' @examples
#' library("CoGAPS")
#' # CR.RNAseq6l3c3t <- CoGAPS(p.RNAseq6l3c3t, params = new("CogapsParams", nPatterns=5))
#' projectR(data=p.ESepiGen4c1l$mRNA.Seq,loadings=CR.RNAseq6l3c3t,
#' dataNames = map.ESepiGen4c1l[["GeneSymbols"]])
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",loadings="LinearEmbeddingMatrix"),function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames = NULL, # a vector with names of data rows
  loadingsNames = NULL, # a vector with names of loadings rows
  NP=NA, # vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
  full=FALSE, # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  model=NA, # optional arguements to choose method for projection
  bootstrapPval=FALSE, # logical to indicate whether to generate p-values using bootstrap
  bootIter=1e3 # No of bootstrap iterations
  ){

  loadings<-loadings@featureLoadings
  ifelse(!is.na(NP),loadings<-loadings[,NP],loadings<-loadings)
  return(projectR(data,loadings = loadings,dataNames = dataNames, loadingsNames = loadingsNames,NP,full,bootstrapPval=bootstrapPval,bootIter=bootIter))

})

#######################################################################################################################################

#' @import limma
#' @importFrom stats var
#' @param center_by_loadings Logical, default \code{FALSE}. If \code{TRUE},
#'   centers the new data using the gene means stored in \code{loadings$center}
#'   (computed from the original training data) before projecting onto the
#'   loading vectors. If \code{FALSE}, each gene is centered by its mean in
#'   the new data prior to projection. Only valid when \code{loadings$center}
#'   is not \code{FALSE}, i.e. when \code{prcomp} was called with
#'   \code{center = TRUE}.
#' @examples
#' pca.RNAseq6l3c3t<-prcomp(t(p.RNAseq6l3c3t))
#' pca.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq,
#' loadings=pca.RNAseq6l3c3t, dataNames = map.ESepiGen4c1l[["GeneSymbols"]])
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",loadings="prcomp"),function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames = NULL, # a vector with names of data rows
  loadingsNames = NULL, # a vector with names of loadings rows
  NP=NA, # vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
  center_by_loadings = FALSE, # If true, use loadings$center to center the `data` prior to projecting onto loadings
  full=FALSE # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ){

  # if `center_by_loadings` is TRUE, verify that the slot exists and is
  # not FALSE (which is true if the `loadings` are calculated with
  # center=FALSE)
  # In addition, require that the names exist. It should as long as the
  # rownames of the loadings matrix exists
  if (center_by_loadings) {
    loadings_center <- loadings$center

    if (is.null(loadings_center) || identical(loadings_center, FALSE)) {
      stop(
        "center_by_loadings = TRUE but loadings$center is missing or FALSE. ",
        "Rerun prcomp() with center = TRUE, or set center_by_loadings = FALSE (default) ",
        "to center each new sample by its own gene means."
      )
    }

    if (is.null(names(loadings_center))) {
        stop(
            "loadings$center is unnamed. Cannot gene-match the centering vector. ",
            "Ensure prcomp() was called with a matrix whose columns are named by gene."
        )
    }
  }

  loadings<-loadings$rotation
  ifelse(!is.na(NP),loadings<-loadings[,NP],loadings<-loadings)

  #match genes in data sets
  if(is.null(dataNames)){
    dataNames <- rownames(data)
  }
  if(is.null(loadingsNames)){
    loadingsNames <- rownames(loadings)
  }
  dataM<-geneMatchR(data1=data, data2=loadings, data1Names=dataNames, data2Names=loadingsNames, merge=FALSE)
  print(paste(as.character(dim(dataM[[2]])[1]),'row names matched between data and loadings'))
  print(paste('Updated dimension of data:',as.character(paste(dim(dataM[[2]]), collapse = ' '))))

  # do projection
  dat2P <- if (center_by_loadings) {
      # enforce geneMatchR ordering, and validate that there are no
      # missing genes. Shouldn't be possible, might as well check anyway
      loadings_center <- loadings_center[rownames(dataM[[2]])]
      if (any(is.na(loadings_center))) {
          missing <- rownames(dataM[[2]])[is.na(loadings_center)]
          stop(
              length(missing), " gene(s) in the matched data are absent from ",
              "loadings$center. This should not happen if loadings is a valid ",
              "prcomp object"
          )
      }
      # transpose to samples x genes
      t(sweep(dataM[[2]], 1, loadings_center, "-"))
  } else {
      # result is samples x genes
      apply(dataM[[2]], 1, function(x) x - mean(x))
  }
  projectionPatterns<- dat2P %*% dataM[[1]] #head(X %*% PCA$rotation)

  if(full==TRUE){
  #calculate percent varience accounted for by each PC in newdata
  #Eigenvalues<-eigen(cov(projectionPatterns))$values
  #PercentVariance<-round(Eigenvalues/sum(Eigenvalues) * 100, digits = 2)

  PercentVariance<-apply(projectionPatterns,2, function(x) 100*var(x)/sum(apply(projectionPatterns,2,var)))

    projectionFit <- list(t(projectionPatterns), PercentVariance) #also need to change this to transpose
    return(projectionFit)
  }
  else{return(t(projectionPatterns))}

})
#######################################################################################################################################

#' @examples
#' pca.RNAseq6l3c3t<-prcomp(t(p.RNAseq6l3c3t))
#' r.RNAseq6l3c3t<-rotatoR(1,1,-1,-1,pca.RNAseq6l3c3t$rotation[,1:2])
#' pca.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq,
#' loadings=r.RNAseq6l3c3t, dataNames = map.ESepiGen4c1l[["GeneSymbols"]])
#'
#' @rdname projectR-methods
#' @aliases projectR

setMethod("projectR",signature(data="matrix",loadings="rotatoR"),function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames = NULL, # a vector with names of data rows
  loadingsNames = NULL, # a vector with names of loadings rows
  NP=NA, # vector of integers indicating which columns of loadings object to use. The default of NP=NA will use entire matrix.
  full=FALSE # logical indicating whether to return the full model solution. By default only the new pattern object is returned.
  ){

  loadings <- loadings@rotatedM
  ifelse(!is.na(NP),loadings<-loadings[,NP],loadings<-loadings)

  #match genes in data sets
  if(is.null(dataNames)){
    dataNames <- rownames(data)
  }
  if(is.null(loadingsNames)){
    loadingsNames <- rownames(loadings)
  }
  dataM<-geneMatchR(data1=data, data2=loadings, data1Names=dataNames, data2Names=loadingsNames, merge=FALSE)
  print(paste(as.character(dim(dataM[[2]])[1]),'row names matched between data and loadings'))
  print(paste('Updated dimension of data:',as.character(paste(dim(dataM[[2]]), collapse = ' '))))
  # do projection
  dat2P<-apply(dataM[[2]],1,function(x) x-mean(x))
  projectionPatterns<- dat2P %*% dataM[[1]] #head(X %*% PCA$rotation)

  if(full==TRUE){
  #calculate percent varience accounted for by each PC in newdata
  Eigenvalues<-eigen(cov(projectionPatterns))$values
  PercentVariance<-round(Eigenvalues/sum(Eigenvalues) * 100, digits = 2)

  projectionFit <- list(t(projectionPatterns), PercentVariance)
  return(projectionFit)
  }
  else{return(t(projectionPatterns))}

})

#######################################################################################################################################

#' @import limma
#' @examples
#' c.RNAseq6l3c3t<-correlateR(genes="T", dat=p.RNAseq6l3c3t, threshtype="N",
#' threshold=10, absR=TRUE)
#' cor.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq, loadings=c.RNAseq6l3c3t,
#' NP="PositiveCOR", dataNames = map.ESepiGen4c1l[["GeneSymbols"]])
#'
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR",signature(data="matrix",loadings="correlateR"),function(
  data, # a dataset to be projected onto
  loadings, # a matrix of continous values to be projected with unique rownames
  dataNames = NULL, # a vector with names of data rows
  loadingsNames = NULL, # a vector with names of loadings rows
  NP=NA, #can be used to select for "NegativeCOR" or "PositiveCOR" list from correlateR class obj containing both. By default is NA
  full=FALSE, # logical indicating whether to return the percent variance accounted for by each projected PC. By default only the new pattern object is returned.
  bootstrapPval=FALSE, # logical to indicate whether to generate p-values using bootstrap
  bootIter=1e3 # No of bootstrap iterations
  ){

  patterns <- loadings@corM
  if(!is.na(NP)){
    patterns<-as.matrix(patterns[[NP]])
    colnames(patterns) <- NP
  }
  else {
  patterns<-loadings
}
  #check length of patterns "PositiveCOR" and "NegativeCOR" or just positive
  if(length(patterns)==2){
    patterns <- do.call(rbind,patterns)
  }
  else{
    patterns <- as.matrix(patterns)
}
  return(projectR(data = data, loadings = patterns,dataNames = dataNames, loadingsNames = loadingsNames,  full = full,
    bootstrapPval = bootstrapPval, bootIter = bootIter))

})

#######################################################################################################################################

#' @param targetNumPatterns desired number of patterns with hclust
#' @param sourceData data used to create cluster object
#' @import limma
#' @import cluster
#' @importFrom stats cutree
#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR", signature(data="matrix", loadings="hclust"),
function(data, loadings, dataNames=NULL, loadingsNames=NULL, full=FALSE,
targetNumPatterns, sourceData,bootstrapPval=FALSE,bootIter=1000)
{
  cut <- cutree(loadings, k=targetNumPatterns)
  patterns <- matrix(0, nrow=nrow(sourceData), ncol=targetNumPatterns)
  rownames(patterns) <- rownames(sourceData)
    for(x in 1:targetNumPatterns)
    {
      patterns[cut==x,x] <- apply(sourceData[cut==x,], 1, cor, y=colMeans(sourceData[cut==x,]))
    }
  return(projectR(data, loadings=patterns, dataNames, loadingsNames, full = full))
})

#' @rdname projectR-methods
#' @aliases projectR
setMethod("projectR", signature(data="matrix", loadings="kmeans"),
function(data, loadings, dataNames=NULL, loadingsNames=NULL, full=FALSE, sourceData,bootstrapPval=FALSE,bootIter=1000)
{
  patterns <- matrix(0, nrow=nrow(sourceData), ncol=length(loadings$size))
  rownames(patterns) <- rownames(sourceData)
  for(x in 1:length(loadings$size))
  {
    patterns[loadings$cluster==x,x] <- apply(sourceData[loadings$cluster==x,], 1, cor, y=colMeans(sourceData[loadings$cluster==x,]))
  }
  return(projectR(data, loadings=patterns, dataNames= dataNames, full = full,bootstrapPval=bootstrapPval,bootIter=bootIter))
})

#########################################################################

#' @examples
#' library("projectR")
#' data(p.RNAseq6l3c3t)
#' nP<-3
#' kClust<-kmeans(t(p.RNAseq6l3c3t),centers=nP)
#' kpattern<-cluster2pattern(clusters = kClust, NP = nP, data = p.RNAseq6l3c3t)
#' p<-as.matrix(p.RNAseq6l3c3t)
#' projectR(p,kpattern)
#'
#' @rdname projectR-methods
#' @aliases projectR

setMethod("projectR", signature(data="matrix", loadings="cluster2pattern"),
function(data, loadings, dataNames=NULL, loadingsNames=NULL, full=FALSE, sourceData,bootstrapPval=FALSE,bootIter=1000)
{
  loadings = loadings@clusterMatrix
  # NA results from cor when sd is zero in some of the groups
  loadings[is.na(loadings)] <- 0
  return(projectR(data, loadings=loadings, dataNames= dataNames, full = full,bootstrapPval=bootstrapPval,bootIter=bootIter))
})

#########################################################################

compareBoots <- function(projection,boots){
mat <- sapply(1:nrow(projection),function(i){
  sapply(1:ncol(projection),function(j){
    val <- sapply(1:length(boots),function(x){
      return(boots[[x]][i,j])
    })
    valD <- ecdf(val)
    qt0 <- valD(0)
    if(qt0 < 0.5){
        return(2*qt0)
      } else {
        return(2*(1-qt0))
      }
  })
})
return(mat)
}
