

####################################################################
## fixOutliers function
####################################################################

fixOutliers <- function(mat, idx, gene.ids) {

  if (class(mat) != "matrix") {
    stop("'mat' needs to be an 'matrix'")
  }

  ##  browser()
  ## Find gene names for outliers
  gene.out <- gene.ids[idx]
  ## Column index (could be multiple)
  col.out <- apply(mat[idx, ], 1, function(x) which(x > mean(x) + 2.665*sd(x)))

  ## pull out values corresponding to that gene, for that sample
  ## omit the outlying value in the average 
  mean.vals <- vector("list", length(gene.out))
  for (i in 1:length(gene.out)) {
    mvals <- numeric(length(col.out[[i]]))
    for (j in seq_along(col.out[[i]])) {
      values <- mat[which(gene.ids==gene.out[i]), col.out[[i]][j] ]
      ## here can set outlying values to missing and then take mean ... 
      mvals[j] <- mean(values[-which.max(values)])  ## ok, assumes 1 outlier among 4 reps
    }
    mean.vals[[i]] <- mvals
  }

  idxs <- cbind(rep(idx, sapply(col.out, length)), unlist(col.out))
  mat[idxs] <- unlist(mean.vals)
  return(mat)
}



####################################################################
## fixMVs function
####################################################################\

fixMVs <- function(mat, idx, gene.ids) {

  if (class(mat) != "matrix") {
    stop("'mat' needs to be an 'matrix'")
  }
  
  ## Find gene names for MVs
  gene.out <- gene.ids[idx]
  ## Column index (could be multiple)
  col.out <- apply(mat[idx, ], 1, function(x) which(is.na(x)))
  ## pull out values corresponding to that gene, for that sample
  ## omit the MV in the average 
  mean.vals <- vector("list", length(gene.out))
  for (i in 1:length(gene.out)) {
    mvals <- numeric(length(col.out[[i]]))
    for (j in seq_along(col.out[[i]])) {
      values <- mat[which(gene.ids==gene.out[i]), col.out[[i]][j]]
      mvals[j] <- mean(values[-which(is.na(values))])
    }
    mean.vals[[i]] <- mvals
  }

  idxs <- cbind(rep(idx, sapply(col.out, length)), unlist(col.out))
  mat[idxs] <- unlist(mean.vals)
  return(mat)
}



####################################################################
## clustPlot function
####################################################################\

clustPlot <- function(cl, mat, nrow, ncol) {
  par(mfrow=c(nrow, ncol))
  for(i in 1:length(table(cl))){
    int <- mat[cl==i,]
    if(is.vector(int)) {
      plot(1:ncol(mat), int, ylim=c(min(int), max(int)+.5), type="l",
           ylab="Log-Ratios",
           main=paste("Cluster",i),
           col=ifelse(is.vector(int),"red","grey"), bty="n", xaxt="n", xlab="")
    } else {
      plot(1:ncol(mat), int[1,], ylim=c(min(int),max(int)+.5),type="l",
           ylab="Log-Ratios",
           main=paste("Cluster",i),
           col="grey", bty="n", xaxt="n", xlab="")
    }
    if (is.null(colnames(mat))) {
      labs <- 1:ncol(mat)
    } else {
      labs <- colnames(mat)
    }
    axis(1, 1:ncol(mat), labels=labs)
    if(!is.vector(int))
      for(j in 2:nrow(int)){	
        lines(1:ncol(mat), int[j,], col="grey")
        lines(1:ncol(mat), apply(int,2,mean), col="red")
      }
    legend("topright", paste(ifelse(is.vector(int),1,nrow(int)),
                             ifelse(is.vector(int),"miRNA","miRNAs")),
           text.col="blue", bty="n")
  }
}


####################################################################
## imputeKNN function
####################################################################

imputeKNN <- function (data, k = 10, distance = "euclidean",
                       rm.na = TRUE, rm.nan = TRUE, rm.inf = TRUE) {

  if (!(is.matrix(data))) {
    stop(message = paste(deparse(substitute(data)),
           " is not a matrix.", sep = ""))
  }

  distance <- match.arg(distance, c("euclidean","correlation"))
  
  nr <- dim(data)[1]
  if (k < 1 | k > nr) {
    stop(message = "k should be between 1 and the number of rows")
  }

  
  if (distance=="correlation"){
    genemeans<-rowMeans(data,na.rm=TRUE)
    genesd<-sd(t(data),na.rm=TRUE)
    data<-(data-genemeans)/genesd
  }
  
  imp.knn <- data
  imp.knn[is.finite(data) == FALSE] <- NA
  t.data<-t(data)
  
  mv.ind <- which(is.na(imp.knn), arr.ind = TRUE)
  arrays <- unique(mv.ind[, 2])
  array.ind <- match(arrays, mv.ind[, 2])
  ngenes <- 1:nr
  for (i in 1:length(arrays)) {
    set <- array.ind[i]:min((array.ind[(i + 1)] - 1), dim(mv.ind)[1],
                            na.rm = TRUE)
    cand.genes <- ngenes[-unique(mv.ind[set, 1])]
    cand.vectors <- t.data[,cand.genes]
    exp.num<- arrays[i]
    for (j in set) {
      gene.num <- mv.ind[j, 1]
      tar.vector <- data[gene.num,]
      if(distance=="correlation")
        r <- cor(cand.vectors,tar.vector, use = "pairwise.complete.obs")
      dist <- switch(distance,
                     euclidean = sqrt(colMeans((tar.vector-cand.vectors)^2, na.rm = TRUE)),
                     correlation = 1 - r)
      dist[is.nan(dist) | is.na(dist)]<-Inf
      dist[dist==0]<-ifelse(is.finite(min(dist[dist>0])), min(dist[dist>0])/2, 1)
      
      if (sum(is.finite(dist)) < k) {
        stop(message = "Fewer than K finite distances found")
      }
      k.genes.ind <- order(dist)[1:k]
      k.genes <- cand.genes[k.genes.ind]
      wghts <- 1/dist[k.genes.ind]/sum(1/dist[k.genes.ind])
      imp.knn[gene.num, exp.num] <- wghts %*% data[k.genes, exp.num]
    }
  }


  if (distance=="correlation") {
    imp.knn <- (imp.knn * genesd) + genemeans
  }
  if (!rm.na) {
    imp.knn[is.na(data) == TRUE & is.nan(data) == FALSE] <- NA
  }
  if (!rm.inf) {
    index <- is.finite(data) == FALSE & is.na(data) == FALSE &
    is.nan(data) == FALSE
    imp.knn[index] <- data[index]
  }
  if (!rm.nan) {
    imp.knn[is.nan(data) == TRUE] <- NaN
  }
  return(imp.knn)
}

