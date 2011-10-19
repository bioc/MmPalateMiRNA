####################################################################
## checkOutliers Method
####################################################################

setGeneric("checkOutliers",
            function(obj) standardGeneric("checkOutliers"))

setMethod("checkOutliers", signature(obj="RGList"), 
          function(obj) {

            ## Detect outliers outside range of mean +/- two std deviations
            ## Returns indexes of outlying observations in each channel (R,Rb and G,Gb)

            sdR <- apply(obj$R, 1, sd)
            sdG <- apply(obj$G, 1, sd)
            meR <- rowMeans(obj$R)
            meG <- rowMeans(obj$G)

            sdRb <- apply(obj$Rb, 1, sd)
            sdGb <- apply(obj$Gb, 1, sd)
            meRb <- rowMeans(obj$Rb)
            meGb <- rowMeans(obj$Gb)
            
            ## RED CHANNEL
            Rout <- which(meR + 2.665*sdR < apply(obj$R, 1, max))
            Rbout <- which(meRb + 2.665*sdRb < apply(obj$Rb, 1, max))
            ## GREEN CHANNEL
            Gout <- which(meG + 2.665*sdG < apply(obj$G, 1, max))
            Gbout <- which(meGb + 2.665*sdGb < apply(obj$Gb, 1, max))
            
            return(list(Rout=Rout, Rbout=Rbout, Gout=Gout, Gbout=Gbout))
          })

####################################################################
## checkMVs Method
####################################################################\

setGeneric("checkMVs",
            function(obj) standardGeneric("checkMVs"))

setMethod("checkMVs", signature(obj="RGList"), 
          function(obj) {

            R.na <- which(apply(obj$R, 1, function(y) any(is.na(y))))
            Rb.na <- which(apply(obj$Rb, 1, function(y) any(is.na(y))))
            G.na <- which(apply(obj$G, 1, function(y) any(is.na(y))))
            Gb.na <- which(apply(obj$Gb, 1, function(y) any(is.na(y))))
            
            return(list(R.na=R.na, Rb.na=Rb.na, G.na=G.na, Gb.na=Gb.na))
            
          })


####################################################################
## filterArray Method
####################################################################\

setGeneric("filterArray",
            function(obj, ... ) standardGeneric("filterArray"))

setMethod("filterArray", signature(obj="RGList"),
          function(obj, keep, frac, number, reps) {

            ## Keeping everything in toKeep
            ## 'valid' values have text in 'keep' AND
            ## have fg levels > frac*bg in at least 'number' samples 
            
            ## 1. Text Filter
            ## NOTE - need to document that text filter is based on 'genes$Name' ... 
            toKeep.txt <- unique(unlist(sapply(keep, function(x) grep(x, obj$genes$Name))))

            ## 2. Background Filter
            pass <- matrix(obj$R > frac*obj$Rb & obj$G > frac*obj$Gb, ncol=ncol(obj))
            toKeep.bg <- which(rowSums(pass) >= number) 

            toKeep <- intersect(toKeep.txt, toKeep.bg)
            reducedSet <- obj[toKeep,] 

            ## Third filter step (keep only probes with at least 'reps' replicates)
            replicates <- table(reducedSet$genes$Gene)
            gene.list <- as.numeric(attr(replicates, "dimnames")[[1]])[replicates>=reps]
            reducedSet <- reducedSet[reducedSet$genes$Gene %in% gene.list,]
            return(reducedSet)
          })



####################################################################
## Methods for lattice plots
## Need for densityplot and levelplot
## One for 'RGList', another for 'list' (with each item in list as 'RGList' object)
####################################################################


####################################################################
## densityplot 
####################################################################

setGeneric("densityplot", function(x, data, ...)
           standardGeneric("densityplot"))

####################################################################
## densityplot - class 'RGList'
####################################################################

setMethod("densityplot", signature(x="RGList", data = "missing"),
          function (x, channel=c("G", "R"), group=NULL, subset=NULL, ... ) {

            channel <- match.arg(channel)
            log2 <- as.data.frame(log2(x[[channel]]))
            ## need to add condition in case 'group' is NULL
            if (!is.null(group)) {
              log2$group <- x$genes[[group]]
            }
            if (!is.null(subset)) {
              log2 <- log2[log2$group%in%subset,]
              log2$group <- factor(log2$group)
            }
            if (is.null(group)) {
              form <- as.formula(paste(" ~ `", paste(names(log2)[-ncol(log2)], collapse = "` + `"), "`", sep=""))
            } else {
              form <- as.formula(paste(" ~ `", paste(names(log2)[-ncol(log2)], collapse = "` + `"),            
                                       "` | ", "group", sep=""))              
            }
            
            xlabtext <- ifelse(channel == "G",
                               "log2 Expression of Green (Control) Channel",
                               "log2 Expression of Red (Experimental) Channel")

            densityplot(form,
                        data=log2, 
                        plot.points=FALSE, 
                        allow.multiple=TRUE,
                        ylab = "Estimated Density",
                        xlab = xlabtext, ... )
          })


####################################################################
## densityplot - class 'list'
## Method for list of either MAList or NChannelSet objects
####################################################################

setMethod("densityplot", signature(x="list", data = "missing"),
          function (x, channel=c("G", "R"), group=NULL, subset=NULL, ... ) {

            channel <- match.arg(channel)
            classes <- sapply(x, class)
            if (any(!classes%in%c("MAList", "NChannelSet"))) {
              stop("Items in list need to be either class 'MAList' or 'NChannelSet'")
            }

            idx1 <- which(classes=="MAList")
            idx2 <- which(classes=="NChannelSet")
            res <- vector("list", length(x))
            res[idx1] <- lapply(x[idx1], function(x) log2(RG.MA(x)[[channel]]))
            res[idx2] <- lapply(x[idx2], function(x) assayData(x)[[channel]])  
            ## assayData() already stored as log2 values
            nlog2 <- as.data.frame(do.call(rbind, res))

            ## Normalization
            if(!is.null(names(x))) {
              nlog2$normalization <- rep(names(x), sapply(ndata.all, nrow))
              nlog2$normalization <- factor(nlog2$normalization, levels=names(x))
            } else {
              nnames <- paste("Normalization", 1:length(x), sep="")
              nlog2$normalization <- factor(rep(nnames, sapply(ndata.all, nrow)))
            }

            
            ## need to add condition in case 'group' is NULL
            if (!is.null(group)) {
              res <- vector("list", length(x))              
              res[idx1] <- lapply(x[idx1], function(x) x$genes[[group]])
              res[idx2] <- lapply(x[idx2], function(x) pData(featureData(x))[[group]])
              nlog2$group <- unlist(res)
            }
            if (!is.null(subset)) {
              nlog2 <- nlog2[nlog2$group%in%subset,]
              nlog2$group <- factor(nlog2$group)
            }
            if (is.null(group)) {
              form <- as.formula(paste(" ~ `", paste(names(nlog2)[1:ncol(x[[1]])], collapse = "` + `"),
                                       "` | ", "normalization", sep=""))
            } else {
              form <- as.formula(paste(" ~ `", paste(names(nlog2)[1:ncol(x[[1]])], collapse = "` + `"),            
                                       "` | ", "group + normalization", sep=""))
            }
          
            xlabtext <- ifelse(channel == "G",
                               "log2 Expression of Green (Control) Channel",
                               "log2 Expression of Red (Experimental) Channel")

            res <- densityplot(form,
                               data=nlog2, 
                               plot.points=FALSE, 
                               allow.multiple=TRUE,
                               ylab = "Estimated Density",
                               xlab = xlabtext, ... )

##            if (!is.null(group)) 
##              res <- useOuterStrips(res)              

            return(res)

          })




####################################################################
## levelplot
####################################################################

setGeneric("levelplot", function(x, data, ...)
           standardGeneric("levelplot"))

####################################################################
## levelplot - class 'RGList'
####################################################################

## Use 'levelplot' in lattice with separation by type of probe 
## Input can be a multi-dimensional array w/3rd dimension giving 
## conditioning variable (here probe type)

setMethod("levelplot", signature(x="RGList", data = "missing"),
          function (x, channel=c("G", "R"), group=NULL, subset=NULL, ... ) {

            channel <- match.arg(channel)
            nc <- ncol(x)
            log2 <- as.data.frame(log2(x[[channel]]))
            if (!is.null(group)) {
              log2$group <- x$genes[[group]]
            }
            if (!is.null(subset)) {
              log2 <- log2[log2$group%in%subset,]
              log2$group <- factor(log2$group)
            }            
            
            if(!is.null(group)) {
              ngrps <- nlevels(log2$group)
              dist.array <- array(NA, dim = c(nc, nc, ngrps))
              dimnames(dist.array) <- list(colnames(x), colnames(x), levels(log2$group))
            } else {
              dist.array <- array(NA, dim = c(nc, nc))
              dimnames(dist.array) <- list(colnames(x), colnames(x))
            }

            if (!is.null(group)) {
              for (i in 1:nc) {
                for (j in 1:nc) {
                  for (k in 1:ngrps) {
                    idx <- which(log2$group == levels(log2$group)[k])
                    dist.array[i,j,k] <- median(abs(log2[idx, i] - log2[idx, j]))
                    diag(dist.array[,,k]) <- NA ## no color for diagonals 
                  }
                }
              }
            } else {
              for (i in 1:nc) {
                for (j in 1:nc) {
                  dist.array[i,j] <- median(abs(log2[, i] - log2[, j]))
                  diag(dist.array) <- NA ## no color for diagonals 
                }
              }
            }

            res <- levelplot(dist.array, 
                             xlab = "Median of absolute differences in log2 expression", 
                             ylab="", ... )


          })


####################################################################
## levelplot - class 'list'
## Method for list of either MAList or NChannelSet objects
####################################################################

## take out 'group' and 'subset' arguments
setMethod("levelplot", signature(x="list", data = "missing"),
          function (x, channel=c("G", "R"), order=NULL, ... ) {

            channel <- match.arg(channel)
            classes <- sapply(x, class)
            if (any(!classes%in%c("MAList", "NChannelSet"))) {
              stop("Items in list need to be either class 'MAList' or 'NChannelSet'")
            }

            idx1 <- which(classes=="MAList")
            idx2 <- which(classes=="NChannelSet")
            log2 <- vector("list", length(x))
            ## condition on whether reordered
            if(is.null(order)) {
              log2[idx1] <- lapply(x[idx1], function(x) log2(RG.MA(x)[[channel]]))
              log2[idx2] <- lapply(x[idx2], function(x) assayData(x)[[channel]])
            } else {
              log2[idx1] <- lapply(x[idx1], function(x) log2(RG.MA(x)[[channel]][,order]))
              log2[idx2] <- lapply(x[idx2], function(x) assayData(x)[[channel]][,order])
            }
            ## assayData() already stored as log2 values

            nc <- ncol(log2[[1]])
            ngrps <- length(x)
            dist.array <- array(NA, dim = c(nc, nc, ngrps))

            if(is.null(order)) {
              cnames <- colnames(log2[[1]])
            } else {
              cnames <- colnames(log2[[1]])[order]
            }
            if(is.null(names(x))) {
              nnames <- paste("Normalization", 1:length(x))
            } else {
              nnames <- names(x)
            }

            dimnames(dist.array) <- list(cnames, cnames, nnames) 

            for (i in 1:nc) {
              for (j in 1:nc) {
                for (k in 1:ngrps) {
                  dist.array[i,j,k] <- median(abs(log2[[k]][,i] - log2[[k]][, j]), na.rm=TRUE)
                  diag(dist.array[,,k]) <- NA ## no color for diagonals 
                }
              }
            }

            res <- levelplot(dist.array, 
                             xlab = "Median of absolute differences in log2 expression", 
                             ylab="", ... )
            

          })


####################################################################
## MADvsMedianPlot - uses xyplot in lattice
####################################################################

setGeneric("MADvsMedianPlot", function(x, ...)
           standardGeneric("MADvsMedianPlot"))


####################################################################
## MADvsMedianPlot - class 'list'
## Method for list of either MAList or NChannelSet objects
####################################################################

setMethod("MADvsMedianPlot", signature(x="list"),
          function (x, channel=c("G", "R"), group=NULL, subset=NULL, ... ) {

            channel <- match.arg(channel)
            classes <- sapply(x, class)
            if (any(!classes%in%c("MAList", "NChannelSet"))) {
              stop("Items in list need to be either class 'MAList' or 'NChannelSet'")
            }

            idx1 <- which(classes=="MAList")
            idx2 <- which(classes=="NChannelSet")
            res <- vector("list", length(x))
            res[idx1] <- lapply(x[idx1], function(x) log2(RG.MA(x)[[channel]]))
            res[idx2] <- lapply(x[idx2], function(x) assayData(x)[[channel]])  
            ## assayData() already stored as log2 values
            ## nlog2 <- as.data.frame(do.call(rbind, res))

            MADs <- sapply(res, function(x) apply(x, 1, mad))
            medians <- sapply(res, function(x) apply(x, 1, median))
            if (!is.null(names(x))) {
              colnames(MADs) <- colnames(medians) <- names(x)
            } else {
              colnames(MADs) <- paste("Normalization", 1:length(x), sep="")
            }

            ## stack this result into single data frame
            res.df <- stack(as.data.frame(MADs))
            names(res.df) <- c("MAD", "Method")
            res.df$Medians <- stack(as.data.frame(medians))$values

            ## need to add condition in case 'group' is NULL
            if (!is.null(group)) {
              res <- vector("list", length(x))              
              res[idx1] <- lapply(x[idx1], function(x) x$genes[[group]])
              res[idx2] <- lapply(x[idx2], function(x) pData(featureData(x))[[group]])
              res.df$group <- unlist(res)
            }
            if (!is.null(subset)) {
              res.df <- res.df[res.df$group%in%subset,]
              res.df$group <- factor(res.df$group)
            }
            
            if (!is.null(group)) {
              res <- xyplot(MAD ~ Medians | Method, data=res.df, groups = group,  auto.key=TRUE, ...)
            } else {
              res <- xyplot(MAD ~ Medians | Method, data=res.df, ...)
            }
            return(res)
            
          })
            
            
            
####################################################################
## MAplot - uses xyplot in lattice
## Note other packages have versions of 'MAplot'
##   which may be superior (affy, ... )
####################################################################

setGeneric("MAplot", function(x, ...)
           standardGeneric("MAplot"))


####################################################################
## MAPlot - class 'MAList'
####################################################################

setMethod("MAplot", signature(x="MAList"),
          function (x,  ... ) {

            res.df <- stack(as.data.frame(x$A))
            names(res.df) <- c("A", "array")
            res.df$M <- stack(as.data.frame(x$M))$values
            xyplot(M ~ A | array, data=res.df,
                   panel = function(x, y, ...) {
                     panel.xyplot(x, y, col="black", ...)
                     panel.loess(x, y, col="red", lwd=2, ...)
                   }, ... )
          })

####################################################################
## MAPlot - class 'NChannelSet'
####################################################################

setMethod("MAplot", signature(x="NChannelSet"),
          function (x,  ... ) {
            
            ## FINISH HERE ..
            x.Mvals <- (assayData(x)$R - assayData(x)$G) 
            x.Avals <- (assayData(x)$R + assayData(x)$G)

            res.df <- stack(as.data.frame(x.Avals))
            names(res.df) <- c("A", "array")
            res.df$M <- stack(as.data.frame(x.Mvals))$values
            xyplot(M ~ A | array, data=res.df,
                   panel = function(x, y, ...) {
                     panel.xyplot(x, y, col="black", ...)
                     panel.loess(x, y, col="red", lwd=2, ...)
                   }, ... )


          })
