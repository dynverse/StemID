#' SCseq class
#'
#' @param expdata Expression data
#' @slot expdata Expression data
#' @slot ndata cell data?
#' @slot fdata feature data
#' @slot distances distances between cells?
#' @slot tsne tsne dimred
#' @slot cluster clustering
#' @slot background ?
#' @slot out ?
#' @slot cpart ?
#' @slot fcol ?
#' @slot filterpar ?
#' @slot clusterpar ?
#' @slot outlierpar ?
#'
#' @import stats
#' @import graphics
#' @import grDevices
#' @import pheatmap
#' @import RColorBrewer
#'
#' @name SCseq-class
#' @rdname SCseq-class
#' @export SCseq
#' @exportClass SCseq
SCseq <- setClass(
  "SCseq",
  slots = c(
    expdata = "data.frame",
    ndata = "data.frame",
    fdata = "data.frame",
    distances = "matrix",
    tsne = "data.frame",
    cluster = "list",
    background = "list",
    out = "list",
    cpart = "vector",
    fcol = "vector",
    filterpar = "list",
    clusterpar = "list",
    outlierpar ="list"
  )
)

setValidity(
  "SCseq",
  function(object) {
    msg <- NULL
    if ( ! is.data.frame(object@expdata) ){
      msg <- c(msg, "input data must be data.frame")
    } else if ( nrow(object@expdata) < 2 ){
      msg <- c(msg, "input data must have more than one row")
    } else if ( ncol(object@expdata) < 2 ){
      msg <- c(msg, "input data must have more than one column")
    } else if (sum( apply( is.na(object@expdata),1,sum ) ) > 0 ){
      msg <- c(msg, "NAs are not allowed in input data")
    } else if (sum( apply( object@expdata,1,min ) ) < 0 ){
      msg <- c(msg, "negative values are not allowed in input data")
    }
    if (is.null(msg)) TRUE
    else msg
  }
)

#' Constructor method of class SCseq
#'
#' @param SCseq scseq object
#'
#' @name SCseq
#' @rdname SCseq-class
setMethod(
  "initialize",
  signature = "SCseq",
  definition = function(.Object, expdata ){
    .Object@expdata <- expdata
    .Object@ndata <- expdata
    .Object@fdata <- expdata
    validObject(.Object)
    return(.Object)
  }
)

#' Method filterdata
#' @name filterdata
#' @rdname filterdata-methods
#' @exportMethod filterdata
setGeneric(
  "filterdata",
  function(object, mintotal=3000, minexpr=5, minnumber=1, maxexpr=Inf, downsample=TRUE, dsn=1, rseed=17000) standardGeneric("filterdata")
)

#' @rdname filterdata-methods
#' @aliases filterdata,SCseq-method
setMethod(
  "filterdata",
  signature = "SCseq",
  definition = function(object,mintotal,minexpr,minnumber,maxexpr,downsample,dsn,rseed) {
    if ( ! is.numeric(mintotal) ) stop( "mintotal has to be a positive number" ) else if ( mintotal <= 0 ) stop( "mintotal has to be a positive number" )
    if ( ! is.numeric(minexpr) ) stop( "minexpr has to be a non-negative number" ) else if ( minexpr < 0 ) stop( "minexpr has to be a non-negative number" )
    if ( ! is.numeric(minnumber) ) stop( "minnumber has to be a non-negative integer number" ) else if ( round(minnumber) != minnumber | minnumber < 0 ) stop( "minnumber has to be a non-negative integer number" )
    if ( ! ( is.numeric(downsample) | is.logical(downsample) ) ) stop( "downsample has to be logical (TRUE/FALSE)" )
    if ( ! is.numeric(dsn) ) stop( "dsn has to be a positive integer number" ) else if ( round(dsn) != dsn | dsn <= 0 ) stop( "dsn has to be a positive integer number" )
    object@filterpar <- list(mintotal=mintotal, minexpr=minexpr, minnumber=minnumber, maxexpr=maxexpr, downsample=downsample, dsn=dsn)
    object@ndata <- object@expdata[,apply(object@expdata,2,sum,na.rm=TRUE) >= mintotal]
    if ( downsample ){
      # set.seed(rseed) # dyngen: don't set a seed
      object@ndata <- downsample(object@expdata,n=mintotal,dsn=dsn)
    }else{
      x <- object@ndata
      object@ndata <- as.data.frame( t(t(x)/apply(x,2,sum))*median(apply(x,2,sum,na.rm=TRUE)) + .1 )
    }
    x <- object@ndata
    object@fdata <- x[apply(x>=minexpr,1,sum,na.rm=TRUE) >= minnumber,]
    x <- object@fdata
    object@fdata <- x[apply(x,1,max,na.rm=TRUE) < maxexpr,]
    return(object)
  }
)


downsample <- function(x,n,dsn){
  x <- round( x[,apply(x,2,sum,na.rm=TRUE) >= n], 0)
  nn <- min( apply(x,2,sum) )
  for ( j in 1:dsn ){
    z  <- data.frame(GENEID=rownames(x))
    rownames(z) <- rownames(x)
    initv <- rep(0,nrow(z))
    for ( i in 1:dim(x)[2] ){
      y <- aggregate(rep(1,nn),list(sample(rep(rownames(x),x[,i]),nn)),sum)
      na <- names(x)[i]
      names(y) <- c("GENEID",na)
      rownames(y) <- y$GENEID
      z[,na] <- initv
      k <- intersect(rownames(z),y$GENEID)
      z[k,na] <- y[k,na]
      z[is.na(z[,na]),na] <- 0
    }
    rownames(z) <- as.vector(z$GENEID)
    ds <- if ( j == 1 ) z[,-1] else ds + z[,-1]
  }
  ds <- ds/dsn + .1
  return(ds)
}


dist.gen <- function(x,method="euclidean", ...) {
  if ( method %in% c("spearman","pearson","kendall") )
    as.dist( 1 - cor(t(x),method=method,...) )
  else
    dist(x,method=method,...)
}

dist.gen.pairs <- function(x,y,...) {
  dist.gen(t(cbind(x,y)),...)
}


binompval <- function(p,N,n){
  pval   <- pbinom(n,round(N,0),p,lower.tail=TRUE)
  pval[!is.na(pval) & pval > 0.5] <- 1-pval[!is.na(pval) & pval > 0.5]
  return(pval)
}

setGeneric("plotgap", function(object) standardGeneric("plotgap"))


setMethod(
  "plotgap",
  signature = "SCseq",
  definition = function(object){
    if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before plotgap")
    if ( sum(is.na(object@cluster$gap$Tab[,3])) > 0 ) stop("run clustexp with do.gap = TRUE first")
    plot(object@cluster$gap,ylim=c( min(object@cluster$gap$Tab[,3] - object@cluster$gap$Tab[,4]),  max(object@cluster$gap$Tab[,3] + object@cluster$gap$Tab[,4])))
  }
)

setGeneric("plotjaccard", function(object) standardGeneric("plotjaccard"))


setMethod(
  "plotjaccard",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before plotjaccard")
            if ( length(unique(object@cluster$kpart)) < 2 ) stop("only a single cluster: no Jaccard's similarity plot")
            barplot(object@cluster$jaccard,names.arg=1:length(object@cluster$jaccard),ylab="Jaccard's similarity")
          }
)

setGeneric("plotoutlierprobs", function(object) standardGeneric("plotoutlierprobs"))

setMethod(
  "plotoutlierprobs",
  signature = "SCseq",
  definition = function(object){
    if ( length(object@cpart) == 0 ) stop("run findoutliers before plotoutlierprobs")
    p <- object@cluster$kpart[ order(object@cluster$kpart,decreasing=FALSE)]
    x <- object@out$cprobs[names(p)]
    fcol <- object@fcol
    for ( i in 1:max(p) ){
      y <- -log10(x + 2.2e-16)
      y[p != i] <- 0
      if ( i == 1 )
        b <- barplot(y,ylim=c(0,max(-log10(x + 2.2e-16))*1.1),col=fcol[i],border=fcol[i],names.arg=FALSE,ylab="-log10prob")
      else
        barplot(y,add=TRUE,col=fcol[i],border=fcol[i],names.arg=FALSE,axes=FALSE)
    }
    abline(-log10(object@outlierpar$probthr),0,col="black",lty=2)
    d <- b[2,1] - b[1,1]
    y <- 0
    for ( i in 1:max(p) ) y <- append(y,b[sum(p <=i),1] + d/2)
    axis(1,at=(y[1:(length(y)-1)] + y[-1])/2,labels=1:max(p))
    box()
  }
)

setGeneric("plotbackground", function(object) standardGeneric("plotbackground"))

#' @importFrom locfit locfit

setMethod(
  "plotbackground",
  signature = "SCseq",
  definition = function(object){
    if ( length(object@cpart) == 0 ) stop("run findoutliers before plotbackground")
    m <- apply(object@fdata,1,mean)
    v <- apply(object@fdata,1,var)
    fit <- locfit(v~lp(m,nn=.7),family="gamma",maxk=500)
    plot(log2(m),log2(v),pch=20,xlab="log2mean",ylab="log2var",col="grey")
    lines(log2(m[order(m)]),log2(object@background$lvar(m[order(m)],object)),col="red",lwd=2)
    lines(log2(m[order(m)]),log2(fitted(fit)[order(m)]),col="orange",lwd=2,lty=2)
    legend("topleft",legend=substitute(paste("y = ",a,"*x^2 + ",b,"*x + ",d,sep=""),list(a=round(coef(object@background$vfit)[3],2),b=round(coef(object@background$vfit)[2],2),d=round(coef(object@background$vfit)[1],2))),bty="n")
  }
)

setGeneric("plotsensitivity", function(object) standardGeneric("plotsensitivity"))


setMethod(
  "plotsensitivity",
  signature = "SCseq",
  definition = function(object){
    if ( length(object@cpart) == 0 ) stop("run findoutliers before plotsensitivity")
    plot(log10(object@out$thr), object@out$stest, type="l",xlab="log10 Probability cutoff", ylab="Number of outliers")
    abline(v=log10(object@outlierpar$probthr),col="red",lty=2)
  }
)

setGeneric("diffgenes", function(object,cl1,cl2,mincount=5) standardGeneric("diffgenes"))

setMethod(
  "diffgenes",
  signature = "SCseq",
  definition = function(object,cl1,cl2,mincount){
    part <- object@cpart
    cl1 <- c(cl1)
    cl2 <- c(cl2)
    if ( length(part) == 0 ) stop("run findoutliers before diffgenes")
    if ( ! is.numeric(mincount) ) stop("mincount has to be a non-negative number") else if (  mincount < 0 ) stop("mincount has to be a non-negative number")
    if ( length(intersect(cl1, part)) < length(unique(cl1)) ) stop( paste("cl1 has to be a subset of ",paste(sort(unique(part)),collapse=","),"\n",sep="") )
    if ( length(intersect(cl2, part)) < length(unique(cl2)) ) stop( paste("cl2 has to be a subset of ",paste(sort(unique(part)),collapse=","),"\n",sep="") )
    f <- apply(object@ndata[,part %in% c(cl1,cl2)],1,max) > mincount
    x <- object@ndata[f,part %in% cl1]
    y <- object@ndata[f,part %in% cl2]
    if ( sum(part %in% cl1) == 1 ) m1 <- x else m1 <- apply(x,1,mean)
    if ( sum(part %in% cl2) == 1 ) m2 <- y else m2 <- apply(y,1,mean)
    if ( sum(part %in% cl1) == 1 ) s1 <- sqrt(x) else s1 <- sqrt(apply(x,1,var))
    if ( sum(part %in% cl2) == 1 ) s2 <- sqrt(y) else s2 <- sqrt(apply(y,1,var))

    d <- ( m1 - m2 )/ apply( cbind( s1, s2 ),1,mean )
    names(d) <- rownames(object@ndata)[f]
    if ( sum(part %in% cl1) == 1 ){
      names(x) <- names(d)
      x <- x[order(d,decreasing=TRUE)]
    }else{
      x <- x[order(d,decreasing=TRUE),]
    }
    if ( sum(part %in% cl2) == 1 ){
      names(y) <- names(d)
      y <- y[order(d,decreasing=TRUE)]
    }else{
      y <- y[order(d,decreasing=TRUE),]
    }
    return(list(z=d[order(d,decreasing=TRUE)],cl1=x,cl2=y,cl1n=cl1,cl2n=cl2))
  }
)


plotdiffgenes <- function(z,gene){
  if ( ! is.list(z) ) stop("first arguments needs to be output of function diffgenes")
  if ( length(z) < 5 ) stop("first arguments needs to be output of function diffgenes")
  if ( sum(names(z) == c("z","cl1","cl2","cl1n","cl2n")) < 5 ) stop("first arguments needs to be output of function diffgenes")
  if ( length(gene) > 1 ) stop("only single value allowed for argument gene")
  if ( !is.numeric(gene) & !(gene %in% names(z$z)) ) stop("argument gene needs to be within rownames of first argument or a positive integer number")
  if ( is.numeric(gene) ){ if ( gene < 0 | round(gene) != gene ) stop("argument gene needs to be within rownames of first argument or a positive integer number") }
  x <- if ( is.null(dim(z$cl1)) ) z$cl1[gene] else t(z$cl1[gene,])
  y <- if ( is.null(dim(z$cl2)) ) z$cl2[gene] else t(z$cl2[gene,])
  plot(1:length(c(x,y)),c(x,y),ylim=c(0,max(c(x,y))),xlab="",ylab="Expression",main=gene,cex=0,axes=FALSE)
  axis(2)
  box()
  u <- 1:length(x)
  rect(u - .5,0,u + .5,x,col="red")
  v <- c(min(u) - .5,max(u) + .5)
  axis(1,at=mean(v),labels = paste(z$cl1n,collapse=","))
  lines(v,rep(mean(x),length(v)))
  lines(v,rep(mean(x)-sqrt(var(x)),length(v)),lty=2)
  lines(v,rep(mean(x)+sqrt(var(x)),length(v)),lty=2)

  u <- ( length(x) + 1 ):length(c(x,y))
  v <- c(min(u) - .5,max(u) + .5)
  rect(u - .5,0,u + .5,y,col="blue")
  axis(1,at=mean(v),labels = paste(z$cl2n,collapse=","))
  lines(v,rep(mean(y),length(v)))
  lines(v,rep(mean(y)-sqrt(var(y)),length(v)),lty=2)
  lines(v,rep(mean(y)+sqrt(var(y)),length(v)),lty=2)
  abline(v=length(x) + .5)
}

setGeneric("plottsne", function(object,final=TRUE) standardGeneric("plottsne"))


setMethod(
  "plottsne",
  signature = "SCseq",
  definition = function(object,final){
    if ( length(object@tsne) == 0 ) stop("run comptsne before plottsne")
    if ( final & length(object@cpart) == 0 ) stop("run findoutliers before plottsne")
    if ( !final & length(object@cluster$kpart) == 0 ) stop("run clustexp before plottsne")
    part <- if ( final ) object@cpart else object@cluster$kpart
    plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,cex=1.5,col="lightgrey")
    for ( i in 1:max(part) ){
      if ( sum(part == i) > 0 )
        text(object@tsne[part == i,1],object@tsne[part == i,2],i,col=object@fcol[i],cex=.75,font=4)
    }
  }
)

setGeneric("plotlabelstsne", function(object,labels=NULL) standardGeneric("plotlabelstsne"))


setMethod(
  "plotlabelstsne",
  signature = "SCseq",
  definition = function(object,labels){
    if ( is.null(labels ) ) labels <- names(object@ndata)
    if ( length(object@tsne) == 0 ) stop("run comptsne before plotlabelstsne")
    plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,cex=1.5,col="lightgrey")
    text(object@tsne[,1],object@tsne[,2],labels,cex=.5)
  }
)

setGeneric("plotsymbolstsne", function(object,types=NULL) standardGeneric("plotsymbolstsne"))



setMethod(
  "plotsymbolstsne",
  signature = "SCseq",
  definition = function(object,types){
    if ( is.null(types) ) types <- names(object@fdata)
    if ( length(object@tsne) == 0 ) stop("run comptsne before plotsymbolstsne")
    if ( length(types) != ncol(object@fdata) ) stop("types argument has wrong length. Length has to equal to the column number of object@ndata")
    coloc <- rainbow(length(unique(types)))
    syms <- c()
    plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,col="grey")
    for ( i in 1:length(unique(types)) ){
      f <- types == sort(unique(types))[i]
      syms <- append( syms, ( (i-1) %% 25 ) + 1 )
      points(object@tsne[f,1],object@tsne[f,2],col=coloc[i],pch=( (i-1) %% 25 ) + 1,cex=1)
    }
    legend("topleft", legend=sort(unique(types)), col=coloc, pch=syms)
  }
)

setGeneric("plotexptsne", function(object,g,n="",logsc=FALSE) standardGeneric("plotexptsne"))



setMethod(
  "plotexptsne",
  signature = "SCseq",
  definition = function(object,g,n="",logsc=FALSE){
    if ( length(object@tsne) == 0 ) stop("run comptsne before plottsne")
    if ( length(intersect(g,rownames(object@ndata))) < length(unique(g)) ) stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
    if ( !is.numeric(logsc) & !is.logical(logsc) ) stop("argument logsc has to be logical (TRUE/FALSE)")
    if ( n == "" ) n <- g[1]
    l <- apply(object@ndata[g,] - .1,2,sum) + .1
    if (logsc) {
      f <- l == 0
      l <- log2(l)
      l[f] <- NA
    }
    mi <- min(l,na.rm=TRUE)
    ma <- max(l,na.rm=TRUE)
    ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    v <- round((l - mi)/(ma - mi)*99 + 1,0)
    layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
    par(mar = c(3,5,2.5,2))
    plot(object@tsne,xlab="Dim 1",ylab="Dim 2",main=n,pch=20,cex=0,col="grey")
    for ( k in 1:length(v) ){
      points(object@tsne[k,1],object@tsne[k,2],col=ColorRamp[v[k]],pch=20,cex=1.5)
    }
    par(mar = c(3,2.5,2.5,2))
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
    layout(1)
  }
)

plot.err.bars.y <- function(x, y, y.err, col="black", lwd=1, lty=1, h=0.1){
  arrows(x,y-y.err,x,y+y.err,code=0, col=col, lwd=lwd, lty=lty)
  arrows(x-h,y-y.err,x+h,y-y.err,code=0, col=col, lwd=lwd, lty=lty)
  arrows(x-h,y+y.err,x+h,y+y.err,code=0, col=col, lwd=lwd, lty=lty)
}

clusGapExt <-function (x, FUNcluster, K.max, B = 100, verbose = interactive(), method="euclidean",random=TRUE, ...) {
  stopifnot(is.function(FUNcluster), length(dim(x)) == 2, K.max >= 2, (n <- nrow(x)) >= 1, (p <- ncol(x)) >= 1)
  if (B != (B. <- as.integer(B)) || (B <- B.) <= 0)
    stop("'B' has to be a positive integer")
  if (is.data.frame(x))
    x <- as.matrix(x)
  ii <- seq_len(n)
  W.k <- function(X, kk) {
    clus <- if (kk > 1)
      FUNcluster(X, kk, ...)$cluster
    else rep.int(1L, nrow(X))
    0.5 * sum(vapply(split(ii, clus), function(I) {
      xs <- X[I, , drop = FALSE]
      sum(dist.gen(xs,method=method)/nrow(xs))
    }, 0))
  }
  logW <- E.logW <- SE.sim <- numeric(K.max)
  if (verbose)
    cat("Clustering k = 1,2,..., K.max (= ", K.max, "): .. ",
        sep = "")
  for (k in 1:K.max) logW[k] <- log(W.k(x, k))
  if (verbose)
    cat("done\n")
  xs <- scale(x, center = TRUE, scale = FALSE)
  m.x <- rep(attr(xs, "scaled:center"), each = n)
  V.sx <- svd(xs)$v
  rng.x1 <- apply(xs %*% V.sx, 2, range)
  logWks <- matrix(0, B, K.max)
  if (random){
    if (verbose)
      cat("Bootstrapping, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n",
          sep = "")
    for (b in 1:B) {
      z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1],
                                                   max = M[2]), nn = n)
      z <- tcrossprod(z1, V.sx) + m.x
      ##z <- apply(x,2,function(m) runif(length(m),min=min(m),max=max(m)))
      ##z <- apply(x,2,function(m) sample(m))
      for (k in 1:K.max) {
        logWks[b, k] <- log(W.k(z, k))
      }
      if (verbose)
        cat(".", if (b%%50 == 0)
          paste(b, "\n"))
    }
    if (verbose && (B%%50 != 0))
      cat("", B, "\n")
    E.logW <- colMeans(logWks)
    SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
  }else{
    E.logW <- rep(NA,K.max)
    SE.sim <- rep(NA,K.max)
  }
  structure(
    class = "clusGap",
    list(
      Tab = cbind(logW, E.logW, gap = E.logW - logW, SE.sim),
      n = n,
      B = B,
      FUNcluster = FUNcluster
    )
  )
}


#' @importFrom cluster pam maxSE
#' @importFrom fpc kmeansCBI clusterboot hclustCBI pamkCBI
clustfun <- function(x,clustnr=20,bootnr=50,metric="pearson",do.gap=FALSE,sat=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000,FUNcluster="kmedoids",distances=NULL,link="single")
{
  if ( clustnr < 2) stop("Choose clustnr > 1")
  di <- if ( FUNcluster == "kmedoids" ) t(x) else dist.gen(t(x),method=metric)
  if ( nrow(di) - 1 < clustnr ) clustnr <-  nrow(di) - 1
  if ( do.gap | sat | cln > 0 ){
    gpr <- NULL
    f <- if ( cln == 0 ) TRUE else FALSE
    if ( do.gap ){
      # set.seed(rseed) # dyneval: don't set a seed
      if ( FUNcluster == "kmeans" )   gpr <- clusGapExt(as.matrix(di), FUNcluster = kmeans, K.max = clustnr, B = B.gap, iter.max=100)
      if ( FUNcluster == "kmedoids" ) gpr <- clusGapExt(as.matrix(di), FUNcluster = function(x,k) pam(dist.gen(x,method=metric),k), K.max = clustnr, B = B.gap, method=metric)
      if ( FUNcluster == "hclust" )   gpr <- clusGapExt(as.matrix(di), FUNcluster = function(x,k){ y <- hclustCBI(x,k,link=link,scaling=FALSE); y$cluster <- y$partition; y }, K.max = clustnr, B = B.gap)
      if ( f ) cln <- maxSE(gpr$Tab[,3],gpr$Tab[,4],method=SE.method,SE.factor)
    }
    if ( sat ){
      if ( ! do.gap ){
        if ( FUNcluster == "kmeans" )   gpr <- clusGapExt(as.matrix(di), FUNcluster = kmeans, K.max = clustnr, B = B.gap, iter.max=100, random=FALSE)
        if ( FUNcluster == "kmedoids" ) gpr <- clusGapExt(as.matrix(di), FUNcluster = function(x,k) pam(dist.gen(x,method=metric),k), K.max = clustnr, B = B.gap, random=FALSE, method=metric)
        if ( FUNcluster == "hclust" )   gpr <- clusGapExt(as.matrix(di), FUNcluster = function(x,k){ y <- hclustCBI(x,k,link=link,scaling=FALSE); y$cluster <- y$partition; y }, K.max = clustnr, B = B.gap, random=FALSE)
      }
      g <- gpr$Tab[,1]
      y <- g[-length(g)] - g[-1]
      mm <- numeric(length(y))
      nn <- numeric(length(y))
      for ( i in 1:length(y)){
        mm[i] <- mean(y[i:length(y)])
        nn[i] <- sqrt(var(y[i:length(y)]))
      }
      if ( f ) cln <- max(min(which( y - (mm + nn) < 0 )),1)
    }
    if ( cln <= 1 ) {
      clb <- list(result=list(partition=rep(1,dim(x)[2])),bootmean=1)
      names(clb$result$partition) <- names(x)
      return(list(x=x,clb=clb,gpr=gpr,di=if ( FUNcluster == "kmedoids" ) dist.gen(di,method=metric) else di))
    }
    if ( FUNcluster == "kmeans" ) clb <- clusterboot(di,B=bootnr,distances=FALSE,bootmethod="boot",clustermethod=kmeansCBI,krange=cln,scaling=FALSE,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
    if ( FUNcluster == "kmedoids" ) clb <- clusterboot(dist.gen(di,method=metric),B=bootnr,bootmethod="boot",clustermethod=pamkCBI,k=cln,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
    if ( FUNcluster == "hclust" ) clb <- clusterboot(di,B=bootnr,distances=FALSE,bootmethod="boot",clustermethod=hclustCBI,k=cln,link=link,scaling=FALSE,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
    return(list(x=x,clb=clb,gpr=gpr,di=if ( FUNcluster == "kmedoids" ) dist.gen(di,method=metric) else di))
  }
}

#' Method clustexp
#' @name clustexp
#' @rdname clustexp-methods
#' @exportMethod clustexp
setGeneric(
  "clustexp",
  function(object,clustnr=20,bootnr=50,metric="pearson",do.gap=FALSE,sat=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000,FUNcluster="kmedoids")
    standardGeneric("clustexp"))

#' @rdname clustexp-methods
#' @aliases clustexp,SCseq-method
setMethod(
  "clustexp",
  signature = "SCseq",
  definition = function(object,clustnr,bootnr,metric,do.gap,sat,SE.method,SE.factor,B.gap,cln,rseed,FUNcluster) {
    if ( ! is.numeric(clustnr) ) stop("clustnr has to be a positive integer") else if ( round(clustnr) != clustnr | clustnr <= 0 ) stop("clustnr has to be a positive integer")
    if ( ! is.numeric(bootnr) ) stop("bootnr has to be a positive integer") else if ( round(bootnr) != bootnr | bootnr <= 0 ) stop("bootnr has to be a positive integer")
    if ( ! ( metric %in% c( "spearman","pearson","kendall","euclidean","maximum","manhattan","canberra","binary","minkowski") ) ) stop("metric has to be one of the following: spearman, pearson, kendall, euclidean, maximum, manhattan, canberra, binary, minkowski")
    if ( ! ( SE.method %in% c( "firstSEmax","Tibs2001SEmax","globalSEmax","firstmax","globalmax") ) ) stop("SE.method has to be one of the following: firstSEmax, Tibs2001SEmax, globalSEmax, firstmax, globalmax")
    if ( ! is.numeric(SE.factor) ) stop("SE.factor has to be a non-negative integer") else if  ( SE.factor < 0 )  stop("SE.factor has to be a non-negative integer")
    if ( ! ( is.numeric(do.gap) | is.logical(do.gap) ) ) stop( "do.gap has to be logical (TRUE/FALSE)" )
    if ( ! ( is.numeric(sat) | is.logical(sat) ) ) stop( "sat has to be logical (TRUE/FALSE)" )
    if ( ! is.numeric(B.gap) ) stop("B.gap has to be a positive integer") else if ( round(B.gap) != B.gap | B.gap <= 0 ) stop("B.gap has to be a positive integer")
    if ( ! is.numeric(cln) ) stop("cln has to be a non-negative integer") else if ( round(cln) != cln | cln < 0 ) stop("cln has to be a non-negative integer")
    if ( ! is.numeric(rseed) ) stop("rseed has to be numeric")
    if ( !do.gap & !sat & cln == 0 ) stop("cln has to be a positive integer or either do.gap or sat has to be TRUE")
    if ( ! ( FUNcluster %in% c("kmeans","hclust","kmedoids") ) ) stop("FUNcluster has to be one of the following: kmeans, hclust,kmedoids")
    object@clusterpar <- list(clustnr=clustnr,bootnr=bootnr,metric=metric,do.gap=do.gap,sat=sat,SE.method=SE.method,SE.factor=SE.factor,B.gap=B.gap,cln=cln,rseed=rseed,FUNcluster=FUNcluster)
    y <- clustfun(object@fdata,clustnr,bootnr,metric,do.gap,sat,SE.method,SE.factor,B.gap,cln,rseed,FUNcluster)
    object@cluster   <- list(kpart=y$clb$result$partition, jaccard=y$clb$bootmean, gap=y$gpr, clb=y$clb)
    object@distances <- as.matrix( y$di )
    # set.seed(111111) # dyneval: don't set a seed
    object@fcol <- sample(rainbow(max(y$clb$result$partition)))
    return(object)
  }
)

#' Method findoutliers
#' @name findoutliers
#' @rdname findoutliers-methods
#' @exportMethod findoutliers
setGeneric(
  "findoutliers",
  function(object,outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.95)
    standardGeneric("findoutliers")
)

#' @rdname findoutliers-methods
#' @aliases findoutliers,SCseq-method
setMethod(
  "findoutliers",
  signature = "SCseq",
  definition = function(object,outminc,outlg,probthr,thr,outdistquant) {
    if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before findoutliers")
    if ( ! is.numeric(outminc) ) stop("outminc has to be a non-negative integer") else if ( round(outminc) != outminc | outminc < 0 ) stop("outminc has to be a non-negative integer")
    if ( ! is.numeric(outlg) ) stop("outlg has to be a non-negative integer") else if ( round(outlg) != outlg | outlg < 0 ) stop("outlg has to be a non-negative integer")
    if ( ! is.numeric(probthr) ) stop("probthr has to be a number between 0 and 1") else if (  probthr < 0 | probthr > 1 ) stop("probthr has to be a number between 0 and 1")
    if ( ! is.numeric(thr) ) stop("thr hast to be a vector of numbers between 0 and 1") else if ( min(thr) < 0 | max(thr) > 1 ) stop("thr hast to be a vector of numbers between 0 and 1")
    if ( ! is.numeric(outdistquant) ) stop("outdistquant has to be a number between 0 and 1") else if (  outdistquant < 0 | outdistquant > 1 ) stop("outdistquant has to be a number between 0 and 1")

    object@outlierpar <- list( outminc=outminc,outlg=outlg,probthr=probthr,thr=thr,outdistquant=outdistquant )
    ### calibrate background model
    m <- log2(apply(object@fdata,1,mean))
    v <- log2(apply(object@fdata,1,var))
    f <- m > -Inf & v > -Inf
    m <- m[f]
    v <- v[f]
    mm <- -8
    repeat{
      fit <- lm(v ~ m + I(m^2))
      if( coef(fit)[3] >= 0 | mm >= -1){
        break
      }
      mm <- mm + .5
      f <- m > mm
      m <- m[f]
      v <- v[f]
    }
    object@background <- list()
    object@background$vfit <- fit
    object@background$lvar <- function(x,object) 2**(coef(object@background$vfit)[1] + log2(x)*coef(object@background$vfit)[2] + coef(object@background$vfit)[3] * log2(x)**2)
    object@background$lsize <- function(x,object) x**2/(max(x + 1e-6,object@background$lvar(x,object)) - x)

    ### identify outliers
    out   <- c()
    stest <- rep(0,length(thr))
    cprobs <- c()
    for ( n in 1:max(object@cluster$kpart) ){
      if ( sum(object@cluster$kpart == n) == 1 ){
        cprobs <- append(cprobs,.5)
        names(cprobs)[length(cprobs)] <- names(object@cluster$kpart)[object@cluster$kpart == n]
        next
      }
      x <- object@fdata[,object@cluster$kpart == n]
      x <- x[apply(x,1,max) > outminc,]
      z <- t( apply(x,1,function(x){ apply( cbind( pnbinom(round(x,0),mu=mean(x),size=object@background$lsize(mean(x),object)) , 1 - pnbinom(round(x,0),mu=mean(x),size=object@background$lsize(mean(x),object)) ),1, min) } ) )
      cp <- apply(z,2,function(x){ y <- p.adjust(x,method="BH"); y <- y[order(y,decreasing=FALSE)]; return(y[outlg]);})
      f <- cp < probthr
      cprobs <- append(cprobs,cp)
      if ( sum(f) > 0 ) out <- append(out,names(x)[f])
      for ( j in 1:length(thr) )  stest[j] <-  stest[j] + sum( cp < thr[j] )
    }
    object@out <-list(out=out,stest=stest,thr=thr,cprobs=cprobs)

    ### cluster outliers
    clp2p.cl <- c()
    cols <- names(object@fdata)
    cpart <- object@cluster$kpart
    di   <- as.data.frame(object@distances)
    for ( i in 1:max(cpart) ) {
      tcol <- cols[cpart == i]
      if ( sum(!(tcol %in% out)) > 1 ) clp2p.cl <- append(clp2p.cl,as.vector(t(di[tcol[!(tcol %in% out)],tcol[!(tcol %in% out)]])))
    }
    clp2p.cl <- clp2p.cl[clp2p.cl>0]

    cadd  <- list()
    if ( length(out) > 0 ){
      if (length(out) == 1){
        cadd <- list(out)
      }else{
        n <- out
        m <- as.data.frame(di[out,out])

        for ( i in 1:length(out) ){
          if ( length(n) > 1 ){
            o   <- order(apply(cbind(m,1:dim(m)[1]),1,function(x)  min(x[1:(length(x)-1)][-x[length(x)]])),decreasing=FALSE)
            m <- m[o,o]
            n <- n[o]
            f <- m[,1] < quantile(clp2p.cl,outdistquant) | m[,1] == min(clp2p.cl)
            ind <- 1
            if ( sum(f) > 1 ) for ( j in 2:sum(f) ) if ( apply(m[f,f][j,c(ind,j)] > quantile(clp2p.cl,outdistquant) ,1,sum) == 0 ) ind <- append(ind,j)
            cadd[[i]] <- n[f][ind]
            g <- ! n %in% n[f][ind]
            n <- n[g]
            m <- m[g,g]
            if ( sum(g) == 0 ) break

          }else if (length(n) == 1){
            cadd[[i]] <- n
            break
          }
        }
      }

      for ( i in 1:length(cadd) ){
        cpart[cols %in% cadd[[i]]] <- max(cpart) + 1
      }
    }

    ### determine final clusters
    for ( i in 1:max(cpart) ){
      if ( sum(cpart == i) == 0 ) next
      f <- cols[cpart == i]
      d <- object@fdata
      if ( length(f) == 1 ){
        cent <- d[,f]
      }else{
        if ( object@clusterpar$FUNcluster == "kmedoids" ){
          x <- apply(as.matrix(dist.gen(t(d[,f]),method=object@clusterpar$metric)),2,mean)
          cent <- d[,f[which(x == min(x))[1]]]
        }else{
          cent <- apply(d[,f],1,mean)
        }
      }
      if ( i == 1 ) dcent <- data.frame(cent) else dcent <- cbind(dcent,cent)
      if ( i == 1 ) tmp <- data.frame(apply(d,2,dist.gen.pairs,y=cent,method=object@clusterpar$metric)) else tmp <- cbind(tmp,apply(d,2,dist.gen.pairs,y=cent,method=object@clusterpar$metric))
    }
    cpart <- apply(tmp,1,function(x) order(x,decreasing=FALSE)[1])

    for  ( i in max(cpart):1){if (sum(cpart==i)==0) cpart[cpart>i] <- cpart[cpart>i] - 1 }

    object@cpart <- cpart

    # set.seed(111111) # dyneval: don't set a seed
    object@fcol <- sample(rainbow(max(cpart)))
    return(object)
  }
)

#' Method comptsne
#' @name comptsne
#' @rdname comptsne-methods
#' @exportMethod comptsne
setGeneric(
  "comptsne",
  function(object,rseed=15555,sammonmap=FALSE,initial_cmd=TRUE,...) standardGeneric("comptsne")
)

#' @rdname comptsne-methods
#' @aliases comptsne,SCseq-method
#' @importFrom MASS sammon
#' @importFrom tsne tsne
setMethod(
  "comptsne",
  signature = "SCseq",
  definition = function(object,rseed,sammonmap,...){
    if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before comptsne")
    # set.seed(rseed) # dyneval: don't set a seed
    di <- if ( object@clusterpar$FUNcluster == "kmedoids") as.dist(object@distances) else dist.gen(as.matrix(object@distances))
    if ( sammonmap ){
      object@tsne <- as.data.frame(sammon(di,k=2)$points)
    }else{
      ts <- if ( initial_cmd ) tsne(di,initial_config=cmdscale(di,k=2),...) else tsne(di,k=2,...)
      object@tsne <- as.data.frame(ts)
    }
    return(object)
  }
)

setGeneric("clustdiffgenes", function(object,pvalue=.01) standardGeneric("clustdiffgenes"))

setMethod(
  "clustdiffgenes",
  signature = "SCseq",
  definition = function(object,pvalue){
    if ( length(object@cpart) == 0 ) stop("run findoutliers before clustdiffgenes")
    if ( ! is.numeric(pvalue) ) stop("pvalue has to be a number between 0 and 1") else if (  pvalue < 0 | pvalue > 1 ) stop("pvalue has to be a number between 0 and 1")
    cdiff <- list()
    x     <- object@ndata
    y     <- object@expdata[,names(object@ndata)]
    part  <- object@cpart
    for ( i in 1:max(part) ){
      if ( sum(part == i) == 0 ) next
      m <-  if ( sum(part != i) > 1 ) apply(x[,part != i],1,mean) else x[,part != i]
      n <-  if ( sum(part == i) > 1 ) apply(x[,part == i],1,mean) else x[,part == i]
      no <- if ( sum(part == i) > 1 ) median(apply(y[,part == i],2,sum))/median(apply(x[,part == i],2,sum)) else sum(y[,part == i])/sum(x[,part == i])
      m <- m*no
      n <- n*no
      pv <- binompval(m/sum(m),sum(n),n)
      d <- data.frame(mean.ncl=m,mean.cl=n,fc=n/m,pv=pv)[order(pv,decreasing=FALSE),]
      cdiff[[paste("cl",i,sep=".")]] <- d[d$pv < pvalue,]
    }
    return(cdiff)
  }
)

setGeneric("plotsaturation", function(object,disp=FALSE) standardGeneric("plotsaturation"))

setMethod(
  "plotsaturation",
  signature = "SCseq",
  definition = function(object,disp){
    if ( length(object@cluster$gap) == 0 ) stop("run clustexp before plotsaturation")
    g <- object@cluster$gap$Tab[,1]
    y <- g[-length(g)] - g[-1]
    mm <- numeric(length(y))
    nn <- numeric(length(y))
    for ( i in 1:length(y)){
      mm[i] <- mean(y[i:length(y)])
      nn[i] <- sqrt(var(y[i:length(y)]))
    }
    cln <- max(min(which( y - (mm + nn) < 0 )),1)
    x <- 1:length(y)
    if (disp){
      x <- 1:length(g)
      plot(x,g,pch=20,col="grey",xlab="k",ylab="log within cluster dispersion")
      f <- x == cln
      points(x[f],g[f],col="blue")
    }else{
      plot(x,y,pch=20,col="grey",xlab="k",ylab="Change in log within cluster dispersion")
      points(x,mm,col="red",pch=20)
      plot.err.bars.y(x,mm,nn,col="red")
      points(x,y,col="grey",pch=20)
      f <- x == cln
      points(x[f],y[f],col="blue")
    }
  }
)

setGeneric("plotsilhouette", function(object) standardGeneric("plotsilhouette"))

#' @importFrom cluster silhouette
setMethod(
  "plotsilhouette",
  signature = "SCseq",
  definition = function(object){
    if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before plotsilhouette")
    if ( length(unique(object@cluster$kpart)) < 2 ) stop("only a single cluster: no silhouette plot")
    kpart <- object@cluster$kpart
    distances  <- if ( object@clusterpar$FUNcluster == "kmedoids" ) as.dist(object@distances) else dist.gen(object@distances)
    si <- silhouette(kpart,distances)
    plot(si)
  }
)

compmedoids <- function(x,part,metric="pearson"){
  m <- c()
  for ( i in sort(unique(part)) ){
    f <- names(x)[part == i]
    if ( length(f) == 1 ){
      m <- append(m,f)
    }else{
      y <- apply(as.matrix(dist.gen(t(x[,f]),method=metric)),2,mean)
      m <- append(m,f[which(y == min(y))[1]])
    }
  }
  m
}

setGeneric("clustheatmap", function(object,final=FALSE,hmethod="single") standardGeneric("clustheatmap"))

setMethod(
  "clustheatmap",
  signature = "SCseq",
  definition = function(object,final,hmethod){
    if ( final & length(object@cpart) == 0 ) stop("run findoutliers before clustheatmap")
    if ( !final & length(object@cluster$kpart) == 0 ) stop("run clustexp before clustheatmap")
    x <- object@fdata
    part <- if ( final ) object@cpart else object@cluster$kpart
    na <- c()
    j <- 0
    for ( i in 1:max(part) ){
      if ( sum(part == i) == 0 ) next
      j <- j + 1
      na <- append(na,i)
      d <- x[,part == i]
      if ( sum(part == i) == 1 ) cent <- d else cent <- apply(d,1,mean)
      if ( j == 1 ) tmp <- data.frame(cent) else tmp <- cbind(tmp,cent)
    }
    names(tmp) <- paste("cl",na,sep=".")
    ld <- if ( object@clusterpar$FUNcluster == "kmedoids" ) dist.gen(t(tmp),method=object@clusterpar$metric) else dist.gen(as.matrix(dist.gen(t(tmp),method=object@clusterpar$metric)))
    if ( max(part) > 1 )  cclmo <- hclust(ld,method=hmethod)$order else cclmo <- 1
    q <- part
    for ( i in 1:max(part) ){
      q[part == na[cclmo[i]]] <- i
    }
    part <- q
    di <-  if ( object@clusterpar$FUNcluster == "kmedoids" ) object@distances else as.data.frame( as.matrix( dist.gen(t(object@distances)) ) )
    pto <- part[order(part,decreasing=FALSE)]
    ptn <- c()
    for ( i in 1:max(pto) ){ pt <- names(pto)[pto == i]; z <- if ( length(pt) == 1 ) pt else pt[hclust(as.dist(t(di[pt,pt])),method=hmethod)$order]; ptn <- append(ptn,z) }
    col <- object@fcol
    mi  <- min(di,na.rm=TRUE)
    ma  <- max(di,na.rm=TRUE)
    layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
    ColorRamp   <- colorRampPalette(brewer.pal(n = 7,name = "RdYlBu"))(100)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    if ( mi == ma ){
      ColorLevels <- seq(0.99*mi, 1.01*ma, length=length(ColorRamp))
    }
    par(mar = c(3,5,2.5,2))
    image(as.matrix(di[ptn,ptn]),col=ColorRamp,axes=FALSE)
    abline(0,1)
    box()

    tmp <- c()
    for ( u in 1:max(part) ){
      ol <- (0:(length(part) - 1)/(length(part) - 1))[ptn %in% names(x)[part == u]]
      points(rep(0,length(ol)),ol,col=col[cclmo[u]],pch=15,cex=.75)
      points(ol,rep(0,length(ol)),col=col[cclmo[u]],pch=15,cex=.75)
      tmp <- append(tmp,mean(ol))
    }
    axis(1,at=tmp,labels = cclmo)
    axis(2,at=tmp,labels = cclmo)
    par(mar = c(3,2.5,2.5,2))
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
    layout(1)
    return(cclmo)
  }
)






