#' @export Ltree
#' @exportClass Ltree
#' @include stemid.R
Ltree <- setClass(
  "Ltree",
  slots = c(
    sc = "SCseq",
    ldata = "list",
    entropy = "vector",
    trproj = "list",
    par = "list",
    prback = "data.frame",
    prbacka = "data.frame",
    ltcoord = "matrix",
    prtree = "list",
    sigcell = "vector",
    cdata = "list",
    trl = "list"
  )
)

setValidity(
  "Ltree",
  function(object) {
    msg <- NULL
    if ( class(object@sc)[1] != "SCseq" ){
      msg <- c(msg, "input data must be of class SCseq")
    }
    if (is.null(msg)) TRUE
    else msg
  }
)

#' @importFrom methods validObject new
setMethod(
  "initialize",
  signature = "Ltree",
  definition = function(.Object, sc ){
    .Object@sc <- sc
    validObject(.Object)
    return(.Object)
  }
)

#' @export
setGeneric(
  "compentropy",
  function(object)
    standardGeneric("compentropy")
)

setMethod(
  "compentropy",
  signature = "Ltree",
  definition = function(object){
    probs   <- t(t(object@sc@ndata)/apply(object@sc@ndata,2,sum))
    object@entropy <- -apply(probs*log(probs)/log(nrow(object@sc@ndata)),2,sum)
    return(object)
  }
)

compproj <- function(pdiloc,lploc,cnloc,mloc,d=NULL){
  pd    <- data.frame(pdiloc)
  k     <- paste("X",sort(rep(1:nrow(pdiloc),length(mloc))),sep="")
  pd$k  <- paste("X",1:nrow(pdiloc),sep="")
  pd    <- merge(data.frame(k=k),pd,by="k")

  if ( is.null(d) ){
    cnv   <- t(matrix(rep(t(cnloc),nrow(pdiloc)),nrow=ncol(pdiloc)))
    pdcl  <- paste("X",lploc[as.numeric(sub("X","",pd$k))],sep="")
    rownames(cnloc) <- paste("X",mloc,sep="")
    pdcn  <- cnloc[pdcl,]
    v     <- cnv - pdcn
  }else{
    v    <- d$v
    pdcn <- d$pdcn
  }
  w <- pd[,names(pd) != "k"] - pdcn

  h <- apply(cbind(v,w),1,function(x){
    x1 <- x[1:(length(x)/2)];
    x2 <- x[(length(x)/2 + 1):length(x)];
    x1s <- sqrt(sum(x1**2)); x2s <- sqrt(sum(x2**2)); y <- sum(x1*x2)/x1s/x2s; return( if (x1s == 0 | x2s == 0 ) NA else y ) })

  rma <- as.data.frame(matrix(h,ncol=nrow(pdiloc)))
  names(rma) <- unique(pd$k)
  pdclu  <- lploc[as.numeric(sub("X","",names(rma)))]
  pdclp  <- apply(t(rma),1,function(x) if (sum(!is.na(x)) == 0 ) NA else mloc[which(abs(x) == max(abs(x),na.rm=TRUE))[1]])
  pdclh  <- apply(t(rma),1,function(x) if (sum(!is.na(x)) == 0 ) NA else x[which(abs(x) == max(abs(x),na.rm=TRUE))[1]])
  pdcln  <-  names(lploc)[as.numeric(sub("X","",names(rma)))]
  names(rma) <- pdcln
  rownames(rma) <- paste("X",mloc,sep="")
  res    <- data.frame(o=pdclu,l=pdclp,h=pdclh)
  rownames(res) <- pdcln
  return(list(res=res[names(lploc),],rma=as.data.frame(t(rma[,names(lploc)])),d=list(v=v,pdcn=pdcn)))
}

pdishuffle <- function(pdi,lp,cn,m,all=FALSE){
  if ( all ){
    d <- as.data.frame(pdi)
    y <- t(apply(pdi,1,function(x) runif(length(x),min=min(x),max=max(x))))
    names(y)    <- names(d)
    rownames(y) <- rownames(d)
    return(y)
  }else{
    fl <- TRUE
    for ( i in unique(lp)){
      if ( sum(lp == i) > 1 ){
        y <-data.frame( t(apply(as.data.frame(pdi[,lp == i]),1,sample)) )
      }else{
        y <- as.data.frame(pdi[,lp == i])
      }
      names(y) <- names(lp)[lp == i]
      rownames(y) <- names(lp)
      z <- if (fl) y else cbind(z,y)
      fl <- FALSE
    }
    z <- t(z[,names(lp)])
    return(z)
  }
}

#' @export
setGeneric(
  "projcells",
  function(object,cthr=0,nmode=FALSE)
    standardGeneric("projcells")
)

setMethod(
  "projcells",
  signature = "Ltree",
  definition = function(object,cthr,nmode){
    if ( ! is.numeric(cthr) ) stop( "cthr has to be a non-negative number" ) else if ( cthr < 0 ) stop( "cthr has to be a non-negative number" )
    if ( ! length(object@sc@cpart == 0) ) stop( "please run findoutliers on the SCseq input object before initializing Ltree" )
    if ( !is.numeric(nmode) & !is.logical(nmode) ) stop("argument nmode has to be logical (TRUE/FALSE)")

    object@par$cthr  <- cthr
    object@par$nmode <- nmode

    lp <- object@sc@cpart
    ld <- object@sc@distances
    n  <- aggregate(rep(1,length(lp)),list(lp),sum)
    n  <- as.vector(n[order(n[,1],decreasing=FALSE),-1])
    m  <- (1:length(n))[n>cthr]
    f  <- lp %in% m
    lp <- lp[f]
    ld <- ld[f,f]

    pdil <- object@sc@tsne[f,]
    cnl  <- aggregate(pdil,by=list(lp),median)
    cnl  <- cnl[order(cnl[,1],decreasing=FALSE),-1]

    pdi <- suppressWarnings( cmdscale(as.matrix(ld),k=ncol(ld)-1) )
    cn <- as.data.frame(pdi[compmedoids(object@sc@fdata[,names(lp)],lp),])
    rownames(cn) <- 1:nrow(cn)

    x <- compproj(pdi,lp,cn,m)
    res <- x$res

    if ( nmode ){
      rma <- x$rma
      z <- paste("X",t(as.vector(apply(cbind(lp,ld),1,function(x){ f <- lp != x[1]; lp[f][which(x[-1][f] == min(x[-1][f]))[1]] }))),sep="")
      k <- apply(cbind(z,rma),1,function(x) (x[-1])[names(rma) == x[1]])
      rn <- res
      rn$l <- as.numeric(sub("X","",z))
      rn$h <- as.numeric(k)
      res <- rn
      x$res <- res
    }

    object@ldata  <- list(lp=lp,ld=ld,m=m,pdi=pdi,pdil=pdil,cn=cn,cnl=cnl)
    object@trproj <- x
    return(object)
  }
)




#' @export
setGeneric(
  "projback",
  function(object,pdishuf=2000,nmode=FALSE,rseed=17000)
    standardGeneric("projback")
)

setMethod(
  "projback",
  signature = "Ltree",
  definition = function(object,pdishuf,nmode,rseed){
    if ( !is.numeric(nmode) & !is.logical(nmode) ) stop("argument nmode has to be logical (TRUE/FALSE)")
    if ( ! is.numeric(pdishuf) ) stop( "pdishuf has to be a non-negative integer number" ) else if ( round(pdishuf) != pdishuf | pdishuf < 0 ) stop( "pdishuf has to be a non-negative integer number" )
    if ( length(object@trproj) == 0 ) stop("run projcells before projback")
    object@par$pdishuf  <- pdishuf
    object@par$rseed    <- rseed

    if ( ! nmode ){
      # set.seed(rseed) # dyneval: don't set a seed
      for ( i in 1:pdishuf ){
        cat("pdishuffle:",i,"\n")
        x <- compproj(pdishuffle(object@ldata$pdi,object@ldata$lp,object@ldata$cn,object@ldata$m,all=TRUE),object@ldata$lp,object@ldata$cn,object@ldata$m,d=object@trproj$d)$res
        y <- if ( i == 1 ) t(x) else cbind(y,t(x))
      }
      ##important
      object@prback <- as.data.frame(t(y))

      x <- object@prback
      x$n <- as.vector(t(matrix(rep(1:(nrow(x)/nrow(object@ldata$pdi)),nrow(object@ldata$pdi)),ncol=nrow(object@ldata$pdi))))
      object@prbacka <- aggregate(data.frame(count=rep(1,nrow(x))),by=list(n=x$n,o=x$o,l=x$l),sum)
    }
    return( object )
  }
)

#' @export
setGeneric(
  "lineagetree",
  function(object,pthr=0.01,nmode=FALSE)
    standardGeneric("lineagetree")
)

setMethod(
  "lineagetree",
  signature = "Ltree",
  definition = function(object,pthr,nmode){
    if ( !is.numeric(nmode) & !is.logical(nmode) ) stop("argument nmode has to be logical (TRUE/FALSE)")
    if ( length(object@trproj) == 0 ) stop("run projcells before lineagetree")
    if ( max(dim(object@prback)) == 0 & ! nmode ) stop("run projback before lineagetree")
    if ( ! is.numeric(pthr) ) stop( "pthr has to be a non-negative number" ) else if ( pthr < 0 ) stop( "pthr has to be a non-negative number" )
    object@par$pthr <- pthr
    cnl    <- object@ldata$cnl
    pdil   <- object@ldata$pdil
    cn    <- object@ldata$cn
    pdi   <- object@ldata$pdi
    m      <- object@ldata$m
    lp     <- object@ldata$lp
    res    <- object@trproj$res
    rma    <- object@trproj$rma
    prback <- object@prback

    cm <- as.matrix(dist(cnl))*0
    linl <- list()
    linn <- list()
    for ( i in 1:length(m) ){
      for ( j in i:length(m) ){
        linl[[paste(m[i],m[j],sep=".")]] <- c()
        linn[[paste(m[i],m[j],sep=".")]] <- c()
      }
    }
    sigcell <- c()
    for ( i in 1:nrow(res) ){
      a <- which( m == res$o[i])
      if ( sum( lp == m[a] ) == 1 ){
        k <- t(cnl)[,a]
        k <- NA
        sigcell <- append(sigcell, FALSE)
      }else{
        b <- which(m == res$l[i] )
        h <- res$h[i]
        if ( nmode ){
          sigcell <- append(sigcell, FALSE)
        }else{
          f <- prback$o == m[a] & prback$l == m[b]
          if ( sum(f) > 0 ){
            ql <- quantile(prback[f,"h"],pthr,na.rm=TRUE)
            qh <- quantile(prback[f,"h"],1 - pthr,na.rm=TRUE)
          }else{
            ql <- 0
            qh <- 0
          }
          sigcell <- if (is.na(h) ) append(sigcell, NA) else if ( h > qh |  h < min(0,ql) ) append(sigcell, TRUE) else append(sigcell, FALSE)
        }
        if ( !is.na(res$h[i]) ){
          w <- t(pdil)[,i] - t(cnl)[,a]
          v <- t(cnl)[,b] - t(cnl)[,a]

          wo <- t(pdi)[,i] - t(cn)[,a]
          vo <-  t(cn)[,b] - t(cn)[,a]
          df <-( h*sqrt(sum(wo*wo)) )/sqrt(sum(vo*vo))*v
          k <- df + t(cnl)[,a]
          cm[a,b] <- cm[a,b] + 1
          so <- m[sort(c(a,b))]
          dfl <-  sign(h)*sqrt( sum( df*df ) )/sqrt(sum(v*v))
          if ( a > b ) dfl <-  1 - dfl
          linn[[paste(so[1],so[2],sep=".")]] <- append( linn[[paste(so[1],so[2],sep=".")]], rownames(pdi)[i] )
          linl[[paste(so[1],so[2],sep=".")]] <- append( linl[[paste(so[1],so[2],sep=".")]], dfl )
        }else{
          k <- t(cnl)[,a]
          for ( j in unique(lp[lp != m[a]]) ){
            b <- which(j == m)
            so <- m[sort(c(a,b))]
            dfl <- 0
            if ( a > b ) dfl <-  1 - dfl
            linn[[paste(so[1],so[2],sep=".")]] <- append( linn[[paste(so[1],so[2],sep=".")]], rownames(pdi)[i] )
            linl[[paste(so[1],so[2],sep=".")]] <- append( linl[[paste(so[1],so[2],sep=".")]], dfl )
          }
        }
      }
      lt <- if ( i == 1 ) data.frame(k) else cbind(lt,k)
    }
    lt <- t(lt)
    cm <- as.data.frame(cm)
    names(cm) <- paste("cl",m,sep=".")
    rownames(cm) <- paste("cl",m,sep=".")
    lt <- as.data.frame(lt)
    rownames(lt) <- rownames(res)
    object@ltcoord <- as.matrix(lt)
    object@prtree  <- list(n=linn,l=linl)
    object@cdata$counts <- cm
    names(sigcell) <- rownames(res)
    object@sigcell <- sigcell

    return( object )
  }
)


setGeneric("comppvalue", function(object,pethr=0.01,nmode=FALSE) standardGeneric("comppvalue"))

setMethod(
  "comppvalue",
  signature = "Ltree",
  definition = function(object,pethr,nmode){
    if ( !is.numeric(nmode) & !is.logical(nmode) ) stop("argument nmode has to be logical (TRUE/FALSE)")
    if ( length(object@prtree) == 0 ) stop("run lineagetree before comppvalue")
    if ( ! is.numeric(pethr) ) stop( "pethr has to be a non-negative number" ) else if ( pethr < 0 ) stop( "pethr has to be a non-negative number" )
    object@par$pethr <- pethr
    cm <- object@cdata$counts
    cmpv   <- cm*NA
    cmpvd  <- cm*NA
    cmbr   <- cm*NA
    cmpvn  <- cm*NA
    cmpvnd <- cm*NA
    cmfr   <- cm/apply(cm,1,sum)
    if ( nmode ){
      N <- apply(cm,1,sum) + 1
      N0 <- sum(N) - N
      n0 <- t(matrix(rep(N,length(N)),ncol=length(N)))
      p <- n0/N0
      n <- cm
      k <- cbind(N,p,n)
      cmpv   <- apply(k,1,function(x){N <- x[1]; p <- x[2:( ncol(cm) + 1 )];  n <- x[( ncol(cm) + 2 ):( 2*ncol(cm) + 1)]; apply(cbind(n,p),1,function(x,N) binom.test(x[1],N,min(1,x[2]),alternative="g")$p.value,N=N)})
      cmpvd   <- apply(k,1,function(x){N <- x[1]; p <- x[2:( ncol(cm) + 1 )];  n <- x[( ncol(cm) + 2 ):( 2*ncol(cm) + 1)]; apply(cbind(n,p),1,function(x,N) binom.test(x[1],N,min(1,x[2]),alternative="l")$p.value,N=N)})
      cmpvn  <- cmpv
      cmpvnd <- cmpvd
      cmbr   <- as.data.frame(n0/N0*N)
      names(cmbr)    <- names(cm)
      rownames(cmbr) <- rownames(cm)
    }else{
      for ( i in 1:nrow(cm) ){
        for ( j in 1:ncol(cm) ){
          f <- object@prbacka$o == object@ldata$m[i] & object@prbacka$l == object@ldata$m[j]
          x <- object@prbacka$count[f]
          if ( sum(f) < object@par$pdishuf ) x <- append(x,rep(0, object@par$pdishuf - sum(f)))
          cmbr[i,j]   <- if ( sum(f) > 0 ) mean(x) else 0
          cmpv[i,j]   <- if ( quantile(x,1 - pethr) < cm[i,j] ) 0 else 0.5
          cmpvd[i,j]  <- if ( quantile(x,pethr) > cm[i,j] ) 0 else 0.5
          cmpvn[i,j]  <- sum( x >= cm[i,j])/length(x)
          cmpvnd[i,j] <- sum( x <= cm[i,j])/length(x)
        }
      }
    }

    diag(cmpv)   <- .5
    diag(cmpvd)  <- .5
    diag(cmpvn)  <- NA
    diag(cmpvnd) <- NA

    object@cdata$counts.br <- cmbr
    object@cdata$pv.e <- cmpv
    object@cdata$pv.d <- cmpvd
    object@cdata$pvn.e <- cmpvn
    object@cdata$pvn.d <- cmpvnd

    m    <- object@ldata$m
    linl <- object@prtree$l
    ls   <- as.data.frame(matrix(rep(NA,length(m)**2),ncol=length(m)))
    names(ls) <- rownames(ls) <- paste("cl",m,sep=".")
    for ( i in 1:( length(m) - 1 )){
      for ( j in (i + 1):length(m) ){
        na <- paste(m[i],m[j],sep=".")
        if ( na %in% names(linl) &  min(cmpv[i,j],cmpv[j,i],na.rm=TRUE) < pethr ){
          y <- sort(linl[[na]])
          nn <- ( 1 - max(y[-1] - y[-length(y)]) )
        }else{
          nn <- 0
        }
        ls[i,j] <- nn
      }
    }
    object@cdata$linkscore <- ls

    return(object)
  }
)

setGeneric("plotlinkpv", function(object) standardGeneric("plotlinkpv"))

setMethod(
  "plotlinkpv",
  signature = "Ltree",
  definition = function(object){
    if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotlinkpv")
    pheatmap(-log2(object@cdata$pvn.e + 1/object@par$pdishuf/2))
  }
)

setGeneric("plotlinkscore", function(object) standardGeneric("plotlinkscore"))

setMethod(
  "plotlinkscore",
  signature = "Ltree",
  definition = function(object){
    if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotlinkscore")
    pheatmap(object@cdata$linkscore,cluster_rows=FALSE,cluster_cols=FALSE)
  }
)

setGeneric("plotmapprojections", function(object) standardGeneric("plotmapprojections"))

#' @importFrom vegan spantree
setMethod(
  "plotmapprojections",
  signature = "Ltree",
  definition = function(object){
    if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotmapprojections")

    cent <- object@sc@fdata[,compmedoids(object@sc@fdata,object@sc@cpart)]
    dc <- as.data.frame(1 - cor(cent))
    names(dc) <- sort(unique(object@sc@cpart))
    rownames(dc) <- sort(unique(object@sc@cpart))
    trl <- spantree(dc[object@ldata$m,object@ldata$m])

    u <- object@ltcoord[,1]
    v <- object@ltcoord[,2]
    cnl <- object@ldata$cnl
    plot(u,v,cex=1.5,col="grey",pch=20,xlab="Dim 1",ylab="Dim 2")
    for ( i in unique(object@ldata$lp) ){ f <- object@ldata$lp == i; text(u[f],v[f],i,cex=.75,font=4,col=object@sc@fcol[i]) }
    points(cnl[,1],cnl[,2])
    text(cnl[,1],cnl[,2],object@ldata$m,cex=2)
    for ( i in 1:length(trl$kid) ){
      lines(c(cnl[i+1,1],cnl[trl$kid[i],1]),c(cnl[i+1,2],cnl[trl$kid[i],2]),col="black")
    }
  }
)

#' @export
setGeneric(
  "compspantree",
  function(object)
    standardGeneric("compspantree")
)

#' @importFrom vegan spantree
setMethod(
  "compspantree",
  signature = "Ltree",
  definition = function(object) {
    medoids <- compmedoids(object@sc@fdata, object@sc@cpart)
    cent <- object@sc@fdata[,medoids]
    dc <- as.data.frame(1 - cor(cent))
    names(dc) <- sort(unique(object@sc@cpart))
    rownames(dc) <- sort(unique(object@sc@cpart))
    trl <- vegan::spantree(dc[object@ldata$m,object@ldata$m])
    object@trl <- list(
      dc = dc,
      cent = cent,
      trl = trl
    )
    object
  }
)

setGeneric(
  "plotmap",
  function(object)
    standardGeneric("plotmap")
)

#' @importFrom vegan spantree
setMethod(
  "plotmap",
  signature = "Ltree",
  definition = function(object){
    if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotmap")

    cent <- object@sc@fdata[,compmedoids(object@sc@fdata,object@sc@cpart)]
    dc <- as.data.frame(1 - cor(cent))
    names(dc) <- sort(unique(object@sc@cpart))
    rownames(dc) <- sort(unique(object@sc@cpart))
    trl <- spantree(dc[object@ldata$m,object@ldata$m])


    u <- object@ldata$pdil[,1]
    v <- object@ldata$pdil[,2]
    cnl <- object@ldata$cnl
    plot(u,v,cex=1.5,col="grey",pch=20,xlab="Dim 1",ylab="Dim 2")
    for ( i in unique(object@ldata$lp) ){ f <- object@ldata$lp == i; text(u[f],v[f],i,cex=.75,font=4,col=object@sc@fcol[i]) }
    points(cnl[,1],cnl[,2])
    text(cnl[,1],cnl[,2],object@ldata$m,cex=2)
    for ( i in 1:length(trl$kid) ){
      lines(c(cnl[i+1,1],cnl[trl$kid[i],1]),c(cnl[i+1,2],cnl[trl$kid[i],2]),col="black")
    }
  }
)

setGeneric(
  "plottree",
  function(object,showCells=TRUE,nmode=FALSE,scthr=0)
    standardGeneric("plottree")
)

setMethod(
  "plottree",
  signature = "Ltree",
  definition = function(object,showCells,nmode,scthr){
    if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotmap")
    if ( !is.numeric(nmode) & !is.logical(nmode) ) stop("argument nmode has to be logical (TRUE/FALSE)")
    if ( !is.numeric(showCells) & !is.logical(showCells) ) stop("argument showCells has to be logical (TRUE/FALSE)")
    if ( ! is.numeric(scthr) ) stop( "scthr has to be a non-negative number" ) else if ( scthr < 0 | scthr > 1 ) stop( "scthr has to be a number between 0 and 1" )


    ramp <- colorRamp(c("darkgreen", "yellow", "brown"))
    mcol <- rgb( ramp(seq(0, 1, length = 101)), max = 255)
    co <- object@cdata
    fc <- (co$counts/( co$counts.br + .5))*(co$pv.e < object@par$pethr)
    fc <- fc*(fc > t(fc)) + t(fc)*(t(fc) >= fc)
    fc <- log2(fc + (fc == 0))

    k <- -log10(sort(unique(as.vector(t(co$pvn.e))[as.vector(t(co$pv.e))<object@par$pethr])) + 1/object@par$pdishuf)
    if (length(k) == 1) k <- c(k - k/100,k)
    mlpv <- -log10(co$pvn.e + 1/object@par$pdishuf)
    diag(mlpv) <- min(mlpv,na.rm=TRUE)
    dcc <- t(apply(round(100*(mlpv - min(k))/(max(k) - min(k)),0) + 1,1,function(x){y <- c(); for ( n in x ) y <- append(y,if ( n < 1 ) NA else mcol[n]); y }))


    cx <- c()
    cy <- c()
    va <- c()
    m <- object@ldata$m
    for ( i in 1:( length(m) - 1 ) ){
      for ( j in ( i + 1 ):length(m) ){
        if ( min(co$pv.e[i,j],co$pv.e[j,i],na.rm=TRUE) < object@par$pethr ){
          if ( mlpv[i,j] > mlpv[j,i] ){
            va <- append(va,dcc[i,j])
          }else{
            va <- append(va,dcc[j,i])
          }
          cx <- append(cx,i)
          cy <- append(cy,j)
        }
      }
    }

    cnl <- object@ldata$cnl
    u <- object@ltcoord[,1]
    v <- object@ltcoord[,2]
    layout( cbind(c(1, 1), c(2, 3)),widths=c(5,1,1),height=c(5,5,1))
    par(mar = c(12,5,1,1))

    if ( showCells ){
      plot(u,v,cex=1.5,col="grey",pch=20,xlab="Dim 1",ylab="Dim 2")
      if ( !nmode ) points(u[object@sigcell],v[object@sigcell],col="black")
    }else{
      plot(u,v,cex=0,xlab="Dim 1",ylab="Dim 2")
    }

    if ( length(va) > 0 ){
      f <- order(va,decreasing=TRUE)
      for ( i in 1:length(va) ){
        if ( object@cdata$linkscore[cx[i],cy[i]] > scthr ){
          if ( showCells ){
            lines(cnl[c(cx[i],cy[i]),1],cnl[c(cx[i],cy[i]),2],col=va[i],lwd=2)
          }else{
            ##nn <- min(10,fc[cx[i],cy[i]])
            lines(cnl[c(cx[i],cy[i]),1],cnl[c(cx[i],cy[i]),2],col=va[i],lwd=5*object@cdata$linkscore[cx[i],cy[i]])
          }
        }
      }
    }

    en <- aggregate(object@entropy,list(object@sc@cpart),median)
    en <- en[en$Group.1 %in% m,]

    mi <- min(en[,2],na.rm=TRUE)
    ma <- max(en[,2],na.rm=TRUE)
    w <- round((en[,2] - mi)/(ma - mi)*99 + 1,0)
    ramp <- colorRamp(c("red","orange", "pink","purple", "blue"))
    ColorRamp <- rgb( ramp(seq(0, 1, length = 101)), max = 255)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    if ( mi == ma ){
      ColorLevels <- seq(0.99*mi, 1.01*ma, length=length(ColorRamp))
    }
    for ( i in m ){
      f <- en[,1] == m
      points(cnl[f,1],cnl[f,2],cex=5,col=ColorRamp[w[f]],pch=20)
    }
    text(cnl[,1],cnl[,2],m,cex=1.25,font=4,col="white")
    par(mar = c(5, 4, 1, 2))
    image(1, ColorLevels,
          matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
          col=ColorRamp,
          xlab="",ylab="",
          xaxt="n")
    coll <- seq(min(k), max(k), length=length(mcol))
    image(1, coll,
          matrix(data=coll, ncol=length(mcol),nrow=1),
          col=mcol,
          xlab="",ylab="",
          xaxt="n")
    layout(1)
  }
)


setGeneric("plotdistanceratio", function(object) standardGeneric("plotdistanceratio"))

setMethod(
  "plotdistanceratio",
  signature = "Ltree",
  definition = function(object){
    if ( length(object@ldata) <= 0 ) stop("run projcells before plotdistanceratio")
    l <- as.matrix(dist(object@ldata$pdi))
    z <- (l/object@ldata$ld)
    hist(log2(z),breaks=100,xlab=" log2 emb. distance/distance",main="")
  }
)


setGeneric("getproj", function(object,i) standardGeneric("getproj"))

setMethod(
  "getproj",
  signature = "Ltree",
  definition = function(object,i){
    if ( length(object@ldata) <= 0 ) stop("run projcells before plotdistanceratio")
    if ( ! i %in% object@ldata$m )  stop(paste("argument i has to be one of",paste(object@ldata$m,collapse=",")))
    x <- object@trproj$rma[names(object@ldata$lp)[object@ldata$lp == i],]
    x <- x[,names(x) != paste("X",i,sep="")]
    f <- !is.na(x[,1])
    x <- x[f,]
    if ( nrow(x) > 1 ){
      y <- x
      y <- as.data.frame(t(apply(y,1,function(x) (x - mean(x))/sqrt(var(x)))))
    }
    names(x) = sub("X","cl.",names(x))
    names(y) = sub("X","cl.",names(y))
    return(list(pr=x,prz=y))
  }
)


setGeneric("projenrichment", function(object) standardGeneric("projenrichment"))


setMethod(
  "projenrichment",
  signature = "Ltree",
  definition = function(object){
    if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotmap")

    ze <- ( object@cdata$pv.e < object@par$pethr | object@cdata$pv.d < object@par$pethr) * (object@cdata$counts + .1)/( object@cdata$counts.br + .1 )
    pheatmap(log2(ze + ( ze == 0 ) ),cluster_rows=FALSE,cluster_cols=FALSE)
  }
)





setGeneric("compscore", function(object,nn=1) standardGeneric("compscore"))

setMethod(
  "compscore",
  signature = "Ltree",
  definition = function(object,nn){
    if ( length(object@cdata) <= 0 ) stop("run comppvalue before compscore")
    if ( ! is.numeric(nn) ) stop( "nn has to be a non-negative integer number" ) else if ( round(nn) != nn | nn < 0 ) stop( "nn has to be a non-negative integer number" )
    x <- object@cdata$counts*(object@cdata$pv.e < object@par$pethr)>0
    y <- x | t(x)

    if ( max(y) > 0 ){
      z <- apply(y,1,sum)
      nl <- list()
      n <- list()
      for ( i in 1:nn ){
        if ( i == 1 ){
          n[[i]] <- as.list(apply(y,1,function(x) grep(TRUE,x)))
          nl <- data.frame( apply(y,1,sum) )
        }
        if ( i > 1 ){
          v <- rep(0,nrow(nl))
          n[[i]] <- list()
          for ( j in 1:length(n[[i-1]]) ){
            cl <- n[[i-1]][[j]]
            if ( length(cl) == 0 ){
              n[[i]][[paste("d",j,sep="")]] <- NA
              v[j] <- 0
            }else{
              k  <- if ( length(cl) > 1 ) apply(y[cl,],2,sum) > 0 else if ( length(cl) == 1 ) y[cl,]
              n[[i]][[paste("d",j,sep="")]] <- sort(unique(c(cl,grep(TRUE,k))))
              v[j] <- length(n[[i]][[paste("d",j,sep="")]])
            }
          }
          names(n[[i]]) <- names(z)
          nl <- cbind(nl,v)

        }
      }
      x <- nl[,i]
      names(x) <- rownames(nl)
    }else{
      x <- rep(0,length(object@ldata$m))
      names(x) <- paste("cl",object@ldata$m,sep=".")
    }

    v <- aggregate(object@entropy,list(object@sc@cpart),median)
    v <- v[v$Group.1 %in% object@ldata$m,]
    w <- as.vector(v[,-1])
    names(w) <- paste("cl.",v$Group.1,sep="")
    w <- w - min(w)

    return(list(links=x,entropy=w,StemIDscore=x*w))
  }
)




setGeneric("plotscore", function(object,nn=1) standardGeneric("plotscore"))

setMethod(
  "plotscore",
  signature = "Ltree",
  definition = function(object,nn){
    if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotscore")
    x <- compscore(object,nn)
    layout(1:3)
    barplot(x$links,names.arg=sub("cl\\.","",object@ldata$m),xlab="Cluster",ylab="Number of links",cex.names=1)
    barplot(x$entropy,names.arg=sub("cl\\.","",object@ldata$m),xlab="Cluster",ylab="Delta-Entropy",cex.names=1)
    barplot(x$StemIDscore,names.arg=sub("cl\\.","",object@ldata$m),xlab="Cluster",ylab="Number of links * Delta-Entropy",cex.names=1)
    layout(1)
  }
)


setGeneric("branchcells", function(object,br) standardGeneric("branchcells"))

setMethod(
  "branchcells",
  signature = "Ltree",
  definition = function(object,br){
    if ( length(object@ldata) <= 0 ) stop("run projcells before branchcells")
    msg <- paste("br needs to be list of length two containing two branches, where each has to be one of", paste(names(object@prtree$n),collapse=","))
    if ( !is.list(br) ) stop(msg) else if ( length(br) != 2 ) stop(msg) else if ( ! br[[1]] %in% names(object@prtree$n) | ! br[[2]] %in% names(object@prtree$n) ) stop(msg)


    n <- list()
    scl <- object@sc
    k <- c()
    cl <- intersect( as.numeric(strsplit(br[[1]],"\\.")[[1]]), as.numeric(strsplit(br[[2]],"\\.")[[1]]))
    if ( length(cl) == 0 ) stop("the two branches in br need to have one cluster in common.")

    for ( i in 1:length(br) ){
      f <- object@sc@cpart[ object@prtree$n[[br[[i]]]] ] %in% cl
      if ( sum(f) > 0 ){
        n[[i]] <- names(object@sc@cpart[ object@prtree$n[[br[[i]]]] ])[f]
        k <- append(k, max( scl@cpart ) + 1)
        scl@cpart[n[[i]]] <- max( scl@cpart ) + 1
      }else{
        stop(paste("no cells on branch",br[[i]],"fall into clusters",cl))
      }
    }
    # set.seed(111111)
    scl@fcol <- sample(rainbow(max(scl@cpart)))
    z <- diffgenes(scl,k[1],k[2])
    return( list(n=n,scl=scl,k=k,diffgenes=z) )
  }
)
