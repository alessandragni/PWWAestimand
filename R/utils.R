
summary.WA <- function(object,type="p",augtype=NULL,...) {# {{{
  
  rmst <- estimate(object$RAW$rmst)
  rmst.test <- estimate(rmst,contrast=rbind(c(1,-1)))
  rmst.log <- estimate(rmst,function(p) log(p))
  rmst.test.log <- estimate(rmst.log,contrast=rbind(c(1,-1)))
  
  meanNtD <- estimate(object$RAW$meanN)
  meanNtD.test <- estimate(meanNtD,contrast=rbind(c(1,-1)))
  meanNtD.log <- estimate(meanNtD,function(p) log(p))
  meanNtD.test.log <- estimate(meanNtD.log,contrast=rbind(c(1,-1)))
  
  eer <- estimate(object$RAW$ratio.means)
  eedr <- estimate(eer,contrast=rbind(c(1,-1)))
  
  if (is.null(augtype)) {
    augtype <- "riskDR"
    if (!is.null(object$ET$riskDRC)) augtype <- "riskDRC" 
  }
  if (augtype=="riskDR")
    ee <- estimate(coef=object$ET$riskDR$riskDR,vcov=object$ET$riskDR$var.riskDR)
  if (augtype=="riskDRC")
    ee <- estimate(coef=object$ET$riskDRC$coef,vcov=object$ET$riskDRC$var)
  eed <- estimate(ee,contrast=rbind(c(1,-1)))
  eer.log <- object$RAW$ratio.means.log
  eedr.log <- estimate(eer.log,contrast=rbind(c(1,-1)))
  eelog <-  estimate(ee,function(p) log(p))
  eedlog <- estimate(eelog,contrast=rbind(c(1,-1)))
  
  res <- list(rmst=rmst,rmst.test=rmst.test,meanNtD=meanNtD,meanNtD.test=meanNtD.test,
              ratio=eer,test.ratio=eedr,meanpt=ee,test.meanpt=eed)
  reslog <- list(rmst=rmst.log,rmst.test=rmst.test.log,
                 meanNtD=meanNtD.log,meanNtD.test=meanNtD.test.log,
                 ratio=eer.log,test.ratio=eedr.log,meanpt=eelog,test.meanpt=eedlog)
  class(res) <- "summary.WA"
  class(reslog) <- "summary.WA"
  attr(res,"log") <- (type!="p")
  attr(reslog,"log") <- (type!="p")
  
  if (type=="p") return(res) else return(reslog)
}# }}}

onerunN <- function(i,n,beta11=-0.3,betad1=-0.3,time=2000,cm=~1,var.z=1,cens=1/4,Yr="EpT",trans=NULL,scale1=1,scaled=1,dep=1,...) { ## {{{
  
  set.seed(seeds[i+1])
  
  nid <- n
  if (i%%500==0)	print(i)
  beta1 <- c(beta11,0.3); betad <- c(betad1,0.3); betac <- 0*c(0.3,-0.3); 
  gamma <- c(0.3,-0.3); 
  x2norm <- FALSE; 
  ###
  X <- matrix(rbinom(n*2,1,0.5),n,2)
  colnames(X) <- paste("X",1:2,sep="")
  if (x2norm)  { X2 <- rnorm(n); X[,2] <- X2; gX2 <- dcut(X2,breaks=2) } else gX2 <- X[,2]
  ###
  r1 <- exp( X %*% beta1)
  rd <- exp( X %*% betad)
  rc <- exp( X %*% betac)
  ###
  rr <- simRecurrent(nid,scalecumhaz(Lam1,scale1),death.cumhaz=scalecumhaz(LamD,scaled),
                     cens=cens,dependence=dep,var.z=var.z,r1=r1,rd=rd,rc=rc)
  rr$Z <- attr(rr,"z")[rr$id]
  rr <- cbind(rr,X[rr$id,])
  ###
  rr <- count.history(rr)
  ###
  nid <- max(rr$id)
  rr$revnr <- revcumsumstrata(rep(1,nrow(rr)),rr$id-1,nid)
  rr$statusD <- rr$status
  rr <- dtransform(rr,statusD=3,death==1 & revnr==1)
  rr$Nt_ <- rr$Count1
  ###
  rr$A <- rr$X1
  rr$x <- rr$X2
  dfactor(rr) <- A.f ~A
  ###
  ## constructing response at 2000 days
  dtable(rr,~statusD)
  rr$starty <- rr$start/365.25
  rr$stopy <- rr$stop/365.25
  rr$starty <- rr$start
  rr$stopy <- rr$stop
  ###
  
  # check simulation
  # cc1 <-  phreg(Event(starty,stopy,statusD==1)~X1+X2+cluster(id),rr)
  # ccd <-  phreg(Event(starty,stopy,statusD==3)~X1+X2+cluster(id),rr)
  
  # km1 <- phreg(Event(starty,stopy,statusD==3)~strata(A.f)+cluster(id),rr)
  # resmean.phreg(km1,times=3)
  
  outS <- WA_recurrent(Event(starty,stopy,statusD)~A.f+cluster(id),rr,trans=trans,time=time,death.code=3)
  doutS <- summary.WA(outS)
  ##
  ldoutS <- summary.WA(outS,type="log")
  
  if (dep!=0) 
  {
    outae <- WA_recurrent(Event(starty,stopy,statusD)~A.f+cluster(id),rr,augmentR=~x+Z,augmentC=~Count1+x+Z,trans=trans,time=time,death.code=3) 
  } else { 
    outae <- WA_recurrent(Event(starty,stopy,statusD)~A.f+cluster(id),rr,augmentR=~x,augmentC=~Count1+x,trans=trans,time=time,death.code=3) }
  
  dd <- outS$RAW$ratio.means$coefmat # Ratio of means E(N(min(D,t)))/E(min(D,t)) 
  alpha = 0.05
  doutae <- summary.WA(outae)
  
  #doutae$test.ratio$coefmat
  #doutae$ratio$coefmat # Ratio of means E(N(min(D,t)))/E(min(D,t)) 
  #doutae$meanpt$coefmat # Mean of Events per time-unit E(N(min(D,t))/min(D,t)) 
  
  # Here I am taking the [treat0] - [treat1] first for Mean of Events per time-unit (2 cases depending on augment R and C, i.e. outS and outae)
  # and then for the Ratio of means (only for outae)
  coef=c(doutS$test.meanpt$coefmat[1,1],doutae$test.meanpt$coefmat[1,1],doutS$test.ratio$coefmat[1,1],ldoutS$test.ratio$coefmat[1,1])
  se.coef=c( doutS$test.meanpt$coefmat[1,2],doutae$test.meanpt$coefmat[1,2],doutS$test.ratio$coefmat[1,2],ldoutS$test.ratio$coefmat[1,2]) # and se
  
  # now i am taking coefficients (treat0 and treat 1) first for  Mean of Events per time-unit then for ratio of means, 
  # first from outae then from outS
  coefe <- c(doutae$meanpt$coefmat[,1],doutae$ratio$coefmat[,1],doutS$meanpt$coefmat[,1],doutS$ratio$coefmat[,1],ldoutS$meanpt$coefmat[,1],ldoutS$ratio$coefmat[,1])
  names(coefe) <-c( rep(c("aug-mean-ratio","ratio"),each=2), rep(c("mean-ratio","ratio"),each=2), rep(c("log-mean-ratio","log-ratio"),each=2))
  se.coefe <- c(doutae$meanpt$coefmat[,2],doutae$ratio$coefmat[,2], doutS$meanpt$coefmat[,2],doutS$ratio$coefmat[,2], ldoutS$meanpt$coefmat[,2],ldoutS$ratio$coefmat[,2])
  
  dS <- summary(outS$ET$riskDR)
  dae <- summary(outae$ET$riskDR)
  coefbate <- c(outS$ET$riskDR$coef,dS$ateDR[,1],outae$ET$riskDR$coef,dae$ateDR[,1])
  se.coefbate <- c(outS$ET$riskDR$se.coef,dS$ateDR[,2],outae$ET$riskDR$se.coef,dae$ateDR[,2])
  
  
  # here i take denominators (treat0 and treat 1) and numerators
  coefb <- c(doutae$rmst$coefmat[,1],doutae$meanNtD$coefmat[,1])
  names(coefb) <- rep(c("rmst","meanN"),each=2)
  se.coefb <- c(doutae$rmst$coefmat[,2],doutae$meanNtD$coefmat[,2])
  names(se.coefb) <- rep(c("rmst","meanN"),each=2)
  
  # and putting all together 
  coefe <- c(coefe,coefb,coefbate)
  se.coefe <- c(se.coefe,se.coefb,se.coefbate)
  
  pvals = c(doutS$test.meanpt$coefmat[1,5],doutae$test.meanpt$coefmat[1,5],doutae$test.ratio$coefmat[1,5],ldoutS$test.ratio$coefmat[1,5]) < alpha
  
  names(coef) <- c("mean","Aug-mean","ratio","log-ratio")
  names(se.coef) <- names(coef)
  names(pvals) = names(coef)
  ###
  res <- list(seed=seeds[i+1],coef=coef,se.coef=se.coef, pvals = pvals,coefe=coefe,se.coefe=se.coefe)
  return(res) 
}  ## }}}

onerunNS <- function(i,n,beta11=-0.3,betad1=-0.3,time=2000,cm=~1,var.z=1,cens=1/4,Yr="EpT",trans=NULL,scale1=1,scaled=1,dep=1,augmentR=NULL,augmentC=NULL,...) { ## {{{
  
  set.seed(seeds[i+1])
  
  nid <- n
  if (i%%500==0)	print(i)
  beta1 <- c(beta11,0.3); betad <- c(betad1,0.3); betac <- 0*c(0.3,-0.3); 
  gamma <- c(0.3,-0.3); 
  x2norm <- FALSE; 
  ###
  X <- matrix(rbinom(n*2,1,0.5),n,2)
  colnames(X) <- paste("X",1:2,sep="")
  if (x2norm)  { X2 <- rnorm(n); X[,2] <- X2; gX2 <- dcut(X2,breaks=2) } else gX2 <- X[,2]
  ###
  r1 <- exp( X %*% beta1)
  rd <- exp( X %*% betad)
  rc <- exp( X %*% betac)
  ###
  rr <- simRecurrent(nid,scalecumhaz(Lam1,scale1),death.cumhaz=scalecumhaz(LamD,scaled),
                     cens=cens,dependence=dep,var.z=var.z,r1=r1,rd=rd,rc=rc)
  rr$Z <- attr(rr,"z")[rr$id]
  rr <- cbind(rr,X[rr$id,])
  ###
  rr <- count.history(rr)
  ###
  nid <- max(rr$id)
  rr$revnr <- revcumsumstrata(rep(1,nrow(rr)),rr$id-1,nid)
  rr$statusD <- rr$status
  rr <- dtransform(rr,statusD=3,death==1 & revnr==1)
  rr$Nt_ <- rr$Count1
  ###
  rr$A <- rr$X1
  rr$x <- rr$X2
  dfactor(rr) <- A.f ~A
  ###
  ## constructing response at 2000 days
  dtable(rr,~statusD)
  rr$starty <- rr$start/365.25
  rr$stopy <- rr$stop/365.25
  rr$starty <- rr$start
  rr$stopy <- rr$stop
  ###
  
  # check simulation
  # cc1 <-  phreg(Event(starty,stopy,statusD==1)~X1+X2+cluster(id),rr)
  # ccd <-  phreg(Event(starty,stopy,statusD==3)~X1+X2+cluster(id),rr)
  
  # km1 <- phreg(Event(starty,stopy,statusD==3)~strata(A.f)+cluster(id),rr)
  # resmean.phreg(km1,times=3)
  
  outS <- WA_recurrent(Event(starty,stopy,statusD)~A.f+cluster(id),rr,trans=trans,time=time,death.code=3)
  doutS <- summary.WA(outS)
  ##
  ldoutS <- summary.WA(outS,type="log")
  
  
  if (is.null(augmentR)) augmentR <-  ~x
  if (is.null(augmentC)) augmentC <-  ~Count1+x
  
  if (dep!=0) 
  {
    outae <- WA_recurrent(Event(starty,stopy,statusD)~A.f+cluster(id),rr,augmentR=augmentR,
                          augmentC=augmentC,trans=trans,time=time,death.code=3) 
  } else { 
    outae <- WA_recurrent(Event(starty,stopy,statusD)~A.f+cluster(id),rr,
                          augmentR=augmentR,augmentC=augmentC,trans=trans,time=time,death.code=3) }
  
  dd <- outS$RAW$ratio.means$coefmat # Ratio of means E(N(min(D,t)))/E(min(D,t)) 
  alpha = 0.05
  doutae <- summary.WA(outae)
  
  #doutae$test.ratio$coefmat
  #doutae$ratio$coefmat # Ratio of means E(N(min(D,t)))/E(min(D,t)) 
  #doutae$meanpt$coefmat # Mean of Events per time-unit E(N(min(D,t))/min(D,t)) 
  
  # Here I am taking the [treat0] - [treat1] first for Mean of Events per time-unit (2 cases depending on augment R and C, i.e. outS and outae)
  # and then for the Ratio of means (only for outae)
  coef=c(doutS$test.meanpt$coefmat[1,1],doutae$test.meanpt$coefmat[1,1],doutS$test.ratio$coefmat[1,1],ldoutS$test.ratio$coefmat[1,1])
  se.coef=c( doutS$test.meanpt$coefmat[1,2],doutae$test.meanpt$coefmat[1,2],doutS$test.ratio$coefmat[1,2],ldoutS$test.ratio$coefmat[1,2]) # and se
  
  # now i am taking coefficients (treat0 and treat 1) first for  Mean of Events per time-unit then for ratio of means, 
  # first from outae then from outS
  coefe <- c(doutae$meanpt$coefmat[,1],doutae$ratio$coefmat[,1],doutS$meanpt$coefmat[,1],doutS$ratio$coefmat[,1],ldoutS$meanpt$coefmat[,1],ldoutS$ratio$coefmat[,1])
  names(coefe) <-c( rep(c("aug-mean-ratio","ratio"),each=2), rep(c("mean-ratio","ratio"),each=2), rep(c("log-mean-ratio","log-ratio"),each=2))
  se.coefe <- c(doutae$meanpt$coefmat[,2],doutae$ratio$coefmat[,2], doutS$meanpt$coefmat[,2],doutS$ratio$coefmat[,2], ldoutS$meanpt$coefmat[,2],ldoutS$ratio$coefmat[,2])
  
  dS <- summary(outS$ET$riskDR)
  dae <- summary(outae$ET$riskDR)
  coefbate <- c(outS$ET$riskDR$coef,dS$ateDR[,1],outae$ET$riskDR$coef,dae$ateDR[,1])
  se.coefbate <- c(outS$ET$riskDR$se.coef,dS$ateDR[,2],outae$ET$riskDR$se.coef,dae$ateDR[,2])
  
  
  # here i take denominators (treat0 and treat 1) and numerators
  coefb <- c(doutae$rmst$coefmat[,1],doutae$meanNtD$coefmat[,1])
  names(coefb) <- rep(c("rmst","meanN"),each=2)
  se.coefb <- c(doutae$rmst$coefmat[,2],doutae$meanNtD$coefmat[,2])
  names(se.coefb) <- rep(c("rmst","meanN"),each=2)
  
  # and putting all together 
  coefe <- c(coefe,coefb,coefbate)
  se.coefe <- c(se.coefe,se.coefb,se.coefbate)
  
  pvals = c(doutS$test.meanpt$coefmat[1,5],doutae$test.meanpt$coefmat[1,5],doutae$test.ratio$coefmat[1,5],ldoutS$test.ratio$coefmat[1,5]) < alpha
  
  names(coef) <- c("mean","Aug-mean","ratio","log-ratio")
  names(se.coef) <- names(coef)
  names(pvals) = names(coef)
  ###
  res <- list(seed=seeds[i+1],coef=coef,se.coef=se.coef, pvals = pvals,coefe=coefe,se.coefe=se.coefe)
  return(res) 
}  ## }}}

onerunNR <- function(i,n,beta11=-0.3,betad1=-0.3,time=2000,cm=~1,var.z=1,cens=1/4,Yr="EpT",trans=NULL,scale1=1,scaled=1,dep=1,...) { ## {{{
  
  set.seed(seeds[i+1])
  
  nid <- n
  if (i%%500==0)	print(i)
  beta1 <- c(beta11,0.3); betad <- c(betad1,0.3); betac <- 0*c(0.3,-0.3); 
  gamma <- c(0.3,-0.3); 
  x2norm <- FALSE; 
  ###
  X <- matrix(rbinom(n*2,1,0.5),n,2)
  colnames(X) <- paste("X",1:2,sep="")
  if (x2norm)  { X2 <- rnorm(n); X[,2] <- X2; gX2 <- dcut(X2,breaks=2) } else gX2 <- X[,2]
  ###
  r1 <- exp( X %*% beta1)
  rd <- exp( X %*% betad)
  rc <- exp( X %*% betac)
  ###
  rr <- simRecurrent(nid,scalecumhaz(Lam1,scale1),death.cumhaz=scalecumhaz(LamD,scaled),
                     cens=cens,dependence=dep,var.z=var.z,r1=r1,rd=rd,rc=rc)
  rr$Z <- attr(rr,"z")[rr$id]
  rr <- cbind(rr,X[rr$id,])
  ###
  rr <- count.history(rr)
  ###
  nid <- max(rr$id)
  rr$revnr <- revcumsumstrata(rep(1,nrow(rr)),rr$id-1,nid)
  rr$statusD <- rr$status
  rr <- dtransform(rr,statusD=3,death==1 & revnr==1)
  rr$Nt_ <- rr$Count1
  ###
  rr$A <- rr$X1
  rr$x <- rr$X2
  dfactor(rr) <- A.f ~A
  ###
  ## constructing response at 2000 days
  dtable(rr,~statusD)
  rr$starty <- rr$start/365.25
  rr$stopy <- rr$stop/365.25
  rr$starty <- rr$start
  rr$stopy <- rr$stop
  ###
  
  outS <- WA_recurrent(Event(starty,stopy,statusD)~A.f+cluster(id),rr,trans=trans,time=time,death.code=3)
  doutS <- summary(outS)
  ##
  ldoutS <- summary.WA(outS,type="log")
  
  if (dep!=0) 
  {
    outae <- WA_recurrent(Event(starty,stopy,statusD)~A.f+cluster(id),rr,augmentR=~x+Z,augmentC=~Count1+x+Z,trans=trans,time=time,death.code=3) 
  } else { 
    outae <- WA_recurrent(Event(starty,stopy,statusD)~A.f+cluster(id),rr,augmentR=~x,augmentC=~Count1+x,trans=trans,time=time,death.code=3) }
  
  doutae <- summary.WA(outae)
  
  coef=c(doutS$test.meanpt$coefmat[1,1],doutae$test.meanpt$coefmat[1,1],doutS$test.ratio$coefmat[1,1],ldoutS$test.ratio$coefmat[1,1])
  se.coef=c( doutS$test.meanpt$coefmat[1,2],doutae$test.meanpt$coefmat[1,2],doutS$test.ratio$coefmat[1,2],ldoutS$test.ratio$coefmat[1,2]) # and se
  
  alpha <- 0.05
  pvals = c(doutS$test.meanpt$coefmat[1,5],
            doutae$test.meanpt$coefmat[1,5],doutae$test.ratio$coefmat[1,5],ldoutS$test.ratio$coefmat[1,5]) < alpha
  
  names(coef) <- c("mean","aug-mean","ratio","log-ratio")
  names(se.coef) <- names(coef)
  names(pvals) = names(coef)
  ###
  res <- list(seed=seeds[i+1],coef=coef, se.coef=se.coef, pvals = pvals)
  return(res) 
}  ## }}}

ana <- function(res,true=NULL) { ## {{{
  
  coef <-do.call("rbind",lapply(res,function(x) x$coef) )
  scoef <-do.call("rbind",lapply(res,function(x) x$se.coef) )
  pvals <-do.call("rbind",lapply(res,function(x) x$pvals) )
  # print(pvals)
  m <- cbind( apply(coef,2,mean), apply(coef,2,sd), apply(scoef,2,mean), apply(pvals,2,mean) )
  colnames(m) <- c("mean","sd","mean-se","power")
  
  coefe <-do.call("rbind",lapply(res,function(x) x$coefe) )
  scoefe <-do.call("rbind",lapply(res,function(x) x$se.coefe) )
  me <- cbind( apply(coefe,2,mean), apply(coefe,2,sd), apply(scoefe,2,mean))
  colnames(me) <- c("mean","sd","mean-se")
  
  if (!is.null(true)) {# {{{
    truedif <- do.call("rbind",lapply(true,function(x) x$coef) )
    truedif <- apply(truedif,2,mean)
    #covdif <- apply(t(t(ee - 1.96*se.ee) < truedif & t(ee + 1.96*se.ee) > truedif),2,mean,na.rm=TRUE)
    covdif <- apply(t(t(coef - 1.96*scoef) < truedif & t(coef + 1.96*scoef) > truedif),2,mean,na.rm=TRUE)
    
    tcoefe <- do.call("rbind",lapply(true,function(x) x$coefe) )
    tcoefe <- apply(tcoefe,2,mean)
    #cove <- apply(t(t(ee - 1.96*se.ee) < tcoefe & t(ee + 1.96*se.ee) > tcoefe),2,mean,na.rm=TRUE)
    cove <- apply(t(t(coefe - 1.96*scoefe) < tcoefe & t(coefe + 1.96*scoefe) > tcoefe),2,mean,na.rm=TRUE)
    
    m <- cbind(m,covdif)
    me <- cbind(me,cove)
    colnames(m)[5] <- "coverage"
    colnames(me)[4] <- "coverage"
  } # }}}
  
  result <- list(summary = round(m, 4), coef=round(me,4)) 
  return(result)
} ## }}}

anaR <- function(res,true=NULL) { ## {{{
  
  coef <-do.call("rbind",lapply(res,function(x) x$coef) )
  scoef <-do.call("rbind",lapply(res,function(x) x$se.coef) )
  pvals <-do.call("rbind",lapply(res,function(x) x$pvals) )
  # print(pvals)
  m <- cbind( apply(coef,2,mean), apply(coef,2,sd), apply(scoef,2,mean), apply(pvals,2,mean) )
  colnames(m) <- c("mean","sd","mean-se","power")
  
  if (!is.null(true)) {# {{{
    truedif <- rep(0,3) 
    #covdif <- apply(t(t(ee - 1.96*se.ee) < truedif & t(ee + 1.96*se.ee) > truedif),2,mean,na.rm=TRUE)
    covdif <- apply(t(t(coef - 1.96*scoef) < truedif & t(coef + 1.96*scoef) > truedif),2,mean,na.rm=TRUE)
    
    m <- cbind(m,covdif)
    colnames(m)[5] <- "coverage"
  } # }}}
  
  result <- list(summary = round(m, 4)) 
  return(result)
} ## }}}