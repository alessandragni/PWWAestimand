##### Demo Script for the Simulation Study in Section 5 #####

# remotes::install_github("kkholst/mets",ref="develop")
library(mets)
library(doMC)

# We base the simulation study on estimated data for patients of Copenhagen.
# In the paper, these simulations are based on the HF-Action clinical trial data

data(base1cumhaz)
data(base4cumhaz)
data(drcumhaz)
dr <- drcumhaz
base1 <- base1cumhaz
base4 <- base4cumhaz

Lam1 <- base1cumhaz
Lam2 <- base4cumhaz
LamD <- drcumhaz


##### Code for one run #####

i = 0

seeds <- sample(1:10^6, size = 5001 , replace = FALSE)

onerunN <- function(i,n,beta11=-0.3,betad1=-0.3,
                    time=1000,cm=~1,var.z=1,cens=1/4,Yr="EpT",
                    trans=NULL,scale1=1,scaled=1,dep=1,...) { ## {{{
  
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
  rr <- simRecurrentII(nid,scalecumhaz(Lam1,scale1),cumhaz2=Lam2,
                       death.cumhaz=scalecumhaz(LamD,scaled),cens=cens,dep=dep,var.z=var.z,r1=r1,rd=rd,rc=rc)
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
  
  outS <- WA_recurrent(Event(starty,stopy,statusD)~A.f+cluster(id),rr,trans=trans,time=time,death.code=3)
  doutS <- summary(outS)
  
  if (dep!=0) 
  {outae <- WA_recurrent(Event(starty,stopy,statusD)~A.f+cluster(id),rr,augmentR=~x+Z,augmentC=~Count1+x+Z,trans=trans,time=time,death.code=3) } else { 
    outae <- WA_recurrent(Event(starty,stopy,statusD)~A.f+cluster(id),rr,augmentR=~x,augmentC=~Count1+x,trans=trans,time=time,death.code=3) }
  
  dd <- outS$RAW$ratio.means$coefmat # Ratio of means E(N(min(D,t)))/E(min(D,t)) 
  
  alpha = 0.05
  doutae <- summary(outae)
  
  #doutae$test.ratio$coefmat
  doutae$ratio$coefmat # Ratio of means E(N(min(D,t)))/E(min(D,t)) 
  doutae$meanpt$coefmat # Mean of Events per time-unit E(N(min(D,t))/min(D,t)) 
  
  # Here I am taking the [treat0] - [treat1] first for Mean of Events per time-unit (2 cases depending on augment R and C, i.e. outS and outae)
  # and then for the Ratio of means (only for outae)
  coef=c(doutS$test.meanpt$coefmat[1,1],doutae$test.meanpt$coefmat[1,1], doutae$test.ratio$coefmat[1,1])
  se.coef=c(doutS$test.meanpt$coefmat[1,2],doutae$test.meanpt$coefmat[1,2], doutae$test.ratio$coefmat[1,2]) # and se
  
  # now i am taking coefficients (treat0 and treat 1) first for  Mean of Events per time-unit then for ratio of means, 
  # first from outae then from outS
  coefe <- c(doutae$meanpt$coefmat[,1],doutae$ratio$coefmat[,1],doutS$meanpt$coefmat[,1],doutS$ratio$coefmat[,1])
  names(coefe) <-c( rep(c("aug-mean-ratio","ratio"),each=2), rep(c("mean-ratio","ratio"),each=2))
  se.coefe <- c(doutae$meanpt$coefmat[,2],doutae$ratio$coefmat[,2],doutS$meanpt$coefmat[,2],doutS$ratio$coefmat[,2])
  
  # here i take denominators (treat0 and treat 1) and numerators
  coefb <- c(doutae$rmst$coefmat[,1],doutae$meanNtD$coefmat[,1])
  names(coefb) <- rep(c("rmst","meanN"),each=2)
  se.coefb <- c(doutae$rmst$coefmat[,2],doutae$meanNtD$coefmat[,2])
  names(se.coefb) <- rep(c("rmst","meanN"),each=2)
  
  # and putting all together 
  coefe <- c(coefe,coefb)
  se.coefe <- c(se.coefe,se.coefb)
  
  pvals = c( doutS$test.meanpt$coefmat[1,5],doutae$test.meanpt$coefmat[1,5],doutae$test.ratio$coefmat[1,5]) < alpha
  
  names(coef) <- c("mean","Aug-mean","ratio")
  names(se.coef) <- names(coef)
  names(pvals) = names(coef)
  ###
  res <- list(seed=seeds[i+1],coef=coef, se.coef=se.coef, pvals = pvals,coefe=coefe, se.coefe=se.coefe)
  return(res) 
}  ## }}}


##### Small test for one run to distinguish the cases #####

###### fralty shared : dep = 1 ######
onerunN(i,n=1000,var.z=2,time=3000,cens=NULL,dep=1,scale1=2,trans=NULL)
# var.z controls for the variance of the frailty

# trans is the transformation one may want to apply within the expectation;
# if set to null, the function will use no transformation
# if set to 0.333, for instance, a cube root transformation will be applied

###### independence : dep = 0 ######
onerunN(i,n=1000,var.z=2,time=3000,cens=NULL,dep=0,scale1=1)

###### frailty only recurrent : dep = 4 ######
onerunN(i,n=1000,var.z=2,time=3000,cens=NULL,dep=4,scale1=1)
onerunN(i,n=1000,var.z=2,time=3000,cens=NULL,dep=4,scale1=2,beta11=0,betad1=0)


##### Function for the analysis of multiple simulation #####

ana <- function(res,true=NULL) { ## {{{
  
  coef <-do.call("rbind",lapply(res,function(x) x$coef) )
  scoef <-do.call("rbind",lapply(res,function(x) x$se.coef) )
  pvals <-do.call("rbind",lapply(res,function(x) x$pvals) )

  m <- cbind( apply(coef,2,mean), apply(coef,2,sd), apply(scoef,2,mean), apply(pvals,2,mean) )
  colnames(m) <- c("mean","sd","mean-se","power")
  
  coefe <-do.call("rbind",lapply(res,function(x) x$coefe) )
  scoefe <-do.call("rbind",lapply(res,function(x) x$se.coefe) )
  me <- cbind( apply(coefe,2,mean), apply(coefe,2,sd), apply(scoefe,2,mean))
  colnames(me) <- c("mean","sd","mean-se")
  
  
  if (!is.null(true)) {# {{{
    truedif <- do.call("rbind",lapply(true,function(x) x$coef) )
    truedif <- apply(truedif,2,mean)
    covdif <- apply(t(t(coef - 1.96*scoef) < truedif & t(coef + 1.96*scoef) > truedif),2,mean,na.rm=TRUE)
    
    m <- cbind(m,covdif)
    colnames(m)[5] <- "coverage"
    
    tcoefe <- do.call("rbind",lapply(true,function(x) x$coefe) )
    tcoefe <- apply(tcoefe,2,mean)
    cove <- apply(t(t(coefe - 1.96*scoefe) < tcoefe & t(coefe + 1.96*scoefe) > tcoefe),2,mean,na.rm=TRUE)
    
    me <- cbind(me,cove)
    colnames(me)[4] <- "coverage"
  } # }}}
  
  result <- list(summary = round(m, 4), coef=round(me,4)) 
  return(result)
} ## }}}


##### Run some simulations #####

cc=detectCores()
registerDoMC(cc)
###

n <- 1000 # in the paper set to 1000
nsim <- 10 # in the paper set to 5000

######  first case : beta1 = betad = -0.3 ######
beta11 <- -0.3
betad1 <- -0.3

# set values (changed with respect to the paper)
depvals <- c(1, 4)
timevals <- c(4000)
cens_vals <- c(1/16)
theta_vals <- c(1)
scale1_vals <- c(2)
scaled_vals <- c(4)


# computation of the "mean","sd","mean-se",("power") 
resl <- list()
outtot <- c()
k <- 0
for (dep in depvals) 
  for (time in timevals) 
    for (cens in cens_vals) 
      for (varz in theta_vals) 
        for (scale1 in scale1_vals) 
          for (scaled in scaled_vals) 
          {
            if (dep==0 & varz>=1) break;
            print(c(dep,time,varz,scale1,scaled))
            k <- k+1
            ###
            resl[[k]] <- 
              foreach (i=0:nsim) %dopar% onerunN(i,n,var.z=varz,time=time,dep=dep,scale1=scale1,scaled=scaled,cens=cens,
                                                 trans=0.33,beta11=beta11,betad1=betad1)
            mm <- ana(resl[[k]])
            print(mm)
            outtot <- rbind(outtot,cbind(n,cens,dep,varz,scale1,mm$summary))
          }

out <- list(resl=resl,outtot=outtot)

save(out,file="simulations/DemoSim.rda")


# computation of the true values to obtain the "coverage"
n <- 100000
nsim <- 10
restrue <- list()
outtot2 <- c()
k <- 0
for (dep in depvals) 
  for (time in timevals) 
    for (cens in cens_vals) 
      for (varz in theta_vals) 
        for (scale1 in scale1_vals) 
          for (scaled in scaled_vals) 
          {
            if (dep==0 & varz>=1) break;
            print(c(dep,time,varz,scale1,scaled))
            k <- k+1
            ###
            restrue[[k]] <- 
              foreach (i=0:nsim) %dopar% onerunN(i,n,var.z=varz,time=time,dep=dep,
                                                 scale1=scale1,scaled=scaled,trans=0.33,cens=NULL)
            mm <- ana(resl[[k]],true=restrue[[k]])
            print(mm)
            outtot2 <- rbind(outtot2,cbind(n,cens,dep,varz,scale1,mm$summary))
          }

out <- list(resl=resl,outtot=outtot,true=restrue,outtot2=outtot2)
ana(out$resl[[1]],true=restrue[[1]])
save(out,file="simulations/DemoSim_withCov.rda")



###### one may want also to set beta11 = betad1 = 0 and repeat the analysis ######



##### Results analysis with beta1 = betad = -0.3 #####

load('simulations/DemoSim_withCov.rda')
resl = out[[1]]
outtot = out[[2]]
true = out[[3]]
outtot2 = out[[4]]


outtotSUMM <- c()
outtotCOEF <- c()
k <- 0
for (dep in depvals) 
  for (time in timevals) 
    for (cens in cens_vals) 
      for (varz in theta_vals) 
        for (scale1 in scale1_vals) 
          for (scaled in scaled_vals) 
          {
            if (dep==0 & varz>=1) break;
            print(c(dep,varz,cens,scale1))
            k <- k+1
            ###
            mm <- ana(resl[[k]], true = true[[k]])
            print(mm)
            outtotSUMM <- rbind(outtotSUMM,cbind(dep,varz,cens,scale1,mm$summary))
            outtotCOEF <- rbind(outtotCOEF,cbind(dep,varz,cens,scale1,mm$coef))
          }

outtotSUMM = as.data.frame(outtotSUMM)
outtotCOEF = as.data.frame(outtotCOEF)


###### Filling Table 1 ######
# For the first estimator in Eq. (7)  we employ aug.mean.ratio (for 0 and 1)
# for the second estimator in Eq. (8) we employ mean.ratio (for 0 and 1)

# e.g. v = 1 and theta = first value and kc = first value
# two cases 0 and 1
outtotCOEF[(outtotCOEF$dep==depvals[1]) & (outtotCOEF$varz==theta_vals[1]) & 
             (outtotCOEF$cens==cens_vals[1]) & (outtotCOEF$scale1==scale1_vals[1]), ] 
# 0 - 1
outtotSUMM[(outtotSUMM$dep==depvals[1]) & (outtotSUMM$varz==theta_vals[1]) & 
             (outtotSUMM$cens==cens_vals[1]) & (outtotSUMM$scale1==scale1_vals[1]), -c(8)] 


###### Filling Table 2 (first part) ######
# we employ aug.mean and ratio

# e.g. v = 1 and theta = first value and kc = first value and scale1 = first value
outtotSUMM[(outtotSUMM$dep==depvals[1]) & (outtotSUMM$varz==theta_vals[1]) & 
             (outtotSUMM$cens==cens_vals[1]) & (outtotSUMM$scale1==scale1_vals[1]),   -c(9)][c(2,3),]





###### CODE for producing tables ######

library(xtable)

####### For Table 1 #######
# depvals <- c(1) #, 4)
# theta_vals <- c(0.5, 1, 2)
# cens_vals <- c(1/4, 2/4)
# scale1_vals <- c(1)

temp = NULL

for (dep in depvals) {
  for (theta in theta_vals) {
    for (cens in cens_vals) {
      for (s1 in scale1_vals) {
        df_row1 = outtotCOEF[(outtotCOEF$dep == dep) & 
                               (outtotCOEF$varz == theta) & 
                               (outtotCOEF$cens == cens) & 
                               (outtotCOEF$scale1 == s1), ][c(1,2,5,6),]
        rownames(df_row1) = NULL
        
        df_row <- outtotSUMM[(outtotSUMM$dep == dep) &
                               (outtotSUMM$varz == theta) &
                               (outtotSUMM$cens == cens) &
                               (outtotSUMM$scale1 == s1), -c(8)]
        rownames(df_row) = NULL
        
        temp = rbind(temp,
                     cbind(df_row1[1:2,-c(4)], df_row1[3:4,c(5:8)]),
                     cbind(df_row[2,-c(4)], df_row[1,c(5:8)]))
        rownames(temp) = NULL
      }
    }
  }
}

print(xtable(temp, digits = 4))




###### For Table 2 (first part) #####
# depvals <- c(1, 4)
# theta_vals <- c(0.5, 1, 2)
# cens_vals <- c(1/4)
# scale1_vals <- c(0.5, 1, 2, 3)

temp = NULL

for (dep in depvals) {
  for (theta in theta_vals) {
    for (cens in cens_vals) {
      for (s1 in scale1_vals) {
        
        df_row <- outtotSUMM[(outtotSUMM$dep == dep) &
                               (outtotSUMM$varz == theta) &
                               (outtotSUMM$cens == cens) &
                               (outtotSUMM$scale1 == s1), -c(9)]
        rownames(df_row) = NULL
        
        temp = rbind(temp,
                     cbind(df_row[2,-c(4)], df_row[3,c(5:8)]))
        rownames(temp) = NULL
      }
    }
  }
}

print(xtable(temp, digits = 4))

