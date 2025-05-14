##### Demo Script for the Case Study in Web Appendix D #####
#           based on Colorectal cancer data                #

# Similar code is employed for the case study in Section 6
# based on the HF-Action clinical trial data

library(mets) 
?WA_recurrent

###### Load data ######
load("data/colorectal.rda")

###### Preprocessing ######
rr <- colorectal
rrid <- countID(rr)
rr <- cbind(rr,rrid)
rr$treat.f <- factor(rr$treatment)
rr <- count.history(rr,status="new.lesions",types=1)
rr$treattime <- (rr$lbnr__id==1)*1
###
rr$status <- rr$new.lesions
rr <- dtransform(rr,status=2,state==1)
dtable(rr,~status)
rr$rid <- rr$reverseCountid
rr$trt <- rr$treatment

# 
res <- resR <- c()
for (tt in seq(0.5,2.9,by=0.1)) {
  dd0 <- WA_recurrent(Event(time0,time1,status)~treat.f+cluster(id),
                      rr,time=tt,death.code=2,
                      augmentR=~age+who.PS+prev.resection,
                      augmentC=~age+who.PS+prev.resection) #, trans = 1/3)
  ee <- estimate(coef=dd0$ET$riskDR$riskDR,vcov=dd0$ET$riskDR$var.riskDR)$coefmat
  eeR <- dd0$RAW$ratio.means$coefmat
  tdif <- abs(diff(ee[,1])/sum(ee[,2]^2)^.5)
  pval <- 2*(1-pnorm(tdif))
  tdifR <- abs(diff(eeR[,1])/sum(eeR[,2]^2)^.5)
  pvalR <- 2*(1-pnorm(tdifR))
  res <- rbind(res,c(tt,ee[,1],ee[1,3:4],ee[2,3:4],pval))
  resR <- rbind(resR,c(tt,eeR[,1],eeR[1,3:4],eeR[2,3:4],pvalR))
}


plotres <- function(res, ylab) {
  matplot(res[,1],res[,2:3],type="l",lwd=3,ylim=c(0.2,1.2),xlab="Time (years)",ylab=ylab)
  plotConfRegion(res[,1],res[,4:5],col=1)
  plotConfRegion(res[,1],res[,6:7],col=2)
  if (ncol(res)==8) {
    sigp <- (res[,8]<0.05)
    points(res[sigp,1],res[sigp,2],pch="*",cex = 2.5, col=1)
    points(res[sigp,1],res[sigp,3],pch="*",cex = 2.5, col=2)
  }
  legend("bottomright",c("treatment:S","treatment:C"),lty=1:2,col=1:2,lwd=2.5)
}


plotres(res, ylab = "PWWA estimand")
plotres(resR, ylab = "EWWA estimand")



