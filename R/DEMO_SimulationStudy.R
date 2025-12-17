##### Demo Script for the Simulation Study in Section 5 #####

# We base the simulation study on estimated data for patients of Copenhagen.
# In the paper, these simulations are based on the HF-Action clinical trial data,
# that however cannot be shared

# remotes::install_github("kkholst/mets",ref="develop")

library(mets)
library(doMC)

# Import functions
source("R/utils.R")


# In the functions onerunN, onerunNS, onerunR:

# - trans is the transformation one may want to apply within the expectation;
#   if set to null, the function will use no transformation
#   if set to 0.333, for instance, a cube root transformation will be applied

# - dep sets the frailty: if frailty shared, set dep = 1 (v = 1 in the paper); 
#   if frailty only recurrent events, set dep = 4 (v = 0 in the paper)

# - var.z controls for the variance of the frailty

# We import Copenhagen data
data(base1cumhaz)
data(base4cumhaz)
data(drcumhaz)
dr <- drcumhaz
base1 <- base1cumhaz
base4 <- base4cumhaz


# We create the the three cases a), b) and c) as in the paper
Lam11 <- cbind(c(0,1,4.31),c(0,2.5,3.45)) # case (c)
Lam12 <- base1cumhaz                      # case (a) - based on other data
Lam13 <- cbind(c(0,1,4.31),c(0,0.5,3.45)) # case (b)
Lams <- list(Lam11,Lam12,Lam13)

Lam2 <- base4cumhaz
LamD <- drcumhaz

seeds <- sample(1:10^6, size = 5001 , replace = FALSE)


##### Run some simulations #####

cc=detectCores()
registerDoMC(cc)

###


######  ONE EXAMPLE FOR TAB 1 : beta1 = betad = -0.3, case (i) ######
beta11 <- -0.3
betad1 <- -0.3

# set values (changed with respect to the paper)
depvals <- c(1, 4)
timevals <- c(4000)
cens_vals <- c(1/16)
theta_vals <- c(1)
scale1_vals <- c(1)
scaled_vals <- c(2)
type_vals <- c(1)


n <- 1000 # in the paper set to 1000
nsim <- 10 # in the paper set to 5000

# computation of the "mean","sd","mean-se",("power") 
resl <- list()
outtot <- c()
k <- 0
for (dep in depvals) 
  for (time in timevals) 
    for (cens in cens_vals) 
      for (type in type_vals)
        for (scale1 in scale1_vals) 
          for (scaled in scaled_vals) 
          {
            Lam1 <- Lams[[type]]
            for (varz in theta_vals)  {
              if (dep==0 & varz>=1) break;
              print(c(dep,type,time,varz,scale1,scaled))
              k <- k+1
              ###
              resl[[k]] <- 
                foreach (i=0:nsim) %dopar% 
                onerunN(i,n,var.z=varz,time=time,dep=dep,scale1=scale1,scaled=scaled,cens=cens,trans=0.33,beta11=beta11,betad1=betad1)
              # for case (ii) change this way:
              # resl[[k]] <- 
              #   foreach (i=0:nsim) %dopar% 
              #   onerunNS(i,n,var.z=varz,time=time,dep=dep,scale1=scale1,scaled=scaled,cens=cens,trans=0.33,beta11=beta11,betad1=betad1)
              mm <- ana(resl[[k]])
              print(mm)
              outtot <- rbind(outtot,cbind(n,cens,type,dep,varz,scale1,scaled,mm$summary))
              }
          }


# computation of the true values 
# as the mean of simulations to obtain also the "bias", "mse" and the "coverage"

n <- 10000
nsim <- 100
restrue <- list()
outtot2 <- c()
k <- 0
for (dep in depvals) 
  for (time in timevals) 
    for (cens in cens_vals) 
      for (type in 2)
        for (scale1 in scale1_vals) 
          for (scaled in scaled_vals) 
          { 
            Lam1 <- Lams[[type]]
            for (varz  in theta_vals)  {
            if (dep==0 & varz>=1) break;
            print(c(dep,type,time,varz,scale1,scaled))
            k <- k+1
            ###
            restrue[[k]] <- 
              foreach (i=0:nsim) %dopar% 
              onerunN(i,n,var.z=varz,time=time,dep=dep,scale1=scale1,scaled=scaled,cens=NULL,trans=0.33,beta11=beta11,betad1=betad1)
            # for case (ii) change this way:
            # restrue[[k]] <- 
            #   foreach (i=0:nsim) %dopar% 
            #   onerunNS(i,n,var.z=varz,time=time,dep=dep,scale1=scale1,scaled=scaled,cens=cens,trans=0.33,beta11=beta11,betad1=betad1)
            }
          }


outtab1 <- list(resl=resl,outtot=outtot,true=restrue,beta11=beta11,betad1=betad1)
# save(outtab1, file="WAhfTab1.rda")





######  ONE EXAMPLE FOR TAB 2 : beta1 = betad = -0.3 ######

beta11 <- -0.3
betad1 <- -0.3 
# set both to 0 for the WEB TABLE 3 (NULL CASE) 
# and save through save(outtab2,file="WAhfTab2_NULL.rda")


n <- 1000 # in the paper set to 1000
nsim <- 10 # in the paper set to 5000


# computation of the "mean","sd","mean-se",("power") 
resl <- list()
outtot <- c()
k <- 0
for (dep in depvals) 
  for (time in timevals) 
    for (cens in cens_vals) 
      for (type in type_vals)
        for (scale1 in scale1_vals) 
          for (scaled in scaled_vals) 
          {
            Lam1 <- Lams[[type]]
            for (varz in theta_vals)  {
              if (dep==0 & varz>=1) break;
              print(c(dep,type,time,varz,scale1,scaled))
              k <- k+1
              ###
              resl[[k]] <- 
                foreach (i=0:nsim) %dopar% 
                onerunNR(i,n,var.z=varz,time=time,dep=dep,scale1=scale1,scaled=scaled,cens=cens,trans=0.33,beta11=beta11,betad1=betad1)
              mm <- anaR(resl[[k]])
              print(mm)
              outtot <- rbind(outtot,cbind(n,cens,type,dep,varz,scale1,scaled,mm$summary))
              print(outtot)
            }
          }

# computation of the true values 
# as the mean of simulations

n <- 10000 
nsim <- 100

restrue <- list()
k <- 0
for (dep in depvals) 
  for (time in timevals) 
    for (cens in cens_vals) 
      for (type in type_vals)
        for (scale1 in scale1_vals) 
          for (scaled in scaled_vals) 
          {
            Lam1 <- Lams[[type]]
            for (varz in theta_vals)  {
              if (dep==0 & varz>=1) break;
              print(c(dep,type,time,varz,scale1,scaled))
              k <- k+1
              ###
              restrue[[k]] <- 
                foreach (i=0:nsim) %dopar% 
                onerunNR(i,n,var.z=varz,time=time,dep=dep,scale1=scale1,scaled=scaled,cens=NULL,trans=0.33,beta11=beta11,betad1=betad1)
            }
          }

outtab2 <- list(resl=resl,outtot=outtot,true=restrue,beta11=beta11,betad1=betad1)
# save(outtab2,file="WAhfTab2.rda")





