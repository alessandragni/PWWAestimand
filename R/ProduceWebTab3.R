
library(xtable)

# Import functions
source("R/utils.R")


##### Load simulations for Table 2 #####
load("WAhfTab2_NULL.rda")


outtotSUMM <- c()
k <- 0
for (dep in c(1,4)) 
  for (time in c(3)) 
    for (cens in c(1/4)) 
      for (type in 1:3)
        for (scale1 in c(1)) 
          for (scaled in c(4,1)) 
          {
            # Lam1 <- Lams[[type]]
            for (varz  in c(2,1))  {
              if (dep==0 & varz>=1) break;
              k <- k+1
              ###
              mm <- anaR(outtab2$resl[[k]],true=outtab2$true[[k]])
              print(mm)
              outtotSUMM <- rbind(outtotSUMM,cbind(type,dep,varz,scale1,scaled,mm$summary))
            }
          }

outtotSUMM = as.data.frame(outtotSUMM) # for 0 and 1


###### Create Web Table 3 (beta1=betad=0) #####

temp = NULL

for (type in c(2,3,1))
  for (dep in c(1,4)) 
    for (varz  in c(1,2))
      for (scale1 in c(1)) 
        for (scaled in c(1,4)) 
          if(dep != 0) {
            df_row <- outtotSUMM[(outtotSUMM$type == type) &
                                   (outtotSUMM$dep == dep) &
                                   (outtotSUMM$varz == varz) &
                                   (outtotSUMM$scale1 == scale1) &
                                   (outtotSUMM$scaled == scaled), ]
            rownames(df_row) = NULL
            
            temp = rbind(temp,
                         cbind(df_row[2,-c(10)], df_row[3,c(6:9)]))
            rownames(temp) = NULL
            
          }


print(xtable(temp, digits = 3))

