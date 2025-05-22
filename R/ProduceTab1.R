
library(xtable)


# Import functions
source("R/utils.R")


##### Load simulations for Table 1, case (i) #####
load("WAhfTab1.rda")


outtotSUMM <- outtotCOEF <- c()
k <- 0
for (dep in c(1,4,0)) 
  for (time in c(3)) 
    for (cens in c(1/4,1/2)) 
      for (type in 2)
        for (scale1 in c(1)) 
          for (scaled in c(1)) 
          {
            # Lam1 <- Lams[[type]]
            for (varz  in c(2,1,0.5))  {
              if (dep==0 & varz>=1) break;
              k <- k+1
              ###
              mm <- ana(outtab1$resl[[k]],true=outtab1$true[[k]])
              print(mm)
              outtotSUMM <- rbind(outtotSUMM,cbind(type,dep,varz,cens,scale1,scaled,mm$summary))
              outtotCOEF <- rbind(outtotCOEF,cbind(type,dep,varz,cens,scale1,scaled,mm$coef))
            }
          }

outtotSUMM = as.data.frame(outtotSUMM) # for 0 and 1
outtotCOEF = as.data.frame(outtotCOEF) # for  (0-1)

###

temp = NULL

for (dep in c(1,4,0)) 
  for (time in c(3)) 
    for (varz  in c(0.5,1,2))
      for (cens in c(1/4,1/2)) 
        for (type in 2)
          for (scale1 in c(1)) 
            for (scaled in c(1))   {
              if(dep != 0) {
        df_row1 = outtotCOEF[(outtotCOEF$dep == dep) & 
                               (outtotCOEF$varz == varz) & 
                               (outtotCOEF$cens == cens) & 
                               (outtotCOEF$scale1 == scale1) &
                               (outtotCOEF$scaled == scaled), ][c(1,2,5,6),]
        print(df_row1)
        rownames(df_row1) = NULL
        
        df_row <- outtotSUMM[(outtotSUMM$dep == dep) &
                               (outtotSUMM$varz == varz) &
                               (outtotSUMM$cens == cens) &
                               (outtotSUMM$scale1 == scale1) &
                               (outtotSUMM$scaled == scaled), -c(10)] # remove power from here
        print(df_row)
        rownames(df_row) = NULL
        
        
        temp = rbind(temp,
                     cbind(df_row1[1:2,], df_row1[3:4,c(7:10)]),
                     cbind(df_row[2,], df_row[1,c(7:10)]))
        rownames(temp) = NULL
              }
            }


print(xtable(temp, digits = 4))



##### Load simulations for Table 1, case (ii) #####

# Simulations augmenting with only L for the RCT part and 
# L+Count1 for the censoring augmentation

###### low censoring 1/4 ######

load("WAhfTab1onlyL.rda") # low censoring 1/4

outtotSUMML <- outtotCOEFL <- c()
k <- 0
for (dep in c(1)) 
  for (time in c(3)) 
    for (cens in c(1/4)) 
      for (type in 2)
        for (scale1 in c(1)) 
          for (scaled in c(1)) 
          {
            # Lam1 <- Lams[[type]]
            for (varz  in c(0.5,1,2))  {
              if (dep==0 & varz>=1) break;
              k <- k+1
              ###
              mm <- ana(outtabL$resl[[k]],true=outtabL$true[[k]])
              print(mm)
              #outtotSUMM <- rbind(outtotSUMM,cbind(dep,varz,cens,scale1,mm$summary))
              outtotSUMML <- rbind(outtotSUMML,cbind(type,dep,varz,cens,scale1,scaled,mm$summary))
              #outtotCOEF <- rbind(outtotCOEF,cbind(dep,varz,cens,scale1,mm$coef))
              outtotCOEFL <- rbind(outtotCOEFL,cbind(type,dep,varz,cens,scale1,scaled,mm$coef))
            }
          }

outtotSUMML = as.data.frame(outtotSUMML) # for 0 and 1
outtotCOEFL = as.data.frame(outtotCOEFL) # for  (0-1)

###

###### high censoring 2/4 ######
load("WAhfTab1onlyH.rda")  # high censoring 2/4
outtotSUMMH <- outtotCOEFH <- c()
k <- 0
for (dep in c(1)) 
  for (time in c(3)) 
    for (cens in c(2/4)) 
      for (type in 2)
        for (scale1 in c(1)) 
          for (scaled in c(1)) 
          {
            # Lam1 <- Lams[[type]]
            for (varz  in c(0.5,1,2))  {
              if (dep==0 & varz>=1) break;
              k <- k+1
              ###
              mm <- ana(outtabH$resl[[k]],true=outtabH$true[[k]])
              print(mm)
              outtotSUMMH <- rbind(outtotSUMMH,cbind(type,dep,varz,cens,scale1,scaled,mm$summary))
              outtotCOEFH <- rbind(outtotCOEFH,cbind(type,dep,varz,cens,scale1,scaled,mm$coef))
            }
          }

outtotSUMMH = as.data.frame(outtotSUMMH) # for 0 and 1
outtotCOEFH = as.data.frame(outtotCOEFH) # for  (0-1)

outtotSUMM = rbind(outtotSUMML, outtotSUMMH)
outtotCOEF = rbind(outtotCOEFL, outtotCOEFH)


##### Table 1: all together #####

temp2 = NULL

for (dep in c(1)) 
  for (time in c(3)) 
    for (varz  in c(0.5,1,2))
      for (cens in c(1/4,1/2)) 
        for (type in 2)
          for (scale1 in c(1)) 
            for (scaled in c(1))   {
              if(dep != 0) {
                df_row1 = outtotCOEF[(outtotCOEF$dep == dep) & 
                                       (outtotCOEF$varz == varz) & 
                                       (outtotCOEF$cens == cens) & 
                                       (outtotCOEF$scale1 == scale1) &
                                       (outtotCOEF$scaled == scaled), ][c(1,2,5,6),]
                print(df_row1)
                rownames(df_row1) = NULL
                
                df_row <- outtotSUMM[(outtotSUMM$dep == dep) &
                                       (outtotSUMM$varz == varz) &
                                       (outtotSUMM$cens == cens) &
                                       (outtotSUMM$scale1 == scale1) &
                                       (outtotSUMM$scaled == scaled), -c(10)] # remove power from here
                print(df_row)
                rownames(df_row) = NULL
                
                temp2 = rbind(temp2,
                             cbind(df_row1[1:2,], df_row1[3:4,c(7:10)]),
                             cbind(df_row[2,], df_row[1,c(7:10)]))
                rownames(temp2) = NULL
              }
            }


print(xtable(temp2, digits = 4))


temp = round(temp,3)
temp = temp[temp$dep==1,]
temp2 = round(temp2,3)
temp$type = temp2$type = NULL
temp$dep = temp2$dep = NULL
temp2 = temp2[,1:8]

x = cbind(temp[,1:8], temp2[,5:8], temp[,9:12])

print(xtable(x,digits=3))


