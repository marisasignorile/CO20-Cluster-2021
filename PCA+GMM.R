library(readxl)
library(writexl)
suppressMessages(library(dplyr))
library(tidyr)
library(mclust)
library(reportRx)
library(survminer)
require("survival")
library(lattice)
require(ggplot2)
library(ggpubr)
library(ggthemes)

wtox<-read_excel("/Users/marisasignorile/Desktop/Clustering2/wtoxco20.xlsx")
wtox$id<-as.factor(wtox$id)
#trt arm 1 (b+c)
dfbc<-wtox[which(wtox$trt==1),]
dfbc<-dfbc[,-2]
#trt arm 2 (p+c)
dfpc<-wtox[which(wtox$trt==2),]
dfpc<-dfpc[,-2]

#remove AST/ALT, hyperglycemia, hypomagnesemia, and hyponatremia from analysis
dfbc<-dfbc[,-c(2,3,4,5,30:41)]
dfpc<-dfpc[,-c(2,3,4,5,30:41)]

#trt arm 1 (b+c)

#week 0
bc0.pr <- prcomp(dfbc[c(2,6,10,14,18,22,26,30,34,38)], center = FALSE, scale = FALSE)
bc0.pr1<-bc0.pr$x[,1]#PC1
#bc0.pr[["rotation"]] #weights for each toxicity

#week 2
bc2.pr <- prcomp(dfbc[c(3,7,11,15,19,23,27,31,35,39)], center = FALSE, scale = FALSE)
bc2.pr1<-bc2.pr$x[,1]#PC1
#bc2.pr[["rotation"]] #weights for each toxicity

#week 4
bc4.pr <- prcomp(dfbc[c(4,8,12,16,20,24,28,32,36,40)], center = FALSE, scale = FALSE)
bc4.pr1<-bc4.pr$x[,1]#PC1
#bc4.pr[["rotation"]] #weights for each toxicity

#week 8
bc8.pr <- prcomp(dfbc[c(5,9,13,17,21,25,29,33,37,41)], center = FALSE, scale = FALSE)
bc8.pr1<-bc8.pr$x[,1]#PC1
#bc8.pr[["rotation"]] #weights for each toxicity

#trt arm 2 (p+c)

#week 0
pc0.pr <- prcomp(dfpc[c(2,6,10,14,18,22,26,30,34,38)], center = FALSE, scale = FALSE)
pc0.pr1<-pc0.pr$x[,1]#PC1
#pc0.pr[["rotation"]] #weights for each toxicity

#week 2
pc2.pr <- prcomp(dfpc[c(3,7,11,15,19,23,27,31,35,39)], center = FALSE, scale = FALSE)
pc2.pr1<-pc2.pr$x[,1]#PC1
#pc2.pr[["rotation"]] #weights for each toxicity

#week 4
pc4.pr <- prcomp(dfpc[c(4,8,12,16,20,24,28,32,36,40)], center = FALSE, scale = FALSE)
pc4.pr1<-pc4.pr$x[,1]#PC1
#pc4.pr[["rotation"]] #weights for each toxicity

#week 8
pc8.pr <- prcomp(dfpc[c(5,9,13,17,21,25,29,33,37,41)], center = FALSE, scale = FALSE)
pc8.pr1<-pc8.pr$x[,1]#PC1
#pc8.pr[["rotation"]] #weights for each toxicity


#GMM trt1 and trt2

PCbc<-cbind(bc0.pr1,bc2.pr1,bc4.pr1,bc8.pr1)
PCpc<-cbind(pc0.pr1,pc2.pr1,pc4.pr1,pc8.pr1)



#GMM trt 1

gmm2t1 = Mclust(PCbc,2)
gmm3t1 = Mclust(PCbc,3)
gmm4t1 = Mclust(PCbc,4)

dfbc$cluster<-gmm2t1$classification

#2 clusters is the most optimal
#plot(gmm2t1, what=c("classification"))




#GMM trt 2

gmm2t2 = Mclust(PCpc,2)
gmm3t2 = Mclust(PCpc,3)
gmm4t2 = Mclust(PCpc,4)

dfpc$cluster<-gmm2t2$classification

#2 clusters is the most optimal
#plot(gmm2t2, what=c("classification"))

#plot trt 1 by symptom


#convert to long form
CVDlong<-reshape(as.data.frame(dfbc[,c(1,2:5,42)]) ,direction="long", idvar="id",varying=c("CVD 0","CVD 2","CVD 4","CVD 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

CVDlong$cluster<-as.factor(CVDlong$cluster) 

s2 <- ggplot(data = CVDlong, aes(x = week, y = Prevalence, group = id))
s2 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("CVD (Treatment arm 1)") 

#convert to long form
Diarrlong<-reshape(as.data.frame(dfbc[,c(1,6:9,42)]) ,direction="long", idvar="id",varying=c("Diarrhea 0","Diarrhea 2","Diarrhea 4","Diarrhea 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

Diarrlong$cluster<-as.factor(Diarrlong$cluster) 

s3 <- ggplot(data = Diarrlong, aes(x = week, y = Prevalence, group = id))
s3 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Diarrhea (Treatment arm 1)")

#convert to long form
dysplong<-reshape(as.data.frame(dfbc[,c(1,10:13,42)]) ,direction="long", idvar="id",varying=c("Dyspnea 0","Dyspnea 2","Dyspnea 4","Dyspnea 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

dysplong$cluster<-as.factor(dysplong$cluster) 

s4 <- ggplot(data = dysplong, aes(x = week, y = Prevalence, group = id))
s4 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Dyspnea (Treatment arm 1)")



#convert to long form
Fatilong<-reshape(as.data.frame(dfbc[,c(1,14:17,42)]) ,direction="long", idvar="id",varying=c("Fatigue 0","Fatigue 2","Fatigue 4","Fatigue 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

Fatilong$cluster<-as.factor(Fatilong$cluster) 

s5 <- ggplot(data = Fatilong, aes(x = week, y = Prevalence, group = id))
s5 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Fatigue (Treatment arm 1)")





#convert to long form
GIlong<-reshape(as.data.frame(dfbc[,c(1,18:21,42)]) ,direction="long", idvar="id",varying=c("GI 0","GI 2","GI 4","GI 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

GIlong$cluster<-as.factor(GIlong$cluster) 

s6 <- ggplot(data = GIlong, aes(x = week, y = Prevalence, group = id))
s6 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("GI (Treatment arm 1)")

#convert to long form
headlong<-reshape(as.data.frame(dfbc[,c(1,22:25,42)]) ,direction="long", idvar="id",varying=c("Head 0","Head 2","Head 4","Head 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

headlong$cluster<-as.factor(headlong$cluster) 

s7 <- ggplot(data = headlong, aes(x = week, y = Prevalence, group = id))
s7 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Head Pain (Treatment arm 1)")


#convert to long form
msklong<-reshape(as.data.frame(dfbc[,c(1,26:29,42)]) ,direction="long", idvar="id",varying=c("MSK 0","MSK 2","MSK 4","MSK 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

msklong$cluster<-as.factor(msklong$cluster) 

s11 <- ggplot(data = msklong, aes(x = week, y = Prevalence, group = id))
s11 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("MSK Pain (Treatment arm 1)")

#convert to long form
rashlong<-reshape(as.data.frame(dfbc[,c(1,30:33,42)]) ,direction="long", idvar="id",varying=c("Rash 0","Rash 2","Rash 4","Rash 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

rashlong$cluster<-as.factor(rashlong$cluster) 

s12 <- ggplot(data = rashlong, aes(x = week, y = Prevalence, group = id))
s12 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Rash (Treatment arm 1)")

#convert to long form
Skinlong<-reshape(as.data.frame(dfbc[,c(1,34:37,42)]) ,direction="long", idvar="id",varying=c("Skin 0","Skin 2","Skin 4","Skin 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

Skinlong$cluster<-as.factor(Skinlong$cluster) 

s13 <- ggplot(data = Skinlong, aes(x = week, y = Prevalence, group = id))
s13 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Skin (Treatment arm 1)")

#convert to long form
wllong<-reshape(as.data.frame(dfbc[,c(1,38:41,42)]) ,direction="long", idvar="id",varying=c("Weight loss 0","Weight loss 2","Weight loss 4","Weight loss 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

wllong$cluster<-as.factor(wllong$cluster) 

s14 <- ggplot(data = wllong, aes(x = week, y = Prevalence, group = id))
s14 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Weight loss (Treatment arm 1)")


b<-s2 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("CVD")+theme_minimal()
c<-s3 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Diarrhea") +theme_minimal()
e<-s4 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Dyspnea") +theme_minimal()
f<-s5 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Fatigue") +theme_minimal()
g<-s6 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("GI ") +theme_minimal()
h<-s7 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Head Pain") +theme_minimal()
l<-s11 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("MSK Pain") +theme_minimal()
m<-s12 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Rash") +theme_minimal()
n<-s13 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Skin") +theme_minimal()
o<-s14 + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Weight loss") +theme_minimal()


ggarrange(b,c,e,f,g,h,l,m,n,o,ncol=5,nrow=2)+theme_minimal()




#plot trt 2 by each symptom


#convert to long form
CVD2long<-reshape(as.data.frame(dfpc[,c(1,2:5,42)]) ,direction="long", idvar="id",varying=c("CVD 0","CVD 2","CVD 4","CVD 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

CVD2long$cluster<-as.factor(CVD2long$cluster) 

s2b <- ggplot(data = CVD2long, aes(x = week, y = Prevalence, group = id))
s2b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("CVD ") 




#convert to long form
Diarr2long<-reshape(as.data.frame(dfpc[,c(1,6:9,42)]) ,direction="long", idvar="id",varying=c("Diarrhea 0","Diarrhea 2","Diarrhea 4","Diarrhea 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

Diarr2long$cluster<-as.factor(Diarr2long$cluster) 

s3b <- ggplot(data = Diarr2long, aes(x = week, y = Prevalence, group = id))
s3b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Diarrhea")




#convert to long form
dysp2long<-reshape(as.data.frame(dfpc[,c(1,10:13,42)]) ,direction="long", idvar="id",varying=c("Dyspnea 0","Dyspnea 2","Dyspnea 4","Dyspnea 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

dysp2long$cluster<-as.factor(dysp2long$cluster) 

s4b <- ggplot(data = dysp2long, aes(x = week, y = Prevalence, group = id))
s4b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Dyspnea")



#convert to long form
Fati2long<-reshape(as.data.frame(dfpc[,c(1,14:17,42)]) ,direction="long", idvar="id",varying=c("Fatigue 0","Fatigue 2","Fatigue 4","Fatigue 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

Fati2long$cluster<-as.factor(Fati2long$cluster) 

s5b <- ggplot(data = Fati2long, aes(x = week, y = Prevalence, group = id))
s5b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Fatigue ")




#convert to long form
GI2long<-reshape(as.data.frame(dfpc[,c(1,18:21,42)]) ,direction="long", idvar="id",varying=c("GI 0","GI 2","GI 4","GI 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

GI2long$cluster<-as.factor(GI2long$cluster) 

s6b <- ggplot(data = GI2long, aes(x = week, y = Prevalence, group = id))
s6b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("GI ")




#convert to long form
head2long<-reshape(as.data.frame(dfpc[,c(1,22:25,42)]) ,direction="long", idvar="id",varying=c("Head 0","Head 2","Head 4","Head 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

head2long$cluster<-as.factor(head2long$cluster) 

s7b <- ggplot(data = head2long, aes(x = week, y = Prevalence, group = id))
s7b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Head Pain ")





#convert to long form
msk2long<-reshape(as.data.frame(dfpc[,c(1,26:29,42)]) ,direction="long", idvar="id",varying=c("MSK 0","MSK 2","MSK 4","MSK 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

msk2long$cluster<-as.factor(msk2long$cluster) 

s11b <- ggplot(data = msk2long, aes(x = week, y = Prevalence, group = id))
s11b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("MSK Pain")



#convert to long form
rash2long<-reshape(as.data.frame(dfpc[,c(1,30:33,42)]) ,direction="long", idvar="id",varying=c("Rash 0","Rash 2","Rash 4","Rash 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

rash2long$cluster<-as.factor(rash2long$cluster) 

s12b <- ggplot(data = rash2long, aes(x = week, y = Prevalence, group = id))
s12b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Rash ")




#convert to long form
Skin2long<-reshape(as.data.frame(dfpc[,c(1,34:37,42)]) ,direction="long", idvar="id",varying=c("Skin 0","Skin 2","Skin 4","Skin 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

Skin2long$cluster<-as.factor(Skin2long$cluster) 

s13b <- ggplot(data = Skin2long, aes(x = week, y = Prevalence, group = id))
s13b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Skin ")




#convert to long form
wl2long<-reshape(as.data.frame(dfpc[,c(1,38:41,42)]) ,direction="long", idvar="id",varying=c("Weight loss 0","Weight loss 2","Weight loss 4","Weight loss 8"), v.names="Prevalence", timevar="week",times=c(0,2,4,8))

wl2long$cluster<-as.factor(wl2long$cluster) 

s14b <- ggplot(data = wl2long, aes(x = week, y = Prevalence, group = id))
s14b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Weight loss ")

b2<-s2b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("CVD") +theme_minimal()
c2<-s3b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Diarrhea") +theme_minimal()
e2<-s4b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Dyspnea") +theme_minimal()
f2<-s5b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Fatigue") +theme_minimal()
g2<-s6b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("GI ") +theme_minimal()
h2<-s7b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Head Pain") +theme_minimal()
l2<-s11b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("MSK Pain") +theme_minimal()
m2<-s12b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Rash") +theme_minimal()
n2<-s13b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Skin") +theme_minimal()
o2<-s14b + stat_summary(aes(group = cluster,col=cluster), geom = "line", fun = mean, size = 1)+ggtitle("Weight loss") +theme_minimal()


ggarrange(b2,c2,e2,f2,g2,h2,l2,m2,n2,o2,ncol=5,nrow=2)+theme_minimal()



#merge demographic and cluster data frames and read in
gmmbc<-read_xlsx("/Users/marisasignorile/Desktop/Clustering2/gmmbc2.xlsx")
gmmpc<-read_xlsx("/Users/marisasignorile/Desktop/Clustering2/gmmpc2.xlsx")
gmmbc<-as.data.frame(gmmbc)
gmmpc<-as.data.frame(gmmpc)
gmmbc$id<-as.factor(gmmbc$id)
gmmbc$cluster<-as.character(gmmbc$cluster)
gmmpc$id<-as.factor(gmmpc$id)
gmmpc$cluster<-as.character(gmmpc$cluster)


#demographics
#rm_covsum(data=gmmbc,covs=c("Age","Sex","Number of metastatic sites","Previous chemo drug classes","Primary site","Presence of liver metastasis","ECOG"),digits=2,maincov = "cluster",all.stats=FALSE,include_missing=T,percentage='column', pvalue=TRUE)
#rm_covsum(data=gmmpc,covs=c("Age","Sex","Number of metastatic sites","Previous chemo drug classes","Primary site","Presence of liver metastasis","ECOG"),digits=2,maincov = "cluster",all.stats=FALSE,include_missing=T, percentage='column',pvalue=TRUE)



#survival for trt1 and trt2

#convert months to years
gmmbc$os<-gmmbc$os/12
gmmpc$os<-gmmpc$os/12
fitb <- survfit(Surv(os,os_event) ~ cluster, data = gmmbc)
ggsurvplot(fitb, data = gmmbc,xlab = "Years",pval=TRUE,xlim=c(0,3),risk.table = TRUE,legend.labs = c("Cluster 1", "Cluster 2"),risk.table.height = 0.35) +labs(title= "KM plot for treatment arm 1 (brivanib and cetuximab)")


fit1b <- survfit(Surv(os,os_event) ~ cluster, data = gmmpc)
ggsurvplot(fit1b, data = gmmpc,xlab = "Years",pval=TRUE,xlim=c(0,3),risk.table = TRUE,legend.labs = c("Cluster 1", "Cluster 2"),risk.table.height = 0.35) +labs(title= "KM plot for treatment arm 2 (placebo and cetuximab)")


