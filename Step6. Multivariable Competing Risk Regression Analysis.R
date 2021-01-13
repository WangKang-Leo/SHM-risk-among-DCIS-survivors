#####################################
################regression###########
####################################
library(tidyverse);library(SEERaBomb);library(ggsci);library(survival);library(survminer);library(bbmle);library(data.table)
library(cmprsk);library(tableone)
thyall=as.data.frame(fread("E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/brca_sec.csv"))
thyall$er=3
thyall$er[thyall$erstatus==1]=1
thyall$er[thyall$erstatus==2]=0
thyall$er=as.factor(thyall$er)

thyall$pr=3
thyall$pr[thyall$prstatus==1]=1
thyall$pr[thyall$prstatus==2]=0
thyall$pr=as.factor(thyall$pr)

thyall$Race[thyall$race=="White"]=1
thyall$Race[thyall$race=="Black"]=2
thyall$Race[thyall$race=="Other"]=3
thyall$Race=as.factor(thyall$Race)

thyall$Grade=4
thyall$Grade[thyall$grade=="1"]=1
thyall$Grade[thyall$grade=="2"]=2
thyall$Grade[thyall$grade%in%c(3,4)]=3
thyall$Grade=as.factor(thyall$Grade)

thyall$surgery=3
#thyall$surgery[(thyall$surgprif == 0) | (thyall$sssurg==0)]=3 #保乳=1
thyall$surgery[(thyall$surgprif %in% c(20:24)) | (thyall$sssurg %in% c(10,20))]=1 #保乳=1
thyall$surgery[(thyall$surgprif %in% c(30:80)) | (thyall$sssurg %in% c(30,80))]=2
thyall$surgery=as.factor(thyall$surgery)

thyall$tumsz[is.na(thyall$tumsz)]=4
thyall$tumsz=as.factor(thyall$tumsz)

thyall$rad[thyall$trt=="eb"]=1
thyall$rad[thyall$trt=="nr"]=0

thyall$agescore=NA
thyall$agescore[thyall$agedx>=60]=0
thyall$agescore[thyall$agedx<60 & thyall$agedx>=40]=1
thyall$agescore[thyall$agedx<40]=2

thyall$gradescore=NA
thyall$gradescore[thyall$grade==1]=0
thyall$gradescore[thyall$grade==2]=1
thyall$gradescore[thyall$grade==3|thyall$grade==4]=2

thyall$riskscore=NA
thyall$riskscore=thyall$size_score+thyall$gradescore+thyall$agescore
###################
thyall$tumszII=NA
thyall$tumszII[thyall$trt=="nr"]=0
thyall$tumszII[thyall$trt=="eb"]=1
thyall$tumszII[thyall$trt=="eb" & (thyall$agedx<50 | thyall$grade%in%c(3,4))]=2
table(thyall$tumszII)

#secs=c("ALL","AML","CLL","CML","MM","HL","NHL")
thyall$AML.crr.status=NA
thyall$AML.crr.time=NA
thyall$AML.crr.status=0
thyall$AML.crr.status[thyall$c2occ==1]=2
thyall$AML.crr.status[thyall$CODS!="alive"]=2
thyall$AML.crr.status[thyall$cancer2=="CML"&thyall$surv>=12]=1
thyall$AML.crr.status[thyall$cancer2%in%c("ALL","AML","CLL","NHL","MM","HL","CMMML","MDS","MPN")]=0
thyall$AML.crr.status[thyall$surv>240]=0
#thyall$AML.crr.status[thyall$CODS=="alive" & thyall$status==0]=0
thyall$AML.crr.time=thyall$surv
thyall$AML.crr.time[thyall$surv>240]=240
thyall=thyall%>%filter(Grade%in%c(3)|agedx<50)
table(is.na(thyall$surv))
submydata=thyall[,c("AML.crr.status","AML.crr.time","yrdx","agedx","Race","Grade","tumsz","er","pr","surgery","rad")]
#submydata=submydata%>%filter(!is.na(submydata$er),!is.na(submydata$surgery))
cov<-model.matrix(~as.factor(rad), data= submydata)[,-1]
summary(crr(submydata$AML.crr.time, submydata$AML.crr.status,cov1=cov))
