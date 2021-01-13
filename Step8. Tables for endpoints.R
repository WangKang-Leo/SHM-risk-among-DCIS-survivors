#########Tables for endpoint###
####HM###############
load("E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER/mrgd/popsae.RData")
load("E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER/mrgd/cancTumSzCRT.RData")
papFolICD=c(8201,8230,8500:8504,8507,8523) 
canc$cancer=as.character(canc$cancer)
canc$sex=as.character(canc$sex)
canc$cancer[canc$seqnum%in%c("0","1") & !is.na(canc$CODS) & canc$trt%in%c("nr","eb") & canc$sex=="Female" & canc$histo3%in%papFolICD & canc$cancer=="breastCIS" & canc$chemo=="0"& (is.na(canc$cstumsiz)|(canc$cstumsiz!="990"&canc$cstumsiz!="997"&canc$cstumsiz!="998"))&(is.na(canc$eod10ex)|canc$eod10ex!=5)& (is.na(canc$eod10sz)|canc$eod10sz!=997)&(is.na(canc$eod10sz)|canc$eod10sz!=998)& (is.na(canc$tumsizs)|canc$tumsizs!=998)]="DCIS" 
table(canc$cancer)
secs=c("AML","ALL","CML","CLL","HL","NHL","MM")
a=p2s(canc,firstS="DCIS",secondS=secs,yrcut=1970) 
a$cancer2
write.csv(a,file="E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/brca_sec.csv")

mydata=as.data.frame(fread("E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/sec.csv"))
secs=c("ALL","AML","CLL","CML","MM","HL","NHL")
mydata=subset(mydata,cancer2!=c("ALL","AML","CLL","CML","MM","HL","NHL","breast","breastCIS","MDS","MPN","CMML"))
mydata=mydata%>%filter(cancer2%in%c("breast","breastCIS"))
library(CBCgrps)
mydatasub<-as.data.frame(mydata[,c("trt","surv")])
mydatasub$rad=as.factor(mydatasub$trt)
twogrps(mydatasub,"trt")

mydata=as.data.frame(fread("E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/brca_sec.csv"))

mydata$yrg=cut(mydata$yrdx,breaks=c(1975,1980,1990,2000,2010,2017),right=F,dig.lab=4)
mydata$ageg=cut(mydata$age86,breaks=c(0,30,40,50,60,70,120),right=F)

mydata$er=NA
mydata$er[mydata$erstatus==1]=1
mydata$er[mydata$erstatus==2]=0

mydata$pr=NA
mydata$pr[mydata$prstatus==1]=1
mydata$pr[mydata$prstatus==2]=0

mydata$Race=NA
mydata$Race[mydata$race=="White"]=1
mydata$Race[mydata$race=="Black"]=2
mydata$Race[mydata$race=="Other"]=3

mydata$Grade=NA
mydata$Grade[mydata$grade=="1"]=1
mydata$Grade[mydata$grade=="2"]=2
mydata$Grade[mydata$grade%in%c(3,4)]=3

mydata$surgery=NA
mydata$surgery[(mydata$surgprif == 0) | (mydata$sssurg==0)]=0 #保乳=1
mydata$surgery[(mydata$surgprif %in% c(20:24)) | (mydata$sssurg %in% c(10,20))]=1 #保乳=1
mydata$surgery[(mydata$surgprif %in% c(30:80)) | (mydata$sssurg %in% c(30,80))]=2

mydata$rad[mydata$trt=="eb"]=1
mydata$rad[mydata$trt=="nr"]=0

mydata$group=NA
mydata$group[mydata$trt=="nr"&mydata$cancer2=="HNL"]=1
mydata$group[mydata$trt=="eb"&mydata$cancer2=="NHL"]=2
mydata$group[mydata$trt=="nr"&(mydata$cancer2!="NHL"|is.na(mydata$cancer2))]=3
mydata$group[mydata$trt=="eb"&(mydata$cancer2!="NHL"|is.na(mydata$cancer2))]=4
table(mydata$group)

library(tableone)
mydatasub=subset(mydata,group%in%c(2,4))
vars=c("ageg","yrg","Race","tumsz","Grade","er","pr","surgery","surv")
catvars=c("ageg","yrg","Race","tumsz","Grade","er","pr","surgery")
tableOne <- CreateTableOne(vars=vars,strata = c("group"),factorVars = catvars,includeNA=T,data = mydatasub)
tab_out<-print(tableOne)
write.csv(tab_out,file="E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/NHL_table.csv")
library(CBCgrps)
head(mydata)
mydatasub<-as.data.frame(mydata[,c("group","surv")])
mydatasub=subset(mydatasub,group==2|group==4)
mydatasub$surv=mydatasub$surv/12
mydatasub$group=as.factor(mydatasub$group)
twogrps(mydatasub,"group")
