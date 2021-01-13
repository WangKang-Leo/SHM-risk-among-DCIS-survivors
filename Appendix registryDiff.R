#############appendix-figure1-registryDiff###########
library(SEERaBomb);library(ggplot2);library(dplyr)
load("E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER/mrgd/popsae.RData")
load("E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER/mrgd/cancTumSzCRT.RData")
canc$sex
(d=canc%>%filter(sex=="Female",cancer=="breastCIS",agedx<40))
(d=d%>%filter(histo3=="8201"|histo3=="8230"|histo3=="8500"|histo3=="8501"|
                histo3=="8502"|histo3=="8503"|histo3=="8504"|histo3=="8507"|
                histo3=="8523"))
d=d%>%filter(is.na(cstumsiz)|(cstumsiz!="990"&cstumsiz!="997"&cstumsiz!="998"))
d=d%>%filter(is.na(eod10ex)|eod10ex!=5,is.na(eod10sz)|eod10sz!=997,is.na(eod10sz)|eod10sz!=998,is.na(tumsizs)|tumsizs!=998)
d=d%>%filter(seqnum==0|seqnum==1) 
d=d%>%filter(chemo=="0")   #191983
d=d%>%filter(trt%in%c("nr","eb")) #no and external beam  184988
####follow up#####  !!!!!!184363table(is.na(d$CODS))
d=d%>%filter(!is.na(d$CODS))

b=canc%>%filter(seqnum==2,sex=="Female")

df<-bind_rows(d,b)

papFolICD=c(8201,8230,8500:8504,8507,8523) 
canc$cancer=as.character(canc$cancer)
canc$sex=as.character(canc$sex)
canc$cancer[canc$seqnum%in%c("0","1") & !is.na(canc$CODS) & canc$trt%in%c("nr","eb") & canc$sex=="Female" & canc$histo3%in%papFolICD & canc$cancer=="breastCIS" & canc$chemo=="0"& (is.na(canc$cstumsiz)|(canc$cstumsiz!="990"&canc$cstumsiz!="997"&canc$cstumsiz!="998"))&(is.na(canc$eod10ex)|canc$eod10ex!=5)& (is.na(canc$eod10sz)|canc$eod10sz!=997)&(is.na(canc$eod10sz)|canc$eod10sz!=998)& (is.na(canc$tumsizs)|canc$tumsizs!=998)]="DCIS" 
table(canc$cancer)
canc$cancer=factor(canc$cancer)

secs=c("AML","ALL","CML","CLL","HL","NHL","MM")
df=p2s(canc,firstS="DCIS",secondS=secs,yrcut=1975) 
table(df$status) #####secondary tumor after DCIS#########
(O=df%>%group_by(trt)%>%summarize(sum(status)))
df=as.data.frame(df)
write.csv(df,file ="E:/Ph.D projects/3.SEER/1.SEERaBomb/figure/Appendix Figure 1/SEERaBomb.csv")
library(plyr);library(data.table)
setwd("E:/Ph.D projects/3.SEER/1.SEERaBomb/figure/Appendix Figure 1")
mydata=as.data.frame(fread("appendix figure1.csv"))
mydata=mydata%>%filter(HM=="NHL")
df <- ddply(mydata,.(registry),transform,len=length(as.factor(YD)))
ggplot(df,aes(x=YD,color=registry)) + geom_step(aes(len=len,y=..y.. * len),stat="ecdf",size=2)+theme_classic()+
  labs(x="Year of diagnosis",y="Number of NHL cases")+ggtitle("HL cases") +theme(legend.position=c(.2, .7)) ##5*3#
library(SEERaBomb)
