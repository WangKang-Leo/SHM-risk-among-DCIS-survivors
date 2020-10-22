R scripts for extracting and generating data from SEER files.
library(tidyverse);library(SEERaBomb);library(ggsci)#load packages         
# Functions will be explained in comments at first occurrence in this script (so, duplicate functions in later figures will not be explained again).
# SEERaBomb is maintained on github. You can install the current github version using:
install.packages("devtools") #if it isn't already installed
library(devtools)
install_github("radivot/SEERaBomb",subdir="SEERaBomb", force=TRUE)

# As a backup to this github/devtools approach, every ~6 months Windows SEERaBomb binary is produced via:
install.packages(c("LaF","RSQLite","dplyr","XLConnect","Rcpp","rgl","reshape2","mgcv","DBI","bbmle")) # first get dependencies from CRAN
install.packages("SEERaBomb",repos="http://epbi-radivot.cwru.edu") #need line above since only SEERaBomb is in the repository

# Finally, once each year, after testing it against the most recent SEER data release (typically in April-May), a novel version of SEERaBomb is uploaded to CRAN. The CRAN version should be stable. It can be installed via:
install.packages("SEERaBomb")
# If important bugs are reported I will fix and update the CRAN version, otherwise, to install new versions with minor bug fixes and/or new features, please use one of the other approaches. I will try to assure that the CRAN version works on indows, mac and ubuntu.
# After installation SEERaBomb help pages can be reached by: library(SEERaBomb);help(pack="SEERaBomb")

# First make a database named cancTumSz.RData using SEERaBomb's function "mkSEER".
# cancTumSz differs from cancDef (SEERaBomb's pickFields() default) in that in includes tumor size and stage data.
rm(list=ls())
library(dplyr)
library(SEERaBomb)
library(dbplyr)
(df=getFields(seerHome="E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER"))
head(df,60)

c=c("casenum","firstprm","reg","sex","histo3","beho3","histo2",
    "eod2","eod4","eod10sz","eodcode", "cstumsiz",           ##"adjajccstg","adjtm6value","adjnm6value","adjm6value",
    "surgprif","sssurg","numnodes","nosurg","seqnum", "cstseval",
   "eod10ex","radiatn","radsurg","race","marstat","agedx",
    "yrbrth","modx","yrdx","lateral","tumsizs","grade","dxconf","reptsrc","eod10pn","erstatus",
    "prstatus","her2","ICD9","ICD10","statrec","COD","surv","chemo")
(rdf=pickFields(df,picks=c))   ##   #"surgscof","adjtm6value","adjnm6value","adjm6value","adjajccstg",
mkSEER(rdf,seerHome="E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER") # This merges all cancer binaries in SEER data folder mrgd.

# Now, we will categorize tumor sizes. eod10, eod4 and cstumsiz definitions were derived from the SEER data 
library(survival);library(survminer);library(bbmle)
load("E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER/mrgd/popsae.RData")
load("E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER/mrgd/cancDef.RData") #loads data.frame canc into memory 
##########pre-processing############
canc$obs <- 1:nrow(canc) # we will need this below.
thy=canc%>%filter(cancer=="breastCIS")
thy$eod4sz=substr(thy$eod4,1,nchar(thy$eod4)-2)
#thy$eod4sz=substr(thy$eod4,1,nchar(thy$eod4)-2) # In eod4, the first 2 digits are the size, so discard the rest.
thy_cs <- thy[is.na(thy$eodcode),] # Make separate data frames to makes things easier below, we will join (rbind) them later.
thy_eod10=thy%>%filter(eodcode==4)
thy_eod4=thy%>%filter(eodcode==3)

unique(thy_cs$cstumsiz)
thy_cs=thy_cs%>%filter(!(thy_cs$cstumsiz>=201 & thy_cs$cstumsiz<=989),thy_cs$cstumsiz!=0,thy_cs$cstumsiz!=998,thy_cs$cstumsiz!=997,thy_cs$cstumsiz!=999) #201-989 mm is too large and probably reflect data entry errors. 000 = no tumor. 999 = missing tumour size. We discard those. 
thy_cs$size_score=NA
thy_cs$size_score[thy_cs$cstumsiz%in%c(001:015,990,991)]=0 # <=15 mm.
thy_cs$size_score[thy_cs$cstumsiz%in%c(016:040,993,994)]=1 # 16-40mm.
thy_cs$size_score[thy_cs$cstumsiz%in%c(041:200,995,996)]=2 # >=41

thy_cs$cstumsiz[thy_cs$cstumsiz%in%c(001:009,990,991)]=0 # <=9 mm.
thy_cs$cstumsiz[thy_cs$cstumsiz%in%c(010:019,992)]=1 # 10-19mm.
thy_cs$cstumsiz[thy_cs$cstumsiz%in%c(020:049,993,994,995)]=2 # 20-49
thy_cs$cstumsiz[thy_cs$cstumsiz%in%c(050:200,996)]=3

unique(thy_eod10$eod10sz)
thy_eod10=thy_eod10%>%filter(!(thy_eod10$eod10sz>=201 & thy_eod10$eod10sz<=999)) # Same as with cstumsiz. We discard those as well.
thy_eod10$size_score=NA
thy_eod10$size_score[thy_eod10$eod10sz%in%c(001:015)]=0
thy_eod10$size_score[thy_eod10$eod10sz%in%c(016:040)]=1
thy_eod10$size_score[thy_eod10$eod10sz%in%c(041:200)]=2

thy_eod10$eod10sz[thy_eod10$eod10sz%in%c(001:009)]=0
thy_eod10$eod10sz[thy_eod10$eod10sz%in%c(010:019)]=1
thy_eod10$eod10sz[thy_eod10$eod10sz%in%c(020:049)]=2
thy_eod10$eod10sz[thy_eod10$eod10sz%in%c(050:200)]=3

thy_eod4=thy_eod4%>%filter(thy_eod4$eod4sz!=99,thy_eod4$eod4sz!=98,thy_eod4$eod4sz!=00) # 00 = no tumour. 99 = missing size. We discard those as well.
thy_eod4=thy_eod4%>%filter(!(is.na(thy_eod4$eod4sz) | thy_eod4$eod4sz==""))
thy_eod4$size_score=NA
thy_eod4$size_score[thy_eod4$eod4sz%in%c(01:15)]=0
thy_eod4$size_score[thy_eod4$eod4sz%in%c(16:40)]=1
thy_eod4$size_score[thy_eod4$eod4sz%in%c(41:97)]=2

thy_eod4$eod4sz[thy_eod4$eod4sz%in%c(01:10)]=0
thy_eod4$eod4sz[thy_eod4$eod4sz%in%c(10:19)]=1
thy_eod4$eod4sz[thy_eod4$eod4sz%in%c(20:49)]=2
thy_eod4$eod4sz[thy_eod4$eod4sz%in%c(49:97)]=3

thy_cs$tumsz=thy_cs$cstumsiz # Rename the new size columns to the new column "tumsz" before rbinding.
thy_eod10$tumsz=thy_eod10$eod10sz
thy_eod4$tumsz=thy_eod4$eod4sz

thySz=rbind(thy_cs,thy_eod10,thy_eod4) # Now, we have 1 large database with all tumours with known tumour sizes again.
unique(thySz$size_score)
thySz=thySz%>%select(-eod4sz) # Get rid of columns we don't need anymore.

canc$tumsz=NA
canc$size_score=NA
obsunique=intersect(canc$obs,thySz$obs) # We use column "obs" to prevent duplicates below. 
cancnoSz=subset(canc,!(canc$obs%in%obsunique))
cancTumSz=rbind(cancnoSz,thySz)
canc=cancTumSz%>%select(-obs)
unique(canc$size_score)

# We also make some columns that we will need throughout the script.
canc$trt="nr" # nr = no radiation. Will be left as 0 (no radiation) and 7 (radiation recommended but refused). Do it this way to initialize the vector.
canc$trt[canc$radiatn%in%c(8)]="uk" # uk = unknown.
canc$trt[canc$radiatn%in%c(1,4)]="eb" # eb = external beam radiotherapy.
canc$trt[canc$radiatn%in%c(2,3,5,6)]="ii" # ii = radioactive isotopes. In the case of WDTC, these all concern radioactive iodine (RAI) but most will be "3".
canc$trt=factor(canc$trt,levels=c("nr","ii","eb","uk"))

# To conform to the WHO 2016 classification, we will change the ICD-O-3 codes that are in canc$cancer==MPN.
canc$cancer=as.character(canc$cancer)
canc$cancer[canc$cancer=="MPN"]="CMML" # CMML is a placeholder for now, we will put all correct ICD-O-3 codes in MPN below.
MPNWHO2016=c(9740:9742,9876,9945,9946,9950,9960:9964,9975) # These are the ICD-O-3 codes for MPN according to the WHO 2016 classification.
canc$cancer[canc$cancer=="CMML" & canc$histo3%in%MPNWHO2016]="MPN" # Other codes are left as CMML but will not be used.

canc$yrmodx=(round(canc$yrdx+((canc$modx-0.5)/12),3)) # Make yrmodx to improve resolution over yrdx: 2014.0 is January 2014, 2014.5 is July 2014.
canc$survy=round(((canc$surv+0.5)/12),3)
canc$cancer=fct_collapse(canc$cancer,AML=c("AML","AMLti","APL"))
canc$cancNo="cg2" # Do this to initialize the vector, cg2 = greater than 2.
canc$cancNo[canc$seqnum==2]="c2"
canc$cancNo[canc$seqnum<2]="c1"
canc$cancNo=factor(canc$cancNo,levels=c("c1","c2","cg2"))
##########################################
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
######################################################
#################Starting Analyses####################
load("E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER/mrgd/popsae.RData")
load("E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER/mrgd/cancTumSzCRT.RData")

d=canc%>%filter(sex=="Female",cancer=="breastCIS")
####DCIS#############241877
(d=d%>%filter(histo3=="8201"|histo3=="8230"|histo3=="8500"|histo3=="8501"|
                histo3=="8502"|histo3=="8503"|histo3=="8504"|histo3=="8507"|
                histo3=="8523"))
####extention######241523
d=d%>%filter(is.na(cstumsiz)|(cstumsiz!="990"&cstumsiz!="997"&cstumsiz!="998"))
d=d%>%filter(is.na(eod10ex)|eod10ex!=5,is.na(eod10sz)|eod10sz!=997,is.na(eod10sz)|eod10sz!=998,is.na(tumsizs)|tumsizs!=998)
#####bilateral###

##d=d%>%filter(lateral!=4)####
####BCS######### 154117
#d=d%>%filter((surgprif>19&surgprif<25)|(sssurg==10|sssurg==20))
table(d$sssurg)
table(d$surgprif)
####first primary### 193280
d=d%>%filter(seqnum==0|seqnum==1)   #0：One primary only in the patient’s lifetime,1：First of two or more primaries
####radiation#####
unique(d$trt)
d=d%>%filter(chemo=="0")   #191983
d=d%>%filter(trt%in%c("nr","eb")) #no and external beam  184988
####follow up#####  !!!!!!184363
table(is.na(d$CODS))
d=d%>%filter(!is.na(d$CODS))

mydata=d
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
library(tableone)
vars=c("agedx","yrdx","Race","tumsz","Grade","er","pr","surgery","surv")
catvars=c("Race","tumsz","Grade","er","pr","surgery")
tableOne <- CreateTableOne(vars=vars,strata = c("trt"),factorVars = catvars,includeNA=T,data = mydata)
tab_out<-print(tableOne, catDigits = 1, contDigits = 2, pDigits = 3,quote = T, missing = T, explain = TRUE, printToggle = TRUE,test = TRUE, smd = T, noSpaces = FALSE, padColnames = FALSE,varLabels = FALSE, format = c("fp", "f", "p", "pf")[1],showAllLevels = T, cramVars = NULL, dropEqual = FALSE,exact = NULL, nonnormal = NULL, minMax = FALSE)

library(CBCgrps)
head(mydata)
mydatasub<-as.data.frame(mydata[,c("rad","agedx","yrdx","Race","tumsz","Grade","er","pr","surgery","surv")])
mydatasub$rad=as.factor(mydatasub$rad)
a<-twogrps(mydatasub,"rad")
write.csv(tab_out,file="E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/table1.csv")

#####low risk####
canc$agescore=NA
canc$agescore[canc$agedx>=60]=0
canc$agescore[canc$agedx<60 & canc$agedx>=40]=1
canc$agescore[canc$agedx<40]=2

canc$gradescore=NA
canc$gradescore[canc$grade==1]=0
canc$gradescore[canc$grade==2]=1
canc$gradescore[canc$grade==3|canc$grade==4]=2

canc$riskscore=NA
canc$riskscore=canc$size_score+canc$gradescore+canc$agescore

table(canc$riskscore,canc$trt)
d=canc%>%filter(riskscore%in%c(0:4)) ##low/intermediate-risk
write.csv(d,file="E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/brca_primarylowrisk.csv")
#################################
######################
load("E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER/mrgd/popsae.RData")
load("E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER/mrgd/cancTumSzCRT.RData")
secs=c("ALL","AML","CLL","CML","MM","HL","NHL")
(p=canc%>%filter(cancer%in%c(secs),sex=="Female"))
#p <-p%>%filter(!is.na(reg),!is.na(db),!is.na(agedx),!is.na(sex),!is.na(race),!is.na(yrdx),!is.na(modx),!is.na(surv))
c<-bind_rows(d,p)
a=p2s(c,firstS="breastCIS",secondS=secs,yrcut=1975)
#########dA including O/E secondary HM#######
papFolICD=c(8201,8230,8500:8504,8507,8523) 
canc$cancer=as.character(canc$cancer)
canc$sex=as.character(canc$sex)
canc$cancer[canc$seqnum%in%c("0","1") & !is.na(canc$CODS) & canc$trt%in%c("nr","eb") & canc$sex=="Female" & canc$histo3%in%papFolICD & canc$cancer=="breastCIS" & canc$chemo=="0"& (is.na(canc$cstumsiz)|(canc$cstumsiz!="990"&canc$cstumsiz!="997"&canc$cstumsiz!="998"))&(is.na(canc$eod10ex)|canc$eod10ex!=5)& (is.na(canc$eod10sz)|canc$eod10sz!=997)&(is.na(canc$eod10sz)|canc$eod10sz!=998)& (is.na(canc$tumsizs)|canc$tumsizs!=998)]="DCIS" 
table(canc$cancer)
canc$cancer=factor(canc$cancer)

##low-risk##
canc$cancer[canc$seqnum%in%c("0","1") & !is.na(canc$CODS) & canc$trt%in%c("nr","eb") & canc$sex=="Female" & canc$histo3%in%papFolICD & canc$cancer=="breastCIS" & canc$chemo=="0"& (is.na(canc$cstumsiz)|(canc$cstumsiz!="990"&canc$cstumsiz!="997"&canc$cstumsiz!="998"))&(is.na(canc$eod10ex)|canc$eod10ex!=5)& (is.na(canc$eod10sz)|canc$eod10sz!=997)&(is.na(canc$eod10sz)|canc$eod10sz!=998)& (is.na(canc$tumsizs)|canc$tumsizs!=998)& canc$agedx>=50]="DCIS" 
canc$cancer[canc$cancer%in%secs]="SHM"
table(canc$cancer)
pf=seerSet(canc,popsae,Sex="Female",ageStart=0,ageEnd=100)#pooled (races) females 
pf=mk2D(pf,secondS=secs)#adds secs (commons.R) background rates to pf
#plot2D(pf)
trts=c("nr","eb")
pf=csd(pf,brkst=c(0,1,3,6,10),trts=trts,firstS="DCIS",exclUnkSurv=FALSE)
dA=pf$DF

write.csv(dA,file="E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/OandE(0,1,3,6,10).csv")
########Figure#######
##############Figure3A and B################
library(tidyverse);library(SEERaBomb);library(ggsci);library(survival);library(survminer);library(bbmle);library(data.table)
dA=as.data.frame(fread("E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/OandE(0,1,3,6,10).csv"))
dA=dA%>%filter(cancer2=="HL")
gp=geom_point();gl=geom_line()
geRR=geom_errorbar(aes(ymin=rrL,ymax=rrU),width=.2)
gh=geom_hline(yintercept=1)
sbb=theme(strip.background=element_blank())
ltb=theme(legend.margin=margin(0,0,0,0),legend.title=element_blank())
ltp=theme(legend.position="top")
lh=theme(legend.direction="horizontal")
sy=scale_y_log10()
jco=scale_color_jco()
tc=function(sz) theme_classic(base_size=sz);
gxi=xlab("Age (Years)")
gyi=ylab(quote(paste("Cases per ",10^5," Person Years")))
gy=ylab("RR of Developing HL")
x=scale_x_continuous(breaks=c(0,1,3,6,10))#surv times
y=scale_y_continuous(breaks=c(0,1,2,3,4,5,6))#surv times
myt=theme(legend.key.height=unit(.25,'lines'),legend.position=c(.5,.95))
cc=coord_cartesian(ylim=c(0,6))#clips high errorbars
gx=xlab("Years Since DCIS Diagnosis (years)")
xlim=xlim(0,10)
dA%>%ggplot(aes(x=t,y=RR,col=trt,width=.15))+gp+gl+gx+gy+gh+geRR+tc(14)+ltp+jco+sbb+ltb+lh+cc+x+y

P1=dA%>%ggplot(aes(x=t,y=RR,col=trt,width=.15))+gp+gl+gx+gy+gh+geRR+tc(14)+ltp+jco+sbb+ltb+lh+x+scale_y_continuous(breaks=c(0,1,2,4,6,8))+coord_cartesian(ylim = c(0,8))
P2=dA%>%ggplot(aes(x=t,y=RR,col=trt,width=.15))+gp+gl+gh+geRR+tc(14)+ltp+jco+sbb+ltb+x+scale_y_continuous(breaks=c(8,12,16))+coord_cartesian(ylim = c(8,16))+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) 
ggarrange(P2,P1,heights=c(1/3,2/3),ncol = 1, nrow = 2) 

thyall=as.data.frame(fread("E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/brca_sec.csv"))
thyall=thyall%>%filter(agedx<=65)
thyall=thyall%>%filter(tumsz<2)
thyall$mds2=0
thyall$mds2[thyall$cancer2=="CLL"]=1
thyall$SurvObj=with(thyall, Surv(surv/12, mds2==1))

#######
secs=c("ALL","AML","CLL","CML","MM","HL","NHL")
thyall$mds2=0
thyall$mds2[thyall$cancer2%in%secs]=1
thyall$SurvObj=with(thyall, Surv(surv/12, mds2==1))
#####
summary(coxph(SurvObj ~ trt, data = thyall))
fit <- survfit(SurvObj ~ trt, data = thyall) 
ggsurvplot(fit,
           log.rank.weights = "n",
           ylim=c(0,0.004),
           xlim=c(0,20),
           xlab="Time After DCIS Diagnosis (Years)",
           ylab="DCIS Patients With CLL(%)",
           break.x.by=5,
           break.y.by=0.001,
           surv.scale="percent",
           conf.int =F,
           cumevents = F,
           linetype = "trt",
           risk.table.y.text=F,
           palette = c("#0073C2FF","#E7C000FF"),
           fun = "cumhaz")+theme_survminer(
             font.main = c(14),
             font.submain = c(14),
             font.caption = c(14),
             font.x = c(14),
             font.y = c(14),
             font.tickslab = c(11)
           )
ggsurvplot(fit,pval = TRUE,pval.size=12,
           log.rank.weights = "n")
###########AGE vs RR##########
load("E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER/mrgd/popsae.RData")
load("E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER/mrgd/cancDef.RData")

d=incidSEER(canc,popsae,"ALL")
d=d%>%filter(age<=85,year>=2000)
d=d%>%mutate(ageG=cut(age,seq(0,85,5)))
d=d%>%group_by(cancer,ageG)%>%
  summarize(age=mean(age),py=sum(py),n=sum(n))%>%
  mutate(incid=n/py,grp="Background")
d=d%>%select(cancer,grp,everything(),-ageG)#reorder columns
#the next 3 lines define a set of non-heme-malignacies
HM=c(secs,"MDS","CMML","MPN","CLL","HCL","OL","NHL","MM","HL","LGL")
sc=canc%>%filter(agedx<50,yrdx<2000)#based on frequency at younger ages
(NHM=setdiff(names(sort(table(sc$cancer),decr=T)[1:8]),HM))#non-hemes
brksa=c(0,40,50,60,70,80)#broad 1st interval avoids 0 CML groups
system.time(D<-riskVsAge(canc,firstS="breastCIS",secondS="ALL",brksa=brksa))#~36s
D=D%>%filter(rad!="Unk",chemo!="Unk")
D=D%>%group_by(cancer2,rad,chemo,age)%>%
  summarize(py=sum(py),n=sum(o),incid=n/py)
D=D%>%rename(cancer=cancer2)%>%unite(grp,rad,chemo,sep=", ")
dd=bind_rows(D,d)%>%filter(age>=20)
ord=c("Rad, Chemo","No Rad, Chemo","Rad, No Chemo",
      "No Rad, No Chemo","Background")
dd$grp=factor(dd$grp,levels=ord)
dd$cancer=factor(dd$cancer,levels=secs)
myt=theme(legend.position=c(.7,.83),legend.key.height=unit(.65,'lines'))
dd=dd%>%mutate(LL=qpois(0.025,n)/py,UL=qpois(0.975,n)/py)#make CI
dd=dd%>%mutate(ages=age+(as.numeric(grp)-3))#shift ages to see CI
dd%>%ggplot(aes(x=ages,y=incid,col=grp))+gl+facet_grid(~cancer)+gxi+
  agts+gyi+jco+sy+tc(11)+ltb+sbb+myt+ge
ggsave("~/Results/tutorial/ageTherapyEx4.pdf",width=3.5,height=2.5)

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

##############PSM-DCIS################
library(tidyverse);library(SEERaBomb);library(ggsci);library(data.table);library(tableone)#load packages          
library(survival);library(survminer);library(bbmle);library(sqldf);library(MatchIt);library(nonrandom);library(Matching)
setwd("E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort")
df=as.data.frame(fread("brca_sec.csv"))
df$er=2
df$er[df$erstatus==1]=1
df$er[df$erstatus==2]=0

df$pr=2
df$pr[df$prstatus==1]=1
df$pr[df$prstatus==2]=0

df$Race[df$race=="White"]=1
df$Race[df$race=="Black"]=2
df$Race[df$race=="Other"]=3

df$Grade=4
df$Grade[df$grade=="1"]=1
df$Grade[df$grade=="2"]=2
df$Grade[df$grade%in%c(3,4)]=3

df$surgery=3
df$surgery[(df$surgprif == 0) | (df$sssurg==0)]=0 #保乳=1
df$surgery[(df$surgprif %in% c(20:24)) | (df$sssurg %in% c(10,20))]=1 #保乳=1
df$surgery[(df$surgprif %in% c(30:80)) | (df$sssurg %in% c(30,80))]=2


table(df$status)
df$tumor_size=4
df$tumor_size[df$tumsz%in%c(0,1)]=1
df$tumor_size[df$tumsz%in%c(2)]=2
df$tumor_size[df$tumsz%in%c(3)]=3

df$status=1
df$status[df$CODS=="alive"]=0

df$c2occ=as.factor(df$c2occ)
df$yrdx=as.numeric(df$yrdx)
df$agedx=as.numeric(df$agedx)
df$Race=as.factor(df$Race)
df$tumor_size=as.factor(df$tumor_size)
df$Grade=as.factor(df$Grade)
df$er=as.factor(df$er)
df$pr=as.factor(df$pr)
df$surgery=as.factor(df$surgery)

data=df[,c("trt","c2occ","status","cancer2","yrdx","agedx","Race","tumor_size","Grade","er","pr","surgery","survy")]

ALL_RT=data%>%filter(c2occ=="0"|cancer2=="NHL",trt=="eb")

ALL_nonRT=data%>%filter(c2occ=="0"|cancer2=="NHL",trt=="nr")
sapply(ALL_nonRT, function(x) sum(is.na(x)))

ALL_RT=ALL_RT[,c("c2occ","yrdx","agedx","Race","tumor_size","Grade","er","pr","surgery","status","survy","trt")]
ALL_nonRT=ALL_nonRT[,c("c2occ","yrdx","agedx","Race","tumor_size","Grade","er","pr","surgery","status","survy","trt")]
ALL_nonRT=na.omit(ALL_nonRT)

f=matchit(c2occ~yrdx+agedx+Race+tumor_size+Grade+er+pr+surgery,data=ALL_RT,method="nearest",discard="none",
          caliper=0.5,reestimate=F,
          ratio=5)
summary(f)
newdata <- match.data(f)
stable1 <- CreateTableOne(vars=c("yrdx","agedx","Race","tumor_size","Grade","er","pr","surgery"),strata="c2occ", data=newdata,factorVars=c("Race","tumor_size","Grade","er","pr","surgery"))
print(stable1,showAllLevels = TRUE)

f1=matchit(c2occ~yrdx+agedx+Race+tumor_size+Grade+er+pr+surgery,data=ALL_nonRT,method="nearest",discard="none",
           reestimate=T,caliper=0.01,
           ratio=5)
ALL_nonRT <- match.data(f1)

ALL=bind_rows(newdata,ALL_nonRT)
write.csv(ALL,file="E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/PSM/DCIS_NHL.csv")
################################
###############################
df=as.data.frame(fread("E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/PSM/DCIS_NHL.csv"))
df$group[df$c2occ==1 & df$trt=="nr"]=1
df$group[df$c2occ==0 & df$trt=="nr"]=2
df$group[df$c2occ==1 & df$trt=="eb"]=3
df$group[df$c2occ==0 & df$trt=="eb"]=4

tapply(df$survy, df$group, FUN=quantile)
tapply(df$yrdx, df$group, FUN=quantile)

library(tableone)
vars=c("yrdx","agedx","Race","tumor_size","Grade","er","pr","surgery","survy")
catvars=c("Race","tumor_size","Grade","er","pr","surgery")
tableOne <- CreateTableOne(vars=vars,strata = c("group"),factorVars = catvars,includeNA=T,data = df)
tab_out<-print(tableOne)
write.csv(tab_out,file="E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/PSM/DCIS_NHL_table.csv")
#######PSM-HM##########
library(tidyverse);library(SEERaBomb);library(ggsci);library(data.table);library(tableone)#load packages          
library(survival);library(survminer);library(bbmle);library(sqldf);library(MatchIt);library(nonrandom);library(Matching)
load("E:/Ph.D projects/3.SEER/1.SEERaBomb/data/SEER/mrgd/cancTumSzCRT.RData")
setwd("E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort")
df=as.data.frame(fread("brca_sec.csv"))

mydata=df%>%filter(cancer2%in%c("ALL"))
mydata=mydata[,c("casenum","trt")]

HM=canc%>%filter(cancer%in%c("ALL"),seqnum%in%c("2"),sex=="Female")
HM=HM[,c("casenum","cancer","histo3","agedx","yrdx","race","chemo","CODS","survy")]

ALL=left_join(mydata,HM,by="casenum")
ALL$group=1
ALL_RT=ALL%>%filter(trt=="eb")
ALL_nonRT=ALL%>%filter(trt=="nr")
ALL_nonRT$survy[is.na(ALL_nonRT$survy)]=0.2
#write.csv(ALL_nonRT,file="E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/PSM/aaaaaaaa.csv")
#ALL_nonRT=as.data.frame(fread("E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/PSM/aaaaaaaa.csv"))

DenoveHM=canc%>%filter(cancer%in%c("ALL"),seqnum%in%c("0","1"),sex=="Female")
DenoveHM=DenoveHM[,c("casenum","trt","cancer","histo3","agedx","yrdx","race","chemo","CODS","survy")]
DenoveHM$group=0
DenoveALL_RT=bind_rows(ALL_RT,DenoveHM)
DenoveALL_nonRT=bind_rows(ALL_nonRT,DenoveHM)
DenoveALL_RT=na.omit(DenoveALL_RT)
DenoveALL_nonRT=na.omit(DenoveALL_nonRT)

DenoveALL_nonRT$Race[DenoveALL_nonRT$race=="White"]=1
DenoveALL_nonRT$Race[DenoveALL_nonRT$race=="Black"]=2
DenoveALL_nonRT$Race[DenoveALL_nonRT$race=="Other"]=3

DenoveALL_RT$Race[DenoveALL_RT$race=="White"]=1
DenoveALL_RT$Race[DenoveALL_RT$race=="Black"]=2
DenoveALL_RT$Race[DenoveALL_RT$race=="Other"]=3

DenoveALL_RT$Race=as.factor(DenoveALL_RT$Race)
DenoveALL_nonRT$Race=as.factor(DenoveALL_nonRT$Race)

DenoveALL_RT$chemo=as.factor(DenoveALL_RT$chemo)
DenoveALL_nonRT$chemo=as.factor(DenoveALL_nonRT$chemo)

DenoveALL_RT$histo3=as.factor(DenoveALL_RT$histo3)
DenoveALL_nonRT$histo3=as.factor(DenoveALL_nonRT$histo3)

DenoveALL_RT$RT="RT"
DenoveALL_nonRT$RT="nonRT"
f=matchit(group~yrdx+agedx+Race+chemo+histo3,data=DenoveALL_RT,method="nearest",discard="none",
          caliper=0.2,reestimate=F,
          ratio=5)
summary(f)
newdata <- match.data(f)
stable1 <- CreateTableOne(vars=c("yrdx","agedx","Race","chemo","histo3"),strata="group", data=newdata,factorVars=c("Race","chemo","histo3"))
print(stable1,showAllLevels = TRUE)

table(DenoveALL_nonRT$group)
DenoveALL_nonRT=DenoveALL_nonRT%>%filter(histo3%in%c("9835","9834","9820"))
f=matchit(group~yrdx+agedx+Race+chemo+histo3,data=DenoveALL_nonRT,method="nearest",discard="none",
          caliper=1.5,reestimate=F,
          ratio=5)
summary(f)
newdata2 <- match.data(f)
stable1 <- CreateTableOne(vars=c("yrdx","agedx","Race","chemo","histo3"),strata="group", data=newdata2,factorVars=c("Race","chemo","histo3"))
print(stable1,showAllLevels = TRUE)

ALL=bind_rows(newdata,newdata2)
write.csv(ALL,file="E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/PSM/ALL_denove.csv")
###############################################
df=as.data.frame(fread("E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/PSM/NHL_denove.csv"))
df$group1[df$group==1 & df$RT=="nonRT"]=1
df$group1[df$group==0 & df$RT=="nonRT"]=2
df$group1[df$group==1 & df$RT=="RT"]=3
df$group1[df$group==0 & df$RT=="RT"]=4

tapply(df$survy, df$group1, FUN=sum)
tapply(df$agedx, df$group1, FUN=quantile)

library(tableone)
vars=c("yrdx","agedx","Race","chemo","histo3")
catvars=c("Race","chemo","histo3")
tableOne <- CreateTableOne(vars=vars,strata = c("group1"),factorVars = catvars,includeNA=T,data = df)
tab_out<-print(tableOne)
write.csv(tab_out,file="E:/Ph.D projects/3.SEER/1.SEERaBomb/Analyses/Cohort/PSM/NHL_denove_table.csv")

############################
#############################



fit <- survfit(Surv(survy,status_survival) ~ status+trt, data = ALL)    #KM分析

ggsurvplot(fit,data=ALL, 
           title="NHL",
           ylab="Overall Survival", #""
           xlab="Follow-up time(month)",
           pval = TRUE, 
           pval.method=T,
           conf.int="F",
           risk.table = TRUE,
           legend.title = "Group",
           risk.table.y.text=F,
           xlim = c(0,120),
           break.time.by = 40,
           legend = "right",
           palette = "science")

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

