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
