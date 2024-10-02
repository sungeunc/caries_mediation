library(survival)
library(survminer)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
require(tidyr)
library(naniar)
require(stringr)
require(lubridate)
library(gdata)
library(survival)
library(matrixStats)
library(Hmisc)
library(glmnet) 
library(cvAUC)
library(doMC)
library(dplyr)
library(tableone)
library(mma)

adm.cens=365*5
d<-data
d<-d[!d$tstop>adm.cens,] 

# age grouping (ages 0-5, 6-10, 11-18) # repeat analysis for each age group
prime<-d[d$agecat3==1,]
child<-d[d$agecat3==2,]
perm<-d[d$agecat3==3,]

## survival analysis 
# survival curve
fit.km.prime  <- survfit(Surv(tstart, tstop, newdecay.t) ~ race_re, id=studyID,
                         conf.type="log-log",
                         data = prime)


plot<- ggsurvplot(
  fit.km.prime, 
  fun = "event", # plot cumulative incidence
  conf.int = T, # include confidence intervals
  conf.int.style = "ribbon",
  censor = FALSE, #
  xlab = "Days since the initial visit", 
  legend = c(0.1, 0.8), 
  legend.labs = c("Black","Hispanic","White","Other"), 
  surv.scale = "percent", # show y-axis in %
  ylab = "Cumulative Incidence (%)", # label y-axis
  ylim = c(0,1), break.y.by = .2,# set limits of y-axis
  xlim = c(0,adm.cens), break.x.by = 200,
  legend.title = "",
  font.legend = c(11, "bold"),
  risk.table = TRUE,
  #title="Age 0-5 Years at Baseline",
  ggtheme=custom_theme(),
  risk.table=TRUE,
  palette =c("#440154","#5ec962","#2a788e","#fde725"))

# survival model for recurrent events 
survc.pri<-coxph(Surv(tstart, tstop, newdecay.t) ~age_visit + age_sq +female+ (race3)+ insure_vist+
                   m_musc+m_neuro+m_growth+smoke+
                   proc_clean_yr+proc_fluor_yr+proc_seal_yr+proc_restor_yr+proc_extract_yr+
                   pct_lessHS_i+ADI_nat,method="breslow",robust=TRUE,
                 data =prime, cluster = studyID)

## Mediation analysis 

## Black vs White 
test<-prime # repeat this for child and perm group and by race 
test<-test[test$race3=="White"|test$race3=="Black",]
test$race3<-droplevels(test$race3)
co.list<-c( "race3","female","age_visit","insure_vist", "m_growth","m_musc","m_neuro","smoke","proc_clean_yr","proc_fluor_yr","proc_seal_yr",
            "proc_restor_yr","proc_extract_yr","ADI_nat","pct_lessHS_i","newtime","newdecay.t","tstart","tstop","county")
co.var <- match(co.list, names(test))

co <- test[,co.var]
co=co[complete.cases(co),]

y=Surv(co$tstart,co$tstop, co$newdecay.t)

x=co[,2:15]
pred =co[,1]


data.surv.contx<-data.org(x,y,pred=pred,mediator=c(3,7:ncol(x)),
                            alpha=0.1,alpha2=0.1)
med.black.pri<-boot.med(data=data.surv.contx,n=10, n2=10,nonlinear=TRUE)  

## Hispanic vs White 
test<-prime
test<-test[test$race3=="White"|test$race2=="Hispanic",]
test$race3<-droplevels(test$race3)
co.list<-c( "race3","female","age_visit","insure_vist", "m_growth","m_musc","m_neuro","smoke","proc_clean_yr","proc_fluor_yr","proc_seal_yr",
            "proc_restor_yr","proc_extract_yr","ADI_nat","pct_lessHS_i","newtime","newdecay.t","tstart","tstop","county")
co.var <- match(co.list, names(test))

co <- test[,co.var]
co=co[complete.cases(co),]

y=Surv(co$tstart,co$tstop, co$newdecay.t)

x=co[,2:15]
pred =co[,1]

data.surv.contx<-data.org(x,y,pred=pred,mediator=c(3,7:ncol(x)),
                          alpha=0.1,alpha2=0.1)
med.latino.prime<-boot.med(data=data.surv.contx,n=10, n2=10,nonlinear=TRUE)  

## Other vs White 
test<-prime
test<-test[test$race3=="White"|test$race3=="Other/NATIVE/ISLAND",]
test$race3<-droplevels(test$race3)
co.list<-c( "race3","female","age_visit","insure_vist", "m_growth","m_musc","m_neuro","smoke","proc_clean_yr","proc_fluor_yr","proc_seal_yr",
            "proc_restor_yr","proc_extract_yr","ADI_nat","pct_lessHS_i","newtime","newdecay.t","tstart","tstop","county")
co.var <- match(co.list, names(test))

co <- test[,co.var]
co=co[complete.cases(co),]

y=Surv(co$tstart,co$tstop, co$newdecay.t)

x=co[,2:15]
pred =co[,1]

data.surv.contx<-data.org(x,y,pred=pred,mediator=c(3,7:ncol(x)),
                          alpha=0.1,alpha2=0.1)

med.other.pri<-boot.med(data=data.surv.contx,n=10, n2=10,nonlinear=TRUE)  
