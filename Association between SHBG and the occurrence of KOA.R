#Author: Junjie Wang, Jingyi Yin, Xiaoyue Zhang, Jianqiao Wang, Xing Xing, Jun Tu, Guoqi Cai
#Object: Association between sex hormone-binding globulin (SHBG) and the occurrence of knee osteoarthritis (KOA)
#The association between estradiol, testosterone and the occurrence of KOA is similar, so we have only upload the code for the analysis of SHBG and the occurrence of KOA

#Abbreviations
#PA: Physical activity; HRT: Hormone replacement therapy; OC: Oral contraceptive


library(haven)
library(dplyr)
library(survival)

#Set up working directory
setwd("E:/KOA/SHBG")

#Import the data
a<-read.csv("SHBG_longitudinal.csv")

#Change the variable type
fvars<-c("HRT","OC","ethnicity","smoking","alcohol","hypertension","diabetes","menopause")
a[fvars]<-lapply(a[fvars],as.factor)

#Standardize the exposure
sd<-sd(a$SHBG)
mean<-mean(a$SHBG)
a$SHBG_sd<-round(((a$SHBG-mean)/sd),2)

###Accelerated failure time (AFT) models
#The log-logistic distribution was selected for AFT models based on the minimum Akaike Information criterion among different survival distributions (e.g.: Weibull, logistic, log-logistic, log-normal, exponential, and Gaussian)

mod1<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd,data=a,dist="loglogistic")
summary(mod1)
exp(coef(mod1))
exp(confint(mod1))

mod2<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd+age+BMI,data=a,dist="loglogistic")
summary(mod2)
exp(coef(mod2))
exp(confint(mod2))

mod3<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,data=a,dist="loglogistic")
summary(mod3)
exp(coef(mod3))
exp(confint(mod3))

#The restricted cubic spline analysis of the association between exposure and the occurrence of KOA was conducted in Stata

#Check for interaction between exposure and menopause
fit<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd*menopause+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA,data=a,dist="loglogistic")
summary(fit)

#Covert the exposure into quartile
#View the quartiles
SHBG_m<-quantile(a$SHBG,c(0.25,0.50,0.75))
SHBG_m
a$SHBG_Q4<-ifelse(a$SHBG<40.71,1,
                  ifelse(a$SHBG<56.79,2,
                         ifelse(a$SHBG<76.63,3,4)))
a$SHBG_q4<-factor(a$SHBG_Q4,levels=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))
a$SHBG_q4<-relevel(a$SHBG_q4,"Q1")

#AFT models
mod1<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4,data=a,dist="loglogistic")
summary(mod1)
exp(coef(mod1))
exp(confint(mod1))

mod2<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4+age+BMI,data=a,dist="loglogistic")
summary(mod2)
exp(coef(mod2))
exp(confint(mod2))

mod3<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,data=a,dist="loglogistic")
summary(mod3)
exp(coef(mod3))
exp(confint(mod3))

#Check for interaction between exposure and menopause
fit<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4*menopause+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA,data=a,dist="loglogistic")
summary(fit)
