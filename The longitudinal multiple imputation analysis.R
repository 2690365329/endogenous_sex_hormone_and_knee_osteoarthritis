#Author: Junjie Wang, Jingyi Yin, Xiaoyue Zhang, Jianqiao Wang, Xing Xing, Jun Tu, Guoqi Cai
#Object: The analysis of multiple imputation of the association of sex hormone-binding globulin (SHBG) and the occurrence of knee osteoarthritis (KOA)
#The analysis of multiple imputation of the association between estradiol, testosterone and the occurrence of KOA is similar, so we have only upload the code for the analysis of SHBG and the occurrence of KOA
#Abbreviations
#PA: Physical activity; HRT: Hormone replacement therapy; OC: Oral contraceptive



library(dplyr)
library(haven)
library(survival)
library(mice)

#Set up working directory
setwd("E:/KOA/SHBG")

#Import the data
a<-read.csv("SHBG_longitudinal.csv")

#Change the variable type
fvars<-c("HRT","OC","ethnicity","smoking","alcohol","hypertension","diabetes","menopause","KOA_occurrence")
a[fvars]<-lapply(a[fvars],as.factor)

#Multiple imputations for covariates with missing data
imp_data<-mice(a,seed=2024,m=5,printFlag=T)

###Analyze the imputed dataset 1.
data1<-complete(imp_data,1)

#Standardize the exposure
sd<-sd(data1$SHBG)
mean<-mean(data1$SHBG)
data1$SHBG_sd<-(data1$SHBG-mean)/sd
data1$SHBG_sd<-round(data1$SHBG_sd,2)

#Convert exposure into quartiles
SHBG_m<-quantile(data1$SHBG,c(0.25,0.50,0.75))
SHBG_m
data1$SHBG_Q4<-ifelse(data1$SHBG<40.71,1,
                      ifelse(data1$SHBG<56.79,2,
                             ifelse(data1$SHBG<76.63,3,4)))
data1$SHBG_q4<-factor(data1$SHBG_Q4,levels=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))
data1$SHBG_q4<-relevel(data1$SHBG_q4,"Q1")

###Accelerated failure time models
#Expousre as continutious
#mod1
mod1_1<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd,data=data1,dist="loglogistic")
#mod2
mod1_2<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd+age+BMI,data=data1,dist="loglogistic")
#mod3
mod1_3<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,data=data1,dist="loglogistic")

#Expusore as quartiles
#mod1
fit1_1<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4,data=data1,dist="loglogistic")
#mod2
fit1_2<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4+age+BMI,data=data1,dist="loglogistic")
#mod3
fit1_3<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,data=data1,dist="loglogistic")

###Analyze the imputed dataset 2.
data2<-complete(imp_data,2)
sd<-sd(data2$SHBG)
mean<-mean(data2$SHBG)
data2$SHBG_sd<-(data2$SHBG-mean)/sd
data2$SHBG_sd<-round(data2$SHBG_sd,2)
SHBG_m<-quantile(data2$SHBG,c(0.25,0.50,0.75))
SHBG_m
data2$SHBG_Q4<-ifelse(data2$SHBG<40.71,1,
                      ifelse(data2$SHBG<56.79,2,
                             ifelse(data2$SHBG<76.63,3,4)))
data2$SHBG_q4<-factor(data2$SHBG_Q4,levels=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))
data2$SHBG_q4<-relevel(data2$SHBG_q4,"Q1")
mod2_1<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd,data=data2,dist="loglogistic")
mod2_2<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd+age+BMI,data=data2,dist="loglogistic")
mod2_3<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,data=data2,dist="loglogistic")
fit2_1<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4,data=data2,dist="loglogistic")
fit2_2<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4+age+BMI,data=data2,dist="loglogistic")
fit2_3<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,data=data2,dist="loglogistic")

###Analyze the imputed dataset 3.
data3<-complete(imp_data,3)
sd<-sd(data3$SHBG)
mean<-mean(data3$SHBG)
data3$SHBG_sd<-(data3$SHBG-mean)/sd
data3$SHBG_sd<-round(data3$SHBG_sd,2)
SHBG_m<-quantile(data3$SHBG,c(0.25,0.50,0.75))
SHBG_m
data3$SHBG_Q4<-ifelse(data3$SHBG<40.71,1,
                      ifelse(data3$SHBG<56.79,2,
                             ifelse(data3$SHBG<76.63,3,4)))
data3$SHBG_q4<-factor(data3$SHBG_Q4,levels=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))
data3$SHBG_q4<-relevel(data3$SHBG_q4,"Q1")
mod3_1<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd,data=data3,dist="loglogistic")
mod3_2<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd+age+BMI,data=data3,dist="loglogistic")
mod3_3<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,data=data3,dist="loglogistic")
fit3_1<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4,data=data3,dist="loglogistic")
fit3_2<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4+age+BMI,data=data3,dist="loglogistic")
fit3_3<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,data=data3,dist="loglogistic")

###Analyze the imputed dataset 4.
data4<-complete(imp_data,4)
sd<-sd(data4$SHBG)
mean<-mean(data4$SHBG)
data4$SHBG_sd<-(data4$SHBG-mean)/sd
data4$SHBG_sd<-round(data4$SHBG_sd,2)
SHBG_m<-quantile(data4$SHBG,c(0.25,0.50,0.75))
SHBG_m
data4$SHBG_Q4<-ifelse(data4$SHBG<40.71,1,
                      ifelse(data4$SHBG<56.79,2,
                             ifelse(data4$SHBG<76.63,3,4)))
data4$SHBG_q4<-factor(data4$SHBG_Q4,levels=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))
data4$SHBG_q4<-relevel(data4$SHBG_q4,"Q1")
mod4_1<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd,data=data4,dist="loglogistic")
mod4_2<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd+age+BMI,data=data4,dist="loglogistic")
mod4_3<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,data=data4,dist="loglogistic")
fit4_1<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4,data=data4,dist="loglogistic")
fit4_2<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4+age+BMI,data=data4,dist="loglogistic")
fit4_3<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,data=data4,dist="loglogistic")

###Analyze the imputed dataset 5.
data5<-complete(imp_data,5)
sd<-sd(data5$SHBG)
mean<-mean(data5$SHBG)
data5$SHBG_sd<-(data5$SHBG-mean)/sd
data5$SHBG_sd<-round(data5$SHBG_sd,2)
SHBG_m<-quantile(data5$SHBG,c(0.25,0.50,0.75))
SHBG_m
data5$SHBG_Q4<-ifelse(data5$SHBG<40.71,1,
                      ifelse(data5$SHBG<56.79,2,
                             ifelse(data5$SHBG<76.63,3,4)))
data5$SHBG_q4<-factor(data5$SHBG_Q4,levels=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))
data5$SHBG_q4<-relevel(data5$SHBG_q4,"Q1")
mod5_1<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd,data=data5,dist="loglogistic")
mod5_2<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd+age+BMI,data=data5,dist="loglogistic")
mod5_3<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_sd+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,data=data5,dist="loglogistic")
fit5_1<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4,data=data5,dist="loglogistic")
fit5_2<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4+age+BMI,data=data5,dist="loglogistic")
fit5_3<-survreg(Surv(time/365.25,KOA_occurrence)~SHBG_q4+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,data=data5,dist="loglogistic")

###Rubin rules
##Expousre as continutious
#The analysis of mod1
e<-summary(pool(list(mod1_1,mod2_1,mod3_1,mod4_1,mod5_1)),conf.int=T,conf.level=0.95)
e$or<-exp(e$estimate)
e$dl<-exp(e$estimate-1.96*e$std.error)      
e$ul<-exp(e$estimate+1.96*e$std.error)

#The analysis of mod2
e<-summary(pool(list(mod1_2,mod2_2,mod3_2,mod4_2,mod5_2)),conf.int=T,conf.level=0.95)
e$or<-exp(e$estimate)
e$dl<-exp(e$estimate-1.96*e$std.error)      
e$ul<-exp(e$estimate+1.96*e$std.error)

#The analysis of mod3
e<-summary(pool(list(mod1_3,mod2_3,mod3_3,mod4_3,mod5_3)),conf.int=T,conf.level=0.95)
e$or<-exp(e$estimate)
e$dl<-exp(e$estimate-1.96*e$std.error)      
e$ul<-exp(e$estimate+1.96*e$std.error)

##Expousre as quartiles
#The analysis of mod1
e<-summary(pool(list(fit1_1,fit2_1,fit3_1,fit4_1,fit5_1)),conf.int=T,conf.level=0.95)
e$or<-exp(e$estimate)
e$dl<-exp(e$estimate-1.96*e$std.error)      
e$ul<-exp(e$estimate+1.96*e$std.error)
#The analysis of mod2
e<-summary(pool(list(fit1_2,fit2_2,fit3_2,fit4_2,fit5_2)),conf.int=T,conf.level=0.95)
e$or<-exp(e$estimate)
e$dl<-exp(e$estimate-1.96*e$std.error)      
e$ul<-exp(e$estimate+1.96*e$std.error)

#The analysis of mod3
e<-summary(pool(list(fit1_3,fit2_3,fit3_3,fit4_3,fit5_3)),conf.int=T,conf.level=0.95)
e$or<-exp(e$estimate)
e$dl<-exp(e$estimate-1.96*e$std.error)      
e$ul<-exp(e$estimate+1.96*e$std.error)

###Convergence diagnosis
#Continuous covariates 
densityplot(imp_data,~BMI)
densityplot(imp_data,~TDI)
densityplot(imp_data,~PA)

#Categorical covariates
stripplot(imp_data,ethnicity)
stripplot(imp_data,smoking)
stripplot(imp_data,alcohol)
