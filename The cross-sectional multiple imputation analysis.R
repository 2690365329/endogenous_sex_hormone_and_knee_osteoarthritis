#Author: Junjie Wang, Jingyi Yin, Xiaoyue Zhang, Jianqiao Wang, Xing Xing, Jun Tu, Guoqi Cai
#Object: The analysis of multiple imputation of the association of sex hormone-binding globulin (SHBG) and the prevalence of knee osteoarthritis (KOA)
#The analysis of multiple imputation of the association between estradiol, testosterone and the prevalence of KOA is similar, so we have only upload the code for the analysis of SHBG and the occurrence of KOA

#Abbreviations
#PA: Physical activity; HRT: Hormone replacement therapy; OC: Oral contraceptive

library(haven)
library(dplyr)
library(mice)

#Set up working directory
setwd("E:/KOA/SHBG")

#Import the data
a<-read.csv("SHBG_cross_sectional.csv")

#Change the variable type
fvars<-c("HRT","OC","ethnicity","smoking","alcohol","hypertension","diabetes","menopause","KOA_baseline")
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
SHBG_m #Q1: 40.55, Q2: 56.60, Q3: 76.42
data1$SHBG_Q4<-ifelse(data1$SHBG<40.55,1,
                      ifelse(data1$SHBG<56.60,2,
                             ifelse(data1$SHBG<76.42,3,4)))
data1$SHBG_q4<-factor(data1$SHBG_Q4,levels=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))
data1$SHBG_q4<-relevel(data1$SHBG_q4,"Q1")

###Logistic regression]
#Expousre as continutious
#mod1
mod1_1<-glm(KOA_baseline~SHBG_sd,family=binomial,data=data1)
#mod2
mod1_2<-glm(KOA_baseline~SHBG_sd+age+BMI,family=binomial,data=data1)
#mod3
mod1_3<-glm(KOA_baseline~SHBG_sd+age+BMI+OC+TDI+ethnicity+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,family=binomial,data=data1)

#Expusore as quartiles
fit1_1<-glm(KOA_baseline~SHBG_q4,family=binomial,data=data1)
fit1_2<-glm(KOA_baseline~SHBG_q4+age+BMI,family=binomial,data=data1)
fit1_3<-glm(KOA_baseline~SHBG_q4+age+BMI+OC+TDI+ethnicity+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,family=binomial,data=data1)

###Analyze the imputed dataset 2.
data2<-complete(imp_data,2)
sd<-sd(data2$SHBG)
mean<-mean(data2$SHBG)
data2$SHBG_sd<-(data2$SHBG-mean)/sd
data2$SHBG_sd<-round(data2$SHBG_sd,2)
SHBG_m<-quantile(data2$SHBG,c(0.25,0.50,0.75))
SHBG_m
data2$SHBG_Q4<-ifelse(data2$SHBG<40.55,1,
                      ifelse(data2$SHBG<56.60,2,
                             ifelse(data2$SHBG<76.42,3,4)))
data2$SHBG_q4<-factor(data2$SHBG_Q4,levels=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))
data2$SHBG_q4<-relevel(data2$SHBG_q4,"Q1")
mod2_1<-glm(KOA_baseline~SHBG_sd,family=binomial,data=data2)
mod2_2<-glm(KOA_baseline~SHBG_sd+age+BMI,family=binomial,data=data2)
mod2_3<-glm(KOA_baseline~SHBG_sd+age+BMI+OC+TDI+ethnicity+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,family=binomial,data=data2)
fit2_1<-glm(KOA_baseline~SHBG_q4,family=binomial,data=data2)
fit2_2<-glm(KOA_baseline~SHBG_q4+age+BMI,family=binomial,data=data2)
fit2_3<-glm(KOA_baseline~SHBG_q4+age+BMI+OC+TDI+ethnicity+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,family=binomial,data=data2)

###Analyze the imputed dataset 3.
data3<-complete(imp_data,3)
sd<-sd(data3$SHBG)
mean<-mean(data3$SHBG)
data3$SHBG_sd<-(data3$SHBG-mean)/sd
data3$SHBG_sd<-round(data3$SHBG_sd,2)
SHBG_m<-quantile(data3$SHBG,c(0.25,0.50,0.75))
SHBG_m
data3$SHBG_Q4<-ifelse(data3$SHBG<40.55,1,
                      ifelse(data3$SHBG<56.60,2,
                             ifelse(data3$SHBG<76.42,3,4)))
data3$SHBG_q4<-factor(data3$SHBG_Q4,levels=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))
data3$SHBG_q4<-relevel(data3$SHBG_q4,"Q1")
mod3_1<-glm(KOA_baseline~SHBG_sd,family=binomial,data=data3)
mod3_2<-glm(KOA_baseline~SHBG_sd+age+BMI,family=binomial,data=data3)
mod3_3<-glm(KOA_baseline~SHBG_sd+age+BMI+OC+TDI+ethnicity+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,family=binomial,data=data3)
fit3_1<-glm(KOA_baseline~SHBG_q4,family=binomial,data=data3)
fit3_2<-glm(KOA_baseline~SHBG_q4+age+BMI,family=binomial,data=data3)
fit3_3<-glm(KOA_baseline~SHBG_q4+age+BMI+OC+TDI+ethnicity+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,family=binomial,data=data3)

###Analyze the imputed dataset 4.
data4<-complete(imp_data,4)
sd<-sd(data4$SHBG)
mean<-mean(data4$SHBG)
data4$SHBG_sd<-(data4$SHBG-mean)/sd
data4$SHBG_sd<-round(data4$SHBG_sd,2)
SHBG_m<-quantile(data4$SHBG,c(0.25,0.50,0.75))
SHBG_m
data4$SHBG_Q4<-ifelse(data4$SHBG<40.55,1,
                      ifelse(data4$SHBG<56.60,2,
                             ifelse(data4$SHBG<76.42,3,4)))
data4$SHBG_q4<-factor(data4$SHBG_Q4,levels=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))
data4$SHBG_q4<-relevel(data4$SHBG_q4,"Q1")
mod4_1<-glm(KOA_baseline~SHBG_sd,family=binomial,data=data4)
mod4_2<-glm(KOA_baseline~SHBG_sd+age+BMI,family=binomial,data=data4)
mod4_3<-glm(KOA_baseline~SHBG_sd+age+BMI+OC+TDI+ethnicity+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,family=binomial,data=data4)
fit4_1<-glm(KOA_baseline~SHBG_q4,family=binomial,data=data4)
fit4_2<-glm(KOA_baseline~SHBG_q4+age+BMI,family=binomial,data=data4)
fit4_3<-glm(KOA_baseline~SHBG_q4+age+BMI+OC+TDI+ethnicity+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,family=binomial,data=data4)

###Analyze the imputed dataset 5.
data5<-complete(imp_data,5)
sd<-sd(data5$SHBG)
mean<-mean(data5$SHBG)
data5$SHBG_sd<-(data5$SHBG-mean)/sd
data5$SHBG_sd<-round(data5$SHBG_sd,2)
SHBG_m<-quantile(data5$SHBG,c(0.25,0.50,0.75))
SHBG_m
data5$SHBG_Q4<-ifelse(data5$SHBG<40.55,1,
                      ifelse(data5$SHBG<56.60,2,
                             ifelse(data5$SHBG<76.42,3,4)))
data5$SHBG_q4<-factor(data5$SHBG_Q4,levels=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))
data5$SHBG_q4<-relevel(data5$SHBG_q4,"Q1")
mod5_1<-glm(KOA_baseline~SHBG_sd,family=binomial,data=data5)
mod5_2<-glm(KOA_baseline~SHBG_sd+age+BMI,family=binomial,data=data5)
mod5_3<-glm(KOA_baseline~SHBG_sd+age+BMI+OC+TDI+ethnicity+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,family=binomial,data=data5)
fit5_1<-glm(KOA_baseline~SHBG_q4,family=binomial,data=data5)
fit5_2<-glm(KOA_baseline~SHBG_q4+age+BMI,family=binomial,data=data5)
fit5_3<-glm(KOA_baseline~SHBG_q4+age+BMI+OC+TDI+ethnicity+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,family=binomial,data=data5)

###Rubin rules
##Expousre as continutious
#The analysis of mod1
summary(pool(list(mod1_1,mod2_1,mod3_1,mod4_1,mod5_1)),conf.int=T,conf.level=0.95)
dl<-exp(-0.3179847-1.96*0.01741846)
dl
ul<-exp(-0.3179847+1.96*0.01741846)
ul
or<-exp(-0.3179847)
or

#The analysis of mod2
summary(pool(list(mod1_2,mod2_2,mod3_2,mod4_2,mod5_2)),conf.int=T,conf.level=0.95)
or<-exp(0.008913958)
or
dl<-exp(0.008913958-1.96*0.020275925)
dl
ul<-exp(0.008913958+1.96*0.020275925)
ul

#The analysis of mod3
summary(pool(list(mod1_3,mod2_3,mod3_3,mod4_2,mod5_3)),conf.int=T,conf.level=0.95)
or<-exp(1.997639e-02)
or
ul<-exp(1.997639e-02+1.96*2.133522e-02)
ul
dl<-exp(1.997639e-02-1.96*2.133522e-02)
dl

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
