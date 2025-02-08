#Author: Junjie Wang, Jingyi Yin, Xiaoyue Zhang, Jianqiao Wang, Xing Xing, Jun Tu, Guoqi Cai
#Object: The analysis of the assocaition between sex hormone-binding globulin (SHBG) and the prevalence of knee osteoarthritis (KOA)
#The analysis of the association between estradiol, testosterone and the prevalence of KOA is similar, so we have only upload the code for the analysis of SHBG and the prevalence of KOA

#Abbreviations
#PA: Physical activity; HRT: Hormone replacement therapy; OC: Oral contraceptive

library(dplyr)
library(haven)

#Set up working directory
setwd("E:/KOA/SHBG")

#Import the data
a<-read.csv("SHBG_cross_sectional.csv")

#Change the variable type
fvars<-c("HRT","OC","ethnicity","smoking","alcohol","hypertension","diabetes","menopause")
a[fvars]<-lapply(a[fvars],as.factor)

#Standardize the exposure
sd<-sd(a$SHBG)
mean<-mean(a$SHBG)
a$SHBG_sd<-round(((a$SHBG-mean)/sd),2)

#Logistic regression models
mod1<-glm(KOA_baseline~SHBG_sd,family=binomial,data=a)
summary(mod1)
exp(coef(mod1))
exp(confint(mod1))

mod2<-glm(KOA_baseline~SHBG_sd+age+BMI,family=binomial,data=a)
summary(mod2)
exp(coef(mod2))
exp(confint(mod2))

mod3<-glm(KOA_baseline~SHBG_sd+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,family=binomial,data=a)
summary(mod3)
exp(coef(mod3))
exp(confint(mod3))

#Check for interaction between exposure and menopause
fit<-glm(KOA_baseline~SHBG_sd*menopause+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA,family=binomial,data=a)
summary(fit)


#Restricted cubic splines
library(rms)
library(ggplot2)

#Restricted the exposure to the 2.5%-97.5% range
value<-quantile(a$SHBG,c(0.025,0.975))
value
a1<-a %>%
  dplyr::filter(SHBG>20.18 & SHBG<128.76)
#RCS
dd<-rms::datadist(a1)
options(datadist = "dd")
res.logitics.rcs<-rms::Glm(KOA_baseline~rcs(SHBG,4)+HRT+age+TDI+ethnicity+BMI+smoking+alcohol+OC+hypertension+diabetes+PA+menopause,data=a1,x=TRUE,y=TRUE) 
pre.logistics<-rms::Predict(res.logitics.rcs,SHBG, fun=exp, type = "predictions", ref.zero=TRUE)
pre.logistics
res.logitics.rcs
anova(res.logitics.rcs)
P2<-ggplot()+
  geom_line(data = pre.logistics, aes(SHBG, yhat), linetype = "solid", linewidth = 1, alpha = 1, colour = "red") +
  geom_ribbon(data = pre.logistics, aes(SHBG, ymin = lower, ymax = upper), alpha = 0.2, fill = "orangered") +
  geom_hline(yintercept = 1, linetype = 2, color = "grey")

#View the quartiles
s<-quantile(a$SHBG,c(0.25,0.50,0.75))
s #Q1: 41.38; Q2: 56.60; Q3: 75.16

P2<-P2+theme_classic()+geom_hline(yintercept = 1,linetype=2,size=1)+
  geom_vline(xintercept =c(41.38,56.60,75.16),linetype=2,size=1)+labs(x="SHBG",y="Odds ratio (95% confidence interval)",title="")+geom_text(aes(x=20.19, y=1.01, label=paste0(
    "P for nonlinear = 0.045")), hjust=0,color = "black",size=6)
P2<-P2+theme(axis.text = element_text(color = "black", size = 20),text=element_text(color = "black", size = 20))
P2
#Save P2 as a PDF file

#Covert the exposure into quartile
#View the quartiles
SHBG_m<-quantile(a$SHBG,c(0.25,0.50,0.75))
SHBG_m
a$SHBG_Q4<-ifelse(a$SHBG<40.55,1,
                  ifelse(a$SHBG<56.60,2,
                         ifelse(a$SHBG<76.42,3,4)))
a$SHBG_q4<-factor(a$SHBG_Q4,levels=c(1,2,3,4),labels=c("Q1","Q2","Q3","Q4"))
a$SHBG_q4<-relevel(a$SHBG_q4,"Q1")

#Logistic regressions
mod1<-glm(KOA_baseline~SHBG_q4,family=binomial,data=a)
summary(mod1)
exp(coef(mod1))
exp(confint(mod1))

mod2<-glm(KOA_baseline~SHBG_q4+age+BMI,family=binomial,data=a)
summary(mod2)
exp(coef(mod2))
exp(confint(mod2))

mod3<-glm(KOA_baseline~SHBG_q4+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA+menopause,family=binomial,data=a)
summary(mod3)
exp(coef(mod3))
exp(confint(mod3))

#Check for interaction between exposure and menopause
fit<-glm(KOA_baseline~SHBG_q4*menopause+OC+age+TDI+ethnicity+BMI+smoking+alcohol+HRT+hypertension+diabetes+PA,family=binomial,data=a)
summary(fit)
