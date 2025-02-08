*Author: Junjie Wang, Jingyi Yin, Xiaoyue Zhang, Jianqiao Wang, Xing Xing, Jun Tu, Guoqi Cai
*Object: The restricted spline analysis of the association between sex hormone-binding globulin (SHBG) and the occurrence of knee osteoarthritis (KOA)
*The restricted spline analysis of the association between estradiol, testosterone and the occurrence of KOA is similar, so we have only upload the code for the analysis of SHBG and KOA
*Abbreviations
*PA: Physical activity; HRT: Hormone replacement therapy; OC: Oral contraceptive


*Import the data
import delimited "E:\KOA\SHBG\SHBG_longitudinal.csv",clear
*Variable transformation
replace tdi = "." if tdi =="NA"
replace ethnicity = "." if ethnicity=="NA"
replace bmi ="." if bmi=="NA"
replace smoking="." if smoking=="NA"
replace alcohol="." if alcohol=="NA"
replace pa="." if pa=="NA"
destring tdi ethnicity bmi smoking drinking pa, replace

*Restricted the exposure to the 2.5%-97.5% range
centile shbg, centile(2.5 97.5)
scalar lower = r(c_1)
scalar upper = r(c_2)
keep if shbg > lower & shbg < upper

*RCS
stset time, failure(koa_occurrence)
mkspline shbg_k=shbg, nknots(4) cubic displayknots
summarize shbg, detail
streg shbg_k* ib(0).hrt ib(0).oc ib(1).ethnicity ib(1).smoking ib(1).alcohol  ib(0).hypertension ib(0).diabetes ib(1).menopause age tdi bmi pa, distribution(loglogistic) time tratio allbaselevels
sort shbg
predictnl xb= _b[shbg_k1]*(shbg_k1-56.79) ///
+_b[shbg_k2]*(shbg_k2-4.325534) ///
+_b[shbg_k3]*(shbg_k3-.125627),ci(ll ul)
gen or=exp(xb)
gen or_ll=exp(ll)
gen or_ul=exp(ul)
sum shbg
twoway (rarea or_ll or_ul shbg, sort color(gs20)) (line or_ll or_ul or shbg,yline(1)) if inrange(shbg,r(min),r(max))

*View the p-value for non-linear relationship 
streg shbg ib(0).hrt ib(0).oc ib(1).ethnicity ib(1).smoking ib(1).alcohol  ib(0).hypertension ib(0).diabetes ib(1).menopause age tdi bmi pa, distribution(loglogistic) time tratio allbaselevels
estat ic
estimates store linear
streg shbg_k* ib(0).hrt ib(0).oc ib(1).ethnicity ib(1).smoking ib(1).alcohol  ib(0).hypertension ib(0).diabetes ib(1).menopause age tdi bmi pa, distribution(loglogistic) time tratio allbaselevels
estat ic
estimates store nonlinear
lrtest linear nonlinear

*save the image to PDF