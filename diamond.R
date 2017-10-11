# FISH trial

outcome = c(430,500)
nonevent = c(4570,4500)

Treats = c("Drug","Placebo")
Treats = as.factor(Treats)
Treats = relevel(Treats, ref = "Placebo")
#
library(rstanarm)

fish = data.frame(outcome,nonevent,Treats)

fitx = stan_glm(cbind(outcome,nonevent)~Treats ,data=fish, family = binomial(link = "logit"))

summary(fitx,digits = 3)

cidrug = posterior_interval(fitx,prob = 0.95)
ORDrug = exp(cidrug)[2,]
ORDrug

class(fitx)
#
fitxsample = as.data.frame(fitx)
fitxsample$OR = exp(fitxsample$TreatsDrug)
mean(fitxsample$OR)

Pref = outcome[2]/(nonevent[2]+outcome[2])

fitxsample$RR = fitxsample$OR/((1-Pref)+Pref*fitxsample$OR)
mean(fitxsample$RR)

fitxsample$RRR = 1-fitxsample$RR
mean(fitxsample$RRR)

# RRR>10
mean(fitxsample$RRR>0)
mean(fitxsample$RRR>0.10)
mean(fitxsample$RRR>0.20)


# Prior
uninformativeprior = normal(location = 0,scale=10)
#change prior, mildlyskeptical(0,0.4),moderatelyskeptical(0,0.07),highlyskeptical(0,0.03)

#enthusiastic
enthusiasticlocation = log(0.85) # 15% reduction leading to OR of 0.85
enthusiasticscale = abs(log(0.85)/1.96) # scale always positive




fitxuniform = update(fitx,prior = uninformativeprior)
summary(fitxuniform,digits=3)

fitxuniformsample = as.data.frame(fitxuniform)
fitxuniformsample$OR = exp(fitxuniformsample$TreatsDrug)
mean(fitxuniformsample$OR)



fitxuniformsample$RR = fitxuniformsample$OR/((1-Pref)+Pref*fitxuniformsample$OR)
mean(fitxuniformsample$RR)

fitxuniformsample$RRR = 1-fitxuniformsample$RR
mean(fitxuniformsample$RRR)

# RRR>10
mean(fitxuniformsample$RRR>0)
mean(fitxuniformsample$RRR>0.10)
mean(fitxuniformsample$RRR>0.20)



