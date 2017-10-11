library(survival)
data(colon)
colon
fit <- survfit(Surv(time,status)~rx, data=colon)
ggkm(fit, timeby=500, ystratalabs=c("Obs","Lev","Lev+5FU"))

a = 0- diff(surv$Evolocumb) # diff to calculate consecutive diff of censored/dead from ppopulation at risk, subt from 0 to convert to positive
a

sum(a)

b1 = surv$Time[-1] # removing first element of vector time

b1

b1[length(b1)] # last elemt of b1
surv$Evolocumb[length(surv$Evolocumb)] # last element of Evolocumab

Group1 = c(rep(b1,a),rep(b1[length(b1)],surv$Evolocumb[length(surv$Evolocumb)]))

str(Group1)

Group1.time = c(rep(1,sum(a)),rep(0,surv$Evolocumb[length(surv$Evolocumb)])) # constructing censor/dead as zero , alive as1


c = 0-diff(surv$Placebo)

Group2 = c(rep(b1,c),rep(b1[length(b1)],surv$Placebo[length(surv$Placebo)]))

str(Group2)

Group2.time = c(rep(1,sum(c)),rep(0,surv$Placebo[length(surv$Placebo)])) # constructing censor/dead as zero , alive as1

GroupCode =c(rep("Evolocumab",surv$Evolocumb[1]),rep("Placebo",surv$Placebo[1]))

Weeks = c(Group1,Group2)
Censor =c(Group1.time,Group2.time)

Final_data = data.frame(GroupCode,Weeks,Censor)


sum(c)

library(survival)
library(KMsurv)

Survdata = survfit(Surv(Final_data$Weeks,Final_data$Censor)~1)

summary(Survdata)
plot(Survdata)

Survdata2 = survfit(Surv(Final_data$Weeks,Final_data$Censor)~Final_data$GroupCode)
summary(Survdata2)
plot(Survdata2)

survdiff(Surv(Final_data$Weeks,Final_data$Censor)~Final_data$GroupCode ,rho=0)