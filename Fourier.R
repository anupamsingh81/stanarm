# Fourier Trial analysis.


PrimaryEndPoint = c(1344,1563)
total = c(13780,13784)

SecondaryEndPoint = c(816,1013)



prop.test(PrimaryEndPoint,total)

library(BayesianFirstAid)
fit = bayes.prop.test(PrimaryEndPoint,total)

plot(fit)

summary(fit)


Fourier = as.data.frame(fit)

Fourier$RelativeRisk = Fourier$theta1/Fourier$theta2

Fourier$OR = (1-Fourier$theta2)/(1-Fourier$theta1) * Fourier$theta1 / Fourier$theta2

Fourier$difference =Fourier$theta2 - Fourier$theta1

Fourier$NNT = 1/Fourier$difference

mean(Fourier$RelativeRisk)
mean(Fourier$OR)
quantile(Fourier$OR,c(0.025,0.975))

quantile(Fourier$RelativeRisk,c(0.025,0.975))

# CompVal , what are chances Relative Risk is less than 0.85
mean(Fourier$RelativeRisk< 0.85)


mean(Fourier$difference)

quantile(Fourier$difference,c(0.025,0.975))

mean(Fourier$NNT)

quantile(Fourier$NNT, c(0.025,0.975))






# Stan

# Frequentist glm


Events  = PrimaryEndPoint
NoEvents = total - PrimaryEndPoint

Treatment = c("Evolocumab","Placebo")

str(Treatment)

Treatment = as.factor(Treatment)
levels(Treatment)

Treatment = relevel(Treatment, ref = "Placebo") # relevel so that reference level is placebo ,

fit2 = glm(cbind(Events,NoEvents)~Treatment , family = binomial())

summary(fit2)


library(arm)
display(fit2)







library(brms)

df= data.frame(Events,total,Treatment)

fit3 = brm(Events | trials(total)~Treatment,data = df ,family = binomial())

summary(fit3)
# samples

samples = as.data.frame(fit3$fit)

mean(samples$ORDrug)
mean(Fourier$OR)

library(rstanarm)

df2 = data.frame(Events,NoEvents,Treatment)

fit4 = stan_glm(cbind(Events,NoEvents)~Treatment ,data=df2, family = binomial(link = "logit"))

summary(fit4)
coef(fit4)

ci95 <- posterior_interval(fit4, prob = 0.95)
ci95[2,] # second row of ci95 containing OR

ORCI = exp(ci95[2,])
ORCI


# prior calculation
# kaul function ,spieghalter
# One way to select the mean (m) and standard deviation (s)of the prior log odds ratio is from previously published data.  
#If such data are not available, one can enumerate the limits of one's belief regarding all possible values of the odds
#ratio, convert these values to logarithms, and compute the mean and standard deviation for a 95% confidence interval
#of the associated Gaussian distribution. Thus, if you thinkthe odds ratio must lie between 0.8 (a 20% relative risk
#reduction) and 1.1 (a 10% relative risk increase), thenm = [ ln(0.8) + ln(1.1) ] / 2 = - 0.064
#ands = [ | ln(0.8) | + | ln(1.1) | ] / 3.92 = 0.081

kaul = function(a,b){
  m = (log(a) + log(b))/2 
  sd = (abs(log(a)) + abs(log(b)))/3.92
  neff = 4/(sd*sd)
  param = c(m,sd,neff)
  param
}

kaul(0.8,1.1) # gives mean,location and neffective for normal prior on log oods ratio
exp
exp(0.8)
log(0.8)


# prior location(mean 0, sd = 2 on intercept weak prior gelman)
fit5 = stan_glm(cbind(Events,NoEvents)~Treatment ,data=df2, prior = cauchy(location=0 , scale=2) ,family = binomial(link = "logit"))
coef(fit5)

# So, one could justify constructing a prior based on the investigator's expectation of a 25% benefit by centering the distribution on a RR of 0.75 (equivalent to a relative risk reduction of 25%) and a very small likelihood (say 2.5% probability) of RR>1.0 and RR <0.75. Thus the 95% CI of this prior distribution (characterized as ENTHUSIASTIC) would be 0.56 to 1.0 RR or a mean ln (RR) of -0.288 and sd ln (RR) of 0.147. 
fit6 = stan_glm(cbind(Events,NoEvents)~Treatment ,data=df2, prior = normal(location= -0.288 , scale=0.147) ,family = binomial(link = "logit"))
coef(fit6 , digits =3 )
summary(fit6)

# The other choice of prior would be based on a mean RR of 1.0 (c/w null effect) and a very small probability of benefit <0.75 RR. Thus the 95% CI of this prior distribution (characterized as SKEPTICAL) would be 0.75 to 1.33 RR or a mean ln (RR) of 0 and sd ln (RR) of 0.147. 

fit7 = stan_glm(cbind(Events,NoEvents)~Treatment ,data=df2, prior = normal(location= 0 , scale=0.147) ,family = binomial(link = "logit"))
coef(fit7)
summary(fit7)


#For the uniform prior, I use a log odds mean of zero and sd of 2. What values do you recommend for a "weakly informative" prior, especially when there is no previous information available?

fit8 = stan_glm(cbind(Events,NoEvents)~Treatment ,data=df2, prior = normal(location= 0 , scale= 2) ,family = binomial(link = "logit"))
coef(fit8)
summary(fit8)

# fit 9 , expectation, very enthusiastic expectation , reduction maximum 60% , minimum 0%
kaul(0.6,1)
fit9 = stan_glm(cbind(Events,NoEvents)~Treatment ,data=df2, prior = normal(location= -0.25 , scale= 0.13) ,family = binomial(link = "logit"))
coef(fit9)


  
credint <- posterior_interval(fit9, prob = 0.95)
ci95[2,] # second row of ci95 containing OR

ORCI = exp(credint[2,])
ORCI

# informed priors

# http://andrewgelman.com/2008/01/24/specifying_a_pr/
# http://andrewgelman.com/2008/02/05/specifying_a_pr_1/

 # odds ratio of 

# 


samples2 = as.data.frame(fit4$stanfit)

samples2$ORDrug = 1-samples2$TreatmentPlacebo

mean(samples2$ORDrug)



# BayesFactor


df3 = data.frame(rbind(Events,NoEvents))

library(BayesFactor)

df3 = as.matrix(df3)
bf = contingencyTableBF(df3, sampleType = "indepMulti", fixedMargin = "cols")
bf

# bayes factor of proportions

# https://mvuorre.github.io/post/2017/bayes-factors-with-brms/

df2
df4 = data.frame(Events,total,Treatment)
library(brms)
get_prior(Events | trials(total) ~ 0 + Treatment, 
          family=binomial(link="identity"),
          data = df4)
Prior <- set_prior("beta(1, 1)", class = "b", lb = 0, ub = 1)

m1 <- brm(Events | trials(total) ~ 0 + Treatment, 
          family=binomial(link="identity"),
          prior = Prior,
          sample_prior = TRUE,
          iter = 1e4,
          data = df4,
          cores = 4)

print(m1,digits=3)

h1 <- hypothesis(m1, "TreatmentPlacebo = TreatmentEvolocumab")
print(h1, digits = 3)

BF10 = 1/h1$hypothesis$Evid.Ratio
BF10





kaul(0.73,0.85)
kaul(0.48,1.19)
