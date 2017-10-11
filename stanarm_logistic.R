library(rstanarm)
library(tidyverse)
data("wells")
wells$dist100 <- wells$dist / 100

ggplot(wells, aes(x = dist100, y = ..density.., fill = switch == 1)) +
  geom_histogram() + 
  scale_fill_manual(values = c("gray30", "skyblue"))


t_prior <- student_t(df = 7, location = 0, scale = 2.5)
fit1 <- stan_glm(switch ~ dist100, data = wells, 
                 family = binomial(link = "logit"), 
                 prior = t_prior, prior_intercept = t_prior)  
                 

fit1

summary(fit1)

# Predicted probability as a function of x
pr_switch <- function(x, ests) plogis(ests[1] + ests[2] * x)
# A function to slightly jitter the binary data
jitt <- function(...) {
  geom_point(aes_string(...), position = position_jitter(height = 0.05, width = 0.1), 
             size = 2, shape = 21, stroke = 0.2)
}
ggplot(wells, aes(x = dist100, y = switch, color = switch)) + 
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  jitt(x="dist100") + 
  stat_function(fun = pr_switch, args = list(ests = coef(fit1)), 
                size = 2, color = "gray35")

# Additional predictor
fit2 <- update(fit1, formula = switch ~ dist100 + arsenic) 
fit2
summary(fit2)

# visualise as function of both variable
pr_switch2 <- function(x, y, ests) plogis(ests[1] + ests[2] * x + ests[3] * y)
grid <- expand.grid(dist100 = seq(0, 4, length.out = 100), 
                    arsenic = seq(0, 10, length.out = 100))
grid$prob <- with(grid, pr_switch2(dist100, arsenic, coef(fit2)))
ggplot(grid, aes(x = dist100, y = arsenic)) + 
  geom_tile(aes(fill = prob)) + 
  geom_point(data = wells, aes(color = factor(switch)), size = 2, alpha = 0.85) + 
  scale_fill_gradient() +
  scale_color_manual("switch", values = c("green", "red"), labels = c("No", "Yes"))


# We can see that the red points (switch=1) are predominantly clustered in the upper-left region of the plot where the predicted probability of switching is highest.

# compare_models

loo1 <- loo(fit1)
loo2 = loo(fit2)

# low loocic is better like AIC/BIC

loo1
loo2
compare_models(loo1, loo2)

# These results favor fit2 over fit1, as the estimated difference in elpd (the expected log pointwise predictive density for a new dataset) is so much larger than its standard error.
# LOO penalizes models for adding additional predictors (this helps counter overfitting), but in this case fit2 represents enough of an improvement over fit1 that the penalty 
#for including arsenic is negligible (as it should be if arsenic is an important predictor).

launch_shinystan(fit2)

pp_check(fit2, plotfun = "error_binned")  # ?bayesplot::ppc_error_binned


prior_summary(fit2)

# Programmatic estimate visualisation
bayesplot::color_scheme_set("green")
plot(fit2)
plot(fit2, regex_pars = c("dist100", "arsenic"))

# posterior predictive check of probability of output
pp_check(fit2)

y_rep= posterior_predict(fit2,draws=200)
y = wells$switch
library(bayesplot)

ppc_hist(y, y_rep[1:5, ])

# Test statistics
prop_zero <- function(x) mean(x == 0) # proportion of zeros


ppc_stat(y, y_rep, stat = "prop_zero", binwidth = 0.005)

# Accuracy statistics


# violin plot , intercept
pplot<-plot(fit2, "areas", prob = 0.95, prob_outer = 1)
pplot+ geom_vline(xintercept = 0)

# AUC

linpred <- posterior_linpred(fit2)
preds <- posterior_linpred(fit2, transform=TRUE)

pred <- colMeans(preds)
pr <- as.integer(pred >= 0.5)

library(caret)
# confusion matrix
confusionMatrix(pr, y)[2:3]
# posterior classification accuracy
round(mean(xor(pr,as.integer(y))),3)
# posterior balanced classification accuracy
round((mean(xor(pr[y==0]>0.5,as.integer(y[y==0])))+mean(xor(pr[y==1]>0.5,as.integer(y[y==1]))))/2,3)


# optimism correction with loo

# PSIS-LOO weights
library(loo)
log_lik=log_lik(fit2, parameter_name = "log_lik")
psis=psislw(-log_lik)
#plot(psis$pareto_k)
#plot(psis$lw_smooth[,1],linpred[,1])
# LOO predictive probabilities
ploo=colSums(preds*exp(psis$lw_smooth))
# LOO classification accuracy
round(mean(xor(ploo>0.5,as.integer(y))),3)
# LOO balanced classification accuracy
round((mean(xor(ploo[y==0]>0.5,as.integer(y[y==0])))+mean(xor(ploo[y==1]>0.5,as.integer(y[y==1]))))/2,2)
# calibration plot rstanarm
plot(pred,ploo)

# AUC

library(pROC)
plot.roc(y,pred,percent=TRUE,col="#1c61b6",  print.auc=TRUE)
plot.roc(y,ploo,percent=TRUE,col="#008600",  print.auc=TRUE, print.auc.y=40, add=TRUE)

legend("bottomright", legend=c("Posterior ROC", "LOO ROC"), col=c("#1c61b6", "#008600"), lwd=2)





# rms

library(rms)


a= datadist(wells)
options(datadist="a")

fit3 = lrm(switch~dist100+arsenic,data=wells,x=T,y=T)

fit3

fit4 = lrm(switch~dist100,data=wells,x=T,y=T)

fit4

validate(fit3,B=200)

plot(calibrate(fit3,B=200))



# library brms
library(brms)

fit5 = brm(formula=switch~dist100+arsenic,family = "bernoulli",data=wells)

summary(fit5)

# test model

prob=as.data.frame(predict(fit5,type="response"))
library(e1071)
pred1=ifelse(prob$Estimate> 0.5,1,0)

confusionMatrix(data=pred1,reference=wells$switch)


#fit 5 with priors # restricting odds ratio to between -0.7 (1/2) and 0.7 (twice)

Prior= get_prior(formula=switch~dist100+arsenic,family = "bernoulli",data=wells)

Prior <- set_prior("beta(1, 1)", class = "b", lb = 0, ub = 1)

?set_prior
priorb= c(set_prior("normal(0,0.35)", class = "b", coef = "arsenic"),
          set_prior("normal(0,0.35)", class = "b", coef = "dist100"))

fit6 = brm(formula=switch~dist100+arsenic,prior= priorb, family = "bernoulli",data=wells)

fit6
summary(fit6)

# plot _effects (plot logit vs loess) # very nice logit plots

plot(marginal_effects(fit6), ask = FALSE)

# very

prob=as.data.frame(predict(fit6,type="response"))
library(e1071)
pred1=ifelse(prob$Estimate> 0.5,1,0)

confusionMatrix(data=pred1,reference=wells$switch)

library(brms)
w5= WAIC(fit5)
w6 = WAIC(fit6)

library(loo)
compare(w5,w6)

compare_ic(w5,w6 , ic=c('loo','waic') )# c.f. compare_models for rstanarm
?compare_ic