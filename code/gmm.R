### packages and data
library(tidyverse)
library(parallel)
library(nloptr)
library(fastGHQuad)
library(flatlandr)
library(gmm)
library(Rearrangement)
library(mgcv)
library(neuralnet)
options('nloptr.show.inequality.warning'=FALSE)
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))

iterations <- readRDS(paste(getwd(), "/data/model_output/iterations.rds", sep = "")) %>%
  filter(complete.cases(.))

### fit smoothed function to iteration data and minimize it
inormed = as.data.frame(lapply(iterations, function(x) (x-min(x))/(max(x)-min(x))))
colmins = unlist(lapply(iterations, min))
colrange = unlist(lapply(iterations, function(x) max(x) - min(x)))

modfunc <- function(m){
  form <- as.formula(paste(m, "gamma + beta + R + cbar + hbar + 
                    lambda + sigma + alpha + delta", sep = "~"))
  mod = neuralnet(form,
                  data = inormed, hidden=11, act.fct = "logistic", 
                  linear.output = FALSE, stepmax = 5e5, rep = 5)
  mod = lm(get(m) ~ .^2 + I(.^2))
  return(mod)
}

momentvars <- paste("m", seq(1,11,1), sep = "_")
modlist <- mclapply(momentvars, modfunc, mc.cores = max(11,detectCores()-2))

smoothedobj <- function(t){
  if ((t[2]*colrange[2] + colmins[2])*(t[3]*colrange[3] + colmins[3])>1){
    r = rep(1,11)} else { 
      newdata = data.frame("gamma" = t[1], "beta" = t[2], "R" = t[3],
                           "cbar" = t[4], "hbar" = t[5], "lambda" = t[6],
                           "sigma" = t[7], "alpha" = t[8], "delta" = t[9])
      momlist = lapply(modlist, function(mod) predict(mod, newdata))
      r = unlist(momlist)*colrange[c(10:20)] + colmins[c(10:20)]
    }
  return(r)
}

sol <- lbfgs(x0 = rep(0,9), fn = function(t) sum(smoothedobj(t)^2), 
             lower = rep(0, 9), upper = rep(1,9), control = list(maxeval = 3000))
theta = sol$par*colrange[c(1:9)] + colmins[c(1:9)]
print(theta)

### load data
df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500) 

df.hh <- df.hh %>%
  filter(designation !="none") %>%
  select(hhid, wave, wt_hh, food_consumption, total_income, var_inc, avg_inc, 
         home_value, home_investment, imputed_bufferstock, quake_aid, M_avg) %>% 
  mutate(home_value = home_value/M_avg,
         imputed_bufferstock = imputed_bufferstock/M_avg,
         food_consumption = food_consumption/M_avg,
         home_investment = home_investment/M_avg,
         quake_aid = quake_aid/M_avg,
         total_income = total_income/M_avg) %>%
  group_by(hhid) %>%
  mutate(lag_h = lag(home_value, order_by = wave),
         lag_x = lag(imputed_bufferstock, order_by = wave),
         lag_y = lag(total_income, order_by = wave),
         lag_c = lag(food_consumption, order_by = wave),
         lag_i = lag(home_investment, order_by = wave)) %>%
  filter(avg_inc>0 & food_consumption>0 & lag_h + home_investment>0)
df.hh <- df.hh[complete.cases(df.hh),]

### Grid size
xn = 40; hn = 40

### GMM using theta as starting point
### reasonable penalty value for R*B>1
momentmat <- lapply(c(1:nrow(df.hh)), momentmatcher, vfx = list(function(c,i){return(mean(df.hh$food_consumption))}, 
                                                                function(c,i){return(mean(df.hh$home_investment))}), 
                    t0 = rep(1,9), data = df.hh[,-c(1:2)])
momentmat <- do.call(rbind, momentmat)
gpenalty <- sum(colSums(momentmat)^2)

gmmwrap <- function(t, x){
  if (t[2]*t[3]>1){
    return(matrix(gpenalty*theta[2]*theta[3], nrow = nrow(x), ncol = 11))
  } else{
    momentmat = gmmmomentmatcher(theta = t, df = x)
    print(sum(colSums(momentmat)^2))
    return(momentmat)
  }
}

smoothderivs <- function(theta, x, eps = .01){
  tnorm <- (theta - colmins[c(1:9)])/colrange[c(1:9)]
  m0 = smoothedobj(tnorm)
  tp = c(tnorm[1] + eps, tnorm[c(2:9)])
  m1 = (smoothedobj(tp) - m0)/eps
  tp = c(tnorm[1], tnorm[2] + eps, tnorm[c(3:9)])
  m2 = (smoothedobj(tp) - m0)/eps
  tp = c(tnorm[1:2], tnorm[3] + eps, tnorm[c(4:9)])
  m3 = (smoothedobj(tp) - m0)/eps
  tp = c(tnorm[1:3], tnorm[4] + eps, tnorm[c(5:9)])
  m4 = (smoothedobj(tp) - m0)/eps
  tp = c(tnorm[1:4], tnorm[5] + eps, tnorm[c(6:9)])
  m5 = (smoothedobj(tp) - m0)/eps
  tp = c(tnorm[1:5], tnorm[6] + eps, tnorm[c(7:9)])
  m6 = (smoothedobj(tp) - m0)/eps
  tp = c(tnorm[1:6], tnorm[7] + eps, tnorm[c(8:9)])
  m7 = (smoothedobj(tp) - m0)/eps
  tp = c(tnorm[1:7], tnorm[8] + eps, tnorm[9])
  m8 = (smoothedobj(tp) - m0)/eps
  tp = c(tnorm[1:8], tnorm[9] + eps)
  m9 = (smoothedobj(tp) - m0)/eps
  m = cbind(m1,m2,m3,m4,m5,m6,m7,m8,m9)
  return(m)
}

g <- gmm(g = gmmwrap, x = df.hh, t0 = theta, gradv = smoothderivs,       
         type = "twoStep", optfct = "nlminb", 
         lower = c(1.01, .6, 1, .01, .01, .01, .01, .1, .5),
         upper = c(10, .99, 1.4, .9, .9, 10, 1, .99, 1),
         control = list(x.tol = 1e-2, rel.tol = 1e-2, abs.tol = 1e-6))
summary(g)

theta = g$coefficients
gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
saveRDS(theta, paste(getwd(), "/data/model_output/theta.rds", sep = ""))

### statespace
xn = 40; hn = 40

### Final Value function
statespace = create.statespace(ubm = c(20,20), theta, method = "log")
v0 = cbind(statespace, data.frame("Tw" = rep(0, nrow(statespace)),
                                  "cfx" = rep(0, nrow(statespace)),
                                  "ifx" = rep(0, nrow(statespace))))
V = VFI(v0, theta)
saveRDS(V, paste(getwd(), "/data/model_output/V.rds", sep = ""))
