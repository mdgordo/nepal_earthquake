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

### Define GMM function

gmmmomentmatcher <- function(theta, df) {
  print(theta)
  saveRDS(theta, paste(getwd(), "/data/model_output/theta.rds", sep = ""))
  
  ### initial guess
  statespace = create.statespace(ubm = c(20,20), theta, method = "log")
  v0 = cbind(statespace, data.frame("Tw" = rep(0, nrow(statespace)),
                                    "cfx" = rep(0, nrow(statespace)),
                                    "ifx" = rep(0, nrow(statespace))))
  
  ### VFI
  V = VFI(v0, theta)
  saveRDS(V, paste(getwd(), "/data/model_output/V.rds", sep = ""))
  finalV <- V[[length(V)]]
  policycfx <- interpolater.creater(finalV, theta, var = "cfx", method = "neuralnet")
  policyifx <- interpolater.creater(finalV, theta, var = "ifx", method = "neuralnet")
  
  ### data
  momentmat <- mclapply(c(1:nrow(df)), momentmatcher, vfx = list(policycfx, policyifx), t0 = theta, 
                        data = df[,-c(1:2)], mc.cores = detectCores())
  momentmat <- do.call(rbind, momentmat)
  return(momentmat)
}

iterations <- data.frame("gamma" <- c(), "beta" <- c(), "R" <- c(), "cbar" <- c(),
                         "hbar" <- c(), "lambda" <- c(), "sigma" <- c(), "alpha" <- c(),
                         "delta" <- c(), "fval" = c())
saveRDS(iterations, paste(getwd(), "/data/model_output/iterations.rds", sep = ""))

### reasonable penalty value for R*B>1
momentmat <- lapply(c(1:nrow(df.hh)), momentmatcher, vfx = list(function(c,i){return(mean(df.hh$food_consumption))}, 
                                                                function(c,i){return(mean(df.hh$home_investment))}), 
                    t0 = rep(1,9), data = df.hh[,-c(1:2)])
momentmat <- do.call(rbind, momentmat)
gpenalty <- sum(colMeans(momentmat^2))

### wrapper functions
globalwrap <- function(theta){
  if (theta[2]*theta[3]>1) {
    return(gpenalty*theta[2]*theta[3]) }
  else { 
    momentmat = gmmmomentmatcher(theta, df.hh)
  }
  g = sum(colMeans(momentmat^2)); print(g)
  iterations = readRDS(paste(getwd(), "/data/model_output/iterations.rds", sep = ""))
  iterations = rbind(iterations, c(theta, g))
  saveRDS(iterations, paste(getwd(), "/data/model_output/iterations.rds", sep = ""))
  return(g)
}

gmmwrap <- function(theta, df){
  if (theta[2]*theta[3]>1) {
    return(matrix(gpenalty*theta[2]*theta[3], nrow = nrow(df), ncol = 11)) }
  else {
    momentmat = gmmmomentmatcher(theta, df)
  }
  g = sum(colMeans(momentmat^2)); print(g)
  iterations = readRDS(paste(getwd(), "/data/model_output/iterations.rds", sep = ""))
  iterations = rbind(iterations, c(theta, g))
  saveRDS(iterations, paste(getwd(), "/data/model_output/iterations.rds", sep = ""))
  return(momentmat)
}

## Global optimizer
g0 <- directL(fn = globalwrap, 
            lower = c(1.01, .6, 1, .01, .01, .01, .01, .1, .5), 
            upper = c(10, .99, 1.2, .9, .9, 10, 1.5, .99, 1),
            control = list(ftol_rel = 1e-2, xtol_rel = 1e-2))
g0

print("starting GMM")

g <- gmm(g = gmmwrap, x = df.hh, t0 = g0$par,         
         gradv = momentgradient, type = "twoStep", optfct = "nlminb", 
         lower = c(1.01, .6, 1, .01, .01, .01, .01, .1, .5),
         upper = c(10, .99, 1.2, .9, .9, 10, 1.5, .99, 1),
         control = list(x.tol = 1e-4, rel.tol = 1e-4, abs.tol = 1e-10))
saveRDS(g, paste(getwd(), "/data/model_output/theta.rds", sep = ""))
summary(g)


theta = g$par
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


