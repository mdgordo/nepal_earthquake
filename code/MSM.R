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

### start guess for GMM
#gamma0 <- 8.75; beta0 <- .84; R0 <- 1.02; cbar0 <- .44; hbar0 <- .35
#lambda0 <- 2.46; sigma0 <- .77; alpha0 <- .79; delta0 <- .95
#theta0 <- c(gamma0, beta0, R0, cbar0, hbar0, lambda0, sigma0, alpha0, delta0)

### Grid size
xn = 40; hn = 40

### Define GMM function

gmmmomentmatcher <- function(theta) {
  print(theta) ### get rid of df arg for global
  if (theta[2]*theta[3]>1) return(1*theta[2]*theta[3]) else {  ### use for global
  #if (theta[2]*theta[3]>1) return(matrix(1*theta[2]*theta[3], nrow = nrow(df), ncol = 11)) else {  ### use for GMM
    saveRDS(theta, paste(getwd(), "/data/model_output/theta.rds", sep = ""))
    
    ### initial guess
    statespace = create.statespace(ubm = c(20,20), theta, method = "equal")
    v0 = cbind(statespace, data.frame("Tw" = rep(0, nrow(statespace)),
                                      "cfx" = rep(0, nrow(statespace)),
                                      "ifx" = rep(0, nrow(statespace))))
    
    ### VFI
    V = VFI(v0, theta)
    saveRDS(V, paste(getwd(), "/data/model_output/V.rds", sep = ""))
    finalV <- V[[length(V)]]
    policycfx <- interpolater.creater(finalV, theta, var = "cfx", method = "neuralnet")
    policyifx <- interpolater.creater(finalV, theta, var = "ifx", method = "neuralnet")
    
    ### data - change to df.hh for global
    momentmat <- mclapply(c(1:nrow(df.hh)), momentmatcher, vfx = list(policycfx, policyifx), t0 = theta, 
                          data = df.hh[,-c(1:2)], mc.cores = detectCores())
    momentmat <- do.call(rbind, momentmat)
    print(sum(colMeans(momentmat)^2))
    return(sum(colMeans(momentmat)^2)) ## Use for global
    #return(momentmat) ## Use for GMM
  }
}

## Global optimizer
g <- mlsl(x0 = c(3, .95, 1.02, .5, .5, 1, .5, .9, .9),  ### reasonable numbers from the lit - also can try directL
          fn = gmmmomentmatcher, gr = momentgradient,  ### if using gradient, need to change df args to df.hh
            lower = c(1.01, .6, 1, .01, .01, .01, .01, .1, .5),  ### conservative bounds
            upper = c(10, .99, 1.2, .9, .9, 10, 1.5, .99, 1),
            control = list(ftol_rel = 1e-3, xtol_rel = 1e-3))  ## also set VFI tol to 1e-5
g

#g <- gmm(g = gmmmomentmatcher, x = df.hh, t0 = theta0,         
#         gradv = momentgradient, type = "twoStep", optfct = "nlminb", 
#         lower = c(1.01, .6, 1, .01, .01, .01, .01, .1, .5),
#         upper = c(10, .99, 1.2, .9, .9, 10, 1.5, .99, 1),
#         control = list(x.tol = 1e-4, rel.tol = 1e-4, abs.tol = 1e-10))
#summarize(g)

## change to g$par for GMM
theta <- g$sol
gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
saveRDS(theta, paste(getwd(), "/data/model_output/theta.rds", sep = ""))

### statespace
xn = 40; hn = 40

### Final Value function
statespace = create.statespace(ubm = c(20,20), theta, method = "equal")
v0 = cbind(statespace, data.frame("Tw" = rep(0, nrow(statespace)),
                                  "cfx" = rep(0, nrow(statespace)),
                                  "ifx" = rep(0, nrow(statespace))))
V = VFI(v0, theta)
saveRDS(V, paste(getwd(), "/data/model_output/V.rds", sep = ""))


