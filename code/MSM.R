### packages and data
library(tidyverse)
library(parallel)
library(nloptr)
library(fastGHQuad)
library(flatlandr)
library(gmm)
library(Rearrangement)
library(rdrobust)
options('nloptr.show.inequality.warning'=FALSE)
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))
source(paste(getwd(), "/code/rdhelpers.R", sep = ""))

df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500) 
#df.hh <- read_csv("/home/mdg59/project/WBHRVS/full_panel.csv", guess_max = 7500)

b <- optbw("consumption", b = "dist_2_seg13", df = df.hh, fuzzy = TRUE, dist.exclude = "none")[1,3]

df.hh <- df.hh %>%
  filter(abs(dist_2_seg13) < b & designation !="none") %>%
  select(hhid, wave, wt_hh, food_consumption, total_income, var_inc, avg_inc, 
         home_value, home_investment, imputed_bufferstock, quake_aid) %>% ## add damages?
  group_by(hhid) %>%
  mutate(lag_h = lag(home_value, order_by = wave),
         lag_x = lag(imputed_bufferstock, order_by = wave),
         lag_y = lag(total_income, order_by = wave),
         lag_c = lag(food_consumption, order_by = wave),
         lag_i = lag(home_investment, order_by = wave)) %>%
  filter(avg_inc>0 & food_consumption>0 & lag_h + home_investment>0)
df.hh <- df.hh[complete.cases(df.hh),]

### parameters
gamma0 <- 4.77; beta0 <- .9; R0 <- 1.01; cbar0 <- .1; hbar0 <- .1
lambda0 <- 1.5; sigma0 <- .06; alpha0 <- .8; delta0 <- .9

theta0 <- c(gamma0, beta0, R0, cbar0, hbar0, lambda0, sigma0, alpha0, delta0)

### Grid size
xn = 40; hn = 40; yn = 10
ygrid = as.vector(quantile(df.hh$avg_inc, seq(0,1,length.out = yn), na.rm = TRUE))

### Define GMM function

gmmmomentmatcher <- function(theta, df) {
  print(theta)
  saveRDS(theta, paste(getwd(), "/data/model_output/theta.rds", sep = ""))
  
  ### VFI
  V = VFI(theta)
  saveRDS(V, paste(getwd(), "/data/model_output/V.rds", sep = ""))
  finalV <- V[[length(V)]]
  
  ### data
  momentmat <- mclapply(c(1:nrow(df)), momentmatcher, vfx = finalV, t0 = theta, 
                        data = df[,-c(1:2)], mc.cores = detectCores())
  momentmat <- do.call(rbind, momentmat)
  print(sum(colMeans(momentmat)^2))
  return(momentmat)
}

g <- gmm(g = gmmmomentmatcher, x = df.hh, t0 = theta0,
         gradv = momentgradient, type = "twoStep", onlyCoefficients=TRUE, optfct = "nlminb", 
         lower = c(1.01, .6, .9, .01, .01, 0, 0, .1, .5),
         upper = c(10, .99, 1.65, 1, 1, 10, 1, .99, 1))

summarize(g)

theta <- g$coefficients
gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
saveRDS(theta, paste(getwd(), "/data/model_output/theta.rds", sep = ""))

### statespace
xn = 40; hn = 40; yn = 10
ygrid = as.vector(quantile(df.hh$avg_inc, seq(0,1,length.out = yn), na.rm = TRUE))

### Final Value function
V = VFI(theta)
saveRDS(V, paste(getwd(), "/data/model_output/V.rds", sep = ""))


