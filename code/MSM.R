### packages and data
library(tidyverse)
library(parallel)
library(nloptr)
library(flatlandr)
library(rdrobust)
library(gmm)
options('nloptr.show.inequality.warning'=FALSE)
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))
source(paste(getwd(), "/code/rdhelpers.r", sep = ""))

df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500) 
#df.hh <- read_csv("/home/mdg59/project/WBHRVS/full_panel.csv", guess_max = 7500)

b <- optbw("consumption", b = "dist_2_seg13", df = df.hh, fuzzy = TRUE, dist.exclude = "none")[1,3]

df.hh <- df.hh %>%
  filter(abs(dist_2_seg13) < b & designation !="none") %>%
  select(hhid, wave, wt_hh, food_consumption, total_income, var_inc, avg_inc, home_value, home_investment, imputed_bufferstock, quake_aid) %>%
  group_by(hhid) %>%
  mutate(lag_h = lag(home_value, order_by = wave),
         lag_x = lag(imputed_bufferstock, order_by = wave),
         lag_y = lag(total_income, order_by = wave),
         lag_c = lag(food_consumption, order_by = wave),
         lag_i = lag(home_investment, order_by = wave)) 
df.hh <- df.hh[complete.cases(df.hh),]

### parameters
gamma0 <- 2.988
beta0 <- .9
R0 <- 1.012
cbar0 <- .733
hbar0 <- .6
lambda0 <- 1.2
sigma0 <- .2
alpha0 <- .9
delta0 <- .95

theta0 <- c(gamma0, beta0, R0, cbar0, hbar0, lambda0, sigma0, alpha0, delta0)

### Shocks - Gaussian Quadrature points
gqpts = c(-2.651961, -1.673552, -.8162879, 0, .8162879, 1.673552, 2.651961)
gqwts = c(.0009717812, .0545155828, .4256072526, .8102646176, .4256072526, .0545155828, .0009717812)/sqrt(pi)

### Grid size
xn = hn = 30; yn = 15
ygrid = seq(8, 13, length.out = yn)


### Define GMM function

gmmmomentmatcher <- function(theta, data) {
  print(theta)
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  
  ### statespace
  Blist = lapply(ygrid, function(y) -lambda*exp(y))
  cminlist = lapply(ygrid, function(y) cbar*exp(y))
  hminlist = lapply(ygrid, function(y) hbar*exp(y))
  
  statespace = lapply(c(1:length(ygrid)), create.statespace, Blist = Blist, cminlist = cminlist, 
                      hminlist = hminlist)
  statespace = do.call(rbind, statespace)
  
  ### Initial guess 
  v0 = guesser(statespace, theta)
  
  ### VFI
  V = VFI(v0, statespace, theta = theta)
  V = compact(V)
  
  finalV <- cbind(statespace, V[[length(V)]])
  
  ### data
  momentmat <- mclapply(c(1:nrow(data)), momentmatcher, vfx = finalV, theta = theta, data = data, mc.cores = detectCores())
  momentmat <- do.call(rbind, momentmat)
  return(momentmat)
}

g <- gmm(gmmmomentmatcher, x = df.hh[,-c(1:2)], t0 = theta0, optfct = "nlminb", 
         upper = c(Inf, 1, Inf, 1, 1, Inf, Inf, 1, 1), lower = c(1, 0, 0, 0, 0, 0, 0, 0, 0),
         onlyCoefficients = TRUE)

saveRDS(g, "g.rds")
