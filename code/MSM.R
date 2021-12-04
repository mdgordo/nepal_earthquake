library(tidyverse)
library(parallel)

df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500)

source(paste(getwd(), "/code/VFIfunctions.r", sep = ""))

### parameter vector

### initial guesses
gamma <- 3
beta <- .95
R <- 1.03
alpha = .3
mu = 11
sigma = .35
cmin = 10000
lambda = .2

theta <- c(gamma, beta, R, alpha, mu, sigma, cmin, lambda)

### we want to match on non durable consumption and non transfer income - percapita?
df <- filter(df.hh, !is.na(income_gross) & !is.na(lag_income_gross) & !is.na(age_hh) & !is.na(class5))
cissebarret <- lm(log(income_gross+1) ~ as.factor(caste_recode) + poly(age_hh, 2) + femalehh + class5 + class10 +
                    slope + elevation + poly(log(lag_income_gross+1), 3) + as.factor(land_qtle),
                  data = df, weights = wt_hh, na.action = na.omit)
winc_guess <- as.vector(coef(cissebarret))

df$resid <- cissebarret$residuals^2
cissebarret_resid <- lm(resid ~ as.factor(caste_recode) + poly(age_hh, 2) + femalehh + class5 + class10 +
                    slope + elevation + poly(log(lag_income_gross+1), 3) + as.factor(land_qtle),
                  data = df, weights = wt_hh, na.action = na.omit)
wresid_guess <- as.vector(coef(cissebarret_resid))
### need to use to predict income and variance - what moments to match here

### PDF Moments

empiricalmoments <- c(0, .04, .25, .52, .64, .77, .86)

### utility and production functional forms
u <- function(x, gamma){x^(1-gamma)/(1-gamma)}
f <- function(k, epsilon, alpha) {epsilon*k^alpha}

### Optimization setup
ROI_setup <- function(Vfx, x, theta, epsilon, B, roots = FALSE){
  objective <- function(ck) {
    c <- ck[1]
    k <- ck[2]
    xtplus1 = theta[3]*(x - c - k) + f(k, epsilon, theta[4])
    r = u(c, theta[1]) + theta[2]*mean(Vfx(xtplus1))
    return(r)
  }
  ### has to hold with certainty, so min of production
  constraint <- function(ck){
    c <- ck[1]
    k <- ck[2]
    return(B - theta[3]*(x - c - k) - min(f(k, epsilon, theta[4])))
  }
  nlprob <- OP(F_objective(objective, 2),
               F_constraint(constraint, dir = "<=", rhs = 0),
               maximum = TRUE,
               bounds = V_bound(li = c(1,2), lb = c(1,1)))
  r <- ROI_solve(nlprob, solver = "nloptr.cobyla", start = c(1.1, 1.1))
  if (roots=TRUE) return(solution(r)) else return(objective(solution(r)))
}

bellman_operator <- function(grid, w, theta, epsilon, B){
  Vfx = approxfun(grid, w, rule = 2)
  Tw = rep(NA, length(grid))
  Tw = mclapply(xgrid, ROI_setup, Vfx = Vfx, theta, epsilon, B, mc.cores = 6, mc.preschedule = FALSE)
  return(unlist(Tw))
}

#grid = xgrid; vinit = v0; tol = 1e-13; maxiter = 300

VFI <- function(grid, vinit, tol = 1e-13, maxiter = 300, theta, epsilon, B){
  w = matrix(0, length(grid), 1)
  w[,1] = vinit
  d = 1
  i = 2
  while (d > tol & i < maxiter){
    w = cbind(w, rep(0, length(grid)))
    w[,i] = bellman_operator(grid, w[,i-1], theta, epsilon, B)
    d = sqrt(sum((w[,i] - w[,i-1])^2))
    i = i+1
  }
  return(w)
}

msm_func <- function(theta){
  ### productivity draws and lower bound
  epsilon <- rlnorm(1000, theta[5], theta[6])
  B <- -min(f(1, epsilon, theta[4]))/(theta[3]-1)
  
  ### Initial grid and guess
  xgrid <- seq(B + 3, 5000000, len = 5000)
  v0 <- rep(0, length(xgrid))
  
  ### Solve for value function
  V = VFI(xgrid, v0, theta, epsilon, B)
  Valfunc = approxfun(xgrid, V[,dim(V)[2]])
  
  ### Solve for policy function
  cfunc = mclapply(xgrid, ROI_setup, Vfx = Vfx, theta, epsilon, B, roots = TRUE, mc.cores = 6, mc.preschedule = FALSE)
  
  ### Calculate everyone's buffer stock
  
  
  ### Calculate fraction of hhs under each threshold with and without aid
  
  modeledmoments = c()
  g = modeledmoments - empiricalmoments
}


