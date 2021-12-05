library(tidyverse)
library(parallel)

df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500)

source(paste(getwd(), "/code/VFIfunctions.r", sep = ""))

### parameter vector

### initial guesses
gamma <- 3
beta <- .95
R <- 1.03
cmin = 10000
lambda = .2
sigma = 3.85

theta <- c(gamma, beta, R, cmin, lambda, sigma)

u <- function(x, gamma){x^(1-gamma)/(1-gamma)}

### we want to match on non durable consumption and non transfer income - percapita?
df <- filter(df.hh, !is.na(income_gross) & !is.na(lag_income_gross) & !is.na(age_hh) & !is.na(class5))
cissebarret <- lm(log(income_gross+1) ~ as.factor(caste_recode) + poly(age_hh, 2) + femalehh + class5 + class10 +
                    slope + elevation + as.factor(land_qtle) #+ poly(log(lag_income_gross+1), 3) 
                  , data = df, weights = wt_hh, na.action = na.omit)
winc_guess <- as.vector(coef(cissebarret))

df$resid <- cissebarret$residuals^2
cissebarret_resid <- lm(resid ~ as.factor(caste_recode) + poly(age_hh, 2) + femalehh + class5 + class10 +
                          slope + elevation + as.factor(land_qtle) #+ poly(log(lag_income_gross+1), 3) 
                        , data = df, weights = wt_hh, na.action = na.omit)
wresid_guess <- as.vector(coef(cissebarret_resid))

#ggplot() +geom_histogram(aes(x = df$mu), bins = 100)
#ggplot() +geom_histogram(aes(x = df$sigma), bins = 100)
#ggplot() +geom_point(aes(x = df$mu, y = df$sigma))
#ggplot() +geom_point(aes(x = log(df$income_gross+1), y = cissebarret$residuals^2))

### coarsen for feasibility - distributions should be similar across hh years because most of w vars don't change
df$mu <- round(cissebarret$fitted.values,1)
df$sigma <- round(sqrt(cissebarret_resid$fitted.values),)

musig <- unique(cbind(df$mu, df$sigma))

### Moments

pdfmoments_t <- c(0, .04, .25, .52, .64, .77, .86)
pdfmoments_u <- c(0.05, 0.18, 0.32, 0.57, 0.63, 0.79, 0.77, 0.86, 0.92)
mkt_clear <- 0

moments <- c(pdfmoments_t, pdfmoments_u, mkt_clear)

### Optimization setup
bellman_operator <- function(grid, w, B, beta, R, mu, sigma, gamma, cmin){
  Valfunc = approxfun(grid, w, rule = 2)
  y = rlnorm(1000, meanlog = mu, sd = sigma)
  optimizer <- function(x){
    if (x==B) {return(u(cmin, gamma) + beta*mean(Valfunc(B)))} else{
      objective <- function(c) {
        xtplus1 = R*(x - c) + y
        xtplus1 = if_else(xtplus1<B, B, xtplus1)
        r = u(c, gamma) + beta*mean(Valfunc(xtplus1))
        return(r)
      }
      l <- optimize(objective, interval = c(1, x - B), maximum = TRUE)
      return(l$objective)
    }
  }
  Tw = rep(NA, length(grid))
  Tw = mclapply(grid, optimizer, mc.cores = 6)
  return(unlist(Tw))
}

VFI <- function(grid, vinit, tol = 1e-9, maxiter = 300, B, beta, R, mu, sigma, gamma, cmin){
  w = matrix(0, length(grid), 1)
  w[,1] = vinit
  d = 1
  i = 2
  while (d > tol & i < maxiter){
    w = cbind(w, rep(0, length(grid)))
    w[,i] = bellman_operator(grid, w[,i-1], B, beta, R, mu, sigma, gamma, cmin)
    d = sqrt(sum((w[,i] - w[,i-1])^2))
    i = i+1
  }
  return(w)
}


### MSM

msm_func <- function(theta){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cmin = theta[4]; lambda = theta[5]
  
  ### Solve for all the value functions
  
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


