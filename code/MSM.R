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

### coarsen for feasibility
df$mu <- round(cissebarret$fitted.values,1)

### Moments
thresholds <- seq(40000, 200000, 20000)
pdfmoments_t <- c(0, 0, 0, .04, .25, .52, .64, .77, .86)
pdfmoments_u <- c(0.05, 0.18, 0.32, 0.57, 0.63, 0.79, 0.77, 0.86, 0.92)
mkt_clear <- 0

moments <- c(pdfmoments_t, pdfmoments_u, mkt_clear)

### Optimization setup
bellman_operator <- function(grid, w, B, beta, R, y, gamma, cmin){
  Valfunc = approxfun(grid, w, rule = 2)
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

VFI <- function(grid, vinit, tol = 1e-9, maxiter = 300, B, beta, R, y, gamma, cmin){
  w = matrix(0, length(grid), 1)
  w[,1] = vinit
  d = 1
  i = 2
  while (d > tol & i < maxiter){
    w = cbind(w, rep(0, length(grid)))
    w[,i] = bellman_operator(grid, w[,i-1], B, beta, R, y, gamma, cmin)
    d = sqrt(sum((w[,i] - w[,i-1])^2))
    i = i+1
  }
  return(w)
}

policyfunc <- function(x, Vfx, B, beta, R, y, gamma, cmin){
  if (x==B) {return(cmin)} else{
    objective <- function(c) {
      xtplus1 = R*(x - c) + y
      xtplus1 = if_else(xtplus1<B, B, xtplus1)
      r = u(c, gamma) + beta*mean(Vfx(xtplus1))
      return(r)
    }
    l <- optimize(objective, interval = c(1, x - B), maximum = TRUE)
    return(l$maximum)
  }
}

bufferstock <- function(hhid, Vlist, xgrid){
  ### pull right parameters
  c = df$food_consumption[df$hhid==hhid]
  mu = df$mu[df$hhid==hhid]
  i = which(lapply(Vlist, function(x) x$mu==mu))
  cfx = Vlist[[i]]$cfx
  
  ### find root of policy function
  cfxrt <- function(x){cfx(x) - c}
  r <- uniroot(cfxrt, interval = c(1, max(xgrid)), extendInt = "upX")$root
  
  ### pre and post aid buffer stock
  bsx_pre = r - df$quake_aid[df$hhid==hhid]
  bsx_post = r + 50000
  
  ### Calculate counterfactual consumption
  cdist_pre = cfx(bsx_pre)
  cdist_post = cfx(bsx_post)
  return(c(cdist_pre, cdist_post))
}


### MSM

msm_func <- function(theta){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cmin = theta[4]; lambda = theta[5]; sigma = theta[6]
  
  Vlist = vector(mode = "list", length = length(unique(df$mu)))
  
  ### Solve for all the value functions - might have to do this for pre and post interest rates
  for (i in sort(unique(df$mu))){
    mu = sort(unique(df$mu))[i]
    Vlist[[i]]$mu = mu
    
    ### productivity draws and lower bound
    y = rlnorm(1000, meanlog = mu, sd = sigma/mu)
    B = -lambda*mean(y)
    
    ### Initial grid and guess
    xgrid <- seq(B, 10*max(y), len = 5000)
    v0 <- rep(0, length(xgrid))
    
    ### Solve for value function
    V = VFI(xgrid, v0, B, beta, R, y, gamma, cmin)
    Vlist[[i]]$Vfx = approxfun(xgrid, V[,dim(V)[2]])
    
    ### Solve for policy function
    cfx <- mclapply(xgrid, policyfunc, Vfx = Vfx, mc.cores = 6)
    Vlist[[i]]$cfx <- approxfun(xgrid, cfx, rule = 2)
    
  }
  
  ### Calculate everyone's buffer stock and counterfactual consumption - 50000 right number to use?
  bsx = mclapply(df$hhid, bufferstock, Vlist, xgrid, mc.cores = 6)
  cdist_pre = unlist(lapply(bsx, function(x) x[1]))
  cdist_post = unlist(lapply(bsx, function(x) x[2]))
  
  ### Calculate fraction of hhs under each threshold with and without aid
  pdf_pre = lapply(thresholds, function(x) sum(cdist_pre < x)/length(cdist_pre))
  pdf_post = lapply(thresholds, function(x) sum(cdist_post < x)/length(cdist_post))
  
  ### Does market clear? how does aid influx affect interest rate?
  
  modeledmoments = c()
  g = modeledmoments - moments
}


