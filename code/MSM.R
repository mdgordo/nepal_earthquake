library(tidyverse)
library(parallel)
library(ROI)
library(ROI.plugin.alabama)
library(ROI.plugin.nloptr)
library(ROI.plugin.deoptim)

df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500)
#df.hh <- read_csv("/home/mdg59/project/WBHRVS/full_panel.csv", guess_max = 7500)

### Functional forms
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))
rm(bellman_operator, bellman_operator_BD, bufferstock, dfplot, policyfunc, ROI_setup, VFI, wtp100k, wtpcurve, u)

u <- function(c, h, gamma, alpha){(c^alpha * h^(1-alpha))^(1-gamma)/(1-gamma)}
### parameter vector

### initial guesses
gamma0 <- 2.988
beta0 <- .9
R0 <- 1.012
cbar0 <- .733
lambda0 <- 1.55
sigma0 <- .032
alpha0 <- .9

theta0 <- c(gamma0, beta0, R0, cbar0, lambda0, sigma0, alpha0)

aid_amt = 100000

### State vars

## Permanent income
### we want to match on non durable consumption and non transfer income - percapita? - use machine learning? - doing wave 1 at the moment
df <- filter(df.hh, wave==1 & !is.na(income_gross) & !is.na(age_hh) & !is.na(class5) & !is.na(class10) & !is.na(femalehh) & !is.na(caste_recode))
cissebarret <- lm(log(income_gross+1) ~ as.factor(caste_recode) + poly(age_hh, 2) + femalehh + class5 + class10 +
                    slope + elevation + as.factor(land_qtle) #+ poly(log(lag_income_gross+1), 3) 
                  , data = df, weights = wt_hh, na.action = na.omit)

### coarsen for feasibility
df$mu <- round(cissebarret$fitted.values,1)

## Housing - do I want to do a regression here on housing characteristics?


### Moments
thresholds <- seq(40000, 200000, 20000)
pdfmoments_t <- c(0, 0, 0, .04, .25, .52, .64, .77, .86)
#pdfmoments_u <- c(0.05, 0.18, 0.32, 0.57, 0.63, 0.77, 0.79, 0.86, 0.92)
mkt_clear <- 0

moments <- c(pdfmoments_t, mkt_clear)

### Optimization setup
beta = beta0; R = R0; gamma=gamma0; alpha=alpha0; sigma = sigma0; lambda = lambda0; cbar = cbar0
ygrid = sort(unique(df$mu))
hgrid = sort(unique(df$home_value)) + 1
v0 = array(0, dim = c(length(hgrid), 300, length(ygrid)))

bellman_operator <- function(xgridlist, yveclist, ygrid, hgrid, w, Blist, beta, R, gamma, alpha, cminlist, sigma){
  optimizer <- function(x, h, yvec, Valfunc, B, cmin){
    if (x-B <= cmin) {
      return(u(cmin, h, gamma, alpha) + beta*Valfunc(B, h))
    } else {
      hhprob = function(ci) {
        c = ci[1]
        i = ci[2]
        xtplus1 = R*(x - c - i) + yvec
        xtplus1[xtplus1<B] = B
        hp = h + i
        r = u(c, hp, gamma, alpha) + beta*mean(Valfunc(xtplus1, rep(hp, length(xtplus1))))
        return(r)
      }
      constraint = function(ci){
        c <- ci[1]
        i <- ci[2]
        return(x - c - i - B)
      }
      nlprob = OP(F_objective(hhprob, 2),
                   F_constraint(constraint, dir = ">=", rhs = 0),
                   maximum = TRUE,
                   bounds = V_bound(li = c(1,2), lb = c(cmin,0),
                                    ui = c(1,2), ub = c(x-B, x-B)))
      rsol = ROI_solve(nlprob, solver = "nloptr.cobyla", start = c(cmin + 1, 1))
      return(rsol$objval)
    }
  }
  
  wp = vector(mode = "list", length = length(ygrid))
  
  for (j in c(1:length(ygrid))) {
    yvec = yveclist[[j]]
    B = Blist[[j]]
    cmin = cminlist[[j]]
    xgrid = xgridlist[[j]]
    Valfunc = approxfun2(xgrid, hgrid, w[,,j])
    xhcross = crossing(x = xgrid, h = hgrid)
    xhcross$Tw = mcmapply(optimizer, x = xhcross$x, h = xhcross$h, yvec = list(yvec), 
                  Valfunc = list(Valfunc), B = list(B), cmin = list(cmin), mc.cores = detectCores()-2)
    wp[[j]] = as.matrix(pivot_wider(xhcross, id_cols = x, names_from = h, values_from = Tw) %>% 
            column_to_rownames(var="x"))
    
  }
  wp = abind(wp, along = 3)
  return(wp)
}


VFI <- function(ygrid, hgrid, vinit, tol = 1e-30, maxiter = 50, beta, R, gamma, alpha, sigma){
  yveclist = lapply(ygrid, function(y) rlnorm(1000, meanlog = y, sd = sigma*y))
  Blist = lapply(yveclist, function(yvec) -lambda*mean(yvec))
  cminlist = lapply(yveclist, function(yvec) cbar*mean(yvec))
  xgridlist <- lapply(yveclist, function(yvec) seq(-lambda*mean(yvec), 10*max(yvec), len = 300))
  w = vector(mode = "list", length = maxiter-1)
  w[[1]] = vinit
  d = 1; i = 2
  while (d > tol & i < maxiter){
    w[[i]] = bellman_operator(xgridlist, yveclist, ygrid, hgrid, w[[i-1]], Blist, cminlist, beta, R, gamma, alpha, sigma)
    d = sqrt(sum((w[[i]] - w[[i-1]])^2))
    i = i+1
  }
  return(pylr::compact(w))
}

VFI(ygrid, hgrid, v0, tol = 1e-30, maxiter = 3, beta, gamma, R, alpha, sigma)

policyfunc <- function(x, Vfx, B, beta, R, y, gamma, cmin){
  if (x - B <= cmin) {return(cmin)} else{
    objective <- function(c) {
      xtplus1 = R*(x - c) + y
      xtplus1 = if_else(xtplus1<B, B, xtplus1)
      r = u(c, gamma) + beta*mean(Vfx(xtplus1))
      return(r)
    }
    l <- optimize(objective, interval = c(cmin, x - B), maximum = TRUE)
    return(l$maximum)
  }
}

bufferstock <- function(hhid, Vlist){
  ### pull right parameters
  c = df$food_consumption[df$hhid==hhid]
  mu = df$mu[df$hhid==hhid]
  i = which(unlist(lapply(Vlist, function(x) x$mu))==mu)
  cfx = Vlist[[i]]$cfx
  xgrid = Vlist[[i]]$xgrid
  
  ### find root of policy function
  cfxrt = function(x){cfx(x) - c}
  if (c>cfx(max(xgrid))) {r = max(xgrid)} else if (c<cfx(min(xgrid))) {
    r = min(xgrid)
  } else {
    r = uniroot(cfxrt, interval = c(1, max(xgrid)), extendInt = "upX")$root
  }
  ### pre and post aid buffer stock
  bsx_pre = r - df$quake_aid[df$hhid==hhid]
  bsx_post = r + aid_amt
  bsx_actual = r
  
  ### Calculate counterfactual consumption
  cdist_pre = cfx(bsx_pre)
  cdist_post = cfx(bsx_post)
  return(c(cdist_pre, cdist_post, bsx_pre, bsx_post, bsx_actual))
}

### start with good initial guess

y0 = rlnorm(1000, meanlog = 10.5, sd = sigma0*10.5)
B0 = -lambda0*mean(y0)
cmin0 = cbar0*mean(y0)
xgrid0 <- seq(B0, 10*max(y0), len = 2000)
v0 <- rep(0, length(xgrid0))
V0 = VFI(grid = xgrid0, vinit = v0, B = B0, beta = beta0, R = R0, y = y0, 
        gamma = gamma0, cmin = cmin0)
Vfx0 = approxfun(xgrid0, V0[,dim(V0)[2]], rule = 2)

rm(y0, B0, cmin0, xgrid0, v0, V0)

### MSM

msm_func <- function(theta){
  print(theta)
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; lambda = theta[5]; sigma = theta[6]
  
  Vlist = rep(list(list("mu" = NA, "xgrid" = NA, "Vfx" = NA, "cfx" = NA)), length(unique(df$mu)))
  
  ### Solve for all the value functions - might have to do this for pre and post interest rates
  for (i in c(1:length(unique(df$mu)))){
    mu = sort(unique(df$mu))[i]
    Vlist[[i]]$mu = mu
    
    ### productivity draws and lower bound
    y = rlnorm(1000, meanlog = mu, sd = sigma*mu)
    B = -lambda*mean(y)
    cmin = cbar*mean(y)
    
    ### Initial grid and guess
    xgrid <- seq(B, 10*max(y), len = 2000)
    Vlist[[i]]$xgrid = xgrid
    v0 <- Vfx0(xgrid)
    
    ### Solve for value function
    V = VFI(grid = xgrid, vinit = v0, B = B, beta = beta, R = R, y = y, 
            gamma = gamma, cmin = cmin)
    Vfx = approxfun(xgrid, V[,dim(V)[2]])
    Vlist[[i]]$Vfx = Vfx
    
    ### Solve for policy function
    cfx = mclapply(xgrid, policyfunc, Vfx = Vfx, 
                   B = B, beta = beta, R = R, y = y, gamma = gamma, cmin = cmin, mc.cores = detectCores()-2)
    Vlist[[i]]$cfx = approxfun(xgrid, cfx, rule = 2)
    
    #print(i)
    #saveRDS(Vlist, "Vlist.rds")
  }
  
  ### Calculate everyone's buffer stock and counterfactual consumption - 50000 right number to use?
  bsx = mclapply(df$hhid, bufferstock, Vlist, mc.cores = detectCores()-2)
  cdist_post = unlist(lapply(bsx, function(x) x[2]))
  bsx_actual = unlist(lapply(bsx, function(x) x[5]))
  
  ### Calculate fraction of hhs under each threshold with and without aid
  pdf_post = unlist(lapply(thresholds, function(x) sum(cdist_post < x)/length(cdist_post)))
  
  ### Does market clear? 
  mkt_mod = sum(bsx_actual - df$food_consumption)
  
  modeledmoments = c(pdf_post, mkt_mod)
  g = modeledmoments - moments
  return(g)
}

#print(msm_func(theta0))

objective <- function(theta) {
  g = msm_func(theta)
  sse = sum(g^2)
  return(sse)
}

nlprob <- OP(F_objective(objective, 6),
             bounds = V_bound(li = c(1,2,3,4,5,6), lb = c(1,.5,1,0,0,0),
                              ui = c(2,3,4), ub = c(1,2,1)))
sol <- ROI_solve(nlprob, solver = "nloptr.cobyla", start = theta0)
sol
thetafin <- solution(sol)
saveRDS(thetafin, "theta.rds")

