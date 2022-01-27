library(tidyverse)
library(parallel)
library(nloptr)
library(flatlandr)
options('nloptr.show.inequality.warning'=FALSE)

df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500)
#df.hh <- read_csv("/home/mdg59/project/WBHRVS/full_panel.csv", guess_max = 7500)

### Functional forms
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))
#source("/home/mdg59/project/WBHRVS/VFIfunctions.R")
rm(bellman_operator, bellman_operator_BD, bufferstock, dfplot, policyfunc, ROI_setup, VFI, wtp100k, wtpcurve, u)

u <- function(c, h, gamma, alpha){(c^alpha * h^(1-alpha))^(1-gamma)/(1-gamma)}
### parameter vector

### initial guesses
gamma0 <- 2.988
beta0 <- .9
R0 <- 1.012
cbar0 <- .733
hbar0 <- .5
lambda0 <- 1.55
sigma0 <- .2
alpha0 <- .9
delta0 <- .95

theta0 <- c(gamma0, beta0, R0, cbar0, hbar0, lambda0, sigma0, alpha0, delta0)

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
xn = hn = 30; yn = 20

theta = theta0
gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]; lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]

ygrid = chebnodes(m = yn, lb = 8, ub = max(df$mu))
gqpts = c(-2.651961, -1.673552, -.8162879, 0, .8162879, 1.673552, 2.651961)
gqwts = c(.0009717812, .0545155828, .4256072526, .8102646176, .4256072526, .0545155828, .0009717812)
Blist = lapply(ygrid, function(y) -lambda*exp(y))
cminlist = lapply(ygrid, function(y) cbar*exp(y))
hminlist = lapply(ygrid, function(y) hbar*exp(y))

statespace = lapply(c(1:length(ygrid)), create.statespace)
statespace = do.call(rbind, statespace)
v0 = data.frame("Tw" = unlist(mapply(u, statespace$cmin, statespace$h, gamma, alpha)),
                "cfx" = statespace$cmin,
                "ifx" = rep(0, xn*hn*yn))

bellman_operator <- function(wlast, statespace){
  optimizer <- function(x, h, y, flist, B, cmin, hmin){
    yvec = sqrt(2)*sigma*exp(y)*gqpts + exp(y)
    Valfunc = flist[[which(ygrid==y)]]
    if (x-B <= cmin) {
      vfx = u(cmin, max(h,hmin), gamma, alpha) + beta*Valfunc(B, delta*max(h,hmin))
      return(c(cmin, max(0, hmin-h), vfx))
    } else if (x - B - cmin < hmin - h) {
      vfx = u(cmin, hmin, gamma, alpha) + beta*Valfunc(B, delta*hmin)
      return(c(cmin, hmin-h, vfx))
    } else {
      hhprob = function(ci) {
        c = ci[1]
        i = ci[2]
        xtplus1 = R*(x - c - i) + yvec
        hp = delta*(h + i)
        great.expectations = unlist(mapply(Valfunc, x = xtplus1, h = hp))*gqwts
        great.expectations[xtplus1<B] = Valfunc(B, delta*h)*gqwts[xtplus1<B]
        r = u(c, hp, gamma, alpha) + beta*sum(great.expectations)*pi^(-1/2)
        return(-r)
      }
      constraint = function(ci){
        c <- ci[1]
        i <- ci[2]
        return(x - c - i - B)
      }
      sps = c(mean(c(cmin, x-B)), mean(cmin, mean(c(cmin, x-B))))
      sol = cobyla(x0 = sps, fn = hhprob, lower = c(cmin,0), upper = c(x-B, x-B), hin = constraint,
                   control = list(xtol_rel = 1e-3))
      return(c(sol$par, -sol$value))
    }
  }
  
  ### interpolation functions
  flist = lapply(ygrid, interpolater.creater, statespace = statespace, Tw = wlast)
  
  Tw = mcmapply(optimizer, x = statespace$x, h = statespace$h, y = statespace$y, 
                flist = list(flist), B = statespace$B, cmin = statespace$cmin, hmin = statespace$hmin, 
                mc.cores = detectCores()-2, mc.preschedule = FALSE)
  
  wnext = data.frame("Tw" = as.vector(Tw[3,]),
                 "cfx" = as.vector(Tw[1,]),
                 "ifx" = as.vector(Tw[2,]))
  return(wnext)
}

howard = function(c, i, x, h, y, B, cmin, vk){
  Valfunc = vk[[which(ygrid==y)]]
  set.seed(37); yvec = rlnorm(1000, meanlog = y, sd = sigma*y)
  if (c==cmin){
    hop = u(cmin, h, gamma, alpha) + beta*Valfunc(B, h)
  } else{
    xtplus1 = R*(x - c - i) + yvec
    xtplus1[xtplus1<B] = B
    hp = h + i
    hop = u(c,hp,gamma,alpha) + beta*mean(unlist(sapply(xtplus1, Valfunc, hp)))
  }
  return(hop)
}

maxiter = 100; vinit = v0
VFI <- function(vinit, tol = 1e-4, maxiter = 100, k = 10){
  w = vector(mode = "list", length = maxiter-1)
  w[[1]] = vinit
  d = 1; i = 2
  while (d > tol & i < maxiter){
    vlhs = bellman_operator(wlast = w[[i-1]]$Tw, statespace)
    ## Howard acceleration
    #twk = vector(mode = "list", length = k)
    #twk[[1]] = vlhs$Tw
    #vk = vector(mode = "list", length = k)
    #vk[[1]] = lapply(ygrid, interpolater.creater, statespace = statespace, Tw = twk[[1]])
    #for (j in c(1:(k-1))) {
    #  twkp = unlist(mcmapply(howard, vlhs$cfx, vlhs$ifx, statespace$x, statespace$h, statespace$y, statespace$B, statespace$cmin, vk = list(vk[[j]]),
    #                  mc.cores = detectCores()-2))
      ## McQueen-Porteus Bounds
      blower = beta/(1-beta)*min(vlhs$Tw - w[[i-1]]$Tw)
      bupper = beta/(1-beta)*max(vlhs$Tw - w[[i-1]]$Tw)
      if (blower <0 & bupper <0) {vlhs$Tw = vlhs$Tw + (blower + bupper)/2}
    #  vk[[j+1]] = lapply(ygrid, interpolater.creater, statespace = statespace, Tw = twk[[j+1]])
    #}
    #vlhs$Tw = twk[[k]]
    w[[i]] = vlhs
    ## check tol
    d = mean(((w[[i]]$cfx - w[[i-1]]$cfx)/w[[i]]$cfx)^2)
    i = i+1
    print(i)
    saveRDS(w, "/users/mdgordo/Desktop/w.rds")
  }
  return(plyr::compact(w))
}

w = VFI(v0)

w = compact(readRDS("/users/mdgordo/Desktop/w.rds"))

### diff bt each iterations
for (i in c(2:length(w))) {
  d = mean(((w[[i]]$Tw - w[[i-1]]$Tw)/w[[i]]$Tw)^2)
  cd = mean(((w[[i]]$cfx - w[[i-1]]$cfx)/w[[i]]$cfx)^2)
  print(paste("vtol:", d, "ctol:", cd, sep = " "))
}

## plot iterations for given y and h
hidx = statespace$h[10]; yidx = ygrid[1]
xs = lapply(w, function(x) filter(cbind(x, statespace), y==yidx, h==hidx))
xs <- data.frame(do.call(rbind, xs)) 
xs$iteration <- rep(c(1:length(w)), each = xn)

ggplot(filter(xs, iteration < 60)) +
  geom_line(aes(x = x, y = Tw, color = iteration, group = iteration)) + 
  scale_color_viridis_c()

## plot iterations for given y and x
xidx = statespace$x[200]; yidx = ygrid[1]
xs = lapply(w, function(x) filter(cbind(x, statespace), y==yidx, x==xidx))
xs <- data.frame(do.call(rbind, xs)) 
xs$iteration <- rep(c(1:length(w)), each = hn)

ggplot(filter(xs, iteration<50)) +
  geom_line(aes(x = h, y = Tw, color = iteration, group = iteration)) + 
  scale_color_viridis_c()

## heatmap for given iteration, y
ggplot(filter(cbind(w[[length(w)]], statespace), y==yidx)) +
  geom_tile(aes(x = x, y = h, fill = log(-Tw))) +
  scale_fill_viridis_c()

### Policy heatmaps
ggplot(filter(cbind(w[[length(w)]], statespace), y==yidx)) +
  geom_tile(aes(x = x, y = h, fill = log(cfx))) +
  scale_fill_viridis_c()

ggplot(filter(cbind(w[[length(w)]], statespace), y==yidx)) +
  geom_tile(aes(x = x, y = h, fill = ifx)) +
  scale_fill_viridis_c()

### constraint satisfied?
df.cstrt <- cbind(w[[length(w)]], statespace) %>%
  group_by(y) %>%
  summarize(x = max(x),
            h = max(h), 
            spend = NA,
            maxifx = NA)

for (i in c(1:nrow(df.cstrt))) {
  a <- filter(cbind(w[[length(w)]], statespace), y==df.cstrt$y[i] & x==df.cstrt$x[i]) %>%
    mutate(spend = cfx + ifx)
  df.cstrt$spend[i] = min(a$spend)
  a <- filter(cbind(w[[length(w)]], statespace), y==df.cstrt$y[i] & h==df.cstrt$h[i])
  df.cstrt$maxifx[i] = max(a$ifx)
}

df.cstrt$ymax <- sqrt(2)*sigma*exp(ygrid)*gqpts[7] + exp(ygrid)
df.cstrt <- mutate(df.cstrt, xt1 = R*(x - spend) + ymax, 
                   ht1 = delta*(h + maxifx))


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

#nlprob <- OP(F_objective(objective, 6),
#             bounds = V_bound(li = c(1,2,3,4,5,6), lb = c(1,.5,1,0,0,0),
#                              ui = c(2,3,4), ub = c(1,2,1)))
#sol <- ROI_solve(nlprob, solver = "nloptr.cobyla", start = theta0)
#sol
#thetafin <- solution(sol)
#saveRDS(thetafin, "theta.rds")


V <- VFI(ygrid, hgrid, v0, tol = 1e-30, beta, gamma, R, alpha, sigma)

saveRDS(V, "V.rds")
