library(tidyverse)
library(parallel)
library(akima)
library(abind)
library(nloptr)
library(flatlandr)
options('nloptr.show.inequality.warning'=FALSE)

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
#ygrid = sort(unique(df$mu))
#hgrid = sort(unique(df$home_value)) + 1
xn = 55; yn = 20; hn = 50
hgrid = exp(seq(0, log(max(df$home_value)), length = hn))
ygrid = seq(8, max(df$mu), length = yn)
v0 = data.frame("Tw" = rep(0, xn*hn*yn),
                "cfx" = rep(0, xn*hn*yn),
                "ifx" = rep(0, xn*hn*yn))

theta = theta0
cbar = theta[4]; lambda = theta[5]; sigma = theta[6]
set.seed(37); yveclist = lapply(ygrid, function(y) rlnorm(1000, meanlog = y, sd = sigma*y))
Blist = lapply(yveclist, function(yvec) -lambda*mean(yvec))
cminlist = lapply(yveclist, function(yvec) cbar*mean(yvec))
xgridlist <- lapply(yveclist, function(yvec) c(seq(-lambda*mean(yvec), cbar*mean(yvec)-1, length = 10), exp(seq(log(cbar*mean(yvec)), log(10*max(yvec)), length = 45))))

statespace = lapply(c(1:length(ygrid)), create.statespace)
statespace = do.call(rbind, statespace)

bellman_operator <- function(w, statespace, theta){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; sigma = theta[6]; alpha = theta[7]
  optimizer <- function(x, h, y, Valfunc, B, cmin){
    set.seed(37); yvec = rlnorm(1000, meanlog = y, sd = sigma*y)
    if (x-B <= cmin) {
      vfx = u(cmin, h, gamma, alpha) + beta*Valfunc(B, h)
      return(c(cmin, 0, vfx))
    } else {
      hhprob = function(ci) {
        c = ci[1]
        i = ci[2]
        xtplus1 = R*(x - c - i) + yvec
        xtplus1[xtplus1<B] = B
        hp = h + i
        r = u(c, hp, gamma, alpha) + beta*mean(Valfunc(xtplus1, rep(hp, length(xtplus1))))
        return(-r)
      }
      constraint = function(ci){
        c <- ci[1]
        i <- ci[2]
        return(x - c - i - B)
      }
      sol = cobyla(x0 = c(mean(c(cmin, x-B)),0), fn = hhprob, lower = c(cmin,0), upper = c(x-B, x-B), hin = constraint, 
                   control = list(xtol_rel = 1e-2))
      return(c(sol$par, -sol$value))
    }
  }
  
  Valfunc = approxfun2(xgrid, hgrid, w[[1]][,,j]) ### need to fix
  
  Tw = mcmapply(optimizer, x = statespace$x, h = statespace$h, y = statespace$y, 
                Valfunc = list(Valfunc), B = statespace$B, cmin = statespace$cmin, mc.cores = detectCores()-2, mc.preschedule = FALSE)
  
  wnext = data.frame("Tw" = unlist(Tw)[seq(3,length(Tw),3)],
                 "cfx" = unlist(Tw)[seq(1,length(Tw),3)],
                 "ifx" = unlist(Tw)[seq(2,length(Tw),3)])
  return(wnext)
}

VFI <- function(vinit, tol = .05, maxiter = 100, k = 10){
  w = vector(mode = "list", length = maxiter-1)
  w[[1]] = vinit
  d = 1; i = 2
  while (d > tol & i < maxiter){
    vlhs = bellman_operator(xgridlist, yveclist, ygrid, hgrid, w[[i-1]], Blist, beta, R, gamma, alpha, cminlist, sigma)
    ## McQueen-Porteus Bounds
    blower = beta/(1-beta)*min(vlhs$Tw - w[[i-1]]$Tw)
    bupper = beta/(1-beta)*max(vlhs$Tw - w[[i-1]]$Tw)
    wpp = vlhs$Tw + (blower + bupper)/2
    ## Howard acceleration
    vk = vector(mode = "list", length = k)
    vk[[1]] = lapply(c(1:dim(wpp)[3]), function(x) approxfun2(xgridlist[[x]], hgrid, wpp[,,x]))
    for (i in c(1:k-1)) {
      howard = function(c, h, gamma, alpha, vk){
        hop = u(c,h,gamma,alpha) + beta*vk() ### need to make this of xtplus1
      }
      vk[[i+1]] = mcmapply(howard, , , gamma, alpha, vk) ### apply to all of the policys
    }
    w[[i]] = ### bellman_operator(vk[[k]])?
    ## check tol
    d = mean(((w[[i]]$cfx - w[[i-1]]$cfx)/w[[i]]$cfx)^2)
    i = i+1
  }
  return(pylr::compact(w))
}

### chebyshev, cmax > ymax

## plot iterations for given y and h
hidx = 1; yidx = 10
xs = lapply(w, function(x) x[,hidx,yidx])
xdat <- data.frame(do.call(cbind, xs)) %>%
  rownames_to_column("Xgrid") %>%
  pivot_longer(-Xgrid, names_to = "iteration", values_to = "fx") %>%
  mutate(iteration = as.integer(substr(iteration, 2, length(iteration))))

ggplot(xdat) +
  geom_line(aes(x = as.numeric(Xgrid), y = fx, color = iteration, group = iteration)) + 
  scale_color_viridis_c()

## heatmap for given iteration, y
hmap <- data.frame(w[[19]][,,1]) %>%
  rownames_to_column("Xgrid")%>%
  pivot_longer(-Xgrid, names_to = "hgrid", values_to = "vfx") %>%
  mutate(h = dense_rank(as.numeric(substr(hgrid, 2, length(hgrid)))),
         x = dense_rank(as.numeric(Xgrid)))

ggplot(hmap) +
  geom_tile(aes(x = x, y = h, fill = log(-vfx))) +
  scale_fill_viridis_c()

### Policy Function
pfx = bellman_operator(xgridlist, yveclist, ygrid, hgrid, w[[19]], Blist, beta, R, gamma, alpha, cminlist, sigma, policy = TRUE)

### Policy heatmaps
ggplot(pfx[[1]] %>% mutate(xrank = dense_rank(x), hrank = dense_rank(h))) +
  geom_tile(aes(x = xrank, y = hrank, fill = log(cfx))) +
  scale_fill_viridis_c()

ggplot(pfx[[1]] %>% mutate(xrank = dense_rank(x), hrank = dense_rank(h))) +
  geom_tile(aes(x = xrank, y = hrank, fill = ifx)) +
  scale_fill_viridis_c()

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
