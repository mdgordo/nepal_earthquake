library(tidyverse)
library(parallel)

### parameters from MSM.R
gamma <- 4.49
beta <- .934
R <- 1.026
cbar <- .5
lambda <- .418
sigma <- 7.418

u <- function(x){
  u <- x^(1-gamma)/(1-gamma)
  return(u)
}

#df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500)
df.hh <- read_csv("/home/mdg59/project/WBHRVS/full_panel.csv", guess_max = 7500)

### we want to match on non durable consumption and non transfer income - percapita? - use machine learning? - doing wave 1 at the moment
df <- filter(df.hh, wave==1 & !is.na(income_gross) & !is.na(age_hh) & !is.na(class5) & !is.na(class10) & !is.na(femalehh) & !is.na(caste_recode))
cissebarret <- lm(log(income_gross+1) ~ as.factor(caste_recode) + poly(age_hh, 2) + femalehh + class5 + class10 +
                    slope + elevation + as.factor(land_qtle) #+ poly(log(lag_income_gross+1), 3) 
                  , data = df, weights = wt_hh, na.action = na.omit)

### coarsen for feasibility
df$mu <- round(cissebarret$fitted.values,1)
mus <- seq(8.5,13,.0025)

### VFI functions
source("/home/mdg59/project/WBHRVS/VFIfunctions.R")

Vlist = rep(list(list("mu" = NA, "xgrid" = NA, "Vfx" = NA, "cfx" = NA)), length(mus))

### Solve for all the value functions - lower bounding at exp(8.5 + (7.4/8)^2/2) = 7538
for (i in c(1:length(mus))){
  mu = mus[i]
  Vlist[[i]]$mu = mu
  
  ### productivity draws and lower bound
  y = rlnorm(1000, meanlog = mu, sd = sigma/mu)
  B = -lambda*mean(y)
  cmin = cbar*mean(y)
  
  ### Initial grid and guess
  xgrid <- seq(B, 10*max(y), len = 7000)
  Vlist[[i]]$xgrid = xgrid
  v0 <- rep(0, length(xgrid))
  
  ### Solve for value function
  V = VFI(grid = xgrid, vinit = v0)
  Vfx = approxfun(xgrid, V[,dim(V)[2]])
  Vlist[[i]]$Vfx = Vfx
  
  ### Solve for policy function
  cfx = mclapply(xgrid, policyfunc, Vfx = Vfx, mc.cores = detectCores()-2)
  Vlist[[i]]$cfx = approxfun(xgrid, cfx, rule = 2)
  
  print(i)
  
}

saveRDS(Vlist, "Vlist.rds")

### Solve for each hhs buffer stock

bsx = mclapply(df$hhid, bufferstock, Vlist, mc.cores = detectCores()-2)
cdist_pre = unlist(lapply(bsx, function(x) x[1]))
cdist_post = unlist(lapply(bsx, function(x) x[2]))
bsx_pre = unlist(lapply(bsx, function(x) x[3]))
bsx_post = unlist(lapply(bsx, function(x) x[4]))
bsx_actual = unlist(lapply(bsx, function(x) x[5]))

wtp100k <- function(hhid){
  x = bsx_pre[df$hhid==hhid]
  mu = df$mu[df$hhid==hhid]
  i = which(unlist(lapply(Vlist, function(x) x$mu))==mu)
  Vfx = Vlist[[i]]$Vfx
  V = Vfx(x)
  val4aid = unlist(lapply(Vlist, function(f) f$Vfx(x+100000)))
  vi = which.min(abs(V - val4aid))
  mu0 = exp(mu + (sigma/mu)^2/2)
  muprime = exp(mus[vi] + (sigma/mu)^2/2)
  return((mu0 - muprime)/(1 - beta))
}

wtps <- mclapply(df$hhid, wtp100k, mc.cores = detectCores()-2)
saveRDS(wtps, "wtps.rds")



