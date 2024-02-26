### packages and data
library(tidyverse)
library(nloptr)
library(fastGHQuad)
#devtools::install_github("mdgordo/flatlandr",ref="master")
library(flatlandr)
library(gmm)
library(Rearrangement)
library(mgcv)
library(neuralnet)
library(parallel)
options('nloptr.show.inequality.warning'=FALSE)
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))

lbounds <- c(1.01, .85, 1, .01, .01, .01, .01, .2, .85)
ranges <- c(9.99, .14, .1, .99, .99, 1.99, 1.99, .7, .14)
ubounds <- lbounds + ranges
ui = rbind(c(0,-1,-1,0,0,0,0,0,0),
           diag(9),
           -1*diag(9))
ci = c(-2.03, lbounds, -1*ubounds) ## 2.03 approximates multiplicative constraint on BR < 1

### load data
df.adj <- readRDS(paste(getwd(), "/data/model_output/df.adj.rds", sep = ""))
df.adj <- df.adj[complete.cases(df.adj),]

### Grid size
xn = 40; hn = 40

### GMM 
## first coarse global
globalwrap = function(t){
  if (t[2]*t[3]>1) {
    return(99999)
  } else {
  momentmat = gmmmomentmatcher(t, df.adj)
  gvec = colMeans(momentmat); 
  g = sum(gvec^2); print(g)
  return(g) 
  }
}

g <- directL(fn = globalwrap,
             upper = ubounds, lower = lbounds,
             control = list(maxeval = 500, xtol_rel = 1e-2))

g
t0 = g$par

g <- gmm(g = gmmmomentmatcher, x = df.adj, t0 = t0, gradv = momentgrad,      
         type = "twoStep", optfct = "constrOptim",   
         ui = ui, ci = ci, 
         control = list(maxit = 100, reltol = 1e-2))
summary(g)

saveRDS(g, paste(getwd(), "/data/model_output/g.rds", sep = ""))
#g <- readRDS(paste(getwd(), "/data/model_output/g.rds", sep = ""))

theta = g$coefficients
gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]*theta[8]; hbar = theta[5]*(1-theta[8])*theta[9]/(1-theta[9])
lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
saveRDS(theta, paste(getwd(), "/data/model_output/theta.rds", sep = ""))

### Final Value function
statespace = create.statespace(ubm = c(5,5), theta, method = "equal")
v0 = firstguesser(statespace, theta)
V = VFI(v0, theta)
saveRDS(V, paste(getwd(), "/data/model_output/V.rds", sep = ""))

