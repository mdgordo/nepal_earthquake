### packages and data
library(tidyverse)
library(nloptr)
library(fastGHQuad)
#devtools::install_github("mdgordo/flatlandr",ref="master")
library(flatlandr)
library(gmm)
library(parallel)
options('nloptr.show.inequality.warning'=FALSE)
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))

# Bounds chosen based on preliminary tests
lbounds = c(4, .95, 1, .001, .001, 0, .75, 0, .6, .9)
ranges = c(8, .049, .03, .25, .5, 1.2, .75, 1.5, .35, .1)
ubounds <- lbounds + ranges
ui = rbind(diag(10),
           -1*diag(10))
ci = c(lbounds, -1*ubounds)

### Grid size
xn = 40; hn = 40

### load data
df.adj <- readRDS(paste(getwd(), "/data/model_output/df.adj.rds", sep = ""))

### estimate initial weight matrix based on midpoints - helps with conditioning of objective function
m1 <- gmmmomentmatcher(theta = (lbounds + ranges/2), df.adj)
W = solve(cov(m1))


### GMM 
## first coarse global
globalwrap = function(t){
  if (t[2]*t[3]>1) {
    return(99*t[2]*t[3])
  } else {
  momentmat = gmmmomentmatcher(t, df.adj)
  gvec = colMeans(momentmat); 
  r = t(gvec) %*% W %*% gvec
  print(r)
  return(as.numeric(r)) 
  }
}

g <- directL(fn = globalwrap,
             upper = ubounds, lower = lbounds,
             control = list(maxeval = 500, xtol_rel = 1e-2))

g
t0 = g$par

### Update weight matrix
m1 <- gmmmomentmatcher(t0, df.adj)
W = solve(cov(m1))

g <- gmm(g = gmmmomentmatcher, x = df.adj, t0 = t0, gradv = momentgrad,      
         weightsMatrix = W, optfct = "constrOptim",   
         ui = ui, ci = ci, 
         control = list(maxit = 100, reltol = 1e-2))
summary(g)

saveRDS(g, paste(getwd(), "/data/model_output/g.rds", sep = ""))
#g <- readRDS(paste(getwd(), "/data/model_output/g.rds", sep = ""))

theta = g$coefficients
gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]*theta[8]; hbar = theta[5]*(1-theta[8])*theta[9]/(1-theta[9])
lambda = theta[6]; sigma = theta[7]; sigmame = theta[8]; alpha = theta[9]; delta = theta[10]
saveRDS(theta, paste(getwd(), "/data/model_output/theta.rds", sep = ""))

### Final Value function - extend grid (same density)
xn = hn = 80
statespace = create.statespace(ubm = c(10,10), theta, method = "equal")
v0 = firstguesser(statespace, theta)
V = VFI(v0, theta)
saveRDS(V, paste(getwd(), "/data/model_output/V.rds", sep = ""))

