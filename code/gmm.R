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

lbounds = c(1, .85, 1, .001, .001, 0, .2, 0, .5, .9)
ranges = c(4, .14, .09, .6, 2, .6, 1, .6, .49, .099)

ubounds <- lbounds + ranges
ui = rbind(diag(10),
           -1*diag(10))
ci = c(lbounds, -1*ubounds)

print(lbounds); print(ubounds)

### Grid size
xn = 40; hn = 40

### load data
df.adj <- readRDS(paste(getwd(), "/data/model_output/df.adj.rds", sep = ""))

guesses <- data.frame(
  matrix(ncol = length(lbounds) + 1, nrow = 0)
)

### GMM 
## first coarse global
## Identity weight matrix - number of moments = 17
W = diag(17)

globalwrap = function(t, W){
  if (t[2]*t[3]>1) {
    r <- 99*t[2]*t[3]
  } else {
      momentmat = gmmmomentmatcher(t, df.adj)
      gvec = colMeans(momentmat); 
      r = t(gvec) %*% W %*% gvec
      print(r)
  }
  .GlobalEnv$guesses[nrow(.GlobalEnv$guesses) + 1, ] <- c(t, r)
  return(as.numeric(r)) 
}

g <- direct(fn = function(t) globalwrap(t, W = W),
             upper = ubounds, lower = lbounds,
             control = list(maxeval = 1000, xtol_rel = 1e-2))

g
t0 = g$par

saveRDS(guesses, paste(getwd(), "/data/model_output/guesses.rds", sep = ""))

### Update weight matrix
m1 <- gmmmomentmatcher(t0, df.adj)
W = solve(cov(m1))

### Update bounds
new_lb <- t0 - ranges/2
new_ub <- t0 + ranges/2
new_lb <- pmax(new_lb, lbounds)
new_ub <- pmin(new_ub, ubounds)
print(new_lb); print(new_ub)

### second iteration
g <- direct(fn = function(t) globalwrap(t, W = W),
            upper = new_ub, lower = new_lb,
            control = list(maxeval = 500, xtol_rel = 1e-2))
g
t0 = g$par

saveRDS(guesses, paste(getwd(), "/data/model_output/guesses.rds", sep = ""))

### Update weight matrix
m1 <- gmmmomentmatcher(t0, df.adj)
W = solve(cov(m1))

### final round
g <- gmm(g = gmmmomentmatcher, x = df.adj, t0 = t0, gradv = momentgrad,      
         weightsMatrix = W, optfct = "constrOptim",   
         ui = ui, ci = ci, 
         control = list(maxit = 100, reltol = 1e-2))
summary(g)

saveRDS(g, paste(getwd(), "/data/model_output/g.rds", sep = ""))
saveRDS(guesses, paste(getwd(), "/data/model_output/guesses.rds", sep = ""))
#g <- readRDS(paste(getwd(), "/data/model_output/g.rds", sep = ""))

theta = g$coefficients
gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
lambda = theta[6]; sigma = theta[7]; sigmame = theta[8]; alpha = theta[9]; delta = theta[10]
saveRDS(theta, paste(getwd(), "/data/model_output/theta.rds", sep = ""))

### Final Value function - extend grid (same density)
xn = hn = 80
statespace = create.statespace(ubm = c(10,10), theta, method = "equal")
v0 = firstguesser(statespace, theta)
V = VFI(v0, theta)
saveRDS(V, paste(getwd(), "/data/model_output/V.rds", sep = ""))

