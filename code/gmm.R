### packages and data
library(tidyverse)
library(parallel)
library(nloptr)
library(fastGHQuad)
library(flatlandr)
library(gmm)
library(Rearrangement)
library(mgcv)
library(neuralnet)
options('nloptr.show.inequality.warning'=FALSE)
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))

iterations <- readRDS(paste(getwd(), "/data/model_output/iterations.rds", sep = "")) %>%
  filter(complete.cases(.))
iterations <- mutate(iterations, obj = m_1^2 + m_2^2 + m_3^2 + m_4^2 + m_5^2 + m_6^2 + 
                       m_7^2 + m_8^2 + m_9^2 + m_10^2 + m_11^2)

#ggplot(iterations) + geom_point(aes(x = cbar, y = alpha, color = obj)) + scale_color_viridis_c(trans = "log")

modfunc <- function(m){
  X = iterations[,c(1:9)]
  X = cbind(X, iterations[,m]); colnames(X) = c(colnames(X)[1:9], m)
  f = as.formula(paste(m, 
                       "~ .^2 + I(gamma^2) + I(beta^2) + I(R^2) + I(cbar^2) + I(hbar^2) + I(lambda^2) + I(sigma^2) + I(alpha^2) + I(delta^2)", 
                       sep = ""))
  mod = lm(f, data = X)
  return(mod)
}

#momentvars <- paste("m", seq(1,11,1), sep = "_")
#modlist <- mclapply(momentvars, modfunc, mc.cores = max(11,detectCores()-2))
mod <- modfunc("obj")

smoothedobj <- function(t){
  if (t[2]*t[3] > 1 | t[4] > t[8] | t[5] > 1 - t[8]){
    r = max(iterations$obj)} else {
    #r = rep(max(iterations$obj),11)} else { 
      newdata = data.frame("gamma" = t[1], "beta" = t[2], "R" = t[3],
                           "cbar" = t[4], "hbar" = t[5], "lambda" = t[6],
                           "sigma" = t[7], "alpha" = t[8], "delta" = t[9])
     #momlist = lapply(modlist, function(mod) predict(mod, newdata))
     #r = unlist(momlist)*colrange[c(10:20)] + colmins[c(10:20)]
    r = predict(mod, newdata)
    }
  return(r)
  #return(sum(r^2))
}

sol <- cobyla(x0 = c(1.2, .926, 1.054, .33, .19, .08, .42, .63, .93), 
              fn = smoothedobj, #function(t) sum(smoothedobj(t)^2), 
              lower = c(1.01, .6, 1, .01, .01, .01, .01, .1, .5), 
              upper = c(10, .99, 1.2, .9, .9, 10, 1.5, .99, 1), control = list(maxeval = 3000))
theta = sol$par
print(theta)

### load data
df.adj <- readRDS(paste(getwd(), "/data/model_output/df.adj.rds", sep = ""))
df.adj <- df.adj[complete.cases(df.adj),]

### Grid size
xn = 50; hn = 50

### GMM using theta as starting point
### penalty value for constraints
gpenalty = apply(iterations[,c(10:20)], 2, max)

gmmwrap <- function(t, x){
  if (t[2]*t[3] > 1 | t[4] > t[8] | t[5] > 1 - t[8]){
    return(t(matrix(rep(gpenalty, nrow(x)), ncol = nrow(x))))
  } else{
    momentmat = gmmmomentmatcher(theta = t, df = x)
    print(sum(colSums(momentmat)^2))
    return(momentmat)
  }
}

g <- gmm(g = gmmwrap, x = df.adj, t0 = theta,      
         type = "twoStep", optfct = "nlminb",  
         lower = c(1.01, .6, 1, .01, .01, .01, .01, .1, .5), 
         upper = c(10, .99, 1.2, .9, .9, 10, 1.5, .99, 1),
         control = list(x.tol = 1e-2, rel.tol = 1e-2))
summary(g)

theta = g$coefficients
gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
saveRDS(theta, paste(getwd(), "/data/model_output/theta.rds", sep = ""))

### statespace
xn = 60; hn = 60

### Final Value function
statespace = create.statespace(ubm = c(5,5), theta, method = "log")
v0 = cbind(statespace, data.frame("Tw" = rep(0, nrow(statespace)),
                                  "cfx" = rep(0, nrow(statespace)),
                                  "ifx" = rep(0, nrow(statespace)),
                                  "def" = rep(0, nrow(statespace))))
V = VFI(v0, theta)
saveRDS(V, paste(getwd(), "/data/model_output/V.rds", sep = ""))
