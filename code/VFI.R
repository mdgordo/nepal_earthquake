library(tidyverse)
library(parallel)
library(nloptr)
library(flatlandr)
options('nloptr.show.inequality.warning'=FALSE)
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))

### parameters
gamma <- 2.988
beta <- .9
R <- 1.012
cbar <- .733
hbar <- .5
lambda <- 1.55
sigma <- .2
alpha <- .9
delta <- .95

### utility function
u <- function(c, h, i) {
  cd = c^alpha * (h + i)^(1-alpha)
  crra = cd^(1-gamma)/(1-gamma)
  return(crra)
}

### statespace
xn = hn = yn = 40
ygrid = seq(8,13, length.out = yn)
ygrid = chebnodes(m = yn, lb = 8, ub = 13)
Blist = lapply(ygrid, function(y) -lambda*exp(y))
cminlist = lapply(ygrid, function(y) cbar*exp(y))
hminlist = lapply(ygrid, function(y) hbar*exp(y))

statespace = lapply(c(1:length(ygrid)), create.statespace, cheb = FALSE)
statespace = do.call(rbind, statespace)

### Shocks - Gaussian Quadrature points
gqpts = c(-2.651961, -1.673552, -.8162879, 0, .8162879, 1.673552, 2.651961)
gqwts = c(.0009717812, .0545155828, .4256072526, .8102646176, .4256072526, .0545155828, .0009717812)/sqrt(pi)

### Policy Function initial guess 
ifx0 <- mapply(function(x, h, B, cmin, hmin) if_else(x-B-cmin<0 | x-B-cmin<hmin-h, max(0,hmin-h),
                                                     if_else(.1*(x-B)<hmin-h, hmin-h,
                                                            if_else(x-B-cmin>h/.9, .1*(x-B), 0))), 
                        x = statespace$x, h = statespace$h, B = statespace$B, cmin = statespace$cmin, hmin = statespace$hmin)
cfx0 <- mapply(function(x, B, cmin, i) if_else(x-B-i<cmin, cmin, x-B-i), 
               x = statespace$x, B = statespace$B, cmin = statespace$cmin, i = ifx0)
t0 = ceiling(pmax(mapply(function(h, hmin) (log(hmin) - log(h))/log(delta), statespace$h, statespace$hmin), 0))
vfx0 = mapply(function(cmin, hmin, h, ifx, cfx, t) if_else(t<=1, u(cfx, h, ifx) + beta/(1-beta)*u(cmin, hmin, 0),
                                                       u(cfx, h, ifx) + beta^(t-1)*u(cmin, h*delta^t, 0) + beta^(t)/(1-beta)*u(cmin, hmin, 0)),
              cmin = statespace$cmin, hmin = statespace$hmin, h = statespace$h, i = ifx0, c = cfx0, t = t0)
#v0 = cbind(statespace, data.frame("ifx" = ifx0, "cfx" = cfx0, "vfx" = vfx0))

### Chebyshev Approximation
Vfxlist = lapply(ygrid, interpolater.creater, statespace, vfx0)
df.test = crossing(x = seq(-4500,44983,200), h = seq(140,90000,200))
df.test$proj = mapply(function(x,h) Vfxlist[[1]](x,h), x = df.test$x, h = df.test$h)
ggplot(df.test) + geom_tile(aes(x = x, y = h, fill = log(-proj))) + scale_fill_viridis_c()


### Neural net
library(neuralnet)
nn=neuralnet(vfx~x+h, data=v0[v0$y==ygrid[1],], hidden=c(20,30,20), act.fct = "logistic",
             linear.output = TRUE)
df.test = crossing(x = seq(-4500,44983,200), h = seq(140,90000,200))
df.test$proj = compute(nn, df.test)$net.result
ggplot(df.test) + geom_tile(aes(x = x, y = h, fill = log(-proj))) + scale_fill_viridis_c()


hhprob <- function(ci){
  c = ci[1]
  i = ci[2]
  xtplus1 = R*(x - c - i) + yvec
  great.expectations = unlist(mapply(Valfunc, x = xtplus1, h = hp))*gqwts
  great.expectations[xtplus1<B] = Valfunc(B, delta*h)*gqwts[xtplus1<B]
  r = u(c, i, h) + beta*
    return(-r)
}

