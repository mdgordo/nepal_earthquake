### Run top of MSM.R

theta = theta0
gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
shockpts = 21

### statespace
xn = 45; hn = 45; yn = 10
ygrid = as.vector(quantile(df.hh$avg_inc, seq(0,1,length.out = yn), na.rm = TRUE))

## for figuring out max grid value
#(exp(max(great.expectations(max(ygrid), sigma, gqpts)))*R - exp(max(ygrid)+(max(ygrid)*sigma)^2/2))/(R-1)

statespace = create.statespace(ubm = 3, theta, method = "log")

guesser <- function(statespace, theta){
  beta = theta[2]; alpha = theta[8]; delta = theta[9]
  ifx0 <- mapply(function(x, h, B, cmin, hmin) if_else(x-B<cmin | x-B-cmin<hmin-h, max(0,hmin-h),
                                                       if_else((1-alpha)*(x-B-cmin)<hmin-h, hmin-h,
                                                               if_else((1-alpha)*(x-B)-h>0, (1-alpha)*(x-B)-h, 0))), 
                 x = statespace$x, h = statespace$h, B = statespace$B, cmin = statespace$cmin, hmin = statespace$hmin)
  cfx0 <- mapply(function(x, B, cmin, i) if_else(x-B-i<cmin, cmin, x-B-i), 
                 x = statespace$x, B = statespace$B, cmin = statespace$cmin, i = ifx0)
  ## periods w 0 investment until h hits hmin
  t0 = mapply(function(h, hmin) if_else(round(hmin-h)>=0, 0, ceiling((log(hmin) - log(h))/log(delta))), statespace$h, statespace$hmin)
  vfx0 = mapply(function(cmin, hmin, h, ifx, cfx, t) if_else(is.infinite(t), u(cfx, h, ifx, theta) + beta/(1-beta)*u(cmin, h, 0, theta),
                                                             if_else(t<=1, u(cfx, h, ifx, theta) + beta/(1-beta)*u(cmin, hmin, 0, theta),
                                                                     u(cfx, h, ifx, theta) + beta^(t-1)*u(cmin, h*delta^t, 0, theta) + beta^(t)/(1-beta)*u(cmin, hmin, 0, theta))),
                cmin = statespace$cmin, hmin = statespace$hmin, h = statespace$h, i = ifx0, c = cfx0, t = t0)
  return(data.frame("Tw" = vfx0,
                    "cfx" = cfx0,
                    "ifx" = ifx0))
}

### Initial guess 
v0 = guesser(statespace, theta)
w = cbind(statespace, v0)
w = filter(w, y==ygrid[1])

hhprob_funkfact <- function(x, h, y, B, hmin, Valfunc, theta, gqpts){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  hhprob = function(ci){
    c = ci[1]; i = ci[2]
    if (is.nan(c)|is.nan(i)) {return(Inf)} else {
      xtplus1 = R*(x - c - i) + exp(great.expectations(y, sigma, gqpts))
      defaults = which(xtplus1<B)
      xtplus1[defaults] = B
      htplus1 = rep(delta*(h+i), length(xtplus1))
      htplus1[defaults] = max(htplus1, hmin)
      payoff = u(c, h, i, theta)
      EV = sum(mapply(Valfunc, x = xtplus1, h = htplus1)*gqpts$w)/sqrt(pi)
      r = payoff + beta*EV
      return(-r)}
  }
  constraint <- function(ci){
    c = ci[1]
    i = ci[2]
    c1 = x-B-c-i
    return(c1)
  }
  return(list(hhprob, constraint))
}


Tw0 = w$Tw
gqpts = gaussHermiteData(shockpts)
Valfunc = interpolater.creater(ygrid[1], w, Tw0)
optwrap = function(rowidx){
  x = statespace$x[rowidx]; h = statespace$h[rowidx]; y = statespace$y[rowidx]
  B = statespace$B[rowidx]; cmin = statespace$cmin[rowidx]; hmin = statespace$hmin[rowidx]
  
  if (round(x - B,10) <= round(cmin,10) | round(x - B - cmin,10) <= round(hmin - h,10)) {
    vfx = u(cmin, h = max(h, hmin), i = 0, theta) + beta*Valfunc(x = B, h = max(delta*h, hmin))
    return(c(cmin, 0, vfx))
  } else {
    hhprob = hhprob_funkfact(x, h, y, B, hmin, Valfunc, theta, gqpts)
    lbi = max(round(hmin - h,10), 0.0)
    sps = incenter(cmin, lbi, x, B)
    sol = cobyla(x0 = sps, fn = hhprob[[1]], hin = hhprob[[2]],
                 lower = c(cmin,lbi), upper = c(x-B, x-B), 
                 control = list(ftol_rel = 1e-6, ftol_abs = 0, xtol_rel = 0, maxeval = 1000))
    sol2 = cobyla(x0 = cmin, fn = function(c){hhprob[[1]](c(c, x-B-c))}, 
                  lower = cmin, upper = x-B, 
                  control = list(ftol_rel = 1e-6, ftol_abs = 0, xtol_rel = 0, maxeval = 1000))
    if (sol$value > sol2$value) {
      par = c(sol2$par, x-B-sol2$par)
      sol = sol2
    } else {par = sol$par}
    vfxdefault = u(cmin, h = max(h, hmin), i = 0, theta) + beta*Valfunc(x = B, h = max(delta*h, hmin))
    if (vfxdefault < -sol$value) {
      return(c(par, -sol$value))
    } else {return(c(cmin, 0, vfxdefault))}
  }
}

Tw = mclapply(c(1:nrow(statespace)), optwrap, mc.cores = detectCores())
wnext = cbind(statespace, data.frame("Tw" = unlist(lapply(Tw, function(x) x[3])),
                                     "cfx" = unlist(lapply(Tw, function(x) x[1])),
                                     "ifx" = unlist(lapply(Tw, function(x) x[2]))))

wnext <- mutate(wnext, savings = x - cfx - ifx,
                dist2def = x - B - cfx - ifx)

ggplot(filter(wnext, h==0, x < 1e6, x > -1e5)) + geom_line(aes(x = x, y = cfx)) #+ geom_abline(intercept = 0, slope = 1)

### hunting for local minima
rowidx = 1281
optgrid = crossing(c = exp(seq(log(cmin), log(x-B-lbi), length.out = 100)), 
                   i = exp(seq(log(lbi), log(x-B-cmin), length.out = 100))) %>%
  filter(c+i<x-B)
optgrid$tw = unlist(lapply(c(1:nrow(optgrid)), function(i) hhprob[[1]](c(optgrid$c[i], optgrid$i[i]))))

optgrid = optgrid %>% group_by(c) %>%
  mutate(ht = log(i) - log(lag(i))) %>% ungroup() %>%
  group_by(i) %>% 
  mutate(wd = log(c) - log(lag(c)))
  
ggplot(filter(optgrid, c>1e5, c < 1e8, i < 1e9)) + geom_tile(aes(x = c, y = i, fill = log(tw), height = ht, width = wd)) + 
  scale_x_continuous(trans='log') + scale_y_continuous(trans='log') + 
  scale_fill_viridis_c()

ip = exp(seq(log(lbi), log(3e8), length.out = 100))
ggplot() + geom_line(aes(x = ip, 
                         y = unlist(lapply(ip, function(x) log(hhprob[[1]](c(1e6, x)))))))
vx = function(ci){
  c = ci[1]; i = ci[2]
  htplus1 = delta*(h+i)
  xtplus1 = R*(x - c - i) 
  mindraw = B - xtplus1
  pdef = plnorm(mindraw, meanlog = y, sdlog = sigma*y)
  draws = exp(great.expectations(y, sigma, gqpts))
  wts = gqpts$w[draws>mindraw]/sqrt(pi)
  wts[1] = 1 - pdef - sum(wts[2:length(wts)])
  draws = draws[draws>mindraw]
  EVp = sum(mapply(Valfunc, x = xtplus1 + draws, h = htplus1)*wts)
  EV = pdef*Valfunc(B, max(htplus1, hmin)) + EVp
  return(EV)
}
ggplot() + geom_line(aes(x = ip, 
                         y = unlist(lapply(ip, function(x) vx(c(1e6,x))))))

optgrid$xt1 = R*(x - (optgrid$c + optgrid$i)) + exp(ygrid[1]+ (sigma*ygrid[1])^2/2)
optgrid$ht1 = delta*(h + optgrid$i)
optgrid$vfx = mapply(Valfunc, x = optgrid$xt1, h = optgrid$ht1)

ggplot(filter(optgrid, c>1e5, c < 1e8, i < 1e9)) + geom_tile(aes(x = c, y = i, fill = log(-1*(vfx)), height = ht, width = wd)) + 
  scale_x_continuous(trans='log') + scale_y_continuous(trans='log') + 
  scale_fill_viridis_c()
### do these with quad points instead of just expectation? Look for local minima

### increasing number of h points helped somewhat, and changed upper bound of h to same as x in create.statespace
### redefined hhprob with built in constraint- didn't end up implementing this
### log optgrid
### incenter plus univariate approach is promising
### increasing maxeval helps somewhat with high values of x
### experimenting with global optimizers - way too slow
if (optimizer=="stogo"){
  sol = stogo(x0 = c(cmin,lbi), fn = hhprob, gr = hhgrad, 
              lower = c(cmin,lbi), upper = c(x-B, x-B))
} else if (optimizer=="isres"){
  sol = isres(x0 = c(cmin,lbi), fn = hhprob,
              lower = c(cmin,lbi), upper = c(x-B, x-B))
} else if (optimizer=="direct"){
  sol = direct(x0 = c(cmin,lbi), fn = hhprob, gr = hhgrad, 
               lower = c(cmin,lbi), upper = c(x-B, x-B))
} else if (optimizer=="crs"){
  sol = crs2lm(x0 = c(cmin,lbi), fn = hhprob, 
               lower = c(cmin,lbi), upper = c(x-B, x-B))
} else {
  sol = genoud(fn = hhprob, nvars = 2, Domains = matrix(c(cmin, lbi, rep(x-B,2)), nrow = 2),
               gr = hhgrad, boundary.enforcement = 2, print.level = 0, pop.size = 1000, wait.generations = 2)
}




