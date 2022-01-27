create.statespace = function(i, ubm = c(15, 30), cheb = "TRUE"){
  B = Blist[[i]]
  cmin = cminlist[[i]]
  hmin = hminlist[[i]]
  if (cheb) {
    css = chebstatespace(m = xn, lb = c(B, 0), ub = c(ubm[1]*exp(ygrid[i]), ubm[2]*exp(ygrid[i])))
  } else {
    x = seq(B, ubm[1]*exp(ygrid[i]), length.out = xn)
    y = seq(0, ubm[2]*exp(ygrid[i]), length.out = yn)
    css = crossing(x, y)
  }
  colnames(css) <- c("x", "h")
  css$y = ygrid[i]
  css$B = B
  css$cmin = cmin
  css$hmin = hmin
  return(css)
}

interpolater.creater = function(yi, statespace, Tw){
  statespace$Tw = Tw
  s <- filter(statespace, y == yi)
  Twm = matrix(s$Tw, ncol = xn)
  cmat = chebcoefs(Twm, degree = xn - 1)
  i = which(ygrid==yi)
  xmin = Blist[[i]]
  xmax = 15*exp(ygrid[i])
  hmax = 30*exp(ygrid[i])
  i = function(x,h){
    r = chebpred(x, h, coefmat = cmat, lb = c(xmin,0), ub = c(xmax, hmax))
    return(r)
  }
  return(i)
}


dfplot <- function(V){
  V <- as.data.frame(V)
  V$grid <- xgrid
  colnames(V) <- c(1:(dim(V)[2] - 1), "grid")
  V <- pivot_longer(V, cols = -grid, names_to = "iteration", values_to = "fx")
  V$iteration <- as.integer(V$iteration)
  return(V)
}

bufferstock <- function(hhid, Vlist, df){
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
  ### pre and post aid buffer stock - aid = 300000?
  bsx_pre = r - df$quake_aid[df$hhid==hhid]
  bsx_post = r + 100000
  bsx_actual = r
  
  ### Calculate counterfactual consumption
  cdist_pre = cfx(bsx_pre)
  cdist_post = cfx(bsx_post)
  return(c(cdist_pre, cdist_post, bsx_pre, bsx_post, bsx_actual))
}

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

wtpcurve <- function(x, mu, amt){
  i = which(unlist(lapply(Vlist, function(x) x$mu))==mu)
  Vfx = Vlist[[i]]$Vfx
  V = Vfx(x)
  val4aid = unlist(lapply(Vlist, function(f) f$Vfx(x+amt)))
  vi = which.min(abs(V - val4aid))
  mu0 = exp(mu + (sigma/mu)^2/2)
  muprime = exp(mus[vi] + (sigma/mu)^2/2)
  return((mu0 - muprime)/(1 - beta))
}



