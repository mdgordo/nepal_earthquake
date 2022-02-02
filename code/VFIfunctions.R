### utility function
u <- function(c, h, i, theta) {
  alpha = theta[8]; gamma = theta[1]
  cd = c^alpha * (h + i)^(1-alpha)
  crra = cd^(1-gamma)/(1-gamma)
  return(crra)
}

### Creates points for gaussian quadrature
great.expectations <- function(y, sigma){
  mu = sqrt(2)*sigma*exp(y)*gqpts + exp(y)
}

create.statespace = function(i, ubm = c(20, 20), Blist, cminlist, hminlist){
  B = Blist[[i]]; cmin = cminlist[[i]]; hmin = hminlist[[i]]
  bc = cmin + B
  x = c(B, cmin+B/2, cmin+B, exp(seq(0, log(ubm[1]*exp(ygrid[i])+bc), length.out = xn-3))-bc)
  y = c(0, hmin/2, exp(seq(log(hmin), log(ubm[2]*exp(ygrid[i])), length.out = hn-2)))
  css = crossing(x, y)
  colnames(css) <- c("x", "h")
  css$y = ygrid[i]
  css$B = B
  css$cmin = cmin
  css$hmin = hmin
  return(css)
}

### generate initial guess for VFI
guesser <- function(statespace, theta){
  beta = theta[2]; delta = theta[9]
  ifx0 <- mapply(function(x, h, B, cmin, hmin) if_else(x-B-cmin<0 | x-B-cmin<hmin-h, max(0,hmin-h),
                                                       if_else(.1*(x-B)<hmin-h, hmin-h,
                                                               if_else(x-B-cmin>h/.9, .1*(x-B), 0))), 
                 x = statespace$x, h = statespace$h, B = statespace$B, cmin = statespace$cmin, hmin = statespace$hmin)
  cfx0 <- mapply(function(x, B, cmin, i) if_else(x-B-i<cmin, cmin, x-B-i), 
                 x = statespace$x, B = statespace$B, cmin = statespace$cmin, i = ifx0)
  ## periods w 0 investment until h hits hmin
  t0 = ceiling(pmax(mapply(function(h, hmin) (log(hmin) - log(h))/log(delta), statespace$h, statespace$hmin), 0))
  vfx0 = mapply(function(cmin, hmin, h, ifx, cfx, t) if_else(t<=1, u(cfx, h, ifx, theta) + beta/(1-beta)*u(cmin, hmin, 0, theta),
                                                             u(cfx, h, ifx, theta) + beta^(t-1)*u(cmin, h*delta^t, 0, theta) + beta^(t)/(1-beta)*u(cmin, hmin, 0, theta)),
                cmin = statespace$cmin, hmin = statespace$hmin, h = statespace$h, i = ifx0, c = cfx0, t = t0)
  return(data.frame("Tw" = vfx0,
                    "cfx" = cfx0,
                    "ifx" = ifx0))
}


### 2D interpolate
interpolater.creater = function(yi, statespace, Tw){
  statespace$Tw = Tw
  s <- filter(statespace, y == yi)
  fi = function(x,h){
    r = simplr(x0 = x, y0 = h, x = s$x, y = s$h, vfx = s$Tw, method = "simplical")
    return(r)
  }
  return(fi)
}

### 3D interpolate 
complexr <- function(x0, statespace, vfx){
  ## find closest points in y and interpolate
  yvals = sort(unique(statespace$y))
  yi = findInterval(x0[2], yvals, all.inside = TRUE)
  yp = c(yvals[yi], yvals[yi+1])
  xyhfilt = filter(statespace, y==yp[1])
  
  vy1 = interpolater.creater(yp[1], statespace, Tw = vfx)
  vy2 = interpolater.creater(yp[2], statespace, Tw = vfx)
  
  ## closest in x and h for y1
  xvals1 = sort(unique(xyhfilt$x))
  xi1 = findInterval(x0[1], xvals1, all.inside = TRUE)
  xp1 = c(xvals1[xi1], xvals1[xi1+1])
  
  hvals1 = sort(unique(xyhfilt$h))
  hi1 = findInterval(x0[3], hvals1, all.inside = TRUE)
  hp1 = c(hvals1[hi1], hvals1[hi1+1])
  
  ### values at vertices
  d = crossing(y = yp, x = xp1, h = hp1) %>%
    mutate(xidx = if_else(x==xp1[2],1,0),
           yidx = if_else(y==yp[2],1,0),
           hidx = if_else(h==hp1[2],1,0))
  d$vfx = mapply(vy1, d$x, d$h)
  d$vfx[d$y==yp[2]] = mapply(vy2, d$x[d$y==yp[2]], d$h[d$y==yp[2]])
  
  ### change of variables
  xcv = (x0[1]-xp1[1])/(xp1[2] - xp1[1])
  ycv = (x0[2]-yp[1])/(yp[2] - yp[1])
  hcv = (x0[3]-hp1[1])/(hp1[2] - hp1[1])
  xcv0 = c(xcv, ycv, hcv)
  
  ### algorithm
  p = order(xcv0)
  s = vector(mode = "list", length = 4)
  fs = vector(mode = "list", length = 4)
  s[[1]] = rep(1, 3)
  fs[[1]] = d$vfx[d$xidx==s[[1]][1] & d$yidx==s[[1]][2] & d$hidx==s[[1]][3]]
  for (i in c(2:4)) {
    svec = rep(0, 3)
    svec[p[i-1]] = 1
    s[[i]] = s[[i-1]] - svec
    f = d$vfx[d$xidx==s[[i]][1] & d$yidx==s[[i]][2] & d$hidx==s[[i]][3]]
    f2 = d$vfx[d$xidx==s[[i-1]][1] & d$yidx==s[[i-1]][2] & d$hidx==s[[i-1]][3]]
    fs[[i]] = fs[[i-1]] + (1 - xcv0[p[i-1]])*(f - f2)
  }
  return(fs[[4]])
}

### Function factory for generation household maximization problem

hhprob_funkfact <- function(x, h, y, B, hmin, Valfunc, theta){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  
  hhprob = function(ci){
    c = ci[1]
    i = ci[2]
    xtplus1 = R*(x - c - i) + great.expectations(y, sigma)
    defaults = which(xtplus1<B)
    xtplus1[defaults] = B
    htplus1 = rep(delta*(h+i), length(xtplus1))
    htplus1[defaults] = max(htplus1, hmin)
    payoff = u(c, h, i, theta)
    EV = sum(mapply(Valfunc, x = xtplus1, h = htplus1)*gqwts)
    r = payoff + beta*EV
    return(-r)
  }
  constraint <- function(ci){
    c = ci[1]
    i = ci[2]
    c1 = x - c - i - B
    c2 = h + i - hmin
    return(c(c1, c2))
  }
  return(list(hhprob, constraint))
}

### Bellman Operator

bellman <- function(w, statespace, theta){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  
  Vfxlist = lapply(ygrid, interpolater.creater, statespace, w$Tw)
  optwrap = function(rowidx){
    x = statespace$x[rowidx]; h = statespace$h[rowidx]; y = statespace$y[rowidx]
    B = statespace$B[rowidx]; cmin = statespace$cmin[rowidx]; hmin = statespace$hmin[rowidx]
    Valfunc = Vfxlist[[which(ygrid==y)]]
    if (x - B <= cmin | x - B - cmin <= hmin - h) {
      vfx = u(cmin, h = max(h, hmin), i = 0, theta) + beta*Valfunc(x = B, h = max(delta*h, hmin))
      return(c(cmin, 0, vfx))
    } else {
      hhprob = hhprob_funkfact(x, h, y, B, hmin, Valfunc, theta)
      sps = c(w$cfx[rowidx], w$ifx[rowidx])
      lbi = max(hmin - h, 0)
      sol = cobyla(x0 = sps, fn = hhprob[[1]], lower = c(cmin,lbi), upper = c(x-B, x-B), hin = hhprob[[2]],
                   control = list(xtol_rel = 1e-3))
      vfxdefault = u(cmin, h = max(h, hmin), i = 0, theta) + beta*Valfunc(x = B, h = max(delta*h, hmin))
      if (vfxdefault < -sol$value) {
        return(c(sol$par, -sol$value))
      } else {return(c(cmin, 0, vfxdefault))}
    }
  }
  Tw = mclapply(c(1:nrow(statespace)), optwrap,
                mc.cores = detectCores(), mc.preschedule = FALSE)
  wnext = data.frame("Tw" = unlist(lapply(Tw, function(x) x[3])),
                     "cfx" = unlist(lapply(Tw, function(x) x[1])),
                     "ifx" = unlist(lapply(Tw, function(x) x[2])))
  return(wnext)
}

### Howard Policy Iteration
howard <- function(wlast, statespace, theta){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  
  Vfxlist = lapply(ygrid, interpolater.creater, statespace, wlast$Tw)
  
  accelerator = function(rowidx){
    x = statespace$x[rowidx]; h = statespace$h[rowidx]; y = statespace$y[rowidx]
    B = statespace$B[rowidx]; cmin = statespace$cmin[rowidx]; hmin = statespace$hmin[rowidx];
    cfx = wlast$cfx[rowidx]; ifx = wlast$ifx[rowidx]
    Valfunc = Vfxlist[[which(ygrid==y)]]
    if (x - B <= cmin | x - B - cmin <= hmin - h) {
      v = u(cmin, hmin, 0, theta) + beta*Valfunc(B, max(delta*h, hmin))
    } else {
      hhprob = hhprob_funkfact(x, h, y, B, hmin, Valfunc, theta)
      v1 = -1*hhprob[[1]](c(cfx, ifx))
      v2 = u(cmin, hmin, 0, theta) + beta*Valfunc(B, max(delta*h, hmin))
      return(max(v1, v2))
    }
  }
  Twh = mclapply(1:nrow(statespace), accelerator, mc.cores = detectCores())
  return(unlist(Twh))
}


### VFI
VFI <- function(vinit, statespace, maxiter = 100, tol = 1e-4, howardk = 10, theta){
  beta = theta[2]
  w = vector(mode = "list", length = maxiter-1)
  w[[1]] = vinit
  d = 1; i = 2
  while (i<maxiter & d>tol) {
    wnext = bellman(w = w[[i-1]], statespace, theta)
    w[[i]] = wnext
    ## howard
    k = vector(mode = "list", length = howardk+1)
    k[[1]] = wnext
    for (j in c(1:howardk)) {
      wnext$Tw = howard(wnext, statespace, theta)
      k[[j+1]] = wnext
    }
    ### MQP error bounds
    b = beta/(1-beta)*mean(k[[howardk+1]]$Tw - k[[howardk]]$Tw)
    wnext$Tw = wnext$Tw + b
    ## check tol
    w[[i]] = wnext
    d = mean(((w[[i]]$cfx - w[[i-1]]$cfx)/w[[i]]$cfx)^2)
    i = i+1
    print(i)
  }
  return(w)
}

momentmatcher <- function(i, vfx, theta, data){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  df = data[i,]
  xyh = c(df$imputed_bufferstock, log(df$avg_inc), df$lag_h)
  e1  = complexr(x0 = xyh, statespace, vfx = vfx$cfx) - df$food_consumption
  e2  = complexr(x0 = xyh, statespace, vfx = vfx$ifx) - df$home_investment
  e3 = df$home_value - delta*(df$lag_h + df$home_investment)
  e4 = df$imputed_bufferstock - df$quake_aid - df$lag_y - R*(df$lag_x - df$lag_c - df$lag_i)
  e5 = df$var_inc - sigma*df$avg_inc
  e6 = e1*df$imputed_bufferstock
  e7 = e1*log(df$avg_inc)
  e8 = e1*df$lag_h
  e9 = e2*df$imputed_bufferstock
  e10 = e2*log(df$avg_inc)
  e11 = e2*df$lag_h
  return(c(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11)*df$wt_hh)
}

### check tolerance between iterations

tolchecker = function(w){
  for (i in c(2:length(w))) {
    d = mean((w[[i]]$Tw - w[[i-1]]$Tw)^2) ### total squared error
    cd = mean(((w[[i]]$cfx - w[[i-1]]$cfx)/w[[i]]$cfx)^2) ### mean squared error
    id = mean((w[[i]]$ifx - w[[i-1]]$ifx)^2)/sum(w[[i-1]]$ifx) ### defined as mse over total investment bc of zeros
    print(paste("vtol:", d, "ctol:", cd, "itol", id, sep = " "))
  }
}

### check to make sure grid is not binding

satisficer = function(w, yidx, sigma){
  xmax = max(statespace$x[statespace$y==yidx])
  hmax = max(statespace$h[statespace$y==yidx])
  ymax = max(great.expectations(yidx, sigma))
  xmax = 3000
  xnext = xmax - (xmax-ymax)/R
  hnext = hmax*(1/delta - 1)
  df.cstrt = cbind(w[[length(w)]], statespace) %>%
    filter(y==yidx & x==xmax) %>%
    mutate(spend = cfx + ifx,
           constrainedx = spend<xnext)
  df.cstrt2 = cbind(w[[length(w)]], statespace) %>%
    filter(y==yidx & h==hmax) %>%
    mutate(constrainedh = ifx>hnext)
  return(c(sum(df.cstrt$constrainedx), sum(df.cstrt2$constrainedh)))
}

### Counterfactuals

counterfactualizer = function(i, vfx, statespace, df.hh){
  x1 = as.double(df.hh[i,"x_noaid"]); x2 = as.double(df.hh[i,"x_aid"])
  y = log(as.double(df.hh[i,"avg_inc"])); h = as.double(df.hh[i,"lag_h"])
  c1 = complexr(x0 = c(x1, y, h), statespace, vfx = vfx$cfx)
  c2 = complexr(x0 = c(x2, y, h), statespace, vfx = vfx$cfx)
  i1 = complexr(x0 = c(x1, y, h), statespace, vfx = vfx$ifx)
  i2 = complexr(x0 = c(x2, y, h), statespace, vfx = vfx$ifx)
  b1 = as.double(df.hh[i,"x_noaid"]) - c1 - i1
  b2 = as.double(df.hh[i,"x_aid"]) - c2 - i2
  return(c(c1, c2, i1, i2, b1, b2))
}

wtpeliciter = function(x, y, h, aidamt, vfx, statespace){
  c1 = complexr(x0 = c(x, y, h), statespace, vfx = vfx$Tw)
  x2 = x + aidamt
  groot = function(yp) {
    r = c1 - complexr(x0 = c(x2, yp, h), statespace, vfx = vfx$Tw)
    return(r)
  }
  if (groot(min(ygrid))>0) {yp = min(ygrid)} else {
    r = uniroot(groot, interval = c(min(ygrid), y))
    yp = r$root
  }
  wtp = (exp(y) - exp(yp))/(1-beta)
  return(wtp)
}


### Functions for plotting

iterplotter = function(w, hidx=NULL, xidx=NULL, yidx, ifilt = 99){
  ## plot iterations for given y and h
  if (is.null(xidx)) {
    xs = lapply(w, function(x) filter(cbind(x, statespace), y==yidx, h==hidx))
    xs <- data.frame(do.call(rbind, xs)) 
    xs$iteration <- rep(c(1:length(w)), each = xn)
    ggplot(filter(xs, iteration < ifilt)) +
      geom_line(aes(x = x, y = Tw, color = iteration, group = iteration)) + 
      scale_color_viridis_c()
  } else {
    xs = lapply(w, function(x) filter(cbind(x, statespace), y==yidx, x==xidx))
    xs <- data.frame(do.call(rbind, xs)) 
    xs$iteration <- rep(c(1:length(w)), each = hn)
    ggplot(filter(xs, iteration < ifilt)) +
      geom_line(aes(x = h, y = Tw, color = iteration, group = iteration)) + 
      scale_color_viridis_c()
  }
}

threedplotter = function(w, yidx, fillvar, iter = length(w), d3 = TRUE, aidamt = NULL){
  if (fillvar=="Tw") {
    Vfx = interpolater.creater(ygrid[yidx], statespace, w[[iter]]$Tw)
  } else if (fillvar %in% c("cfx", "aidcfx")) {
    Vfx = interpolater.creater(ygrid[yidx], statespace, w[[iter]]$cfx)
  } else {Vfx = interpolater.creater(ygrid[yidx], statespace, w[[iter]]$ifx)}
  
  df.test = crossing(x = seq(Blist[[yidx]], 20*exp(ygrid[yidx]), 200), h = seq(0, 20*exp(ygrid[yidx]), 200))
  df.test$proj = mapply(Vfx, x = df.test$x, h = df.test$h)
  
  if (fillvar %in% c("aidcfx", "aidifx")) {
    df.test$proj2 = mapply(Vfx, df.test$x + aidamt, df.test$h)
    df.test$efx = df.test$proj2 - df.test$proj
  }
  
  if (fillvar=="Tw") {
    threed = ggplot(df.test) + geom_tile(aes(x = x, y = h, fill = -1*log(-proj))) + scale_fill_viridis_c()
  } else if (fillvar=="cfx") {
    threed = ggplot(df.test) + geom_tile(aes(x = x, y = h, fill = proj)) + scale_fill_viridis_c()
  } else if (fillvar=="ifx") {
    threed = ggplot(df.test) + geom_tile(aes(x = x, y = h, fill = proj)) + scale_fill_viridis_c()
  } else {threed = ggplot(df.test) + geom_tile(aes(x = x, y = h, fill = efx)) + scale_fill_viridis_c()}
  
  if (d3) p = plot_gg(threed, multicore=TRUE) else p = threed
  p
}