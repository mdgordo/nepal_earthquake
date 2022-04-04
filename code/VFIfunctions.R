incenter <- function(cmin, lbi, x, B) {
  return(c((cmin*2 + x-B-lbi)/3, (lbi*2 + x-B-cmin)/3))
}

### utility function
u <- function(c, h, i, theta) {
  alpha = theta[8]; gamma = theta[1]
  cd = c^alpha * (h + i)^(1-alpha)
  crra = cd^(1-gamma)/(1-gamma)
  return(crra)
}

dudc <- function(c, h, i, theta){
  alpha = theta[8]; gamma = theta[1]
  d = alpha*c^(alpha-1)*(c^alpha * (h+i)^(1-alpha))^(-gamma)/(h+i)^(alpha-1)
}

dudi <- function(c, h, i, theta){
  alpha = theta[8]; gamma = theta[1]
  d = (1-alpha)*c^(alpha)*(c^alpha * (h+i)^(1-alpha))^(-gamma)/(h+i)^(alpha)
}

simplrdifferentiator <- function(x0, y0, statespace, deriv){
  x = statespace$x
  y = statespace$h
  vfx = statespace$Tw
  ## find closest points
  xvals = sort(unique(x))
  xi = findInterval(x0, xvals, all.inside = TRUE)
  x1 = xvals[xi]
  x2 = xvals[xi+1]
  yvals = sort(unique(y))
  yi = findInterval(y0, yvals, all.inside = TRUE)
  y1 = yvals[yi]
  y2 = yvals[yi+1]
  if (x0>max(x)) x0 = max(x) else if (x0<min(x)) x0 = min(x)
  if (y0>max(y)) y0 = max(y) else if (y0<min(y)) y0 = min(y)
  
  ## function vals at corners
  v1 = vfx[which(x==x1 & y==y1)]
  v2 = vfx[which(x==x1 & y==y2)]
  v3 = vfx[which(x==x2 & y==y1)]
  v4 = vfx[which(x==x2 & y==y2)]
  
  ## change of variables
  xp = (x0-x1)/(x2-x1)
  yp = (y0-y1)/(y2-y1)
  ## derivatives
  if (xp+yp <= 1) {
    if (deriv == "dx"){
      dv = v3/(x2-x1) - v1/(x2-x1) 
    } else {
      dv = v2/(y2-y1) - v1/(y2-y1)
    }
  } else {
    if (deriv == "dx"){
      dv = v4/(x2-x1) - v2/(x2-x1)
    } else {
      dv = v4/(y2-y1) - v3/(y2-y1)
    }
  }
  return(dv)
}

### Creates points for gaussian quadrature
great.expectations <- function(y, sigma, gqpts){
  mu = sqrt(2)*sigma*y*gqpts$x + y
  return(mu)
}

create.statespace = function(ubm = 1.4, theta, method = "uneven"){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  
  defaultpts = (cbar-lambda)*exp(ygrid+sigma^2/2)
  
  if (method=="chebyshev") {
    dlist = list()
    for (i in c(1:length(defaultpts))) {
      cnodesx = chebnodes(xn, lb = defaultpts[i], ub = exp(ubm*ygrid[i]))
      x = c(-lambda*exp(ygrid[i]+sigma^2/2), cnodesx)
      h = chebnodes(hn, lb = 0, ub = exp(ubm*ygrid[i]))
      dlist[[i]] = crossing(x, h)
      dlist[[i]]$y = ygrid[i]
    }
    css = do.call(rbind, dlist)
  } else if (method=="log") {
    dlist = list()
    for (i in c(1:length(defaultpts))) {
      h = exp(seq(0, log(exp(ubm*ygrid[i])), length.out = hn))
      x = h + defaultpts[i]
      dlist[[i]] = crossing(x, h)
      dlist[[i]]$y = ygrid[i]
    }
    css = do.call(rbind, dlist)
  } else {
    debtpts = seq(max(defaultpts), exp(median(ygrid)+sigma^2/2), length.out = 30)
    pospts = exp(seq(median(ygrid+sigma^2/2), ubm*max(ygrid + sigma^2/2), length.out = xn-yn-30))
    x = sort(c(defaultpts, debtpts, pospts))
    h = c(0, exp(seq(log(hbar*exp(ygrid[1]+sigma^2/2)), log(max(x)), length.out = hn-1)))
    css = crossing(y = ygrid, h, x) 
  }
  css = css %>%
    mutate(B = -lambda*exp(y+sigma^2/2),
           cmin = cbar*exp(y+sigma^2/2),
           hmin = hbar*exp(y+sigma^2/2))
  return(css)
}

### 2D interpolate
interpolater.creater = function(yi, statespace, Tw){
  statespace$Tw = Tw
  s <- filter(statespace, y == yi)
  fi = function(x,h){
    r = simplr(x0 = x, y0 = h, x = s$x, y = s$h, vfx = s$Tw, method = "simplical", extrapolation.warning = FALSE)
    return(r)
  }
  #maxdefpt = max(statespace$B + statespace$cmin)
  #cs = filter(s, x>maxdefpt)
  #tmat = matrix(cs$Tw, ncol = xn)
  #cmat = chebcoefs(tmat, degree = xn-1)
  #fi = function(x, h){
  #  if(x<=maxdefpt){
  #    r = simplr(x0 = x, y0 = h, x = s$x, y = s$h, vfx = s$Tw, method = "simplical", extrapolation.warning = FALSE)
  #  } else {
  #    r = chebpred(x,h, cmat, lb = c(maxdefpt, 0), ub = c(ubm*max(ygrid),max(df.hh$home_value, na.rm = TRUE)))
  #  }
  #}
  return(fi)
}

### 3D interpolate 
complexr <- function(x0, statespace, vfx){
  ## find closest points in y and interpolate
  yvals = sort(unique(statespace$y))
  yi = findInterval(x0[2], yvals, all.inside = TRUE)
  if (x0[2]>max(yvals)) x0[2] <- max(yvals)
  if (x0[2]<min(yvals)) x0[2] <- min(yvals)
  yp = c(yvals[yi], yvals[yi+1])
  xyhfilt = filter(statespace, y==yp[1])
  
  vy1 = interpolater.creater(yp[1], statespace, Tw = vfx)
  vy2 = interpolater.creater(yp[2], statespace, Tw = vfx)
  
  ## closest in x and h for y1
  xvals1 = sort(unique(xyhfilt$x))
  xi1 = findInterval(x0[1], xvals1, all.inside = TRUE)
  if (x0[1]>max(xvals1)) x0[1] <- max(xvals1)
  if (x0[1]<min(xvals1)) x0[1] <- min(xvals1)
  xp1 = c(xvals1[xi1], xvals1[xi1+1])
  
  hvals1 = sort(unique(xyhfilt$h))
  hi1 = findInterval(x0[3], hvals1, all.inside = TRUE)
  if (x0[3]>max(hvals1)) x0[3] <- max(hvals1)
  if (x0[3]<min(hvals1)) x0[3] <- min(hvals1)
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

monotonizer <- function(w, yi){
  w = filter(w, y==yi)
  xvals = sort(unique(w$x))
  vrC = w %>%
    pivot_wider(id_cols = x, names_from = h, values_from = cfx) %>%
    select(!x) %>% as.matrix()
  finalVrC = rearrangement(x = list(sort(unique(w$x)), sort(unique(w$h))), 
                           y = vrC, avg = TRUE) %>% as_tibble()
  finalVrC$x = xvals
  finalVrC <- pivot_longer(finalVrC, !x, names_to = "h", values_to = "cfx")
  finalVrC$h = as.double(finalVrC$h)
  vrI = w %>%
    pivot_wider(id_cols = x, names_from = h, values_from = ifx) %>%
    select(!x) %>% as.matrix()
  finalVrI = rearrangement(x = list(sort(unique(w$x)), -1*sort(unique(w$h))), 
                           y = vrI[,rev(c(1:ncol(vrI)))], avg = TRUE) %>% as_tibble()
  finalVrI$x = xvals
  finalVrI <- pivot_longer(finalVrI, !x, names_to = "h", values_to = "ifx")
  finalVrI$h = as.double(finalVrI$h)
  vrT = w %>%
    pivot_wider(id_cols = x, names_from = h, values_from = Tw) %>%
    select(!x) %>% as.matrix()
  finalVrT = rearrangement(x = list(sort(unique(w$x)), sort(unique(w$h))), 
                           y = vrT, avg = TRUE) %>% as_tibble()
  finalVrT$x = xvals
  finalVrT <- pivot_longer(finalVrT, !x, names_to = "h", values_to = "Tw")
  finalVrT$h = as.double(finalVrT$h)
  finalVr = merge(finalVrT, finalVrC, by = c("x", "h"))
  finalVr = merge(finalVr, finalVrI, by = c("x", "h"))
  finalVr$y = w$y; finalVr$B = w$B; finalVr$cmin = w$cmin; finalVr$hmin = w$hmin
  finalVr = finalVr[order(finalVr$y,finalVr$x,finalVr$h),]
  finalVr = finalVr[,c("x", "h", "y", "B", "cmin", "hmin", "Tw", "cfx", "ifx")]
  return(finalVr)
}

### generate initial guess for VFI
guesser <- function(statespace, theta){
  beta = theta[2]; alpha = theta[8]; delta = theta[9]
  ifx0 <- mapply(function(x, h, B, cmin, hmin) if_else(x-B<cmin | x-B-cmin<hmin-h, max(0,hmin-h),
                                                       (1-alpha)*(x-B-cmin-max(0,hmin-h))), 
                 x = statespace$x, h = statespace$h, B = statespace$B, cmin = statespace$cmin, hmin = statespace$hmin)
  cfx0 <- mapply(function(x, B, cmin, i) if_else(x-B-i<cmin, cmin, x-B-i), 
                 x = statespace$x, B = statespace$B, cmin = statespace$cmin, i = ifx0)
  ## periods w 0 investment until h hits hmin
  t0 = mapply(function(h, hmin) if_else(round(hmin-h)>=0, 0, ceiling((log(hmin) - log(h))/log(delta))), statespace$h, statespace$hmin)
  vfx0 = mapply(function(cmin, hmin, h, ifx, cfx, t) if_else(is.infinite(t), u(cfx, h, ifx, theta) + beta/(1-beta)*u(cmin, h, 0, theta),
                                                             if_else(t<=1, u(cfx, h, ifx, theta) + beta/(1-beta)*u(cmin, hmin, 0, theta),
                                                                     u(cfx, h, ifx, theta) + beta^(round(t/2,0))*u(cmin, mean(c(h,hmin)), 0, theta) + beta^(t)/(1-beta)*u(cmin, hmin, 0, theta))),
                cmin = statespace$cmin, hmin = statespace$hmin, h = statespace$h, i = ifx0, c = cfx0, t = t0)
  return(data.frame("Tw" = vfx0,
                    "cfx" = cfx0,
                    "ifx" = ifx0))
}

### Function factory for generation household maximization problem

hhprob_funkfact <- function(x, h, y, B, hmin, Valfunc, theta, gqpts){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  
  hhprob = function(ci){
    c = ci[1]; i = ci[2]
    if (is.nan(c)|is.na(c)) {return(Inf)} else {
      payoff = u(c, h, i, theta)
      htplus1 = delta*(h+i)
      xtplus1 = R*(x - c - i) 
      mindraw = B - xtplus1
      pdef = plnorm(mindraw, meanlog = y, sdlog = sigma*y)
      draws = exp(great.expectations(y, sigma, gqpts))
      if (any(draws>mindraw)){
        wts = gqpts$w[draws>mindraw]/sqrt(pi)
        wts[1] = 1 - pdef - sum(wts[2:length(wts)])
        draws = draws[draws>mindraw]
        EVp = sum(mapply(Valfunc, x = xtplus1 + draws, h = htplus1)*wts)
        EV = pdef*Valfunc(B, max(htplus1, hmin)) + EVp
      } else {
        EV = pdef*Valfunc(B, max(htplus1, hmin)) + (1-pdef)*Valfunc(xtplus1+mindraw+1, htplus1)
      }
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

hhgradfact <- function(x, h, yi, B, hmin, theta, s, Tw, gqpts){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  s$Tw = Tw
  s = filter(s, y==yi)
  hhgrad <- function(ci){
    c = ci[1]; i = ci[2]
    xtplus1 = R*(x - c - i) + exp(great.expectations(yi, sigma, gqpts))
    defaults = which(xtplus1<B)
    xtplus1[defaults] = B
    htplus1 = rep(delta*(h+i), length(xtplus1))
    htplus1[defaults] = max(htplus1, hmin)
    Edvdc = sum(mapply(simplrdifferentiator, x0 = xtplus1, y0 = htplus1, statespace = list(s), deriv = list("dx"))*gqpts$w)/sqrt(pi)
    Edvdi = sum(mapply(simplrdifferentiator, x0 = xtplus1, y0 = htplus1, statespace = list(s), deriv = list("dy"))*gqpts$w)/sqrt(pi)
    dhdc = dudc(c, h, i, theta) - beta*R*Edvdc
    dhdi = dudi(c, h, i, theta) + beta*delta*Edvdi - beta*R*Edvdc
    return(c(dhdc, dhdi))
  }
  return(hhgrad)
}

constraintgrad <- function(ci){
  return(c(-1, -1))
}

### Bellman Operator

bellman <- function(w, theta, shockpts = 13, m = 10000){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  Tw0 = w$Tw
  statespace = w[,c("x", "h", "y", "B", "cmin", "hmin")]
  gqpts = gaussHermiteData(shockpts)
  Vfxlist = lapply(ygrid, interpolater.creater, statespace, Tw0)
  optwrap = function(rowidx){
    x = statespace$x[rowidx]; h = statespace$h[rowidx]; y = statespace$y[rowidx]
    B = statespace$B[rowidx]; cmin = statespace$cmin[rowidx]; hmin = statespace$hmin[rowidx]
    Valfunc = Vfxlist[[which(ygrid==y)]]
    if (round(x - B,10) <= round(cmin,10) | round(x - B - cmin,10) <= round(hmin - h,10)) {
      vfx = u(cmin, h = max(h, hmin), i = 0, theta) + beta*Valfunc(x = B, h = max(delta*h, hmin))
      return(c(cmin, 0, vfx))
    } else {
      hhprob = hhprob_funkfact(x, h, y, B, hmin, Valfunc, theta, gqpts)
      #hhgrad = hhgradfact(x, h, yi = y, B, hmin, theta, s = statespace, Tw = Tw0, gqpts)
      lbi = max(round(hmin - h,10), 0.0)
      #optgrid = crossing(c = seq(cmin, x-B, length.out = optgridpts), i = seq(lbi, x-B, length.out = optgridpts))
      #optgrid = rbind(optgrid, unique(data.frame("c" = seq(cmin, x-B-lbi, length.out = optgridpts), "i" = x-B - seq(cmin, x-B-lbi, length.out = optgridpts)))) %>%
      #  filter(round(x - c - i - B, 10) >= 0 & i > lbi)
      #optgrid$tw = unlist(lapply(c(1:nrow(optgrid)), function(j) hhprob[[1]](c(optgrid$c[j], optgrid$i[j]))))
      #sps = c(optgrid$c[which.min(optgrid$tw)], optgrid$i[which.min(optgrid$tw)])
      sps = incenter(cmin, lbi, x, B)
      sol = cobyla(x0 = sps, fn = hhprob[[1]], hin = hhprob[[2]],
                   lower = c(cmin,lbi), upper = c(x-B, x-B), 
                   control = list(ftol_rel = 1e-8, ftol_abs = 0, xtol_rel = 0, maxeval = m))
      sol2 = cobyla(x0 = cmin, fn = function(c){hhprob[[1]](c(c, x-B-c))}, 
                    lower = cmin, upper = x-B, 
                    control = list(ftol_rel = 1e-8, ftol_abs = 0, xtol_rel = 0, maxeval = m))
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
  #rlist = lapply(ygrid, monotonizer, w = wnext)
  #wnext = do.call(rbind, rlist)
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
VFI <- function(theta, maxiter = 30, tol = 1e-8){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  
  ### statespace
  statespace = create.statespace(ubm = 3, theta, method = "log")
  ### Initial guess 
  v0 = guesser(statespace, theta)
  
  w = vector(mode = "list", length = maxiter)
  w[[1]] = cbind(statespace, v0)
  d = 1; i = 2
  while (i<maxiter & d>tol) {
    wnext = bellman(w = w[[i-1]], theta)
    #w[[i]] = wnext
    ## howard
    #k = vector(mode = "list", length = howardk+1)
    #k[[1]] = wnext
    #for (j in c(1:howardk)) {
    #  wnext$Tw = howard(wnext, statespace, theta)
    #  k[[j+1]] = wnext
    #}
    ### MQP error bounds
    #b = beta/(1-beta)*mean(k[[howardk+1]]$Tw - k[[howardk]]$Tw)
    #wnext$Tw = wnext$Tw + b
    ## check tol
    w[[i]] = wnext
    #saveRDS(w, "/users/mdgordo/Desktop/w.rds")
    d = mean((w[[i]]$Tw - w[[i-1]]$Tw)^2)/mean((w[[1]]$Tw)^2)
    print(d)
    i = i+1
    print(i)
  }
  w[[i]] = bellman(w = w[[i-1]], statespace, theta, shockpts = 31, m = 10000)
  return(compact(w))
}

momentmatcher <- function(i, vfx, t0, data){
  gamma = t0[1]; beta = t0[2]; R = t0[3]; cbar = t0[4]; hbar = t0[5]
  lambda = t0[6]; sigma = t0[7]; alpha = t0[8]; delta = t0[9]
  df = data[i,]
  wt = df$wt_hh/sum(data$wt_hh)
  xyh = c(df$imputed_bufferstock, df$avg_inc, df$lag_h)
  e1  = log(complexr(x0 = xyh, statespace = vfx[,c(1:6)], vfx = vfx$cfx)) - log(df$food_consumption)
  e2  = log(complexr(x0 = xyh, statespace = vfx[,c(1:6)], vfx = vfx$ifx)+1) - log(df$home_investment+1)
  e3 = log(df$home_value+1) - log(delta*(df$lag_h + df$home_investment)+1)
  e4 = df$imputed_bufferstock - df$quake_aid - df$lag_y - R*(df$lag_x - df$lag_c - df$lag_i)
  e5 = sqrt(df$var_inc) - sigma*df$avg_inc
  e6 = e1*df$imputed_bufferstock
  e7 = e1*df$avg_inc
  e8 = e1*log(df$lag_h+1)
  e9 = e2*df$imputed_bufferstock
  e10 = e2*df$avg_inc
  e11 = e2*log(df$lag_h+1)
  return(c(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11)*wt)
}

hhgradient <- function(i, df, vfx, dvlist){
  hh = df[i,]
  wt = hh$wt_hh/sum(df$wt_hh)
  xyh = c(hh$imputed_bufferstock, hh$avg_inc, hh$lag_h)
  
  cfx0 = complexr(x0 = xyh, statespace = vfx[,c(1:6)], vfx = vfx$cfx)
  ifx0 = complexr(x0 = xyh, statespace = vfx[,c(1:6)], vfx = vfx$ifx)
  
  dcdg = (complexr(x0 = xyh, statespace = dvlist[[1]][,c(1:6)], vfx = dvlist[[1]]$cfx) - cfx0)/.01
  didg = (complexr(x0 = xyh, statespace = dvlist[[1]][,c(1:6)], vfx = dvlist[[1]]$ifx) - ifx0)/.01
  dcdB = (complexr(x0 = xyh, statespace = dvlist[[2]][,c(1:6)], vfx = dvlist[[2]]$cfx) - cfx0)/.01
  didB = (complexr(x0 = xyh, statespace = dvlist[[2]][,c(1:6)], vfx = dvlist[[2]]$ifx) - ifx0)/.01
  dcdR = (complexr(x0 = xyh, statespace = dvlist[[3]][,c(1:6)], vfx = dvlist[[3]]$cfx) - cfx0)/.01
  didR = (complexr(x0 = xyh, statespace = dvlist[[3]][,c(1:6)], vfx = dvlist[[3]]$ifx) - ifx0)/.01
  dcdc = (complexr(x0 = xyh, statespace = dvlist[[4]][,c(1:6)], vfx = dvlist[[4]]$cfx) - cfx0)/.01
  didc = (complexr(x0 = xyh, statespace = dvlist[[4]][,c(1:6)], vfx = dvlist[[4]]$ifx) - ifx0)/.01
  dcdh = (complexr(x0 = xyh, statespace = dvlist[[5]][,c(1:6)], vfx = dvlist[[5]]$cfx) - cfx0)/.01
  didh = (complexr(x0 = xyh, statespace = dvlist[[5]][,c(1:6)], vfx = dvlist[[5]]$ifx) - ifx0)/.01
  dcdl = (complexr(x0 = xyh, statespace = dvlist[[6]][,c(1:6)], vfx = dvlist[[6]]$cfx) - cfx0)/.01
  didl = (complexr(x0 = xyh, statespace = dvlist[[6]][,c(1:6)], vfx = dvlist[[6]]$ifx) - ifx0)/.01
  dcds = (complexr(x0 = xyh, statespace = dvlist[[7]][,c(1:6)], vfx = dvlist[[7]]$cfx) - cfx0)/.001
  dids = (complexr(x0 = xyh, statespace = dvlist[[7]][,c(1:6)], vfx = dvlist[[7]]$ifx) - ifx0)/.001
  dcda = (complexr(x0 = xyh, statespace = dvlist[[8]][,c(1:6)], vfx = dvlist[[8]]$cfx) - cfx0)/.01
  dida = (complexr(x0 = xyh, statespace = dvlist[[8]][,c(1:6)], vfx = dvlist[[8]]$ifx) - ifx0)/.01
  dcdd = (complexr(x0 = xyh, statespace = dvlist[[9]][,c(1:6)], vfx = dvlist[[9]]$cfx) - cfx0)/.01
  didd = (complexr(x0 = xyh, statespace = dvlist[[9]][,c(1:6)], vfx = dvlist[[9]]$ifx) - ifx0)/.01
  
  r1 = wt*c(dcdg, dcdB, dcdR, dcdc, dcdh, dcdl, dcds, dcda, dcdd)/cfx0
  r2 = wt*c(didg, didB, didR, didc, didh, didl, dids, dida, didd)/(ifx0+1)
  r3 = wt*c(rep(0,8), -(hh$lag_h + hh$home_investment)/(delta*(hh$lag_h + hh$home_investment)))
  r4 = wt*c(0,0, -(hh$lag_x - hh$lag_c - hh$lag_i), rep(0,6))
  r5 = wt*c(rep(0,6), hh$avg_inc, 0, 0)
  r6 = wt*hh$imputed_bufferstock*r1
  r7 = wt*hh$avg_inc*r1
  r8 = wt*log(hh$lag_h+1)*r1
  r9 = wt*hh$imputed_bufferstock*r2
  r10 = wt*hh$avg_inc*r2
  r11 = wt*log(hh$lag_h+1)*r2
  return(rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11))
}

momentgradient <- function(theta, df){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  
  V = VFI(theta)
  vfx = V[[length(V)]]
  ## numerical for policy functions
  tp = c(theta[1] + .01, theta[c(2:9)])
  dVdg = VFI(tp); print("D1 done")
  dVdg = dVdg[[length(dVdg)]]
  tp = c(theta[1], theta[2] + .01, theta[c(3:9)])
  dVdB = VFI(tp); print("D2 done")
  dVdB = dVdB[[length(dVdB)]]
  tp = c(theta[c(1:2)], theta[3] + .01, theta[c(4:9)])
  dVdR = VFI(tp); print("D3 done")
  dVdR = dVdR[[length(dVdR)]]
  tp = c(theta[c(1:3)], theta[4] + .01, theta[c(5:9)])
  dVdc = VFI(tp); print("D4 done")
  dVdc = dVdc[[length(dVdc)]]
  tp = c(theta[c(1:4)], theta[5] + .01, theta[c(6:9)])
  dVdh = VFI(tp); print("D5 done")
  dVdh = dVdh[[length(dVdh)]]
  tp = c(theta[c(1:5)], theta[6] + .01, theta[c(7:9)])
  dVdl = VFI(tp); print("D6 done")
  dVdl = dVdl[[length(dVdl)]]
  tp = c(theta[c(1:6)], theta[7] + .001, theta[c(8:9)])
  dVds = VFI(tp); print("D7 done")
  dVds = dVds[[length(dVds)]]
  tp = c(theta[c(1:7)], theta[8] + .01, theta[c(9)])
  dVda = VFI(tp); print("D8 done")
  dVda = dVda[[length(dVda)]]
  tp = c(theta[c(1:8)], theta[9] + .01)
  dVdd = VFI(tp); print("D9 done")
  dVdd = dVdd[[length(dVdd)]]
  
  dvlist = list(dVdg, dVdB, dVdR, dVdc, dVdh, dVdl, dVds, dVda, dVdd)
  gradlist = mclapply(c(1:nrow(df)), hhgradient, df, vfx, dvlist, mc.cores = detectCores())
  grad = apply(simplify2array(gradlist), c(1,2), mean)
  return(grad)
}

### check tolerance between iterations

tolchecker = function(w){
  for (i in c(2:length(w))) {
    d = mean((w[[i]]$Tw - w[[i-1]]$Tw)^2)
    cd = mean((w[[i]]$cfx - w[[i-1]]$cfx)^2) 
    id = mean((w[[i]]$ifx - w[[i-1]]$ifx)^2) 
    print(paste("vtol:", d, "ctol:", cd, "itol", id, sep = " "))
  }
}

### check to make sure grid is not binding

satisficer = function(w, yidx, sigma){
  Vfin = w[[length(w)]]
  xmax = max(Vfin$x[Vfin$y==yidx])
  hmax = max(Vfin$h[Vfin$y==yidx])
  ymax = max(great.expectations(yidx, sigma, gqpts))
  xnext = xmax - (xmax-exp(ymax))/R
  hnext = hmax*(1/delta - 1)
  df.cstrt = Vfin %>%
    filter(y==yidx & x==xmax) %>%
    mutate(spend = cfx + ifx,
           constrainedx = spend<xnext)
  df.cstrt2 = Vfin %>%
    filter(y==yidx & h==hmax) %>%
    mutate(constrainedh = ifx>hnext)
  return(c(sum(df.cstrt$constrainedx), sum(df.cstrt2$constrainedh)))
}

### Counterfactuals

counterfactualizer = function(i, vfx, df.hh){
  x1 = as.double(df.hh[i,"x_noaid"]); x2 = as.double(df.hh[i,"x_aid"])
  y = as.double(df.hh[i,"avg_inc"]); h = as.double(df.hh[i,"lag_h"])
  c1 = complexr(x0 = c(x1, y, h), vfx, vfx = vfx$cfx)
  c2 = complexr(x0 = c(x2, y, h), vfx, vfx = vfx$cfx)
  i1 = complexr(x0 = c(x1, y, h), vfx, vfx = vfx$ifx)
  i2 = complexr(x0 = c(x2, y, h), vfx, vfx = vfx$ifx)
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
  if (groot(min(ygrid))<0) {yp = min(ygrid)} else {
    r = uniroot(groot, interval = c(min(ygrid), y))
    yp = r$root
  }
  wtp = (exp(y + (sigma*y)^2/2) - exp(yp + (sigma*yp)^2/2))/(1-beta)
  return(wtp)
}


### Functions for plotting

iterplotter = function(w, hidx=NULL, xidx=NULL, yidx, ifilt = c(1:99), xfilt = NULL, var = "Tw"){
  statespace = w[[length(w)]]
  ## plot iterations for given y and h
  if (is.null(xidx)) {
    xs = lapply(w, function(x) filter(x, y==yidx, h==hidx))
    xs <- data.frame(do.call(rbind, xs)) 
    xs$iteration <- rep(c(1:length(w)), each = length(unique(xs$x)))
    if (!is.null(xfilt)) xs = filter(xs, x!=xfilt)
    ggplot(filter(xs, iteration %in% ifilt)) +
      geom_line(aes(x = x, y = get(var), color = iteration, group = iteration)) + 
      scale_color_viridis_c()
  } else {
    xs = lapply(w, function(x) filter(x, y==yidx, x==xidx))
    xs <- data.frame(do.call(rbind, xs)) 
    xs$iteration <- rep(c(1:length(w)), each = length(unique(xs$h)))
    ggplot(filter(xs, iteration %in% ifilt)) +
      geom_line(aes(x = h, y = get(var), color = iteration, group = iteration)) + 
      scale_color_viridis_c()
  }
}

threedplotter = function(w, yidx, fillvar, iter = length(w), d3 = TRUE, aidamt = NULL, ubm = 1.1, lbm = NULL, lbh = NULL){
  statespace = w[[iter]][,c(1:6)]
  if (fillvar=="Tw") {
    Vfx = interpolater.creater(ygrid[yidx], statespace, w[[iter]]$Tw)
  } else if (fillvar %in% c("cfx", "aidcfx")) {
    Vfx = interpolater.creater(ygrid[yidx], statespace, w[[iter]]$cfx)
  } else {Vfx = interpolater.creater(ygrid[yidx], statespace, w[[iter]]$ifx)}
  if(is.null(lbm)) lbm = cbar*exp(ygrid[yidx])
  if(is.null(lbh)) lbh = hbar*exp(ygrid[yidx])
  df.test = crossing(x = seq(lbm, exp(ubm*ygrid[yidx]), length.out = 200), h = seq(lbh, exp(ubm*ygrid[yidx]), length.out = 200))
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

yplotter = function(w, hidx, iter = length(w), d3 = TRUE, xvals, yvals){
  statespace = w[[iter]][,c(1:6)]
  xp = seq(xvals[1], xvals[2], length.out = 200)
  yp = seq(yvals[1], yvals[2], length.out = 200)
  df.test = crossing(x = xp, y = yp)
  statespace$Tw = w[[iter]]$Tw
  s <- filter(statespace, h == hidx)
  fi = function(xp,yp){
    r = simplr(x0 = xp, y0 = yp, x = s$x, y = s$y, vfx = s$Tw, method = "simplical", extrapolation.warning = FALSE)
    return(r)
  }
  df.test$proj = mapply(fi, x = df.test$x, y = df.test$y)
  threed = ggplot(df.test) + geom_tile(aes(x = x, y = y, fill = -1*log(-proj))) + scale_fill_viridis_c()
  if (d3) p = plot_gg(threed, multicore=TRUE) else p = threed
  p
}



