incenter <- function(cbar, lbi, x, lambda) {
  return(c((cbar*2 + x+lambda-lbi)/3, (lbi*2 + x+lambda-cbar)/3))
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
great.expectations <- function(sigma, gqpts){
  mu = sqrt(2)*sigma*gqpts$x - sigma^2/2
  return(mu)
}

create.statespace = function(ubm = c(17,12), theta, method = "equal"){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  
  defaultpt = (cbar-lambda)
  
  if (method=="chebyshev") {
    cnodesx = chebnodes(xn, lb = defaultpt, ub = ubm[1])
    x = c(-lambda, cnodesx)
    h = chebnodes(hn, lb = 0, ub = ubm[2])
  } else if (method=="log") {
    h = c(0, hbar/2, exp(seq(log(hbar), log(ubm[2]), length.out = hn-1)))
    x = c(-lambda, exp(seq(0, log(ubm[1] - defaultpt + 1), length.out = hn)) + defaultpt - 1)
  } else if (method=="uneven") {
    x = c(-lambda, seq(defaultpt, ubm[1], length.out = xn))
    h = unique(c(seq(0, ubm[2]/2, length.out = round(2*hn/3,0)), seq(ubm[2]/2, ubm[2], length.out = round(1*hn/3,0)+1)))
  } else {
    x = c(-lambda, seq(defaultpt, ubm[1], length.out = xn))
    h = seq(0, ubm[2], length.out = hn)
  }
  css = crossing(x,h)
  return(css)
}

### 2D interpolate
interpolater.creater = function(w, theta, method = "simplical", var = "Tw"){
  cbar = theta[4]; hbar = theta[5]; lambda = theta[6]
  if(var=="Tw") Tw = w$Tw else if (var %in% c("cfx", "aidcfx")) Tw = w$cfx else if (var %in% c("ifx", "aidifx")) Tw = w$ifx
  if (method=="simplical") {
    fi = function(x,h){
      if (x+h<=hbar+cbar-lambda) {x = cbar-lambda} ## ensures default zone is enforced
      r = simplr(x0 = x, y0 = h, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
      return(r)
    }
  } else {
    cs = filter(w, x>cbar-lambda)
    if(var=="Tw") TwC = cs$Tw else if (var %in% c("cfx", "aidcfx")) TwC = cs$cfx else if (var %in% c("ifx", "aidifx")) TwC = cs$ifx
    if (method=="chebyshev") {
      tmat = matrix(TwC, ncol = xn)
      cmat = chebcoefs(tmat, degree = hn-1)
      fi = function(x, h){
        if (x+h<=hbar+cbar-lambda) {x = cbar-lambda} ## ensures default zone is enforced
        if(x<=cbar-lambda){
          r = simplr(x0 = x, y0 = h, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
        } else {
          if (x>max(w$x)) x = max(w$x); if (h>max(w$h)) h = max(w$h)
          r = chebpred(x,h, cmat, lb = c(cbar-lambda, 0), ub = c(max(w$x), max(w$h)))
          return(r)
        }
      } 
  } else if (method=="spline"){
      obj = (TwC - min(TwC))/(max(TwC) - min(TwC))
      xnorm = (cs$x - min(cs$x))/(max(cs$x) - min(cs$x))
      hnorm = (cs$h - min(cs$h))/(max(cs$h) - min(cs$h))
      mod = gam(obj ~ te(xnorm, hnorm))
      fi = function(x, h){
        if (x+h<=hbar+cbar-lambda) {x = cbar-lambda} ## ensures default zone is enforced
        xn = (x - min(cs$x))/(max(cs$x) - min(cs$x))
        hn = (h - min(cs$h))/(max(cs$h) - min(cs$h))
        if(x<=cbar-lambda){
          r = simplr(x0 = x, y0 = h, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
        } else {
          if (x>max(w$x)) xn = 1; if (h>max(w$h)) hn = 1
          r = predict(mod, newdata = data.frame("x" = xn, "h" = hn))
          return(r*(max(TwC) - min(TwC)) + min(TwC))
        }
      } 
    } else if (method == "neuralnet") {
        obj = (TwC - min(TwC))/(max(TwC) - min(TwC))
        xnorm = (cs$x - min(cs$x))/(max(cs$x) - min(cs$x))
        hnorm = (cs$h - min(cs$h))/(max(cs$h) - min(cs$h))
        mod = neuralnet(obj ~ xnorm + hnorm, data = data.frame(obj = obj, xnorm = xnorm, hnorm = hnorm), 
                        hidden=11, act.fct = "logistic", linear.output = FALSE, stepmax = 5e5)
        fi = function(x, h){
          if (x+h<=hbar+cbar-lambda) {x = cbar-lambda} ## ensures default zone is enforced
          xn = (x - min(cs$x))/(max(cs$x) - min(cs$x))
          hn = (h - min(cs$h))/(max(cs$h) - min(cs$h))
          if(x<=cbar-lambda){
            r = simplr(x0 = x, y0 = h, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
          } else {
            if (x>max(w$x)) xn = 1; if (h>max(w$h)) hn = 1
            r = predict(mod, newdata = data.frame("x" = xn, "h" = hn))
            return(r*(max(TwC) - min(TwC)) + min(TwC))
          }
      } 
    }
  }
  return(fi)
}

### 3D interpolate 
complexr <- function(x0, statespace, vfx, theta){
  ## find closest points in y and interpolate
  yvals = sort(unique(statespace$y))
  yi = findInterval(x0[2], yvals, all.inside = TRUE)
  if (x0[2]>max(yvals)) x0[2] <- max(yvals)
  if (x0[2]<min(yvals)) x0[2] <- min(yvals)
  yp = c(yvals[yi], yvals[yi+1])
  xyhfilt = filter(statespace, y==yp[1])
  
  vy1 = interpolater.creater(yp[1], theta, statespace, Tw = vfx)
  vy2 = interpolater.creater(yp[2], theta, statespace, Tw = vfx)
  
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

monotonizer <- function(w, policy = TRUE){
  xvals = sort(unique(w$x))
  vrT = w %>%
    pivot_wider(id_cols = x, names_from = h, values_from = Tw) %>%
    select(!x) %>% as.matrix()
  finalVrT = rearrangement(x = list(sort(unique(w$x)), sort(unique(w$h))), 
                           y = vrT, avg = TRUE) %>% as_tibble()
  finalVrT$x = xvals
  finalVrT <- pivot_longer(finalVrT, !x, names_to = "h", values_to = "Tw")
  finalVrT$h = as.double(finalVrT$h)
  if (policy){
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
    finalVrI = rearrangement(x = list(sort(unique(w$x)), sort(-1*unique(w$h))), 
                             y = vrI[,rev(c(1:ncol(vrI)))], avg = TRUE) %>% as_tibble()
    finalVrI$x = xvals
    finalVrI <- pivot_longer(finalVrI, !x, names_to = "h", values_to = "ifx")
    finalVrI$h = as.double(finalVrI$h)
    finalVrP = merge(finalVrC, finalVrI, by = c("x", "h"))
  } else {
    finalVrP = w[,c("x", "h", "cfx", "ifx")]
  }
  finalVr = merge(finalVrT, finalVrP, by = c("x", "h"))
  finalVr = finalVr[order(finalVr$x,finalVr$h),]
  finalVr = finalVr[,c("x", "h", "Tw", "cfx", "ifx")]
  rownames(finalVr) = NULL
  return(finalVr)
}

### Function factory for generation household maximization problem

hhprob_funkfact <- function(x, h, Valfunc, theta, gqpts){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]; B = -1*lambda
  
  hhprob = function(ci){
    c = ci[1]; i = ci[2]
    if (is.nan(c)|is.na(c)|c+i>x-B) {return((c+i)*1e10)} else {  ### built in constraint 
      payoff = u(c, h, i, theta)
      htplus1 = delta*(h+i)
      xtplus1 = R*(x - c - i) 
      mindraw = B+cbar+max(hbar-htplus1, 0) - xtplus1 ## minimum draw to not default
      pdef = plnorm(mindraw, meanlog = -sigma^2/2, sdlog = sigma)
      draws = exp(great.expectations(sigma, gqpts))
      if (any(draws>mindraw)){
        wts = gqpts$w[draws>mindraw]/sqrt(pi)
        draws = draws[draws>mindraw]
        if (pdef!=0){  ### adjust wts and draws for probability of default
          extrawt = 1 - pdef - sum(wts)
          if (extrawt>0) wts = c(extrawt, wts) else wts[1] = extrawt + wts[1]
          extradraw = qlnorm(pdef + wts[1]/2, meanlog = -sigma^2/2, sdlog = sigma)
          if (extrawt>0) draws = c(extradraw, draws) else draws[1] = extradraw
        }
        EVp = sum(mapply(Valfunc, x = xtplus1 + draws, h = htplus1)*wts)
        EV = pdef*Valfunc(B, max(htplus1, delta*hbar)) + EVp
      } else {
        EV = pdef*Valfunc(B, max(htplus1, delta*hbar)) + (1-pdef)*Valfunc(xtplus1+mindraw+1.856, htplus1) ### mean excess for lognormal
      }
      r = payoff + beta*EV
      return(-r)}
  }
  return(hhprob)
}

hhgradfact <- function(x, h, theta, w, gqpts){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]; B = -1*lambda
  hhgrad <- function(ci){
    c = ci[1]; i = ci[2]
    htplus1 = delta*(h+i)
    xtplus1 = R*(x - c - i) 
    mindraw = B+cbar+max(hbar-htplus1, 0) - xtplus1
    pdef = plnorm(mindraw, meanlog = -sigma^2/2, sdlog = sigma)
    draws = exp(great.expectations(sigma, gqpts))
    if (any(draws>mindraw)){
      wts = gqpts$w[draws>mindraw]/sqrt(pi)
      draws = draws[draws>mindraw]
      if (pdef!=0){  ### adjust wts and draws for probability of default
        extrawt = 1 - pdef - sum(wts)
        if (extrawt>0) wts = c(extrawt, wts) else wts[1] = extrawt + wts[1]
        extradraw = qlnorm(pdef + wts[1]/2, meanlog = -sigma^2/2, sdlog = sigma)
        if (extrawt>0) draws = c(extradraw, draws) else draws[1] = extradraw
      }
      Edvdc = sum(mapply(simplrdifferentiator, x0 = xtplus1 + draws, y0 = htplus1, statespace = list(w), deriv = list("dx"))*wts)
      Edvdi = sum(mapply(simplrdifferentiator, x0 = xtplus1 + draws, y0 = htplus1, statespace = list(w), deriv = list("dy"))*wts)
      EVc = pdef*simplrdifferentiator(x0 = B, y0 = max(htplus1, delta*hbar), statespace = w, deriv = list("dx")) + Edvdc
      EVi = pdef*simplrdifferentiator(x0 = B, y0 = max(htplus1, delta*hbar), statespace = w, deriv = list("dy")) + Edvdi
    } else {
      EVc = pdef*simplrdifferentiator(x0 = B, y0 = max(htplus1, delta*hbar), statespace = w, deriv = list("dx")) + 
        (1-pdef)*simplrdifferentiator(x0 = xtplus1+mindraw+1.856, y0 = htplus1, statespace = w, deriv = list("dx"))
      EVi = pdef*simplrdifferentiator(x0 = B, y0 = max(htplus1, delta*hbar), statespace = w, deriv = list("dy")) + 
        (1-pdef)*simplrdifferentiator(x0 = xtplus1+mindraw+1.856, y0 = htplus1, statespace = w, deriv = list("dy"))
    }
    dhdc = beta*R*EVc - dudc(c, h, i, theta)
    dhdi = beta*R*EVc - dudi(c, h, i, theta) - beta*delta*EVi
    return(c(dhdc, dhdi))
  }
  return(hhgrad)
}

constraintgrad <- function(ci){
  return(c(-1, -1))
}

### Bellman Operator

bellman <- function(w, theta, shockpts = 13, m = 2000){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  gqpts = gaussHermiteData(shockpts)
  Valfunc = interpolater.creater(w, theta, method = "simplical")
  optwrap = function(rowidx){
    x = w$x[rowidx]; h = w$h[rowidx]
    if (round(x + lambda,10) <= round(cbar,10) | round(x + lambda - cbar,10) <= round(hbar - h,10)) {
      vfx = u(cbar, h = max(h, hbar), i = 0, theta) + beta*Valfunc(x = -1*lambda, h = max(delta*h, hbar))
      return(c(cbar, 0, vfx))
    } else {
      hhprob = hhprob_funkfact(x, h, Valfunc, theta, gqpts)
      lbi = max(round(hbar - h,10), 0.0)
      sol = directL(fn = hhprob, 
                  lower = c(cbar,lbi), upper = c(x+lambda-lbi, x+lambda-cbar), 
                  control = list(ftol_rel = 1e-4, ftol_abs = 0, xtol_rel = 0, maxeval = m))
      r = c(sol$par, -sol$value)
      vfxdefault = u(cbar, h = max(h, hbar), i = 0, theta) + beta*Valfunc(x = -1*lambda, h = delta*max(h, hbar))
      if (vfxdefault < -sol$value) {
        return(r)
      } else {return(c(cbar, 0, vfxdefault))}
    }
  }
  Tw = mclapply(c(1:nrow(w)), optwrap, mc.cores = detectCores())
  wnext = cbind(w[,c("x", "h")], data.frame("Tw" = unlist(lapply(Tw, function(x) x[3])),
                                       "cfx" = unlist(lapply(Tw, function(x) x[1])),
                                       "ifx" = unlist(lapply(Tw, function(x) x[2]))))
  return(wnext)
}

### Howard Policy Iteration
howard <- function(w, theta, shockpts = 13){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  gqpts = gaussHermiteData(shockpts)
  Valfunc = interpolater.creater(w, theta, method = "simplical")
  
  accelerator = function(rowidx){
    x = w$x[rowidx]; h = w$h[rowidx]
    cfx = w$cfx[rowidx]; ifx = w$ifx[rowidx]
    if (round(x + lambda,10) <= round(cbar,10) | round(x + lambda - cbar,10) <= round(hbar - h,10)) {
      v = u(cbar, h = max(h, hbar), i = 0, theta) + beta*Valfunc(x = -1*lambda, h = delta*max(h, hbar))
      return(v)
    } else {
      hhprob = hhprob_funkfact(x, h, Valfunc, theta, gqpts)
      v1 = -1*hhprob(c(cfx, ifx))
      v2 = u(cbar, h = max(h, hbar), i = 0, theta) + beta*Valfunc(x = -1*lambda, h = delta*max(h, hbar))
      return(max(v1, v2))
    }
  }
  Twh = mclapply(1:nrow(w), accelerator, mc.cores = detectCores())
  return(unlist(Twh))
}

### VFI
VFI <- function(v0, theta, maxiter = 30, tol = 1e-5, howardk = 10){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  w = vector(mode = "list", length = maxiter)
  w[[1]] = v0
  tol1 = FALSE; tol2 = FALSE; i = 2; d = 1; s = 13
  while (i<=maxiter & !tol2) {
    wnext = bellman(w = w[[i-1]], theta, shockpts = s)
    ## howard
    k = vector(mode = "list", length = howardk+1)
    k[[1]] = wnext
    for (j in c(1:howardk)) {
      wnext$Tw = howard(wnext, theta, shockpts = s)
      k[[j+1]] = wnext
    }
    ### MQP error bounds
    b = beta/(1-beta)*mean(k[[howardk+1]]$Tw - k[[howardk]]$Tw)
    wnext$Tw = wnext$Tw + b
    w[[i]] = wnext
    ## check tol
    d = mean((w[[i]]$Tw - w[[i-1]]$Tw)^2)/mean((w[[i-1]]$Tw)^2)
    ### if tol is met increase gq points
    if (d < tol) {
      if (tol1) tol2 = TRUE else s = 101; tol1 = TRUE
    }
    print(d)
    i = i+1
    print(i)
  }
  return(compact(w))
}

momentmatcher <- function(i, vfx, t0, data){
  gamma = t0[1]; beta = t0[2]; R = t0[3]; cbar = t0[4]; hbar = t0[5]
  lambda = t0[6]; sigma = t0[7]; alpha = t0[8]; delta = t0[9]
  df = data[i,]
  wt = df$wt_hh/sum(data$wt_hh)
  e1  = log(vfx[[1]](df$imputed_bufferstock, df$lag_h)) - log(df$food_consumption)
  e2  = log(vfx[[2]](df$imputed_bufferstock, df$lag_h)+1) - log(df$home_investment+1)
  e3 = log(df$home_value+1) - log(delta*(df$lag_h + df$home_investment)+1)
  e4 = df$imputed_bufferstock - df$quake_aid - df$lag_y - R*(df$lag_x - df$lag_c - df$lag_i)
  e5 = sqrt(df$var_inc) - sigma
  e6 = e1*df$imputed_bufferstock
  e7 = e1*log(df$M_avg)
  e8 = e1*log(df$lag_h+1)
  e9 = e2*df$imputed_bufferstock
  e10 = e2*log(df$M_avg)
  e11 = e2*log(df$lag_h+1)
  return(c(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11)*wt)
}

hhgradient <- function(i, df, vfx, dvlist, theta){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  hh = df[i,]
  wt = hh$wt_hh/sum(df$wt_hh)
  
  cfx0 = vfx[[2]](hh$imputed_bufferstock, hh$lag_h)
  ifx0 = vfx[[3]](hh$imputed_bufferstock, hh$lag_h)
  
  dcdg = (dvlist[[1]][[2]](hh$imputed_bufferstock, hh$lag_h) - cfx0)/.01
  didg = (dvlist[[1]][[3]](hh$imputed_bufferstock, hh$lag_h) - ifx0)/.01
  dcdB = (dvlist[[2]][[2]](hh$imputed_bufferstock, hh$lag_h) - cfx0)/.01
  didB = (dvlist[[2]][[3]](hh$imputed_bufferstock, hh$lag_h) - ifx0)/.01
  dcdR = (dvlist[[3]][[2]](hh$imputed_bufferstock, hh$lag_h) - cfx0)/.01
  didR = (dvlist[[3]][[3]](hh$imputed_bufferstock, hh$lag_h) - ifx0)/.01
  dcdc = (dvlist[[4]][[2]](hh$imputed_bufferstock, hh$lag_h) - cfx0)/.01
  didc = (dvlist[[4]][[3]](hh$imputed_bufferstock, hh$lag_h) - ifx0)/.01
  dcdh = (dvlist[[5]][[2]](hh$imputed_bufferstock, hh$lag_h) - cfx0)/.01
  didh = (dvlist[[5]][[3]](hh$imputed_bufferstock, hh$lag_h) - ifx0)/.01
  dcdl = (dvlist[[6]][[2]](hh$imputed_bufferstock, hh$lag_h) - cfx0)/.01
  didl = (dvlist[[6]][[3]](hh$imputed_bufferstock, hh$lag_h) - ifx0)/.01
  dcds = (dvlist[[7]][[2]](hh$imputed_bufferstock, hh$lag_h) - cfx0)/.001
  dids = (dvlist[[7]][[3]](hh$imputed_bufferstock, hh$lag_h) - ifx0)/.001
  dcda = (dvlist[[8]][[2]](hh$imputed_bufferstock, hh$lag_h) - cfx0)/.01
  dida = (dvlist[[8]][[3]](hh$imputed_bufferstock, hh$lag_h) - ifx0)/.01
  dcdd = (dvlist[[9]][[2]](hh$imputed_bufferstock, hh$lag_h) - cfx0)/.01
  didd = (dvlist[[9]][[3]](hh$imputed_bufferstock, hh$lag_h) - ifx0)/.01
  
  r1 = wt*c(dcdg, dcdB, dcdR, dcdc, dcdh, dcdl, dcds, dcda, dcdd)/cfx0
  r2 = wt*c(didg, didB, didR, didc, didh, didl, dids, dida, didd)/(ifx0+1)
  r3 = wt*c(rep(0,8), -(hh$lag_h + hh$home_investment)/(delta*(hh$lag_h + hh$home_investment)+1))
  r4 = wt*c(0,0, -(hh$lag_x - hh$lag_c - hh$lag_i), rep(0,6))
  r5 = wt*c(rep(0,6), -1, 0, 0)
  r6 = wt*hh$imputed_bufferstock*r1
  r7 = wt*log(hh$M_avg)*r1
  r8 = wt*log(hh$lag_h+1)*r1
  r9 = wt*hh$imputed_bufferstock*r2
  r10 = wt*log(hh$M_avg)*r2
  r11 = wt*log(hh$lag_h+1)*r2
  return(rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11))
}

momentgradient <- function(theta, df){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  
  statespace = create.statespace(ubm = c(20,20), theta, method = "log")
  v0 = cbind(statespace, data.frame("Tw" = rep(0, nrow(statespace)),
                                    "cfx" = rep(0, nrow(statespace)),
                                    "ifx" = rep(0, nrow(statespace))))
  
  V = VFI(v0, theta)
  V = V[[length(V)]]
  vfx = list(V, interpolater.creater(V, theta, var = "cfx", method = "neuralnet"), interpolater.creater(V, theta, var = "ifx", method = "neuralnet"))
  ## numerical for policy functions
  tp = c(theta[1] + .01, theta[c(2:9)])
  dVdg = VFI(V, tp); print("D1 done")
  dVdg = dVdg[[length(dVdg)]]
  dVdg = list(dVdg, interpolater.creater(dVdg, theta, var = "cfx", method = "neuralnet"), interpolater.creater(dVdg, theta, var = "ifx", method = "neuralnet"))
  tp = c(theta[1], theta[2] + .01, theta[c(3:9)])
  dVdB = VFI(V, tp); print("D2 done")
  dVdB = dVdB[[length(dVdB)]]
  dVdB = list(dVdB, interpolater.creater(dVdB, theta, var = "cfx", method = "neuralnet"), interpolater.creater(dVdB, theta, var = "ifx", method = "neuralnet"))
  tp = c(theta[c(1:2)], theta[3] + .01, theta[c(4:9)])
  dVdR = VFI(V, tp); print("D3 done")
  dVdR = dVdR[[length(dVdR)]]
  dVdR = list(dVdR, interpolater.creater(dVdR, theta, var = "cfx", method = "neuralnet"), interpolater.creater(dVdR, theta, var = "ifx", method = "neuralnet"))
  tp = c(theta[c(1:3)], theta[4] + .01, theta[c(5:9)])
  dVdc = VFI(V, tp); print("D4 done")
  dVdc = dVdc[[length(dVdc)]]
  dVdc = list(dVdc, interpolater.creater(dVdc, theta, var = "cfx", method = "neuralnet"), interpolater.creater(dVdc, theta, var = "ifx", method = "neuralnet"))
  tp = c(theta[c(1:4)], theta[5] + .01, theta[c(6:9)])
  dVdh = VFI(V, tp); print("D5 done")
  dVdh = dVdh[[length(dVdh)]]
  dVdh = list(dVdh, interpolater.creater(dVdh, theta, var = "cfx", method = "neuralnet"), interpolater.creater(dVdh, theta, var = "ifx", method = "neuralnet"))
  tp = c(theta[c(1:5)], theta[6] + .01, theta[c(7:9)])
  dVdl = VFI(V, tp); print("D6 done")
  dVdl = dVdl[[length(dVdl)]]
  dVdl = list(dVdl, interpolater.creater(dVdl, theta, var = "cfx", method = "neuralnet"), interpolater.creater(dVdl, theta, var = "ifx", method = "neuralnet"))
  tp = c(theta[c(1:6)], theta[7] + .001, theta[c(8:9)])
  dVds = VFI(V, tp); print("D7 done")
  dVds = dVds[[length(dVds)]]
  dVds = list(dVds, interpolater.creater(dVds, theta, var = "cfx", method = "neuralnet"), interpolater.creater(dVds, theta, var = "ifx"))
  tp = c(theta[c(1:7)], theta[8] + .01, theta[c(9)])
  dVda = VFI(V, tp); print("D8 done")
  dVda = dVda[[length(dVda)]]
  dVda = list(dVda, interpolater.creater(dVda, theta, var = "cfx", method = "neuralnet"), interpolater.creater(dVda, theta, var = "ifx", method = "neuralnet"))
  tp = c(theta[c(1:8)], theta[9] + .01)
  dVdd = VFI(V, tp); print("D9 done")
  dVdd = dVdd[[length(dVdd)]]
  dVdd = list(dVdd, interpolater.creater(dVdd, theta, var = "cfx", method = "neuralnet"), interpolater.creater(dVdd, theta, var = "ifx", method = "neuralnet"))
  
  dvlist = list(dVdg, dVdB, dVdR, dVdc, dVdh, dVdl, dVds, dVda, dVdd)
  gradlist = mclapply(c(1:nrow(df)), hhgradient, df = df, vfx = vfx, dvlist = dvlist, theta = theta, mc.cores = detectCores())
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

satisficer = function(w, sigma){
  Vfin = w[[length(w)]]
  xmax = max(Vfin$x)
  hmax = max(Vfin$h)
  ymax = max(exp(great.expectations(sigma, gqpts)))
  xnext = xmax - (xmax-ymax)/R
  hnext = hmax*(1/delta - 1)
  df.cstrt = Vfin %>%
    filter(x==xmax) %>%
    mutate(spend = cfx + ifx,
           constrainedx = spend<xnext,
           xnext = xnext)
  df.cstrt2 = Vfin %>%
    filter(h==hmax) %>%
    mutate(constrainedh = ifx>hnext,
           hnext = hnext)
  return(list(df.cstrt, df.cstrt2))
}

### Counterfactuals

counterfactualizer = function(i, vfx, df.hh){
  x1 = as.double(df.hh[i,"x_noaid"]); x2 = as.double(df.hh[i,"x_aid"])
  h = as.double(df.hh[i,"lag_h"])
  c1 = vfx[[1]](x1, h)
  c2 = vfx[[1]](x2, h)
  i1 = vfx[[2]](x1, h)
  i2 = vfx[[2]](x2, h)
  b1 = as.double(df.hh[i,"x_noaid"]) - c1 - i1
  b2 = as.double(df.hh[i,"x_aid"]) - c2 - i2
  return(c(c1, c2, i1, i2, b1, b2))
}

wtpeliciter = function(x, h, y, aidamt, vfx, theta){
  vfx = interpolater.creater(vfx, theta)
  aidfrac = aidamt/y
  tau = 1 - (vfx(x, h)/vfx(x+aidfrac,h))^(1/(1-theta[1]))
  return(tau)
}

utilbooster = function(x, h, y, aidamt, vfx, theta){
  vfx = interpolater.creater(vfx, theta)
  aidfrac = aidamt/y
  baselineutils = vfx(x, h)*y^(1-theta[1])
  boost  = vfx(x+aidfrac, h)*y^(1-theta[1])
  return(boost - baselineutils)
}

conditionalizer = function(x, h, y, aidamt, vfx, theta){
  vfx = interpolater.creater(vfx, theta)
  aidfrac = aidamt/y
  tau = 1 - (vfx(x, h)/vfx(x,h+aidfrac))^(1/(1-theta[1]))
  return(tau)
}

### Functions for plotting

iterplotter = function(w, hidx=NULL, xidx=NULL, ifilt = c(1:99), xfilt = NULL, var = "Tw"){
  statespace = w[[length(w)]]
  ## plot iterations for given y and h
  if (is.null(xidx)) {
    xs = lapply(w, function(x) filter(x, round(h,10)==round(hidx,10)))
    xs <- data.frame(do.call(rbind, xs)) 
    xs$iteration <- rep(c(1:length(w)), each = length(unique(round(xs$x,10))))
    if (!is.null(xfilt)) xs = filter(xs, !(x%in%xfilt))
    ggplot(filter(xs, iteration %in% ifilt)) +
      geom_line(aes(x = x, y = get(var), color = iteration, group = iteration)) + 
      scale_color_viridis_c()
  } else {
    xs = lapply(w, function(x) filter(x, round(x,10)==round(xidx,10)))
    xs <- data.frame(do.call(rbind, xs)) 
    xs$iteration <- rep(c(1:length(w)), each = length(unique(round(xs$h,10))))
    ggplot(filter(xs, iteration %in% ifilt)) +
      geom_line(aes(x = h, y = get(var), color = iteration, group = iteration)) + 
      scale_color_viridis_c()
  }
}

threedplotter = function(w, theta, fillvar, iter = length(w), method = "simplical",
                         d3 = TRUE, aidamt = NULL, ubm = 9, lbm = NULL, lbh = NULL){
  Vfx = interpolater.creater(w[[iter]], theta, var = fillvar, method = method)
  if(is.null(lbm)) lbm = -lambda
  if(is.null(lbh)) lbh = 0
  df.test = crossing(x = seq(lbm, ubm, length.out = 200), h = seq(lbh, ubm, length.out = 200))
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

hhprobtroubleshooter <- function(rowidx, w){
  gqpts = gaussHermiteData(13)
  Valfunc = interpolater.creater(w, theta, method = "simplical")
  x = w$x[rowidx]; h = w$h[rowidx]
  hhprob = hhprob_funkfact(x, h, Valfunc, theta, gqpts)
  lbi = max(round(hbar - h,10), 0.0)
  cs = seq(cbar, x+lambda-lbi, length.out = 40)
  is = seq(lbi, x+lambda-cbar, length.out = 40)
  testgrid = crossing(cs, is) %>%
    filter(cs + is <= x + lambda)
  testgrid$tw = sapply(c(1:nrow(testgrid)), function(i)hhprob(c(testgrid$cs[i], testgrid$is[i])))
  ggplot(testgrid) + geom_tile(aes(x = cs, y = is, fill = tw)) + scale_fill_viridis_c()
}

