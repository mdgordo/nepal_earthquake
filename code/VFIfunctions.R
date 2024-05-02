euclid <- function(a, b) {sqrt(sum((a - b)^2))}

ihs <- function(x){log(x + sqrt(x ^ 2 + 1))}
ihsinv <- function(x){(exp(x) - exp(-x))/2}

annuitycalc <- function(R, T){
  R^T * (R-1)/(R^T - 1)
}

### utility function
u <- function(c, h, i, theta) {
  alpha = theta[9]; gamma = theta[1]
  cd = c^alpha * (h + i)^(1-alpha)
  crra = cd^(1-gamma)/(1-gamma)
  return(crra)
}

### Creates points for gaussian quadrature
great.expectations <- function(sigma, gqpts){
  mu = sqrt(2)*sigma*gqpts$x - sigma^2/2
  return(mu)
}

create.statespace = function(ubm = c(5,5), theta, method = "equal"){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]*theta[9]; hbar = theta[5]*(1-theta[9])*theta[10]/(1-theta[10])
  lambda = theta[6]; sigma = theta[7]; sigmame = theta[8]; alpha = theta[9]; delta = theta[10]
  
  if (method=="chebyshev") {
    x = c(-lambda + .01, chebnodes(xn, lb = -lambda + cbar, ub = ubm[1]))
    h = chebnodes(hn, lb = 0.01, ub = ubm[2])
  } else if (method=="log") { 
    x =  exp(seq(0.01, log(ubm[1] + lambda + 1), length.out = hn)) - (lambda+1)
    h =  exp(seq(0.01, log(ubm[2]+1), length.out = hn))-1
  } else if (method=="uneven") {
    x = c(-lambda + .01, seq(-lambda + cbar, ubm[1], length.out = xn))
    h = unique(c(seq(.01, ubm[2]/2, length.out = round(2*hn/3,0)), seq(ubm[2]/2, ubm[2], length.out = round(1*hn/3,0)+1)))
  } else {
    x = c(-lambda + .01, seq(-lambda + cbar, ubm[1], length.out = xn))
    h = seq(.01, ubm[2], length.out = hn)
  }
  css = crossing(x,h)
  return(css)
}

simplr <- function (x0, y0, x, y, xvals, yvals, vfx, method = "bilinear") {
  xi = findInterval(x0, xvals, all.inside = TRUE)
  x1 = xvals[xi]
  x2 = xvals[xi + 1]
  yi = findInterval(y0, yvals, all.inside = TRUE)
  y1 = yvals[yi]
  y2 = yvals[yi + 1]
  if (x0 > max(x)) 
    x0 = max(x)
  else if (x0 < min(x)) 
    x0 = min(x)
  if (y0 > max(y)) 
    y0 = max(y)
  else if (y0 < min(y)) 
    y0 = min(y)
  v1 = vfx[which(x == x1 & y == y1)]
  v2 = vfx[which(x == x1 & y == y2)]
  v3 = vfx[which(x == x2 & y == y1)]
  v4 = vfx[which(x == x2 & y == y2)]
  if (method == "bilinear") {
    xp = -1 + 2 * (x0 - x1)/(x2 - x1)
    yp = -1 + 2 * (y0 - y1)/(y2 - y1)
    vp = 0.25 * (v1 * (1 - xp) * (1 - yp) + v2 * (1 - xp) * 
                   (1 + yp) + v3 * (1 + xp) * (1 - yp) + v4 * (1 + xp) * (1 + yp))
  }
  else {
    xp = (x0 - x1)/(x2 - x1)
    yp = (y0 - y1)/(y2 - y1)
    if (xp + yp <= 1) {
      vp = v1 * (1 - xp - yp) + v2 * yp + v3 * xp
    }
    else {
      vp = v2 * (1 - xp) + v3 * (1 - yp) + v4 * (xp + yp - 1)
    }
  }
  return(vp)
}

###3D interpolate
extrapolator <- function(w, x, h, xvals, hvals, Tw, var) {
  a = max(xvals); b = max(hvals); 
  a0 = xvals[length(xvals)-1];  b0 = hvals[length(hvals)-1]
  a0p = xvals[length(xvals)-3];  b0p = hvals[length(hvals)-3]
  if (x > a+1e-9 & h > b+1e-9) {  
    bslp = (h-b)/(x-a)
    x0 = a; h0 = b; z0 = simplr(x0 = a, y0 = b, x = w$x, y = w$h, xvals = xvals, yvals = hvals, 
                                vfx = Tw, method = "simplical")
    if (bslp >= b/a) {
      x1 = (b0-b)/bslp + a; h1 = b0; z1 = simplr(x0 = x1, y0 = h1, x = w$x, y = w$h, xvals = xvals, yvals = hvals, 
                                                 vfx = Tw, method = "simplical")
      x2 = (b0p-b)/bslp + a; h2 = b0p; z2 = simplr(x0 = x2, y0 = h2, x = w$x, y = w$h, xvals = xvals, yvals = hvals, 
                                                   vfx = Tw, method = "simplical")
    } else {
      x1 = a0; h1 = (a0-a)*bslp + b; z1 = simplr(x0 = x1, y0 = h1, x = w$x, y = w$h, xvals = xvals, yvals = hvals, 
                                                 vfx = Tw, method = "simplical")
      x2 = a0p; h2 = (a0p-a)*bslp + b; z2 = simplr(x0 = x2, y0 = h2, x = w$x, y = w$h, xvals = xvals, yvals = hvals, 
                                                   vfx = Tw, method = "simplical")
    }
  } else if (x > a+1e-9) {
    x0 = a; h0 = h; z0 = simplr(x0 = x0, y0 = h, x = w$x, y = w$h, xvals = xvals, yvals = hvals, 
                                vfx = Tw, method = "simplical")
    x1 = a0; h1 = h; z1 = simplr(x0 = x1, y0 = h1, x = w$x, y = w$h, xvals = xvals, yvals = hvals, 
                                 vfx = Tw, method = "simplical")
    x2 = a0p; h2 = h; z2 = simplr(x0 = x2, y0 = h2, x = w$x, y = w$h, xvals = xvals, yvals = hvals,
                                  vfx = Tw, method = "simplical")
  } else if (h > b+1e-9) {
    x0 = x; h0 = b; z0 = simplr(x0 = x0, y0 = h, x = w$x, y = w$h, xvals = xvals, yvals = hvals, 
                                vfx = Tw, method = "simplical")
    x1 = x; h1 = b0; z1 = simplr(x0 = x1, y0 = h1, x = w$x, y = w$h, xvals = xvals, yvals = hvals, 
                                 vfx = Tw, method = "simplical")
    x2 = x; h2 = b0p; z2 = simplr(x0 = x2, y0 = h2, x = w$x, y = w$h, xvals = xvals, yvals = hvals, 
                                  vfx = Tw, method = "simplical")
  }
  if (z1-z2 <= 1e-9 | z0-z1 <= 1e-9) rise = 0 else {
    deriv = (z0-z1)/(z1-z2)
    if (deriv >= 1) {
      steps = euclid(c(x,h), c(x0,h0))/euclid(c(x0,h0), c(x2,h2))
      rise = (z0-z2)*steps
    } else {
      steps = euclid(c(x,h), c(x0,h0))/euclid(c(x0,h0), c(x1,h1))
      rise = (z0-z1)*(1-deriv^steps)/(1-deriv) - (z0-z1)*deriv
    }
  }
  r = z0 + rise
  if (var=="Tw") r = min(r, 0)
  return(r)
}

interpolater.creater = function(w, theta, method = "simplical", var = "Tw"){
  
  cbar = theta[4]*theta[9]; hbar = theta[5]*(1-theta[9])*theta[10]/(1-theta[10]); lambda = theta[6]
  if (var %in% c("aidcfx", "aidifx")) {
    Tw = if (var == "aidcfx") w$cfx else w$ifx
  } else {
    Tw = w[[var]]
  }
  xvals = sort(unique(w$x)); hvals = sort(unique(w$h))
  xmax = max(xvals); hmax = max(hvals); Tmax = max(Tw)
  xmin = min(xvals); hmin = min(hvals); Tmin = min(Tw)
    
  if (method == "chebyshev") { ### currently not working
    tmat = matrix(TwC, ncol = length(unique(cs$x)))
    cmat = chebcoefs(tmat, degree = length(unique(cs$h))-1)
  }
  
  if (method %in% c("spline", "neuralnet")) {
    obj = (Tw - min(Tw))/(max(Tw) - min(Tw))
    xnorm = (w$x - min(w))/(max(w$x) - min(w$x))
    hnorm = (w$h - min(w$h))/(max(w$h) - min(w$h))
    if (method=="spline") {
      mod = gam(obj ~ te(xnorm, hnorm))
    } else {
      mod = neuralnet(obj ~ xnorm + hnorm, data = data.frame(obj = obj, xnorm = xnorm, hnorm = hnorm), 
                      hidden=11, act.fct = "logistic", linear.output = FALSE, stepmax = 1e5, rep = 5)
    }
  }
  
  fi = function(x,h){
    if (x < -1*lambda) x = -1*lambda
    if (x > xmax+1e-9 | h > hmax+1e-9) {
      r = extrapolator(w, x, h, xvals, hvals, Tw, var)
      return(r)
    } else if (method=="simplical") { 
      r = simplr(x0 = x, y0 = h, x = w$x, y = w$h, xvals = xvals, yvals = hvals, 
                 vfx = Tw, method = "simplical")
    } else if (method=="chebyshev") {
      r = chebpred(x,h, cmat, lb = c(cbar-lambda, 0), ub = c(xmax, hmax))
    } else if (method=="rollmean") {
      xi = findInterval(x, xvals, all.inside = TRUE)
      x1 = xvals[xi]; x2 = xvals[xi+1]
      hi = findInterval(h, hvals, all.inside = TRUE)
      h1 = hvals[hi]; h2 = hvals[hi+1]
      r = mean(Tw[w$x %in% c(x1, x2) & w$h %in% c(h1, h2)])
    } else {
      xn = (x - xmin)/(xmax - xmin)
      hn = (h - hmin)/(hmax - hmin)
      r = predict(mod, newdata = data.frame("xnorm" = xn, "hnorm" = hn))
      r = r*(Tmax - Tmin) + Tmin
    } 
    return(as.numeric(r))
  }
  return(fi)
}

### Function factory for generation household maximization problem (minimizes the negative - which is positive)

hhprob_funkfact <- function(x, h, Valfunc, theta, gqpts){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]*theta[9]; hbar = theta[5]*(1-theta[9])*theta[10]/(1-theta[10])
  lambda = theta[6]; sigma = theta[7]; sigmame = theta[8]; alpha = theta[9]; delta = theta[10] 
  
  draws = exp(great.expectations(sigma, gqpts))
  wts = gqpts$w/sqrt(pi)
  penalty = -1*u(cbar/2, hbar/2, 0, theta) - beta*Valfunc(-lambda*R, hbar/2)
  
  hhprob = function(ci, default = FALSE){
    c = ci[1]; i = ci[2]
    if (is.nan(c)|is.na(c)) {return((c+i+1)*penalty)} else { 
      if (default) {
        payoff = u(cbar, max(h, hbar), 0, theta)
        htplus1 = delta*max(hbar, h)
        xtplus1 = -lambda*R
      } else if (c+i > x+lambda+1e-9) {return((c+i+1)*penalty)} else { ### built in constraint 
        payoff = u(c, h, i, theta)
        htplus1 = delta*(h+i)
        xtplus1 = R*(x - c - i)
      }
      EV = sum(wts*mapply(Valfunc, x = xtplus1 + draws, h = htplus1))
      r = payoff + beta*EV
      return(-r)
    }
  }
  return(hhprob)
}

firstguesser <- function(w, theta){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]*theta[9]; hbar = theta[5]*(1-theta[9])*theta[10]/(1-theta[10])
  lambda = theta[6]; sigma = theta[7]; sigmame = theta[8]; alpha = theta[9]; delta = theta[10]
  hsteady = (1-alpha)*delta/(1-delta)
  guesswrap = function(rowidx){
    x = w$x[rowidx]; h = w$h[rowidx]
    if (h > hsteady) {i = 0; c = x + lambda} else {
      i = min((1-alpha)*(x+lambda), hsteady - h)
      c = x + lambda - i
    }
    u = u(c, h, i, theta)
    return(c(c, i, u, -lambda, 0))
  }
  Tw = mclapply(c(1:nrow(w)), guesswrap, mc.cores = detectCores())
  wnext = cbind(w[,c("x", "h")], data.frame("Tw" = unlist(lapply(Tw, function(x) x[3])),
                                            "cfx" = unlist(lapply(Tw, function(x) x[1])),
                                            "ifx" = unlist(lapply(Tw, function(x) x[2])),
                                            "savings" = unlist(lapply(Tw, function(x) x[4])),
                                            "def" = unlist(lapply(Tw, function(x) x[5]))))
  return(wnext)
}

### Bellman Operator

bellman <- function(w, theta, shockpts = 30, m = 2000){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]*theta[9]; hbar = theta[5]*(1-theta[9])*theta[10]/(1-theta[10])
  lambda = theta[6]; sigma = theta[7]; sigmame = theta[8]; alpha = theta[9]; delta = theta[10]
  gqpts = gaussHermiteData(shockpts)
  Valfunc = interpolater.creater(w, theta, method = "simplical")
  optwrap = function(rowidx){
    x = w$x[rowidx]; h = w$h[rowidx]
    hhprob = hhprob_funkfact(x, h, Valfunc, theta, gqpts)
    
    ## default enforced at constraint
    Vdef = -1*hhprob(c(cbar, 0), default = TRUE)
    if (x + lambda < 1e-9){
      r = c(cbar, 0, Vdef, -lambda, 1)
    } else {
      ### optimization
      sol = directL(fn = hhprob, 
                    lower = c(0,0), upper = c(x+lambda, x+lambda), 
                    control = list(ftol_rel = 1e-6, ftol_abs = 0, xtol_rel = 0, maxeval = m))
      if (x + lambda - sum(sol$par) < 1e-9) { ### double check solutions near boundary - tend to get stuck in local minima
        c0 = min(w$cfx[rowidx], x + lambda); i0 = min(w$ifx[rowidx], x + lambda)
        sol2 = cobyla(fn = hhprob, x0 = c(c0, i0),
                      lower = c(0,0), upper = c(x+lambda, x+lambda), 
                      control = list(ftol_rel = 1e-6, ftol_abs = 0, xtol_rel = 0, maxeval = m))
        if (sol2$value < sol$value) {sol <- sol2}
      }
      if (Vdef > -sol$value) {
        r = c(c(cbar, 0), Vdef, -lambda, 1)} else {
          r = c(sol$par, -sol$value, x - sum(sol$par), 0)} ### default decision
    }
    return(r)
  }
  Tw = mclapply(c(1:nrow(w)), optwrap, mc.cores = detectCores())
  wnext = cbind(w[,c("x", "h")], data.frame("Tw" = unlist(lapply(Tw, function(x) x[3])),
                                            "cfx" = unlist(lapply(Tw, function(x) x[1])),
                                            "ifx" = unlist(lapply(Tw, function(x) x[2])),
                                            "savings" = unlist(lapply(Tw, function(x) x[4])),
                                            "def" = unlist(lapply(Tw, function(x) x[5]))))
  return(wnext)
}

### Howard Policy Iteration
howard <- function(w, theta, shockpts = 30){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]*theta[9]; hbar = theta[5]*(1-theta[9])*theta[10]/(1-theta[10])
  lambda = theta[6]; sigma = theta[7]; sigmame = theta[8]; alpha = theta[9]; delta = theta[10]
  gqpts = gaussHermiteData(shockpts)
  Valfunc = interpolater.creater(w, theta, method = "simplical")
  
  accelerator = function(rowidx){
      x = w$x[rowidx]; h = w$h[rowidx]
      cfx = w$cfx[rowidx]; ifx = w$ifx[rowidx]; def = w$def[rowidx]==1
      hhprob = hhprob_funkfact(x, h, Valfunc, theta, gqpts)
      v1 = -1*hhprob(c(cfx, ifx), default = def)
      vdef = -1*hhprob(c(cfx, ifx), default = TRUE)
      return(max(v1, vdef))
    }
  Twh = mclapply(1:nrow(w), accelerator, mc.cores = detectCores())
  return(unlist(Twh))
}

### VFI
VFI <- function(v0, theta, maxiter = 30, tol = 2.5e-3, howardk = 0, mqp = FALSE, shockpts = c(8,30)){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]*theta[9]; hbar = theta[5]*(1-theta[9])*theta[10]/(1-theta[10])
  lambda = theta[6]; sigma = theta[7]; sigmame = theta[8]; alpha = theta[9]; delta = theta[10]
  w = vector(mode = "list", length = maxiter)
  w[[1]] = v0
  tol1 = FALSE; i = 2; d = 1
  while (i<=maxiter & d > tol) {
    if (i==2) {s = shockpts[1]} else {s = shockpts[2]}
    wnext = bellman(w = w[[i-1]], theta, shockpts = s)
    ## howard
    if (howardk>0) {
      k = vector(mode = "list", length = howardk+1)
      k[[1]] = wnext
      for (j in c(1:howardk)) {
        wnext$Tw = howard(wnext, theta, shockpts = s)
        k[[j+1]] = wnext
      }      
    }
    if (mqp) {
      ### MQP error bounds
      b = beta/(1-beta)*mean(k[[howardk+1]]$Tw - k[[howardk]]$Tw)
      wnext$Tw = wnext$Tw + b
    }
    w[[i]] = wnext
    ### check tol
    d = mean((w[[i]]$cfx - w[[i-1]]$cfx)^2)
    print(d)
    i = i+1
    print(i)
  }
  return(compact(w))
}

### Define GMM function

gmmmomentmatcher <- function(theta, df, v0 = NULL, returnv = FALSE, shockpts = c(8,30), ubm = c(5,5)) {
  print(theta)
  saveRDS(theta, paste(getwd(), "/data/model_output/theta.rds", sep = ""))
  
  ### initial guess
  if (is.null(v0)) {
    statespace = create.statespace(ubm = ubm, theta, method = "equal")
    v0 = firstguesser(statespace, theta)
  } else {v0 = v0[[length(v0)]]}
  
  ### VFI
  V = VFI(v0, theta, shockpts = shockpts)
  saveRDS(V, paste(getwd(), "/data/model_output/V.rds", sep = ""))
  finalV <- V[[length(V)]]
  policycfx <- interpolater.creater(finalV, theta, var = "cfx", method = "rollmean")
  policyifx <- interpolater.creater(finalV, theta, var = "ifx", method = "rollmean")
  policydef <- interpolater.creater(finalV, theta, var = "def", method = "rollmean")
  
  ### data
  momentmat <- momentmatcher(vfx = list(policycfx, policyifx, policydef), t0 = theta, data = df)
  print(sum(colMeans(momentmat)^2))
  if (returnv) {
    return(list(momentmat, V))
  } else {
    return(momentmat)
  }
}

momentmatcher <- function(vfx, t0, data, savings_measure = "liquid"){
  gamma = t0[1]; beta = t0[2]; R = t0[3]; cbar = t0[4]*t0[9]; hbar = t0[5]*(1-t0[9])*t0[10]/(1-t0[10])
  lambda = t0[6]; sigma = t0[7]; sigmame = t0[8]; alpha = t0[9]; delta = t0[10]
  
  set.seed(1646); meshocks <- exp(rnorm(nrow(data), -sigmame^2/2, sigmame))
  data$M_avg <- data$M_avg*meshocks
  
  data <- mutate(data, liquidity = liquidity_hat/M_avg,
                 liquidity_plus = liquidity_plus_hat/M_avg,  
                 food_consumption = food_consumption_hat/M_avg,
                 home_value = home_value_hat/M_avg,
                 home_investment = home_investment_hat/M_avg,
                 quake_aid = quake_aid/M_avg,
                 total_income = total_income_hat/M_avg) %>%
    group_by(hhid) %>%
    mutate(lag_h = lag(home_value, order_by = wave))
  data <- data[complete.cases(data),]
  
  wt = sqrt(data$wt_hh/sum(data$wt_hh))
  x = if (savings_measure=="liquid") {data$liquidity} else {data$liquidity_plus}
  cfx <- mcmapply(vfx[[1]], x, data$lag_h, mc.cores = detectCores())
  ifx <- mcmapply(vfx[[2]], x, data$lag_h, mc.cores = detectCores())
  dfx <- mcmapply(vfx[[3]], x, data$lag_h, mc.cores = detectCores())
  
  ## consumption and investment
  e1 = (cfx - data$food_consumption)/mean(data$food_consumption)
  e2 = (ifx - data$home_investment)/mean(data$home_investment)
  e3 = (dfx - data$default_hat)/mean(data$default_hat)
  ###savings and variance of income
  e4 = ((R-1)*data$cash_savings - data$cap_gains)/mean(data$cap_gains) ### use current cash rather than lag b/c of how survey question is asked
  e5 = (data$loans_made_1yr*annuitycalc(R, 1) + data$loans_made_2yr*annuitycalc(R, 2) + data$loans_made_3yr*annuitycalc(R, 3) - data$loan_payments_recd_ann)/mean(data$loan_payments_recd_ann)
  e6 = (data$loans_taken_1yr*annuitycalc(R, 1) + data$loans_taken_2yr*annuitycalc(R, 2) + data$loans_taken_3yr*annuitycalc(R, 3) - data$amount_owed)/mean(data$amount_owed)
  ### depreciation
  e7 = data$years_ago_built_resid*(data$home_value_resid - (delta-1)*data$years_ago_built_resid)
  ### interactions
  e8 = e1*ihs(x)/(mean(data$food_consumption)*mean(ihs(data$liquidity)))
  e9 = e1*log(data$M_avg)/(mean(data$food_consumption)*mean(log(data$M_avg)))
  e10 = e1*log(data$lag_h+1)/(mean(data$food_consumption)*mean(log(data$lag_h+1)))
  e11 = e2*ihs(x)/(mean(data$home_investment)*mean(ihs(data$liquidity)))
  e12 = e2*log(data$M_avg)/(mean(data$home_investment)*mean(log(data$M_avg)))
  e13 = e2*log(data$lag_h+1)/(mean(data$home_investment)*mean(log(data$lag_h+1)))
  ### variance of income and measurement error
  e14 = (log(data$total_income) + sigma^2/2 + sigmame^2/2)/mean(log(data$total_income))
  e15 = (log(data$total_income)^2 - (sigma^2/2 + sigmame^2/2)^2 - sigma^2 - sigmame^2)/mean(log(data$total_income)^2)
  return(cbind(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15)*wt)
}

partialderivs <- function(idx, tol, v0, theta){
  thetap <- theta
  thetap[idx] <- theta[idx] + tol
  Vp1 <- VFI(v0, thetap, shockpts = c(30,30), tol = 1e-3)
  return(Vp1)
}

momentgrad <- function(theta, df, tol = .005){
  print(theta)
  
  ### initial value and policy functions
  statespace = create.statespace(ubm = c(5,5), theta, method = "equal")
  v0 = firstguesser(statespace, theta)
  V = VFI(v0, theta)
  finalV <- V[[length(V)]]
  policycfx0 <- interpolater.creater(finalV, theta, var = "cfx", method = "rollmean")
  policyifx0 <- interpolater.creater(finalV, theta, var = "ifx", method = "rollmean")
  policydef0 <- interpolater.creater(finalV, theta, var = "def", method = "rollmean")
  momentmat <- momentmatcher(vfx = list(policycfx0, policyifx0, policydef0), t0 = theta, data = df)
  moments0 = colMeans(momentmat)
  
  ### derivative approximations
  partiallist <- lapply(c(1:length(theta)), partialderivs, tol = tol, v0 = finalV, theta = theta)
  
  ddp = matrix(NA, ncol = length(theta), nrow = length(moments0))
  
  ### data
  for (i in c(1:length(theta))){
    vprime <- partiallist[[i]][[length(partiallist[[i]])]]
    tprime = theta; tprime[i] = theta[i] + tol
    policycfxp <- interpolater.creater(vprime, tprime, var = "cfx", method = "rollmean")
    policyifxp <- interpolater.creater(vprime, tprime, var = "ifx", method = "rollmean")
    policydefp <- interpolater.creater(vprime, tprime, var = "def", method = "rollmean")
    momentmatp <- momentmatcher(vfx = list(policycfxp, policyifxp, policydefp), t0 = tprime, data = df)
    momentsp <- colMeans(momentmatp)
    ddp[,i] <- (momentsp - moments0)/tol
  }
  return(ddp)
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
                         d3 = TRUE, aidamt = NULL, ubm = 9, ubh = 9, lbm = NULL, lbh = NULL){
  if (method=="double") {Vfx = interpolater.creater(w[[iter]], theta, var = fillvar, method = "rollmean")
  } else {
    Vfx = interpolater.creater(w[[iter]], theta, var = fillvar, method = method)
  }
  if(is.null(lbm)) lbm = -lambda
  if(is.null(lbh)) lbh = 0
  df.test = crossing(x = seq(lbm, ubm, length.out = 150), h = seq(lbh, ubh, length.out = 150))
  df.test$proj = mapply(Vfx, x = df.test$x, h = df.test$h)
  
  if (fillvar %in% c("aidcfx", "aidifx")) {
    df.test$proj2 = mapply(Vfx, df.test$x + aidamt, df.test$h)
    df.test$efx = df.test$proj2 - df.test$proj
  }
  
  if (method=="double"){
    if (fillvar %in% c("aidcfx", "aidifx")) {
      Vfx2 = interpolater.creater(df.test, theta, var = "efx", method = "spline")
      df.test$efx = mapply(Vfx2, x = df.test$x, h = df.test$h)
    } else {
      Vfx2 = interpolater.creater(df.test, theta, var = "proj", method = "spline")
      df.test$proj = mapply(Vfx2, x = df.test$x, h = df.test$h)
    }
  }
  
  if (fillvar=="Tw") {
    threed = ggplot(df.test) + geom_tile(aes(x = x, y = h, fill = -1*log(-proj))) + scale_fill_viridis_c()
  } else if (fillvar=="cfx") {
    threed = ggplot(df.test) + geom_tile(aes(x = x, y = h, fill = proj)) + scale_fill_viridis_c()
  } else if (fillvar=="ifx") {
    threed = ggplot(df.test) + geom_tile(aes(x = x, y = h, fill = proj)) + scale_fill_viridis_c()
  } else if (fillvar=="savings") {
    threed = ggplot(df.test) + geom_tile(aes(x = x, y = h, fill = proj)) + scale_fill_viridis_c()
  } else {threed = ggplot(df.test) + geom_tile(aes(x = x, y = h, fill = efx)) + scale_fill_viridis_c()}
  
  if (d3) p = plot_gg(threed, multicore=TRUE) else p = threed
  p
}

hhprobtroubleshooter <- function(rowidx, w, crange = c(0,x+lambda), irange = c(0, x+lambda)){
  gqpts = gaussHermiteData(13)
  Valfunc = interpolater.creater(w, theta, method = "simplical")
  x = w$x[rowidx]; h = w$h[rowidx]
  hhprob = hhprob_funkfact(x, h, Valfunc, theta, gqpts)
  cs = seq(crange[1], crange[2], length.out = 40)
  is = seq(irange[1], irange[2], length.out = 40)
  testgrid = crossing(cs, is) %>%
    filter(cs + is <= x + lambda)
  testgrid$tw = sapply(c(1:nrow(testgrid)), function(i) hhprob(c(testgrid$cs[i], testgrid$is[i])))
  ggplot(testgrid) + geom_tile(aes(x = cs, y = is, fill = -1*log(tw))) + scale_fill_viridis_c()
}

constrainttroubleshooter = function(rowidx, w, crange = c(0,x+lambda)){
  gqpts = gaussHermiteData(13)
  Valfunc = interpolater.creater(w, theta, method = "simplical")
  x = w$x[rowidx]; h = w$h[rowidx]
  hhprob = hhprob_funkfact(x, h, Valfunc, theta, gqpts)
  cs = seq(crange[1], crange[2], length.out = 40)
  hs = sapply(cs, function(c) hhprob(c(c, x+lambda-c)))
  ggplot() + geom_line(aes(x = cs, y = hs))
}

