incenter <- function(cbar, lbi, x, lambda) {
  return(c((cbar*2 + x+lambda-lbi)/3, (lbi*2 + x+lambda-cbar)/3))
}

euclid <- function(a, b) {sqrt(sum((a - b)^2))}

ihs <- function(x){log(x + sqrt(x ^ 2 + 1))}

### utility function
u <- function(c, h, i, theta) {
  alpha = theta[8]; gamma = theta[1]
  cd = c^alpha * (h + i)^(1-alpha)
  crra = cd^(1-gamma)/(1-gamma)
  return(crra)
}

### Creates points for gaussian quadrature
great.expectations <- function(sigma, gqpts){
  mu = sqrt(2)*sigma*gqpts$x - sigma^2/2
  return(mu)
}

create.statespace = function(ubm = c(5,5), theta, method = "log"){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  
  if (method=="chebyshev") {
    x = chebnodes(xn, lb = -lambda, ub = ubm[1])
    h = chebnodes(hn, lb = 0, ub = ubm[2])
  } else if (method=="log") { 
    x =  exp(seq(0, log(ubm[1] + lambda +1), length.out = hn)) - (lambda+1)
    h =  exp(seq(0, log(ubm[2]+1), length.out = hn))-1
  } else if (method=="uneven") {
    x = seq(-lambda, ubm[1], length.out = xn)
    h = unique(c(seq(0, ubm[2]/2, length.out = round(2*hn/3,0)), seq(ubm[2]/2, ubm[2], length.out = round(1*hn/3,0)+1)))
  } else {
    x = seq(-lambda, ubm[1], length.out = xn)
    h = seq(0, ubm[2], length.out = hn)
  }
  css = crossing(x,h)
  return(css)
}

###3D interpolate
extrapolator <- function(w, x, h, Tw, var) {
  xvals = unique(w$x); hvals = unique(w$h)
  a = max(xvals); b = max(hvals); 
  a0 = sort(xvals,partial=length(xvals)-1)[length(xvals)-1];  b0 = sort(hvals,partial=length(hvals)-1)[length(hvals)-1]
  a0p = sort(xvals,partial=length(xvals)-3)[length(xvals)-3];  b0p = sort(hvals,partial=length(hvals)-3)[length(hvals)-3]
  if (x > a+1e-10 & h > b+1e-10) {  
    bslp = (h-b)/(x-a)
    x0 = a; h0 = b; z0 = simplr(x0 = a, y0 = b, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
    if (bslp >= b/a) {
      x1 = (b0-b)/bslp + a; h1 = b0; z1 = simplr(x0 = x1, y0 = h1, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
      x2 = (b0p-b)/bslp + a; h2 = b0p; z2 = simplr(x0 = x2, y0 = h2, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
    } else {
      x1 = a0; h1 = (a0-a)*bslp + b; z1 = simplr(x0 = x1, y0 = h1, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
      x2 = a0p; h2 = (a0p-a)*bslp + b; z2 = simplr(x0 = x2, y0 = h2, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
    }
  } else if (x > a+1e-10) {
    x0 = a; h0 = h; z0 = simplr(x0 = x0, y0 = h, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
    x1 = a0; h1 = h; z1 = simplr(x0 = x1, y0 = h1, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
    x2 = a0p; h2 = h; z2 = simplr(x0 = x2, y0 = h2, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
  } else if (h > b+1e-10) {
    x0 = x; h0 = b; z0 = simplr(x0 = x0, y0 = h, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
    x1 = x; h1 = b0; z1 = simplr(x0 = x1, y0 = h1, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
    x2 = x; h2 = b0p; z2 = simplr(x0 = x2, y0 = h2, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
  }
  if (z1-z2 <= 1e-10 | z0-z1 <= 1e-10) rise = 0 else {
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
  
  cbar = theta[4]; hbar = theta[5]; lambda = theta[6]
  if (var=="proj") Tw = w$proj else if (var=="efx") Tw = w$efx else if (var=="Tw") Tw = w$Tw else if (var %in% c("cfx", "aidcfx")) Tw = w$cfx else if (var %in% c("ifx", "aidifx")) Tw = w$ifx
  
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
    if (x > max(w$x)+1e-10 | h > max(w$h)+1e-10) {
      r = extrapolator(w, x, h, Tw, var)
      return(r)
    } else if (x <= cbar-lambda | method=="simplical") {  ### default zone always simplical
      r = simplr(x0 = x, y0 = h, x = w$x, y = w$h, vfx = Tw, method = "simplical", extrapolation.warning = FALSE)
    } else if (method=="chebyshev") {
      r = chebpred(x,h, cmat, lb = c(cbar-lambda, 0), ub = c(max(w$x), max(w$h)))
    } else if (method=="rollmean") {
      xvals = sort(unique(w$x)); hvals = sort(unique(w$h))
      xi = findInterval(x, xvals, all.inside = TRUE)
      x1 = xvals[xi]; x2 = xvals[xi+1]
      hi = findInterval(h, hvals, all.inside = TRUE)
      h1 = hvals[hi]; h2 = hvals[hi+1]
      r = mean(Tw[w$x %in% c(x1, x2) & w$h %in% c(h1, h2)])
    } else {
      xn = (x - min(w$x))/(max(w$x) - min(w$x))
      hn = (h - min(w$h))/(max(w$h) - min(w$h))
      r = predict(mod, newdata = data.frame("xnorm" = xn, "hnorm" = hn))
      r = r*(max(Tw) - min(Tw)) + min(Tw)
    } 
    return(as.numeric(r))
  }
  return(fi)
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

### Function factory for generation household maximization problem (minimizes the negative - which is positive)

hhprob_funkfact <- function(x, h, Valfunc, theta, gqpts){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]; 
  
  draws = exp(great.expectations(sigma, gqpts))
  wts = gqpts$w/sqrt(pi)
  
  hhprob = function(ci, default = FALSE){
    c = ci[1]; i = ci[2]
    if (is.nan(c)|is.na(c)) {return((c+i)*1e10)} else { 
      if (default) {
        payoff = u(cbar, max(h, hbar), 0, theta)
        htplus1 = delta*max(hbar, h)
        xtplus1 = -1*R*lambda
      } else if (c+i>x+lambda) {return((c+i)*1e10)} else { ### built in constraint 
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

hhprob = function(ci){
  c = ci[1]; i = ci[2]
  if (is.nan(c)|is.na(c)|c+i>x-B) {return((c+i)*1e10)} else {  ### built in constraint 
    payoff = u(c, h, i, theta)
    htplus1 = delta*(h+i)
    xtplus1 = R*(x - c - i) 
    mindraw = B+cbar+max(hbar-htplus1, 0) - xtplus1 ## minimum draw to not default
    pdef = plnorm(mindraw, meanlog = -sigma^2/2, sdlog = sigma)
    if (any(draws>mindraw)){
      wts = gqpts$w[draws>mindraw]/sqrt(pi)
      drawsr = draws[draws>mindraw]
      if (pdef!=0){  ### adjust wts and draws for probability of default
        extrawt = 1 - pdef - sum(wts)
        if (extrawt>0) wts = c(extrawt, wts) else wts[1] = extrawt + wts[1]
        extradraw = qlnorm(pdef + wts[1]/2, meanlog = -sigma^2/2, sdlog = sigma)
        if (extrawt>0) drawsr = c(extradraw, drawsr) else drawsr[1] = extradraw
      }
      EVp = mapply(Valfunc, x = xtplus1 + drawsr, h = htplus1)
      EV = pdef*Valfunc(B, max(htplus1, delta*hbar)) + sum(EVp*wts)
    } else {
      EV = pdef*Valfunc(B, max(htplus1, delta*hbar)) + (1-pdef)*Valfunc(xtplus1+mindraw+1.856, htplus1) ### mean excess for lognormal
    }
    r = payoff + beta*EV
    return(-r)}
}

### Bellman Operator

bellman <- function(w, theta, shockpts = 13, m = 2000){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  gqpts = gaussHermiteData(shockpts)
  Valfunc = interpolater.creater(w, theta, method = "simplical")
  optwrap = function(rowidx){
    x = w$x[rowidx]; h = w$h[rowidx]
    hhprob = hhprob_funkfact(x, h, Valfunc, theta, gqpts)
    
    ## default decision
    Vdef = -1*hhprob(c(cbar, 0), default = TRUE)
    optc = alpha*(x+h+lambda); opti = x + lambda - optc ### derived from FOCs to flow
    if (opti<0) {optc = x+lambda; opti = 0}  ## investment constraint
    umax = -1*hhprob(c(optc, opti))
    if (h+opti <=0 | optc <=0 | umax < Vdef) {
      r = c(cbar, 0, Vdef, 1)
    } else {
      ### optimization
      sol = directL(fn = hhprob, 
                    lower = c(0,0), upper = c(x+lambda, x+lambda), 
                    control = list(ftol_rel = 1e-4, ftol_abs = 0, xtol_rel = 0, maxeval = m))
      r = c(sol$par, -sol$value, 0)
    }
    return(r)
  }
  Tw = mclapply(c(1:nrow(w)), optwrap, mc.cores = detectCores())
  wnext = cbind(w[,c("x", "h")], data.frame("Tw" = unlist(lapply(Tw, function(x) x[3])),
                                            "cfx" = unlist(lapply(Tw, function(x) x[1])),
                                            "ifx" = unlist(lapply(Tw, function(x) x[2])),
                                            "def" = unlist(lapply(Tw, function(x) x[4]))))
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
    cfx = w$cfx[rowidx]; ifx = w$ifx[rowidx]; def = w$def[rowidx]
    hhprob = hhprob_funkfact(x, h, Valfunc, theta, gqpts)
    ## default decision
    if (def==1) {
      v = -1*hhprob(c(cbar, 0), default = TRUE) 
      } else {
      v = -1*hhprob(c(cfx, ifx))
      }
    return(v)
  }
  Twh = mclapply(1:nrow(w), accelerator, mc.cores = detectCores())
  return(unlist(Twh))
}

### VFI
VFI <- function(v0, theta, maxiter = 30, tol = 1e-5, howardk = 10, s2 = 51, monotonizer = FALSE){
  gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
  lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
  w = vector(mode = "list", length = maxiter)
  w[[1]] = v0
  tol1 = FALSE; tol2 = FALSE; i = 2; d = 1; s = 13
  while (i<=maxiter & !tol2) {
    wnext = bellman(w = w[[i-1]], theta, shockpts = s)
    if (monotonizer) wnext = monotonizer(wnext, policy = FALSE)
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
    if (tol1) {
      d = mean((w[[i]]$cfx - w[[i-1]]$cfx)^2)*1e-3 ### effective absolute tol of 1e-2 on consumption
    } else {d = mean(((w[[i]]$Tw - w[[i-1]]$Tw)/w[[i-1]]$Tw)^2)}  ### relative tol
    ### if tol is met increase gq points
    if (d < tol) {
      if (tol1) tol2 = TRUE else s = s2; tol1 = TRUE
    }
    print(d)
    i = i+1
    print(i)
  }
  return(compact(w))
}

### Define GMM function

gmmmomentmatcher <- function(theta, df) {
  print(theta)
  saveRDS(theta, paste(getwd(), "/data/model_output/theta.rds", sep = ""))
  
  ### initial guess
  statespace = create.statespace(ubm = c(5,5), theta, method = "log")
  v0 = cbind(statespace, data.frame("Tw" = rep(0, nrow(statespace)),
                                    "cfx" = rep(0, nrow(statespace)),
                                    "ifx" = rep(0, nrow(statespace)),
                                    "def" = rep(0, nrow(statespace))))
  
  ### VFI
  V = VFI(v0, theta)
  saveRDS(V, paste(getwd(), "/data/model_output/V.rds", sep = ""))
  finalV <- V[[length(V)]]
  policycfx <- interpolater.creater(finalV, theta, var = "cfx", method = "rollmean")
  policyifx <- interpolater.creater(finalV, theta, var = "ifx", method = "rollmean")
  
  ### data
  momentmat <- mclapply(c(1:nrow(df)), momentmatcher, vfx = list(policycfx, policyifx), t0 = theta, 
                        data = df[,-c(1:2)], mc.cores = detectCores())
  momentmat <- do.call(rbind, momentmat)
  return(momentmat)
}

momentmatcher <- function(i, vfx, t0, data){
  gamma = t0[1]; beta = t0[2]; R = t0[3]; cbar = t0[4]; hbar = t0[5]
  lambda = t0[6]; sigma = t0[7]; alpha = t0[8]; delta = t0[9]
  df = data[i,]
  wt = sqrt(df$wt_hh/sum(data$wt_hh))
  e1  = vfx[[1]](df$imputed_bufferstock, df$lag_h) - df$food_consumption
  e2  = vfx[[2]](df$imputed_bufferstock, df$lag_h) - df$home_investment
  e3 = df$home_value - delta*(df$lag_h + df$home_investment)
  e4 = df$imputed_bufferstock - df$quake_aid - df$lag_y - R*(df$lag_x - df$lag_c - df$lag_i)
  e5 = log(df$total_income)^2 - sigma^2
  e6 = e1*ihs(df$imputed_bufferstock)
  e7 = e1*log(df$M_avg)
  e8 = e1*log(df$lag_h+1)
  e9 = e2*ihs(df$imputed_bufferstock)
  e10 = e2*log(df$M_avg)
  e11 = e2*log(df$lag_h+1)
  return(c(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11)*wt)
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
  testgrid$tw = sapply(c(1:nrow(testgrid)), function(i) hhprob(c(testgrid$cs[i], testgrid$is[i])))
  ggplot(testgrid) + geom_tile(aes(x = cs, y = is, fill = tw)) + scale_fill_viridis_c()
}

