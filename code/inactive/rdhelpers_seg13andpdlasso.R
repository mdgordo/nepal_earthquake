### RD functions
robustses <- function(x) {
  sqrt(diag(vcovHC(x, cluster=~strata, type = "HC1")))
}

covmatmaker <- function(b, df, vars){
  fulllist = c("wave2", "wave3", "logloss", "shake_pga", "high_caste", "hhmembers", "gorkha_hh", 
               "class5", "class10", "age_hh", "slope", "elevation", "aspect", "gorkha_loss_ever",
               "time_to_market", "time_to_health", "time_to_bank", "time_to_school")
  dfx <- df[,fulllist]
  form = as.formula(paste("~.+0", sep = ""))
  covs <- model.matrix(form, dfx)
  if (vars=="none") {covs = NULL} else if (vars == "all") {covs <- covs} else {
    covs = covs[,colnames(covs) %in% vars]
  }
  ### add segment
  if (b=="dist_2_14") {
    segvars = model.matrix(~seg14_2 + seg14_3 + seg14_2dist + seg14_3dist + seg14_2distpos + seg14_3distpos + 0, data = df)
  } else if (b=="dist_2_seg13"){
    segvars = model.matrix(~seg13_2 + seg13dist + seg13distpos + 0, data = df)
  } else {segvars = NULL}
  covs = cbind(covs, segvars)
  return(covs)
}

dfprep <- function(df, donut, b, dist.exclude, hpop = Inf){
  df <- filter(df, abs(get(b))>donut, !is.na(hhmembers), !is.na(age_hh), !is.na(high_caste), !is.na(class5), !is.na(age_hh), !is.na(gorkha_hh), 
                   !is.na(time_to_market), !is.na(time_to_health), !is.na(time_to_bank), !is.na(time_to_school))
  if (!is.null(dist.exclude)) df <- filter(df, designation!=dist.exclude & abs(get(b))<hpop)
  return(df)
}

vectorprep <- function(df, v, b, vars, vce, fuzzy, ihs, weights){
  if (fuzzy==TRUE) {f <- df$aid_cumulative_bin} else if (fuzzy==FALSE) {f <- NULL} else if (fuzzy=="inv") {f <- 1-df$aid_cumulative_bin} else {f <- unlist(df[, fuzzy])}
  Y <- unlist(df[, v])
  X <- unlist(df[, b])
  if (ihs==TRUE) Y = log(Y+1)
  if (weights==TRUE) w <- df$wt_hh else w <- NULL
  return(list(f, Y, X, w))
}

pdlvarselect <- function(v, maxiter = 10, tol = 1, df, dist.exclude, donut, k, vce, ihs, fuzzy = TRUE){
  ### only works for dist2seg13 at the moment and triangular kernel
  i = 0
  bandinit <- optbw(v, b = "dist_2_seg13", df = df, fuzzy = fuzzy, donut = donut, dist.exclude = dist.exclude, 
                    k = k, vce = vce, ihs = ihs, vars = "all")
  Y = unlist(df[, v])
  if (ihs==TRUE) df$Y = log(Y+1) else df$Y = Y
  hlast <- bandinit[1,1]
  
  fulllist = c("wave2", "wave3", "logloss", "shake_pga", "high_caste", "hhmembers", "gorkha_hh", 
               "class5", "class10", "age_hh", "slope", "elevation", "aspect",
               "time_to_market", "time_to_health", "time_to_bank", "time_to_school")
  
  while (i < maxiter){
    h0 = hlast[length(hlast)]
    
    df.lasso = filter(df, abs(dist_2_seg13)<h0) %>%
      mutate(kwt = wt_hh*(1-abs(dist_2_seg13)/h0))
    
    reg1 <- lm(Y ~ dist_2_seg13 + dist_pos13 + seg13_2 + seg13dist + seg13distpos,
               data = df.lasso, weights = kwt, subset = abs(dist_2_seg13)<h0)
    
    reg2 <- lm(aid_cumulative_bin ~ dist_2_seg13 + dist_pos13 + seg13_2 + seg13dist + seg13distpos,
               data = df.lasso, weights = kwt, subset = abs(dist_2_seg13)<h0)
    
    reg3 <- lm(dist_14 ~ dist_2_seg13 + dist_pos13 + seg13_2 + seg13dist + seg13distpos,
               data = df.lasso, weights = kwt, subset = abs(dist_2_seg13)<h0)
    
    df.lasso$yresid <-reg1$residuals*sqrt(df.lasso$kwt)
    df.lasso$xresid <-reg2$residuals*sqrt(df.lasso$kwt)
    df.lasso$zresid <-reg3$residuals*sqrt(df.lasso$kwt)
    
    dfx <- df.lasso[,fulllist]
    #c = paste0("I(", colnames(dfx), "^2)", collapse = "+")
    form = as.formula(paste("~.",#^2+", #c, 
                            "+0", sep = ""))
    covmat <- model.matrix(form, dfx)*sqrt(df.lasso$kwt)
    
    pdl1 <- rlasso(y = df.lasso$yresid, x = covmat, post = FALSE, intercept = FALSE)
    pdl2 <- rlasso(y = df.lasso$xresid, x = covmat, post = FALSE, intercept = FALSE)
    pdl3 <- rlasso(y = df.lasso$zresid, x = covmat, post = FALSE, intercept = FALSE)
    
    vars <- colnames(pdl1$model)[pdl1$index | pdl2$index | pdl3$index]; vars <- vars[c(2:length(vars))]
    bandinit <- optbw(v, b = "dist_2_seg13", df = df, fuzzy = fuzzy, dist.exclude = dist.exclude, 
                      k = k, vce = vce, ihs = ihs, vars = vars)
    h0 <- bandinit[1,1]
    if (h0 %in% hlast) i = maxiter else i = i+1; hlast = c(hlast, h0)
  }
  return(vars)
}

optbw <- function(v, b, df, fuzzy = FALSE, k = "triangular", weights = TRUE, donut = 0, 
                  vce = NULL, dist.exclude=NULL, vars = NULL, ihs = FALSE, bwscrit = "mserd") {
  df <- dfprep(df, donut, b, dist.exclude)
  vs = vectorprep(df, v, b, vars, vce, fuzzy, ihs, weights)
  f = vs[[1]]; Y = vs[[2]]; X = vs[[3]]; w = vs[[4]]
  if (vars=="opt") {
    vopt <- pdlvarselect(v, maxiter = 10, tol = 1, df, dist.exclude, donut, k, vce, ihs)
    covs <- covmatmaker(b, df, vopt)
  } else {covs <- covmatmaker(b, df, vars)}
  bwsct <- rdbwselect(y = Y, x = X, c = 0, weights = w, fuzzy = f, covs = covs, bwselect = bwscrit,
                      kernel = k, vce = vce, silent = TRUE)
  bws <- data.frame(bwsct$bws)
  return(bws)
}

regout <- function(v, b, df, h = NULL, b0 = NULL, fuzzy = FALSE, k = "triangular", weights = TRUE, donut = 0, 
                   vce = NULL, dist.exclude=NULL, vars = "none", ihs = FALSE, poly = 1){
  df <- dfprep(df, donut, b, dist.exclude)
  vs = vectorprep(df, v, b, vars, vce, fuzzy, ihs, weights)
  f = vs[[1]]; Y = vs[[2]]; X = vs[[3]]; w = vs[[4]]
  if (vars=="opt") {
    vopt <- pdlvarselect(v, maxiter = 10, tol = 1, df, dist.exclude, donut, k, vce, ihs)
    covs <- covmatmaker(b, df, vopt)
  } else {covs <- covmatmaker(b, df, vars)}
  out <- rdrobust(y = Y, x = X, c = 0, fuzzy = f, h = h, b = b0, weights = w, kernel = k,
                  p = poly, vce = vce, covs = covs)
  return(out)
}

plotvar <- function(v, b, df, h=NULL, ihs=FALSE, span = 1, k = "triangular", weights = TRUE, 
                    vce = NULL, donut = 0, dist.exclude=NULL, vars = "none", p = 1) {
  df <- dfprep(df, donut, b, dist.exclude)
  vs = vectorprep(df, v, b, vars, vce = vce, fuzzy = FALSE, ihs, weights)
  f = vs[[1]]; Y = vs[[2]]; X = vs[[3]]; w = vs[[4]]
  if (vars=="opt") {
    vopt <- pdlvarselect(v, maxiter = 10, tol = 1, df, dist.exclude, donut, k, vce, ihs)
    covs <- covmatmaker(b, df, vopt)
  } else {covs <- covmatmaker(b, df, vars)}
  if (!is.null(covs)){
    resreg1 <- lm(Y ~ covs, subset = X>=0)
    resreg2 <- lm(Y ~ covs, subset = X<0)
    resregres <- rep(NA, length(X))
    resregres[X>=0] <- resreg1$residuals
    resregres[X<0] <- resreg2$residuals
    r <- rdplot(y = resregres, x = X, c = 0, weights = w, binselect = "esmv", kernel = k, 
                poly = TRUE, p = p, title = NULL, x.label = b, y.label = v, hide = TRUE, span = span, method = "lm", h = h)
  } else {
    r <- rdplot(y = Y, x = X, c = 0, weights = w, binselect = "esmv", kernel = k, 
                poly = TRUE, p = p, title = NULL, x.label = b, y.label = v, hide = TRUE, span = span, method = "lm", h = h)
  }
  return(r$rdplot)
}

histfunc <- function(b, df, h=40, dist.exclude=NULL){
  df <- dfprep(df, 0, b, dist.exclude)
  X <- unlist(df[, b])
  W <- df$wt_hh
  freqdf <- aggregate(x = list(Freq = W), by = list(X = X), FUN = sum)
  ggplot() +
    geom_histogram(aes(x = X, weight = W), bins = 400) +
    geom_smooth(data = filter(freqdf, X >= 0), aes(x = X, y = Freq)) + 
    geom_smooth(data = filter(freqdf, X < 0), aes(x = X, y = Freq)) + 
    theme_light() + labs(y = "weighted count", x = paste("distance to ", b))
}

rdgazer <- function(rdlist, dvlabs = NULL, xlines = NULL, se_r = "Robust", type = "text", ...){
  dummymods <- list(); coef <- list(); se <- list(); bw <- c(); nobs <- c(); untreatedmean <- c()
  for (i in 1:length(rdlist)) {
    dummymods[[i]] <- lm(rdlist[[i]]$Y ~ rdlist[[i]]$X)
    coef[[i]] <- c(0, rdlist[[i]]$coef["Conventional",])
    se[[i]] <- c(1, rdlist[[i]]$se[se_r,])
    bw[i] <- round(rdlist[[i]]$bws[1,1],2)
    nobs[i] <- sum(rdlist[[i]]$N_h)
  }
  s <- stargazer(dummymods, type = type, coef = coef, se = se, column.labels = dvlabs,
                 omit.stat = "all", digits = 2, df = FALSE, omit = c("Constant"), 
                 covariate.labels = c("Treatment"),
                 add.lines = c(list(c("N", nobs),
                                    c("Bandwidth", bw)), xlines), ...)
}

cutvars <- function(cutpt, var, df){
  ct = ifelse(df[, var] <= cutpt & df$aid_cumulative_bin==1, 1, 0)
  cu = ifelse(df[, var] <= cutpt & df$aid_cumulative_bin==0, 1, 0)
  c = ifelse(df[, var] <= cutpt, 1, 0)
  df <- cbind(ct, cu, c)
  colnames(df) <- paste(c("t", "u", "c"), cutpt, sep = "")
  return(df)
}

rdquant <- function(Y, x, fuzzy = NULL, grid = quantile(Y, seq(.1,.9,.1), na.rm = TRUE), c = 0, ...) {
  if (is.null(fuzzy)) fuzzy <- ifelse(x>c, 1, 0)
  fuzzy0 <- 1 - fuzzy
  coefs <- function(qinv, f){
    y1d <- as.numeric(Y<=qinv)*f
    rd1out <- rdrobust(y = y1d, x = x, fuzzy = f, c = c, ...)
    return(rd1out)
  }
  mods1 <- lapply(grid, FUN = coefs, f = fuzzy)
  mods0 <- lapply(grid, FUN = coefs, f = fuzzy0)
  coefs1 <- as.numeric(lapply(mods1, function(x) x$coef["Conventional",]))
  coefs0 <- as.numeric(lapply(mods0, function(x) x$coef["Conventional",]))
  ses1 <- as.numeric(lapply(mods1, function(x) x$se["Robust",]))
  ses0 <- as.numeric(lapply(mods0, function(x) x$se["Robust",]))
  r1 <- Rearrangement::rearrangement(x = data.frame(grid), y = coefs1)
  r0 <- Rearrangement::rearrangement(x = data.frame(grid), y = coefs0)
  ci1 <- list("x" = grid, "y" = coefs1, "sortedx" = grid, "Lower" = coefs1 - ses1*1.96, "Upper" = coefs1 + ses1*1.96, "cef" = coefs1)
  attr(ci1, "class") <- "conint"
  rci1 <- Rearrangement::rconint(ci1)
  ci0 <- list("x" = grid, "y" = coefs0, "sortedx" = grid, "Lower" = coefs0 - ses0*1.96, "Upper" = coefs0 + ses0*1.96, "cef" = coefs0)
  attr(ci0, "class") <- "conint"
  rci0 <- Rearrangement::rconint(ci0)
  pdfdf <- data.frame("yvals" = rep(grid, 2),
                      "coefs" = c(coefs1, coefs0),
                      "se" = c(ses1, ses0),
                      "rcoefs" = c(r1, r0),
                      "rlower" = c(rci1$Lower, rci0$Lower),
                      "rupper" = c(rci1$Upper, rci0$Upper),
                      "treat" = c(rep(1, length(grid)), rep(0, length(grid))))
  return(pdfdf)
}


qplot <- function(qvar, df, plot = TRUE, grid = NULL, vars = NULL, k = "triangular", h = NULL, b = NULL, 
                  dist.exclude = "none", donut = 0, poly = 1, weights = TRUE, ihs = FALSE, vce = "hc1", fuzzy = FALSE){
  ### only works for segment 13
  df <- dfprep(df, donut, "dist_2_seg13", dist.exclude)
  vs = vectorprep(df, qvar, "dist_2_seg13", vars, vce, fuzzy, ihs, weights)
  f = vs[[1]]; Y = vs[[2]]; X = vs[[3]]; w = vs[[4]]
  if (vars=="opt") {
    vopt <- pdlvarselect(v, maxiter = 10, tol = 1, df, dist.exclude, donut, k, vce, ihs)
    covs <- covmatmaker("dist_2_seg13", df, vopt)
  } else {covs <- covmatmaker("dist_2_seg13", df, vars)}
  if (is.null(grid)) grid = quantile(Y, seq(.1,.9,.1), na.rm = TRUE)
  qtab <- rdquant(Y = Y, x = X, c = 0, fuzzy = f, grid, weights = w, cluster = c, 
                  vce = vce, covs = covs, kernel = k, h = h, b = b, p = poly)
  qtab = filter(unique(qtab), round(rcoefs,2) >=-.1 & round(rcoefs,2) <=1.1)
  
  if (plot) {
    p = ggplot() +
      geom_line(data = qtab, aes(x = yvals, y = rcoefs, group = as.factor(treat), color = as.factor(treat))) +
      geom_errorbar(data = qtab, aes(x = yvals, y = rcoefs, group = as.factor(treat), color = as.factor(treat), ymin = rlower, ymax = rupper, 
                                   fill = as.factor(treat)), alpha = .3) + 
      scale_color_manual(values = c("grey18", "deepskyblue2")) +
      scale_fill_manual(values = c("grey18", "deepskyblue2")) +
      geom_rug(aes(x = grid), color = "red") + xlim(NA, max(grid)) +
      theme_bw() + labs(x = qvar, y = "P(X < x)", fill = "Received Aid", color = "Received Aid")
    p
  } else {
    qtab
  }
}

earthmovers <- function(data, idx, cutpts){
  df = data[idx,]
  cutvars = paste("c", cutpts, sep = "")
  rdcc = lapply(cutvars, regout, b = "dist_2_seg13", df = df, h = h0, b0 = b0, k = "triangular", vce = "hc1", fuzzy = TRUE)
  dF = unlist(lapply(rdcc, function(x) x$coef[1,1]))
  intfun = approxfun(x = cutpts, y = dF)
  emd = integrate(intfun, min(cutpts), max(cutpts), stop.on.error = FALSE)
  return(emd$value)
}

rearranger <- function(q, y) {
  r <- Rearrangement::rearrangement(x = as.data.frame(q), y)
  r[r<0] <- 0
  r[r>1] <- 1
  return(as.vector(r))
}
