### RD functions
robustses <- function(x) {
  sqrt(diag(vcovHC(x, cluster=~strata, type = "HC1")))
}

covmatmaker <- function(df, vars){
  fulllist = c("wave2", "wave3", "shake_pga", "high_caste", "hhmembers", "gorkha_hh", 
               "class5", "class10", "age_hh", "slope", "elevation", "aspect", "gorkha_loss_ever",
               "time_to_market", "time_to_health", "time_to_bank", "time_to_school")
  dfx <- df[,fulllist]
  covs <- model.matrix(~.+0, dfx)
  if (vars=="none") {covs = NULL} else if (vars == "all") {covs <- covs} else {
    covs = covs[,colnames(covs) %in% vars]
  }
  return(covs)
}

dfprep <- function(df, donut, b, dist.exclude, hpop = Inf){
  df <- filter(df, abs(get(b))>donut, !is.na(hhmembers), !is.na(age_hh), !is.na(high_caste), !is.na(class5), !is.na(age_hh), !is.na(gorkha_hh), 
               !is.na(time_to_market), !is.na(time_to_health), !is.na(time_to_bank), !is.na(time_to_school))
  if (!is.null(dist.exclude)) df <- filter(df, designation!=dist.exclude & abs(get(b))<hpop)
  return(df)
}

vectorprep <- function(df, v, b, vce, fuzzy, ihs, weights){
  if (fuzzy!=FALSE & fuzzy != "inv") {f <- unlist(df[, fuzzy])} else if (fuzzy=="inv") {f <- 1 - df$reconstruction_aid_bin} else {f <- NULL}
  Y <- unlist(df[, v])
  X <- unlist(df[, b])
  if (ihs==TRUE) Y = log(Y+1)
  if (weights==TRUE) w <- df$wt_hh else w <- NULL
  return(list(f, Y, X, w))
}

optbw <- function(v, b, df, fuzzy = FALSE, k = "triangular", weights = TRUE, donut = 0, 
                  vce = NULL, dist.exclude=NULL, vars = NULL, ihs = FALSE) {
  df <- dfprep(df, donut, b, dist.exclude)
  vs = vectorprep(df, v, b, vce, fuzzy, ihs, weights)
  f = vs[[1]]; Y = vs[[2]]; X = vs[[3]]; w = vs[[4]]
  covs <- covmatmaker(df, vars)
  bwsct <- rdbwselect(y = Y, x = X, c = 0, weights = w, fuzzy = f, covs = covs,
                      kernel = k, vce = vce, silent = TRUE)
  bws <- data.frame(bwsct$bws)
  return(bws)
}

regout <- function(v, b, df, h = NULL, b0 = NULL, fuzzy = FALSE, k = "triangular", weights = TRUE, donut = 0, 
                   vce = NULL, dist.exclude=NULL, vars = "none", ihs = FALSE, poly = 1){
  df <- dfprep(df, donut, b, dist.exclude)
  vs = vectorprep(df, v, b, vce, fuzzy, ihs, weights)
  f = vs[[1]]; Y = vs[[2]]; X = vs[[3]]; w = vs[[4]]
  covs <- covmatmaker(df, vars)
  out <- rdrobust(y = Y, x = X, c = 0, fuzzy = f, h = h, b = b0, weights = w, kernel = k,
                  p = poly, vce = vce, covs = covs)
  return(out)
}

plotvar <- function(v, b, df, h=NULL, ihs=FALSE, span = 1, k = "triangular", weights = TRUE, 
                    vce = NULL, donut = 0, dist.exclude=NULL, vars = "none", p = 1, poly = TRUE, method = "lm",
                    residualizer = FALSE, showbw = FALSE) {
  df <- dfprep(df, donut, b, dist.exclude)
  vs = vectorprep(df, v, b, vce = vce, fuzzy = FALSE, ihs, weights)
  Y = vs[[2]]; X = vs[[3]]; w = vs[[4]]
  if (residualizer) {
    df$Y <- Y
    resreg <- lm(Y ~ shake_pga + gorkha_loss_ever + slope + high_caste + time_to_market + time_to_health,
                  data = df)
    Y <- resreg$residuals
  }
  r <- rdplot(y = Y, x = X, c = 0, weights = w, binselect = "esmv", kernel = k, 
              p = p, title = NULL, x.label = b, y.label = v, hide = TRUE, span = span, 
              method = method, h = h, poly = poly)
  r <- r$rdplot
  return(r)
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

