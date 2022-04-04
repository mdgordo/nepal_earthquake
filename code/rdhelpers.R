### RD functions

plotvar <- function(v, b, df, h=20, ihs=FALSE, span = 1, k = "triangular", weights = TRUE, donut = 0, dist.exclude=NULL) {
  if (!is.null(dist.exclude)) df <- filter(df, designation!=dist.exclude)
  df <- df[, c(v, b, "wt_hh", "border14_segment", "border_segment13", "elevation", "shake_pga", "wave", "quake_losses")]
  df <- filter(df, abs(get(b))<h & complete.cases(df) & abs(get(b))>donut)
  Y <- unlist(df[, v])
  X <- unlist(df[, b])
  if (ihs==TRUE) Y <- ihs(Y)
  if (b=="dist_2_14") {
    segcovs <- model.matrix(~as.factor(df$border14_segment)+as.factor(df$border14_segment)*df$dist_2_14 + as.factor(df$wave))
    segcovs <- segcovs[, !(colnames(segcovs) %in% c("(Intercept)", "df$dist_2_14"))]
    resreg <- lm(Y ~ segcovs)
    resregres <- resreg$residuals
  } else if (b=="dist_2_seg13"){
    segcovs <- model.matrix(~as.factor(df$border_segment13)+as.factor(df$border_segment13)*df$dist_2_seg13 + as.factor(df$wave))
    segcovs <- segcovs[, !(colnames(segcovs) %in% c("(Intercept)", "df$dist_2_seg13"))]
    resreg <- lm(Y ~ segcovs)
    resregres <- resreg$residuals
  } else {
    segcovs <- model.matrix(df$shake_pga + df$quake_losses + as.factor(df$wave))
    segcovs <- segcovs[, !(colnames(segcovs) %in% c("(Intercept)"))]
    resreg <- lm(Y ~ segcovs)
    resregres <- resreg$residuals
  }
  if (weights==TRUE) w <- df$wt_hh else w <- NULL
  r <- rdplot(y = resregres, x = X, c = 0, weights = w, binselect = "esmv", kernel = k, 
              poly = TRUE, p = 1, title = NULL, x.label = b, y.label = v, hide = TRUE, span = span, method = "lm", h = h)
  return(r$rdplot)
}

histfunc <- function(b, df, h=40, dist.exclude=NULL){
  df <- filter(df, abs(get(b))<h)
  if (!is.null(dist.exclude)) df <- filter(df, designation!=dist.exclude)
  X <- unlist(df[, b])
  W <- df$wt_hh
  freqdf <- aggregate(x = list(Freq = W), by = list(X = X), FUN = sum)
  ggplot() +
    geom_histogram(aes(x = X, weight = W), bins = 400) +
    geom_smooth(data = filter(freqdf, X >= 0), aes(x = X, y = Freq)) + 
    geom_smooth(data = filter(freqdf, X < 0), aes(x = X, y = Freq)) + 
    theme_light() + labs(y = "weighted count", x = paste("distance to ", b))
}

regout <- function(v, b, df, h = NULL, b0 = NULL, fuzzy = FALSE, ihs = FALSE, k = "epanechnikov", weights = TRUE, donut = 0, vce = "hc1", dist.exclude=NULL){
  if (donut!=0) df <- filter(df, abs(get(b))>donut)
  if (!is.null(dist.exclude)) df <- filter(df, designation!=dist.exclude)
  if (vce=="nn") c <- NULL else c <- df$strata
  if (fuzzy==TRUE) {f <- df$aid_cumulative_bin} else if (fuzzy==FALSE) {f <- NULL} else 
    if (fuzzy=="inv") {f <- 1-df$aid_cumulative_bin} else if (fuzz=="random") f <- rbinom(nrow(df), 1, mean(df$aid_cumulative_bin)) else {f <- unlist(df[, fuzzy])}
  if (!is.null(h) & is.null(b0)) b0 <- 2*h
  if (b=="dist_2_14") {
    covs <- model.matrix(~as.factor(df$border14_segment)+as.factor(df$border14_segment)*df$dist_2_14 + df$shake_pga + as.factor(df$wave))
    covs <- covs[, !(colnames(covs) %in% c("(Intercept)", "df$dist_2_14"))]
  } else if (b=="dist_2_seg13"){
    covs <- model.matrix(~as.factor(df$border_segment13)+as.factor(df$border_segment13)*df$dist_2_seg13 + df$shake_pga + as.factor(df$wave))
    covs <- covs[, !(colnames(covs) %in% c("(Intercept)", "df$dist_2_seg13"))]
  } else {
    covs <- model.matrix(~df$shake_pga + as.factor(df$wave))
    covs <- covs[, !(colnames(covs) %in% c("(Intercept)"))]}
  X <- unlist(df[, b])
  Y <- unlist(df[, v])
  if (ihs==TRUE) Y <- ihs(Y)
  if (weights==TRUE) w <- df$wt_hh else w <- NULL
  out <- rdrobust(y = Y, x = X, c = 0, fuzzy = f, h = h, b = b0, weights = w, cluster = c, kernel = k, vce = vce, covs = covs)
  return(out)
}

optbw <- function(v, b, df, fuzzy = FALSE, k = "epanechnikov", weights = TRUE, donut = 0, vce = "hc1", dist.exclude=NULL){
  if (donut!=0) df <- filter(df, abs(get(b))>donut)
  if (!is.null(dist.exclude)) df <- filter(df, designation!=dist.exclude)
  if (vce=="nn") c <- NULL else c <- df$strata
  Y <- unlist(df[, v])
  X <- unlist(df[, b])
  if (fuzzy==TRUE) f <- df$aid_cumulative_bin else if (fuzzy==FALSE) f <- NULL else if (fuzzy=="inv") f <- 1-df$aid_cumulative_bin else f <- unlist(df[, fuzzy])
  if (b=="dist_2_14") {
    covs <- model.matrix(~as.factor(df$border14_segment)+as.factor(df$border14_segment)*df$dist_2_14 + df$shake_pga + as.factor(df$wave))
    covs <- covs[, !(colnames(covs) %in% c("(Intercept)", "df$dist_2_14"))]
  } else if (b=="dist_2_seg13"){
    covs <- model.matrix(~as.factor(df$border_segment13)+as.factor(df$border_segment13)*df$dist_2_seg13 + df$shake_pga + as.factor(df$wave))
    covs <- covs[, !(colnames(covs) %in% c("(Intercept)", "df$dist_2_seg13"))]
  } else {
    covs <- model.matrix(~df$shake_pga + as.factor(df$wave))
    covs <- covs[, !(colnames(covs) %in% c("(Intercept)"))]}
  if (weights==TRUE) w <- df$wt_hh else w <- NULL
  bwsct <- rdbwselect(y = Y, x = X, c = 0, weights = w, fuzzy = f, covs = covs,
                      cluster = c, kernel = k, vce = vce, silent = TRUE)
  bws <- data.frame(bwsct$bws)
  bws$var <- v
  return(bws)
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


qplot <- function(qvar, df, plot = TRUE, grid = NULL){
  qcovs <- model.matrix(~as.factor(df$border_segment13)+as.factor(df$border_segment13)*df$dist_2_seg13 + df$shake_pga + as.factor(df$wave))[,-c(1,3)]
  y <- unlist(df[,qvar])
  if (is.null(grid)) grid = quantile(y, seq(.1,.9,.1), na.rm = TRUE)
  qtab <- rdquant(Y = y, x = df$dist_2_seg13, c = 0, fuzzy = df$aid_cumulative_bin, grid, weights = df$wt_hh, cluster = df$strata, 
                  vce = "hc1", covs = qcovs, kernel = "triangular", h = h0, b = b0)
  qtab = filter(unique(qtab), round(rcoefs,2) >=-.1 & round(rcoefs,2) <=1.1)
  qtiles <- quantile(y, seq(.1, .9, .1), na.rm = TRUE)
  
  if (plot) {
    p = ggplot() +
      geom_line(data = qtab, aes(x = yvals, y = rcoefs, group = as.factor(treat), color = as.factor(treat))) +
      geom_errorbar(data = qtab, aes(x = yvals, y = rcoefs, group = as.factor(treat), color = as.factor(treat), ymin = rlower, ymax = rupper, 
                                   fill = as.factor(treat)), alpha = .3) + 
      guides(color = "none") +
      scale_color_manual(values = c("grey18", "deepskyblue2")) +
      scale_fill_manual(values = c("grey18", "deepskyblue2")) +
      geom_rug(aes(x = qtiles), color = "red") + xlim(NA, max(grid)) +
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


