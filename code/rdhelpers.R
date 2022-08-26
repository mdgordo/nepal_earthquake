### RD functions
robustses <- function(x) {
  sqrt(diag(vcovHC(x, cluster=~strata, type = "HC1")))
}

covmatmaker <- function(b, df, vars){
    if (vars=="all") {vars = c("shake_pga", "as.factor(wave)2", "as.factor(wave)3", "high_caste", "gorkha_hh", "log(gorkha_loss_ever+1)",
                                "class5", "class10", "age_hh", "I(age_hh^2)", "slope", "elevation", "aspect",
                                "time_to_market", "time_to_health", "time_to_bank", "time_to_school")
    } else if (vars=="none") {vars = NULL}
    if (b=="dist_2_14") {
      covs <- model.matrix(~ as.factor(border14_segment) + as.factor(border14_segment):dist_2_14 +
                             shake_pga + as.factor(wave) + high_caste + gorkha_hh + log(gorkha_loss_ever+1) +
                             class5 + class10 + age_hh + I(age_hh^2) + slope + elevation + aspect +
                             time_to_market + time_to_health + time_to_bank + time_to_school, 
                           data = df)
      covs <- covs[, colnames(covs) %in% c("as.factor(border14_segment)2", "as.factor(border14_segment)3",
                                           "dist_2_14", "as.factor(border14_segment)2:dist_2_14",
                                           "as.factor(border14_segment)3:dist_2_14", vars)]
    } else if (b=="dist_2_seg13"){
      covs <- model.matrix(~ as.factor(border_segment13) + as.factor(border_segment13):dist_2_14 +
                             shake_pga + as.factor(wave) + high_caste + gorkha_hh + log(gorkha_loss_ever+1) +
                             class5 + class10 + age_hh + I(age_hh^2) + slope + elevation + aspect +
                             time_to_market + time_to_health + time_to_bank + time_to_school, 
                           data = df)
      covs <- covs[, colnames(covs) %in% c("as.factor(border_segment13)2", "dist_2_14", 
                                           "as.factor(border_segment13)2:dist_2_14", vars)]
    } else {
      covs <- model.matrix(~ shake_pga + as.factor(wave) + high_caste + gorkha_hh + log(gorkha_loss_ever+1) +
                             class5 + class10 + age_hh + I(age_hh^2) + slope + elevation + aspect +
                             time_to_market + time_to_health + time_to_bank + time_to_school, 
                           data = df)
      covs <- covs[, colnames(covs) %in% vars]
    } 
}

optbw <- function(v, b, df, fuzzy = FALSE, k = "epanechnikov", weights = TRUE, donut = 0, 
                  vce = "hc1", dist.exclude=NULL, vars = NULL, ihs = FALSE) {
  if (donut!=0) df <- filter(df, abs(get(b))>donut)
  if (!is.null(dist.exclude)) df <- filter(df, designation!=dist.exclude)
  if (!is.null(!is.null(vars) | vars=="none")) df <- filter(df, !is.na(high_caste), !is.na(class5), !is.na(age_hh), !is.na(gorkha_hh), !is.na(time_to_market), !is.na(time_to_health))
  if (vce=="nn") c <- NULL else c <- df$strata
  if (fuzzy==TRUE) {f <- df$aid_cumulative_bin} else if (fuzzy==FALSE) {f <- NULL} else if (fuzzy=="inv") {f <- 1-df$aid_cumulative_bin} else if (fuzzy=="random") f <- rbinom(nrow(df), 1, mean(df$aid_cumulative_bin)) else {f <- unlist(df[, fuzzy])}
  covs <- covmatmaker(b, df, vars)
  Y <- unlist(df[, v])
  X <- unlist(df[, b])
  if (ihs==TRUE) Y = log(Y+1)
  if (weights==TRUE) w <- df$wt_hh else w <- NULL
  bwsct <- rdbwselect(y = Y, x = X, c = 0, weights = w, fuzzy = f, covs = covs,
                      cluster = c, kernel = k, vce = vce, silent = TRUE)
  bws <- data.frame(bwsct$bws)
  bws$var <- v
  return(bws)
}

pdlvarselect <- function(y, maxiter = 10, tol = 1){
  vars = "all"; vlast = "none"; i = 0
  bandinit <- optbw(y, b = "dist_2_seg13", df = df.hh, fuzzy = TRUE, dist.exclude = de, k = kern, vce = serobust, ihs = takelogs, vars = vars)
  hlast <- bandinit[1,1]
  while (i < maxiter){
    h0 = hlast[length(hlast)]
    
    df.lasso <- filter(df.hh, designation!="none", !is.na(age_hh), !is.na(class5), 
                       abs(dist_2_seg13)<h0) %>%
      mutate(dist_pos = if_else(dist_2_seg13>0, dist_2_seg13, 0),
             dist_neg = if_else(dist_2_seg13<0, dist_2_seg13, 0),
             seg2 = if_else(border_segment13==2, 1, 0),
             seg2dist = if_else(border_segment13==2, dist_2_seg13, 0),
             kwt = wt_hh*(1-abs(dist_2_seg13)/h0))
    
    reg1 <- lm(log(get(y)+1) ~ dist_pos + dist_neg + seg2 + seg2dist,
               data = df.lasso, weights = kwt, subset = abs(dist_2_seg13)<h0)
    
    reg2 <- lm(aid_cumulative_bin ~ dist_pos + dist_neg + seg2 + seg2dist,
               data = df.lasso, weights = kwt, subset = abs(dist_2_seg13)<h0)
    
    df.lasso$yresid <-reg1$residuals*sqrt(df.lasso$kwt)
    df.lasso$xresid <-reg2$residuals*sqrt(df.lasso$kwt)
    covmat <- model.matrix(~ shake_pga + as.factor(wave) + high_caste + gorkha_hh + log(gorkha_loss_ever+1) +
                             class5 + class10 + age_hh + I(age_hh^2) + slope + elevation + aspect +
                             time_to_market + time_to_health + time_to_bank + time_to_school,
                           data = df.lasso)*sqrt(df.lasso$kwt)
    
    
    pdl1 <- rlasso(y = df.lasso$yresid, x = covmat, post = FALSE, intercept = FALSE)
    pdl2 <- rlasso(y = df.lasso$xresid, x = covmat, post = FALSE, intercept = FALSE)
    
    vars <- colnames(pdl1$model)[pdl1$index | pdl2$index]; vars <- vars[c(2:length(vars))]
    bandinit <- optbw(y, b = "dist_2_seg13", df = df.hh, fuzzy = TRUE, dist.exclude = de, k = kern, vce = serobust, ihs = takelogs, vars = vars)
    h0 <- bandinit[1,1]
    if (h0 %in% hlast) i = maxiter else i = i+1; hlast = c(hlast, h0)
  }
  return(vars)
}

regout <- function(v, b, df, h = NULL, b0 = NULL, fuzzy = FALSE, ihs = FALSE, 
                   k = "triangular", weights = TRUE, donut = 0, poly = 1, vce = "hc1", dist.exclude=NULL, vars = NULL){
  if (donut!=0) df <- filter(df, abs(get(b))>donut)
  if (!is.null(dist.exclude)) df <- filter(df, designation!=dist.exclude)
  if (!is.null(vars) | vars=="none") df <- filter(df, !is.na(high_caste), !is.na(class5), !is.na(age_hh), !is.na(gorkha_hh), !is.na(time_to_market), !is.na(time_to_health))
  if (vce=="nn") c <- NULL else c <- df$strata
  if (fuzzy==TRUE) {f <- df$aid_cumulative_bin} else if (fuzzy==FALSE) {f <- NULL} else if (fuzzy=="inv") {f <- 1-df$aid_cumulative_bin} else if (fuzzy=="random") f <- rbinom(nrow(df), 1, mean(df$aid_cumulative_bin)) else {f <- unlist(df[, fuzzy])}
  if (!is.null(h) & is.null(b0)) b0 <- 2*h
  X <- unlist(df[, b])
  Y <- unlist(df[, v])
  if (ihs==TRUE) Y = log(Y+1)
  if (vars=="opt") {
    vopt <- pdlvarselect(v)
    covs <- covmatmaker(b, df, vopt)
  } else {covs <- covmatmaker(b, df, vars)}
  if (weights==TRUE) w <- df$wt_hh else w <- NULL
  out <- rdrobust(y = Y, x = X, c = 0, fuzzy = f, h = h, b = b0, weights = w, cluster = c, kernel = k,
                  p = poly, vce = vce, covs = covs)
  return(out)
}

plotvar <- function(v, b, df, h=20, ihs=FALSE, span = 1, k = "triangular", weights = TRUE, donut = 0, dist.exclude=NULL, vars = NULL) {
  if (!is.null(dist.exclude)) df <- filter(df, designation!=dist.exclude)
  if (!is.null(!is.null(vars) | vars=="none")) {df <- df[, c(v, b, "wt_hh", "border14_segment", "border_segment13", "elevation", "shake_pga", "wave", 
                         "quake_losses", "high_caste", "class5", "age_hh", "gorkha_hh", "slope", 
                         "time_to_market", "time_to_health")]} 
  else {df <- df[, c(v, b, "wt_hh", "border14_segment", "border_segment13", "elevation", "slope", "shake_pga", "wave", "quake_losses")]}
  df <- filter(df, abs(get(b))<h & complete.cases(df) & abs(get(b))>donut)
  Y <- unlist(df[, v])
  X <- unlist(df[, b])
  if (ihs==TRUE) Y <- ihs(Y)
  segcovs <- covmatmaker(b, df, vars)
  resreg <- lm(Y ~ segcovs)
  resregres <- resreg$residuals
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


qplot <- function(qvar, df, plot = TRUE, grid = NULL, vars = NULL, k = "triangular", h = h0, b = b0, donut = 0, poly = 1){
  if (donut!=0) df <- filter(df, abs(dist_2_seg13)>donut)
  if (!is.null(vars) | vars=="none") df <- filter(df, !is.na(high_caste), !is.na(class5), !is.na(age_hh), !is.na(gorkha_hh), !is.na(time_to_market), !is.na(time_to_health))
  qcovs <- covmatmaker("dist_2_seg13", df, vars)
  y <- unlist(df[,qvar])
  if (is.null(grid)) grid = quantile(y, seq(.1,.9,.1), na.rm = TRUE)
  qtab <- rdquant(Y = y, x = df$dist_2_seg13, c = 0, fuzzy = df$aid_cumulative_bin, grid, weights = df$wt_hh, cluster = df$strata, 
                  vce = "hc1", covs = qcovs, kernel = k, h = h, b = b, p = poly)
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
