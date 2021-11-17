### RD functions

plotvar <- function(v, b, df, h=50, ihs=FALSE, span = 1, k = "epanechnikov", weights = TRUE, donut = 0) {
  df <- df[, c(v, b, "wt_hh", "border14_segment", "border_segment13", "elevation", "shake_pga", "wave")]
  df <- filter(df, abs(get(b))<h & complete.cases(df) & abs(get(b))>donut)
  Y <- unlist(df[, v])
  X <- unlist(df[, b])
  if (ihs==TRUE) Y <- ihs(Y)
  if (b=="dist_2_14") {
    segcovs <- model.matrix(~as.factor(df$border14_segment)+as.factor(df$border14_segment)*df$dist_2_14 + df$shake_pga + as.factor(df$wave))
    segcovs <- segcovs[, !(colnames(segcovs) %in% c("(Intercept)", "df$dist_2_14"))]
    resreg <- lm(Y ~ segcovs)
    resregres <- resreg$residuals
  } else if (b=="dist_2_seg13"){
    segcovs <- model.matrix(~as.factor(df$border_segment13)+as.factor(df$border_segment13)*df$dist_2_seg13 + df$shake_pga + as.factor(df$wave))
    segcovs <- segcovs[, !(colnames(segcovs) %in% c("(Intercept)", "df$dist_2_seg13"))]
    resreg <- lm(Y ~ segcovs)
    resregres <- resreg$residuals
  } else {
    segcovs <- model.matrix(df$shake_pga + as.factor(df$wave))
    segcovs <- segcovs[, !(colnames(segcovs) %in% c("(Intercept)"))]
    resreg <- lm(Y ~ segcovs)
    resregres <- resreg$residuals
  }
  if (weights==TRUE) w <- df$wt_hh else w <- NULL
  r <- rdplot(y = resregres, x = X, c = 0, weights = w, binselect = "qsmv", kernel = k, 
              poly = FALSE, title = NULL, x.label = b, y.label = v, hide = TRUE, span = span)
  return(r$rdplot)
}

histfunc <- function(b, df, h=50){
  df <- filter(df, abs(get(b))<h)
  X <- unlist(df[, b])
  W <- df$wt_hh
  freqdf <- aggregate(x = list(Freq = W), by = list(X = X), FUN = sum)
  ggplot() +
    geom_histogram(aes(x = X, weight = W), bins = 400) +
    geom_smooth(data = filter(freqdf, X >= 0), aes(x = X, y = Freq)) + 
    geom_smooth(data = filter(freqdf, X < 0), aes(x = X, y = Freq)) + 
    theme_light() + labs(y = "weighted count", x = paste("distance to ", b))
}

regout <- function(v, b, df, h = NULL, b0 = NULL, fuzzy = FALSE, ihs = FALSE, k = "epanechnikov", weights = TRUE, donut = 0, vce = "hc1"){
  if (donut!=0) df <- filter(df, abs(get(b))>donut)
  if (vce=="nn") c <- NULL else c <- df$strata
  if (fuzzy==TRUE) f <- df$aid_cumulative_bin else if (fuzzy==FALSE) f <- NULL else if (fuzzy=="inv") f <- 1-df$aid_cumulative_bin else f <- unlist(df[, fuzzy])
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

optbw <- function(v, b, df, fuzzy = FALSE, k = "epanechnikov", weights = TRUE, donut = 0, vce = "hc1"){
  if (donut!=0) df <- filter(df, abs(get(b))>donut)
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

rdgazer <- function(rdlist, dvlabs = NULL, xlines = NULL, se_r = "Conventional", type = "text", ...){
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