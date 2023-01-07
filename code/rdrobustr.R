library(MASS, exclude = "select")

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

rdplot = function(y, x, c=0, p=4, nbins = NULL, binselect = "esmv", scale = NULL, 
                  kernel = "uni", weights = NULL, h = NULL, 
                  covs = NULL,  covs_eval = "mean", covs_drop = TRUE, ginv.tol = 1e-20,
                  support = NULL, subset = NULL, masspoints = "adjust",
                  hide = FALSE, ci = NULL, shade = FALSE, 
                  title = NULL, x.label = NULL, y.label = NULL, x.lim = NULL, y.lim = NULL, 
                  col.dots = NULL, col.lines = NULL) {
  
  ############################################################################################
  #start_time <- Sys.time()
  ############################################################################################
  if (!is.null(subset)) {
    x <- x[subset]
    y <- y[subset]
  }
  na.ok <- complete.cases(x) & complete.cases(y)
  
  if (!is.null(covs)){
    if (!is.null(subset))  covs <- subset(covs,subset)
    na.ok <- na.ok & complete.cases(covs)
  } 
  
  if (!is.null(weights)){
    if (!is.null(subset)) weights <- weights[subset]
    na.ok <- na.ok & complete.cases(weights) & weights>=0
  } 
  
  if (p == "loess") {p = 0; line = "loess"} else {p=p; line = FALSE}
  
  x <- x[na.ok]
  y <- y[na.ok]
  
  if (!is.null(covs))    covs    = as.matrix(covs)[na.ok, , drop = FALSE]
  if (!is.null(weights)) weights = as.matrix(weights[na.ok])  
  
  x_min = min(x);	x_max = max(x)
  ind_l = x<c;  ind_r = x>=c
  x_l = x[ind_l]; x_r = x[ind_r]	
  y_l = y[ind_l];	y_r = y[ind_r]
  
  if (!is.null(support)) {
    support_l = support[1]
    support_r = support[2]
    if (support_l<x_min) x_min = support_l
    if (support_r>x_max) x_max = support_r
  }
  
  range_l = c - x_min
  range_r = x_max - c
  
  n_l = length(x_l)
  n_r = length(x_r)
  n = n_l + n_r
  meth="es"
  
  if (is.null(scale)) {
    scale = scale_l = scale_r = 1  
  } else{
    if (length(scale)==1) scale_l = scale_r = scale
    if (length(scale)==2) {
      scale_l = scale[1]
      scale_r = scale[2]
    }
  }
  
  if (!is.null(nbins)) {
    if (length(nbins)==1) nbins_l = nbins_r = nbins
    if (length(nbins)==2) {
      nbins_l = nbins[1]
      nbins_r = nbins[2]
    }
  }
  
  if (is.null(h)) {
    h_l = range_l
    h_r = range_r
  } else{
    if (length(h)==1) h_l = h_r = h
    if (length(h)==2) {
      h_l = h[1]
      h_r = h[2]
    }
  }
  
  flag_no_ci <- FALSE
  if (is.null(ci)) {
    ci<- 95
    flag_no_ci <- TRUE
  }
  
  kernel_type = "Uniform"
  if (kernel=="epanechnikov" | kernel=="epa") kernel_type = "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel_type = "Triangular"
  
  ### Mass Points
  if (is.null(masspoints)) masspoints=FALSE
  mN = n;  M_l = n_l;  M_r = n_r
  if (masspoints=="check" | masspoints=="adjust") {
    X_uniq_l = sort(unique(x_l), decreasing=TRUE)
    X_uniq_r = unique(x_r)
    M_l = length(X_uniq_l)
    M_r = length(X_uniq_r)
    M = M_l + M_r
    mass_l = 1-M_l/n_l
    mass_r = 1-M_r/n_r				
    if (mass_l>=0.2 | mass_r>=0.2){
      print("Mass points detected in the running variable.")
      if (masspoints=="check") print("Try using option masspoints=adjust.")
      if (masspoints=="adjust") {
        if (binselect=="es")    binselect="espr"
        if (binselect=="esmv")  binselect="esmvpr"
        if (binselect=="qs")    binselect="qspr"
        if (binselect=="qsmv")  binselect="qsmvpr"
      }
    }
  }				
  
  ############## COLLINEARITY
  covs_drop_coll=dZ=0
  if (!is.null(covs)) dZ = ncol(covs)
  if (covs_drop == TRUE) covs_drop_coll = 1 
  
  if (!is.null(covs) & isTRUE(covs_drop)) {
    covs.names = colnames(covs)
    if (is.null(covs.names)) {
      covs.names = paste("z",1:ncol(covs),sep="")
      colnames(covs) = covs.names
    }
    covs = covs[,order(nchar(covs.names))]
    covs = as.matrix(covs)
    dZ = length(covs.names)
    covs.check = covs_drop_fun(covs)
    if (covs.check$ncovs < dZ) {
      covs  <- as.matrix(covs.check$covs)
      dZ    <- covs.check$ncovs
      warning("Multicollinearity issue detected in covs. Redundant covariates dropped.")  
    }
  }
  
  #####  ERRORS
  exit=0
  if (c<=x_min | c>=x_max){
    print("c should be set within the range of x")
    exit = 1
  }
  
  if (kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
    print("kernel incorrectly specified")
    exit = 1
  }
  
  if (p<0 ){
    print("p should be a positive number")
    exit = 1
  }
  
  if (scale<=0 |scale_l<=0 |scale_r<=0){
    print("scale should be a positive number")
    exit = 1
  }
  
  p_ceiling = ceiling(p)/p
  
  if (p_ceiling!=1 & p>0) {
    print("p should be an integer number")
    exit = 1
  }
  
  if (n<20){
    print("Not enough observations to perform bin calculations")
    exit = 1
  }
  
  if (exit>0) stop()
  
  ############################################################################################
  #cat(paste("Stop 1: Preps     -> ",  Sys.time()-start_time,"\n", sep=""))
  #start_time <- Sys.time()
  ############################################################################################
  
  
  ###################################################################
  ### Polynomial curve (order = p) ##################################
  ###################################################################
  W_h_l = rdrobust_kweight(x_l, c, h_l, kernel)
  W_h_r = rdrobust_kweight(x_r, c, h_r, kernel)
  
  n_h_l = sum(W_h_l>0)
  n_h_r = sum(W_h_r>0)
  
  if (!is.null(weights)) {
    fw_l = weights[ind_l];  W_h_l = fw_l*W_h_l
    fw_r = weights[ind_r];	W_h_r = fw_r*W_h_r
  }
  
  R_p_l = outer(x_l-c, Y=0:p, FUN = "^")
  R_p_r = outer(x_r-c, Y=0:p, FUN = "^")	
  invG_p_l  = qrXXinv((sqrt(W_h_l)*R_p_l))	
  invG_p_r  = qrXXinv((sqrt(W_h_r)*R_p_r))
  
  gamma_p = NULL
  if (is.null(covs)) {
    gamma_p1_l = invG_p_l%*%crossprod(R_p_l*W_h_l, y_l)	
    gamma_p1_r = invG_p_r%*%crossprod(R_p_r*W_h_r, y_r)
  } else {
    z_l  = covs[ind_l,]; z_r  = covs[ind_r,]	
    D_l  = cbind(y_l,z_l); D_r = cbind(y_r,z_r)
    U_p_l = crossprod(R_p_l*W_h_l,D_l); U_p_r = crossprod(R_p_r*W_h_r,D_r)
    beta_p_l = invG_p_l%*%crossprod(R_p_l*W_h_l,D_l); beta_p_r = invG_p_r%*%crossprod(R_p_r*W_h_r,D_r); 
    
    ZWD_p_l  = crossprod(z_l*W_h_l,D_l)
    ZWD_p_r  = crossprod(z_r*W_h_r,D_r)
    colsZ = 2:max(c(2+dZ-1,2))
    UiGU_p_l =  crossprod(U_p_l[,colsZ],invG_p_l%*%U_p_l) 
    UiGU_p_r =  crossprod(U_p_r[,colsZ],invG_p_r%*%U_p_r) 
    ZWZ_p_l = ZWD_p_l[,colsZ] - UiGU_p_l[,colsZ] 
    ZWZ_p_r = ZWD_p_r[,colsZ] - UiGU_p_r[,colsZ]     
    ZWY_p_l = ZWD_p_l[,1] - UiGU_p_l[,1] 
    ZWY_p_r = ZWD_p_r[,1] - UiGU_p_r[,1]     
    ZWZ_p = ZWZ_p_r + ZWZ_p_l
    ZWY_p = ZWY_p_r + ZWY_p_l
    if (covs_drop_coll == 0) gamma_p = chol2inv(chol(ZWZ_p))%*%ZWY_p
    if (covs_drop_coll == 1) gamma_p = ginv(ZWZ_p, tol = ginv.tol)%*%ZWY_p
    s_Y = c(1 ,  -gamma_p[,1])
    
    gamma_p1_l = t(s_Y%*%t(beta_p_l))
    gamma_p1_r = t(s_Y%*%t(beta_p_r))
  }
  
  
  ###############################################
  ### Preparte data for polynomial curve plot ###
  ###############################################
  
  nplot = 500
  x_plot_l = seq(c-h_l, c, length.out = nplot)
  x_plot_r = seq(c, c+h_r, length.out = nplot)
  
  rplot_l = outer(x_plot_l-c, Y = 0:p, FUN = "^")
  rplot_r = outer(x_plot_r-c, Y = 0:p, FUN = "^")
  
  y_hat_l = rplot_l%*%gamma_p1_l
  y_hat_r = rplot_r%*%gamma_p1_r
  
  if (!is.null(covs) & covs_eval=="mean" ) {
    gammaZ = colMeans(covs)%*%gamma_p
    y_hat_l = rplot_l%*%gamma_p1_l + c(gammaZ)
    y_hat_r = rplot_r%*%gamma_p1_r + c(gammaZ)
  }
  
  
  ############################################################################################
  #cat(paste("Stop 2: Polynomial fit -> ",  Sys.time()-start_time,"\n", sep=""))
  #start_time <- Sys.time()
  ############################################################################################
  ###############################################
  ### Optimal Bins (using polynomial order k) ###
  ###############################################
  k=4
  
  rk_l = outer(x_l, Y = 0:k, FUN = "^")
  rk_r = outer(x_r, Y = 0:k, FUN = "^")
  
  invG_k_l = try(qrXXinv(rk_l),silent=TRUE)
  invG_k_r = try(qrXXinv(rk_r),silent=TRUE)
  
  if (class(invG_k_l)[1] == "try-error" | class(invG_k_r)[1] == "try-error") {
    k = 3
    rk_l = outer(x_l, Y = 0:k, FUN = "^")
    rk_r = outer(x_r, Y = 0:k, FUN = "^")
    invG_k_l = try(qrXXinv(rk_l),silent=TRUE)
    invG_k_r = try(qrXXinv(rk_r),silent=TRUE)		
  }
  
  if (class(invG_k_l)[1] == "try-error" | class(invG_k_r)[1] == "try-error") {
    k = 2
    rk_l = outer(x_l, Y = 0:k, FUN = "^")
    rk_r = outer(x_r, Y = 0:k, FUN = "^")
    invG_k_l = qrXXinv(rk_l)
    invG_k_r = qrXXinv(rk_r)		
  }
  
  gamma_k1_l = invG_k_l%*%crossprod(rk_l, y_l)  
  gamma_k2_l = invG_k_l%*%crossprod(rk_l, y_l^2)
  gamma_k1_r = invG_k_r%*%crossprod(rk_r, y_r)  
  gamma_k2_r = invG_k_r%*%crossprod(rk_r, y_r^2)
  
  ### Bias w/sample
  mu0_k1_l = rk_l%*%gamma_k1_l
  mu0_k1_r = rk_r%*%gamma_k1_r
  mu0_k2_l = rk_l%*%gamma_k2_l
  mu0_k2_r = rk_r%*%gamma_k2_r
  drk_l = matrix(NA,n_l,k)
  drk_r = matrix(NA,n_r,k)
  for (j in 1:k) {
    drk_l[,j] = j*x_l^(j-1)
    drk_r[,j] = j*x_r^(j-1)
  }
  
  ind_l = order(x_l); ind_r = order(x_r)
  x_i_l = x_l[ind_l]; y_i_l = y_l[ind_l]
  x_i_r = x_r[ind_r]; y_i_r = y_r[ind_r]
  
  dxi_l=(x_i_l[2:n_l]-x_i_l[1:(n_l-1)]); dyi_l=(y_i_l[2:n_l]-y_i_l[1:(n_l-1)])
  dxi_r=(x_i_r[2:n_r]-x_i_r[1:(n_r-1)]); dyi_r=(y_i_r[2:n_r]-y_i_r[1:(n_r-1)])
  
  x_bar_i_l = (x_i_l[2:n_l] + x_i_l[1:(n_l-1)]) /2
  x_bar_i_r = (x_i_r[2:n_r] + x_i_r[1:(n_r-1)]) /2
  
  drk_i_l = matrix(NA,n_l-1,k);	rk_i_l  = matrix(NA,n_l-1,(k+1))
  drk_i_r = matrix(NA,n_r-1,k);	rk_i_r  = matrix(NA,n_r-1,(k+1))
  
  for (j in 1:(k+1)) {
    rk_i_l[,j] = x_bar_i_l^(j-1)    
    rk_i_r[,j] = x_bar_i_r^(j-1)
  }
  
  for (j in 1:k) {
    drk_i_l[,j] = j*x_bar_i_l^(j-1)
    drk_i_r[,j] = j*x_bar_i_r^(j-1)
  }
  
  mu1_i_hat_l = drk_i_l%*%(gamma_k1_l[2:(k+1)])
  mu1_i_hat_r = drk_i_r%*%(gamma_k1_r[2:(k+1)])
  
  mu0_i_hat_l = rk_i_l%*%gamma_k1_l
  mu0_i_hat_r = rk_i_r%*%gamma_k1_r
  mu2_i_hat_l = rk_i_l%*%gamma_k2_l
  mu2_i_hat_r = rk_i_r%*%gamma_k2_r
  
  mu0_hat_l = rk_l%*%gamma_k1_l
  mu0_hat_r = rk_r%*%gamma_k1_r
  mu2_hat_l = rk_l%*%gamma_k2_l
  mu2_hat_r = rk_r%*%gamma_k2_r
  
  mu1_hat_l = drk_l%*%(gamma_k1_l[2:(k+1)])
  mu1_hat_r = drk_r%*%(gamma_k1_r[2:(k+1)])
  
  mu1_i_hat_l = drk_i_l%*%(gamma_k1_l[2:(k+1)])
  mu1_i_hat_r = drk_i_r%*%(gamma_k1_r[2:(k+1)])
  
  var_y_l = var(y_l)
  var_y_r = var(y_r)
  
  sigma2_hat_l_bar = mu2_i_hat_l - mu0_i_hat_l^2
  sigma2_hat_r_bar = mu2_i_hat_r - mu0_i_hat_r^2
  ind_s2_l = sigma2_hat_l_bar<0
  ind_s2_r = sigma2_hat_r_bar<0
  sigma2_hat_l_bar[ind_s2_l] = var_y_l 
  sigma2_hat_r_bar[ind_s2_r] = var_y_r  
  
  sigma2_hat_l = mu2_hat_l - mu0_hat_l^2
  sigma2_hat_r = mu2_hat_r - mu0_hat_r^2
  ind_s2_l = sigma2_hat_l<0
  ind_s2_r = sigma2_hat_r<0
  sigma2_hat_l[ind_s2_l] = var_y_l 
  sigma2_hat_r[ind_s2_r] = var_y_r  
  
  B_es_hat_dw = c( ((c-x_min)^2/(12*n))*sum(mu1_hat_l^2),((x_max-c)^2/(12*n))*sum(mu1_hat_r^2))
  V_es_hat_dw = c((0.5/(c-x_min))*sum(dxi_l*dyi_l^2),(0.5/(x_max-c))*sum(dxi_r*dyi_r^2))
  V_es_chk_dw = c((1/(c-x_min))*sum(dxi_l*sigma2_hat_l_bar),(1/(x_max-c))*sum(dxi_r*sigma2_hat_r_bar))
  J_es_hat_dw = J.fun(B_es_hat_dw, V_es_hat_dw, n)
  J_es_chk_dw = J.fun(B_es_hat_dw, V_es_chk_dw, n)
  
  B_qs_hat_dw = c((n_l^2/(24*n))*sum(dxi_l^2*mu1_i_hat_l^2), (n_r^2/(24*n))*sum(dxi_r^2*mu1_i_hat_r^2))
  V_qs_hat_dw = c((1/(2*n_l))*sum(dyi_l^2),(1/(2*n_r))*sum(dyi_r^2))
  V_qs_chk_dw = c((1/n_l)*sum(sigma2_hat_l), (1/n_r)*sum(sigma2_hat_r))
  J_qs_hat_dw = J.fun(B_qs_hat_dw, V_qs_hat_dw, n)
  J_qs_chk_dw = J.fun(B_qs_hat_dw, V_qs_chk_dw, n)
  
  J_es_hat_mv  = c(ceiling((var_y_l/V_es_hat_dw[1])*(n/log(n)^2)), ceiling((var_y_r/V_es_hat_dw[2])*(n/log(n)^2)))
  J_es_chk_mv  = c(ceiling((var_y_l/V_es_chk_dw[1])*(n/log(n)^2)), ceiling((var_y_r/V_es_chk_dw[2])*(n/log(n)^2)))
  J_qs_hat_mv  = c(ceiling((var_y_l/V_qs_hat_dw[1])*(n/log(n)^2)), ceiling((var_y_r/V_qs_hat_dw[2])*(n/log(n)^2)))
  J_qs_chk_mv  = c(ceiling((var_y_l/V_qs_chk_dw[1])*(n/log(n)^2)), ceiling((var_y_r/V_qs_chk_dw[2])*(n/log(n)^2)))
  
  
  #########################################################
  if (binselect=="es") {
    J_star_orig = J_es_hat_dw
    meth="es"
    binselect_type="IMSE-optimal evenly-spaced method using spacings estimators"
    J_IMSE = J_es_hat_dw
    J_MV   = J_es_hat_mv
  }
  if (binselect=="espr") {
    J_star_orig = J_es_chk_dw
    meth="es"
    binselect_type="IMSE-optimal evenly-spaced method using polynomial regression"
    J_IMSE = J_es_chk_dw
    J_MV   = J_es_chk_mv
  }
  if (binselect=="esmv" ) {
    J_star_orig = J_es_hat_mv
    meth="es"
    binselect_type="mimicking variance evenly-spaced method using spacings estimators"
    J_IMSE = J_es_hat_dw
    J_MV   = J_es_hat_mv
  }
  if (binselect=="esmvpr" ) {
    J_star_orig = J_es_chk_mv
    meth="es"
    binselect_type="mimicking variance evenly-spaced method using polynomial regression"
    J_IMSE = J_es_chk_dw
    J_MV   = J_es_chk_mv
  }
  if (binselect=="qs" ) {
    J_star_orig = J_qs_hat_dw
    meth="qs"
    binselect_type="IMSE-optimal quantile-spaced method using spacings estimators"
    J_IMSE = J_qs_hat_dw
    J_MV   = J_qs_hat_mv
  }
  if (binselect=="qspr" ) {
    J_star_orig = J_qs_chk_dw
    meth="qs"
    binselect_type="IMSE-optimal quantile-spaced method using polynomial regression"
    J_IMSE = J_qs_chk_dw
    J_MV   = J_qs_chk_mv
  }
  if (binselect=="qsmv" ) {
    J_star_orig = J_qs_hat_mv
    meth="qs"
    binselect_type="mimicking variance quantile-spaced method using spacings estimators"
    J_IMSE = J_qs_hat_dw
    J_MV   = J_qs_hat_mv
  }
  if (binselect=="qsmvpr" ) {
    J_star_orig = J_qs_chk_mv
    meth="qs"
    binselect_type="mimicking variance quantile-spaced method using polynomial regression"
    J_IMSE = J_qs_chk_dw
    J_MV   = J_qs_chk_mv
  }
  
  J_star_l = scale_l*J_star_orig[1]
  J_star_r = scale_r*J_star_orig[2]
  
  if (!is.null(nbins)) {
    J_star_l = nbins_l
    J_star_r = nbins_r
    binselect_type="manually evenly spaced"
  }
  
  if (var_y_l==0) {
    J_star_l = J_star_l_orig = 1
    print("Warning: not enough variability in the outcome variable below the threshold")
  }
  
  if (var_y_r==0) {
    J_star_r = J_star_r_orig = 1
    print("Warning: not enough variability in the outcome variable above the threshold")
  }
  
  rscale_l = J_star_l / J_IMSE[1]
  rscale_r = J_star_r / J_IMSE[2]
  
  ############################################################################################
  #cat(paste("Stop 3: Optimal Bins   -> ",  Sys.time()-start_time,"\n", sep=""))
  #start_time <- Sys.time()
  ############################################################################################
  
  jump_l  = range_l/J_star_l;jump_r = range_r/J_star_r;
  
  if (meth=="es") {
    jumps_l=seq(x_min,c,jump_l)
    jumps_r=seq(c,x_max,jump_r)
  }   else if (meth=="qs") {
    jumps_l=quantile(x_l,probs=seq(0,1,1/J_star_l))
    jumps_r=quantile(x_r,probs=seq(0,1,1/J_star_r))
  }
  
  bin_x_l = findInterval(x_l, jumps_l,  rightmost.closed=TRUE) - J_star_l-1
  bin_x_r = findInterval(x_r, jumps_r,  rightmost.closed=TRUE) 
  
  ############################################################################################
  #cat(paste("Stop 4a: Loop bins     -> ",  Sys.time()-start_time,"\n", sep=""))
  #start_time <- Sys.time()
  ############################################################################################
  
  rdplot_l = aggregate(cbind(y_l,x_l), list(bin_x_l), FUN=mean)
  rdplot_r = aggregate(cbind(y_r,x_r), list(bin_x_r), FUN=mean)
  
  rdplot_bin_l       =  rdplot_l[,1]
  rdplot_mean_y_l    =  rdplot_l[,2]
  rdplot_mean_x_l    =  rdplot_l[,3]
  
  rdplot_bin_r       =  rdplot_r[,1]
  rdplot_mean_y_r    =  rdplot_r[,2]
  rdplot_mean_x_r    =  rdplot_r[,3]
  
  if (!is.null(covs) & covs_eval=="mean") {
    covs_model_l = lm(y_l~ z_l + factor(bin_x_l)) 
    covs_model_r = lm(y_r~ z_r + factor(bin_x_r)) 
    yhatZ_l = predict(covs_model_l)
    yhatZ_r = predict(covs_model_r)
    rdplot_mean_y_l = aggregate(yhatZ_l, list(bin_x_l), FUN=mean)[,2]
    rdplot_mean_y_r = aggregate(yhatZ_r, list(bin_x_r), FUN=mean)[,2]
  }
  
  t_ind_l = 1:J_star_l
  t_ind_r = 1:J_star_r
  rdplot_mean_bin_l = rowMeans(cbind(matrix(jumps_l[t_ind_l]),matrix(jumps_l[(t_ind_l+1)])))
  rdplot_mean_bin_r = rowMeans(cbind(matrix(jumps_r[t_ind_r]),matrix(jumps_r[(t_ind_r+1)])))
  
  #rdplot_mean_bin_l = rdplot_mean_bin_l[rev(-rdplot_bin_l)]
  rdplot_mean_bin_l = rdplot_mean_bin_l[(rdplot_bin_l+J_star_l+1)]
  rdplot_mean_bin_r = rdplot_mean_bin_r[rdplot_bin_r]
  
  bin_x = c(bin_x_l,bin_x_r)
  rdplot_mean_bin = c(rdplot_mean_bin_l, rdplot_mean_bin_r)
  rdplot_mean_x   = c(rdplot_mean_x_l,   rdplot_mean_x_r)
  rdplot_mean_y   = c(rdplot_mean_y_l,   rdplot_mean_y_r)
  
  
  ############################################################################################
  #cat(paste("Stop 4b: Loop bins     -> ",  Sys.time()-start_time,"\n", sep=""))
  #start_time <- Sys.time()
  ############################################################################################
  rdplot_N_l    = aggregate(y_l, list(-bin_x_l), FUN=length)[,2]
  rdplot_N_r    = aggregate(y_r, list(bin_x_r),  FUN=length)[,2]
  rdplot_sd_y_l = aggregate(y_l, list(-bin_x_l), FUN=sd)[,2]
  rdplot_sd_y_r = aggregate(y_r, list(bin_x_r),  FUN=sd)[,2]
  
  rdplot_N = c(rev(rdplot_N_l),rdplot_N_r)
  
  quant = -qt((1-(ci/100))/2,pmax(rdplot_N-1,1))
  
  rdplot_sd_y_l[is.na(rdplot_sd_y_l)] =0
  rdplot_sd_y_r[is.na(rdplot_sd_y_r)] =0
  
  rdplot_sd_y = c(rev(rdplot_sd_y_l),rdplot_sd_y_r)
  rdplot_se_y <- rdplot_sd_y/sqrt(rdplot_N)
  
  rdplot_cil_bin = rdplot_mean_y - quant*rdplot_se_y
  rdplot_cir_bin = rdplot_mean_y + quant*rdplot_se_y
  
  temp_plot =NULL
  
  rdplot_min_bin_l = jumps_l[1:J_star_l]
  rdplot_max_bin_l = jumps_l[2:(J_star_l+1)]
  
  rdplot_min_bin_r = jumps_r[1:J_star_r]
  rdplot_max_bin_r = jumps_r[2:(J_star_r+1)]
  
  bin_length_l = rdplot_max_bin_l - rdplot_min_bin_l
  bin_length_r = rdplot_max_bin_r - rdplot_min_bin_r
  
  bin_avg_l = mean(bin_length_l)
  bin_avg_r = mean(bin_length_r)
  
  bin_med_l = median(bin_length_l)
  bin_med_r = median(bin_length_r)
  
  rdplot_min_bin = c(rdplot_min_bin_l[rev(-rdplot_bin_l)], rdplot_min_bin_r[rdplot_bin_r])
  rdplot_max_bin = c(rdplot_max_bin_l[rev(-rdplot_bin_l)], rdplot_max_bin_r[rdplot_bin_r])
  bin_length     = c(bin_length_l, bin_length_r)
  bin_avg        = c(bin_avg_l, bin_avg_r)
  bin_med        = c(bin_med_l, bin_avg_r)
  
  vars_bins = data.frame("rdplot_mean_bin"=rdplot_mean_bin,"rdplot_mean_x"=rdplot_mean_x, "rdplot_mean_y"=rdplot_mean_y, "rdplot_min_bin"=rdplot_min_bin, "rdplot_max_bin"=rdplot_max_bin, "rdplot_se_y"=rdplot_se_y, "rdplot_N"=rdplot_N, "rdplot_ci_l"=rdplot_cil_bin, "rdplot_ci_r"=rdplot_cir_bin)
  vars_poly = data.frame("rdplot_x"= c(x_plot_l, x_plot_r), "rdplot_y"= c(y_hat_l, y_hat_r))
  
  ############################################################################################
  #cat(paste("Stop 4d: Loop bins     -> ",  Sys.time()-start_time,"\n", sep=""))
  #start_time <- Sys.time()
  ############################################################################################
  
  if (is.null(col.lines)) col.lines = "red"
  if (is.null(col.dots))  col.dots  = "darkblue"
  #if (is.null(type.dots)) type.dots = 20
  if (is.null(title)) title="RD Plot"
  if (is.null(x.label)) x.label="X axis"
  if (is.null(y.label)) y.label="Y axis"
  #if (is.null(x.lim)) x.lim=c(min(x_l),max(x_r))
  #if (is.null(y.lim)) y.lim=c(min(c(y_l,y_r)),max(c(y_l,y_r)))
  #if (is.null(y.lim)) y.lim=c(min(rdplot_mean_y),max(rdplot_mean_y))
  
  data_bins <- data.frame(rdplot_mean_bin, rdplot_mean_y, rdplot_cil_bin, rdplot_cir_bin)
  data_poly <- data.frame(x_plot_l, y_hat_l, x_plot_r, y_hat_r)
  
  if (line=="loess") {
    temp_plot <- ggplot() + theme_bw() +
      geom_point(data=data_bins, aes(x=rdplot_mean_bin, y=rdplot_mean_y), col=col.dots, na.rm=TRUE) +
      geom_smooth(aes(x=x[x>=c], y=y[x>=c], weight = weights[x>=c]), 
                  se = FALSE, color = col.lines, method = "loess") +
      geom_smooth(aes(x=x[x<c], y=y[x<c], weight = weights[x<c]), 
                  se = FALSE, color = col.lines, method = "loess")
  } else {
    temp_plot <- ggplot() + theme_bw() +
      geom_point(data=data_bins, aes(x=rdplot_mean_bin, y=rdplot_mean_y), col=col.dots, na.rm=TRUE) +
      geom_line( data=data_poly, aes(x=x_plot_l, y=y_hat_l), col=col.lines, na.rm=TRUE) +
      geom_line( data=data_poly, aes(x=x_plot_r, y=y_hat_r), col=col.lines, na.rm=TRUE) 
    
    if (flag_no_ci==FALSE)
      temp_plot <- temp_plot +
        geom_errorbar(data=data_bins, aes(x=rdplot_mean_bin, ymin=rdplot_cil_bin, ymax=rdplot_cir_bin), linetype = 1) 
    if (shade==TRUE){
      temp_plot <- temp_plot +
        geom_ribbon(data=data_bins, aes(x=rdplot_mean_bin, ymin=rdplot_cil_bin, ymax=rdplot_cir_bin))
    }
  }
  
  temp_plot <- temp_plot + labs(x = x.label, y = y.label) + ggtitle(title)+
    coord_cartesian(xlim = x.lim, ylim = y.lim) +
    theme(legend.position = "None") +
    geom_vline(xintercept = c, size = 0.5) 
  
  if (hide == FALSE) print(temp_plot)
  
  
  
  ############################################################################################
  #cat(paste("Stop 5: Plot      -> ",  Sys.time()-start_time,"\n", sep=""))
  #start_time <- Sys.time()
  ############################################################################################
  coef = cbind(gamma_p1_l,gamma_p1_r)
  colnames(coef)=c("Left","Right")
  
  out=list(coef=coef, rdplot=temp_plot, vars_bins=vars_bins, vars_poly=vars_poly,
           J=c(J_star_l,J_star_r), J_IMSE=J_IMSE, J_MV=J_MV, 
           scale=c(scale_l,scale_r), rscale=c(rscale_l,rscale_r),
           bin_avg=c(bin_avg_l,bin_avg_r), bin_med=c(bin_med_l,bin_med_r),
           p=p, c=c, h=c(h_l,h_r), N=c(n_l,n_r), N_h=c(n_h_l,n_h_r), 
           binselect=binselect_type, kernel=kernel_type, coef_covs=gamma_p)
  
  out$call <- match.call()
  class(out) <- "rdplot"
  return(invisible(out))
}

rdrobust = function(y, x, c = NULL, fuzzy = NULL, deriv = NULL,  
                    p = NULL, q = NULL, h = NULL, b = NULL, rho = NULL, 
                    covs = NULL, covs_drop = TRUE, ginv.tol = 1e-20,
                    kernel = "tri", weights = NULL, bwselect = "mserd",
                    vce = "nn", cluster = NULL, nnmatch = 3, level = 95, 
                    scalepar = 1, scaleregul = 1, sharpbw = FALSE, 
                    all = NULL, subset = NULL, masspoints = "adjust",
                    bwcheck = NULL, bwrestrict=TRUE, stdvars=FALSE) {
  
  #print("Start Code")
  #start_time <- Sys.time()
  
  if (!is.null(subset)) {
    x <- x[subset]
    y <- y[subset]
  }
  
  if (is.null(c)) c <- 0
  if (is.null(p) & !is.null(deriv)) {p = deriv+1}
  
  if (length(p) == 0) {
    flag_no_p <- TRUE
    p <- 1
  } else if ((length(p) != 1) | !(p[1]%in%0:20)) {
    stop("Polynomial order p incorrectly specified.\n")
  } else {
    flag_no_p <- FALSE
  }
  
  if (length(q) == 0) {
    flag_no_q <- TRUE
    q <- p + 1
  } else if ((length(q) > 1) | !(q[1]%in%c(0:20)) | (q[1]<p)) {
    stop("Polynomial order (for bias correction) q incorrectly specified.\n")
  } else {
    flag_no_q <- FALSE
  }
  
  if (length(deriv) == 0) {
    flag_no_deriv <- TRUE
    deriv <- 0
  } else if ((length(deriv) > 1) | !(deriv[1]%in%c(0:20)) | (deriv[1]>p)) {
    stop("Derivative order incorrectly specified.\n")
  } else {
    flag_no_deriv <- FALSE
  }
  
  na.ok <- complete.cases(x) & complete.cases(y)
  
  if (!is.null(cluster)){
    if (!is.null(subset))  cluster <- cluster[subset]
    na.ok <- na.ok & complete.cases(cluster)
  } 
  
  if (!is.null(covs)){
    if (!is.null(subset))  covs <- subset(covs,subset)
    na.ok <- na.ok & complete.cases(covs)
  } 
  
  if (!is.null(fuzzy)){
    if (!is.null(subset)) fuzzy <- fuzzy[subset]
    na.ok <- na.ok & complete.cases(fuzzy)
  } 
  
  if (!is.null(weights)){
    if (!is.null(subset)) weights <- weights[subset]
    na.ok <- na.ok & complete.cases(weights) & weights>=0
  } 
  
  x = as.matrix(x[na.ok])
  y = as.matrix(y[na.ok])
  
  if (!is.null(covs))    covs    = as.matrix(covs)[na.ok, , drop = FALSE]
  if (!is.null(fuzzy))   fuzzy   = as.matrix(  fuzzy[na.ok])
  if (!is.null(cluster)) cluster = as.matrix(cluster[na.ok])
  if (!is.null(weights)) weights = as.matrix(weights[na.ok])
  
  if (is.null(masspoints)) masspoints=FALSE
  
  if (vce=="nn" | masspoints=="check" |masspoints=="adjust") {
    order_x = order(x)
    x = x[order_x,,drop=FALSE]
    y = y[order_x,,drop=FALSE]
    if (!is.null(covs))    covs    =  as.matrix(covs)[order_x,,drop=FALSE]
    if (!is.null(fuzzy))   fuzzy   =   fuzzy[order_x,,drop=FALSE]
    if (!is.null(cluster)) cluster = cluster[order_x,,drop=FALSE]
    if (!is.null(weights)) weights = weights[order_x,,drop=FALSE]
  }
  
  ############## COLLINEARITY
  covs_drop_coll=dZ=0
  if (covs_drop == TRUE) covs_drop_coll = 1 
  if (!is.null(covs)) dZ = ncol(covs)
  
  if (!is.null(covs) & isTRUE(covs_drop)) {
    covs.names = colnames(covs)
    if (is.null(covs.names)) {
      covs.names = paste("z",1:ncol(covs),sep="")
      colnames(covs) = covs.names
    }
    covs = covs[,order(nchar(covs.names))]
    covs = as.matrix(covs)
    dZ = length(covs.names)
    covs.check = covs_drop_fun(covs)
    if (covs.check$ncovs < dZ) {
      covs  <- as.matrix(covs.check$covs)
      dZ    <- covs.check$ncovs
      warning("Multicollinearity issue detected in covs. Redundant covariates dropped.")  
    }
  }
  
  kernel   = tolower(kernel)
  bwselect = tolower(bwselect)
  vce      = tolower(vce)
  
  if (is.null(h)) {
    x_iq = quantile(x,.75,type=2) - quantile(x,.25,type=2)
    BWp = min(c(sd(x),x_iq/1.349))
    x_sd = y_sd = 1
    c_orig = c
    if (isTRUE(stdvars)) { 
      y_sd = sd(y)
      x_sd = sd(x)
      y = y/y_sd
      x = x/x_sd
      c = c/x_sd
      BWp = min(c(1,(x_iq/x_sd)/1.349))
    }
  }
  
  ind_l = x<c;  ind_r = x>=c
  X_l = x[ind_l,,drop=FALSE];  X_r = x[ind_r,,drop=FALSE]
  Y_l = y[ind_l,,drop=FALSE];  Y_r = y[ind_r,,drop=FALSE]
  x_min = min(x);  x_max = max(x)
  range_l = abs(c-x_min);  range_r = abs(c-x_max)
  N_l = length(X_l);   N_r = length(X_r)
  N = N_r + N_l
  quant = -qnorm(abs((1-(level/100))/2))
  
  dT = 0
  T_l = T_r = NULL
  perf_comp = FALSE
  if (!is.null(fuzzy)) {
    dT = 1
    T_l  = fuzzy[ind_l,,drop=FALSE];  T_r  = fuzzy[ind_r,,drop=FALSE]
    
    if (var(T_l)==0 | var(T_r)==0) perf_comp=TRUE
    
    if (perf_comp==TRUE | sharpbw==TRUE) {
      dT = 0
      T_l = T_r = NULL
    }
  }
  
  
  
  
  Z_l = Z_r = NULL
  if (!is.null(covs)) {
    Z_l  = covs[ind_l,,drop=FALSE];   Z_r  = covs[ind_r,,drop=FALSE]
  }
  
  g_l = g_r = 0       
  C_l = C_r = NULL
  if (!is.null(cluster)) {
    C_l = cluster[ind_l,,drop=FALSE]; g_l = length(unique(C_l))
    C_r = cluster[ind_r,,drop=FALSE]; g_r = length(unique(C_r))
  }
  
  fw_l = fw_r = 0 
  if (!is.null(weights)) {
    fw_l = weights[ind_l,,drop=FALSE];  fw_r = weights[ind_r,,drop=FALSE]
  }  	
  
  vce_type = "NN"
  if (vce=="hc0")         vce_type = "HC0"
  if (vce=="hc1")         vce_type = "HC1"
  if (vce=="hc2")         vce_type = "HC2"
  if (vce=="hc3")      	  vce_type = "HC3"
  if (!is.null(cluster))	vce_type = "Cluster"
  
  if (vce=="nn") {
    nn_l = rep(1,N_l)
    nn_r = rep(1,N_r)
    dups_l   = ave(nn_l, X_l, FUN = sum)
    dups_r   = ave(nn_r, X_r, FUN = sum)
    dupsid_l = ave(nn_l, X_l, FUN = cumsum)
    dupsid_r = ave(nn_r, X_r, FUN = cumsum)
  }          
  
  #####################################################   CHECK ERRORS
  exit=0
  if (kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
    warning("kernel incorrectly specified")
    exit = 1
  }
  
  if  (bwselect!="mserd" & bwselect!="msetwo" & bwselect!="msesum" & bwselect!="msecomb1" & bwselect!="msecomb2" & bwselect!="cerrd" & bwselect!="certwo" & bwselect!="cersum" & bwselect!="cercomb1" & bwselect!="cercomb2" & bwselect!=""){
    warning("bwselect incorrectly specified")  
    exit = 1
  }
  
  if (bwselect=="cct" | bwselect=="ik" | bwselect=="cv"){
    warning("bwselect options IK, CCT and CV have been depricated. Please see help for new options")  
    exit = 1
  }
  
  if (vce!="nn" & vce!="" & vce!="hc1" & vce!="hc2" & vce!="hc3" & vce!="hc0"){ 
    warning("vce incorrectly specified")
    exit = 1
  }
  
  if (c<=x_min | c>=x_max){
    warning("c should be set within the range of x")
    exit = 1
  }
  
  if (level>100 | level<=0){
    warning("level should be set between 0 and 100")
    exit = 1
  }
  
  if (!is.null(rho)){  
    if (rho<0){
      warning("rho should be greater than 0")
      exit = 1
    }
  }
  
  if (exit>0) stop()
  if (!is.null(h)) bwselect = "Manual"
  if (!is.null(h) & is.null(rho) & is.null(b)) {
    rho = 1
    b = h
  }
  if (!is.null(h) & !is.null(rho) ) b = h/rho
  
  
  if (N<20){
    warning("Not enough observations to perform bandwidth calculations. Estimates computed using entire sample")
    h = b = max(range_l,range_r)
    bwselect = "Manual"
  }
  
  if (kernel=="epanechnikov" | kernel=="epa") {
    kernel_type = "Epanechnikov"
    C_c = 2.34
  }  else if (kernel=="uniform" | kernel=="uni") {
    kernel_type = "Uniform"
    C_c = 1.843
  }   else  {
    kernel_type = "Triangular"
    C_c = 2.576
  }
  
  vce_type = "NN"
  if (vce=="hc0")     		vce_type = "HC0"
  if (vce=="hc1")      	  vce_type = "HC1"
  if (vce=="hc2")      	  vce_type = "HC2"
  if (vce=="hc3")      	  vce_type = "HC3"
  if (vce=="cluster")  	  vce_type = "Cluster"
  if (vce=="nncluster") 	vce_type = "NNcluster"
  
  
  ############################################################################################
  #cat(paste("Preparing data -> ",  Sys.time()-start_time,"\n", sep=""))
  #start_time <- Sys.time()
  ############################################################################################
  mN = N;  M_l = N_l;  M_r = N_r
  
  if (is.null(h)) {
    
    if (masspoints=="check" | masspoints=="adjust") {
      X_uniq_l = sort(unique(X_l), decreasing=TRUE)
      X_uniq_r = unique(X_r)
      M_l = length(X_uniq_l)
      M_r = length(X_uniq_r)
      M = M_l + M_r
      mass_l = 1-M_l/N_l
      mass_r = 1-M_r/N_r				
      if (mass_l>=0.2 | mass_r>=0.2){
        warning("Mass points detected in the running variable.")
        if (masspoints=="check") warning("Try using option masspoints=adjust.")
        if (is.null(bwcheck) & masspoints=="adjust") bwcheck <- 10
      }				
    }
    
    
    c_bw = C_c*BWp*N^(-1/5)
    if (masspoints=="adjust") c_bw = C_c*BWp*M^(-1/5)
    
    if (isTRUE(bwrestrict)) {
      bw_max_l = abs(c-x_min)
      bw_max_r = abs(c-x_max)
      bw_max = max(bw_max_l, bw_max_r)
      c_bw <- min(c_bw, bw_max)
    }
    
    bw.adj <- 0
    if (!is.null(bwcheck)) {
      bwcheck_l = min(bwcheck, M_l)
      bwcheck_r = min(bwcheck, M_r)
      bw_min_l = abs(X_uniq_l-c)[bwcheck_l] + 1e-8
      bw_min_r = abs(X_uniq_r-c)[bwcheck_r] + 1e-8
      c_bw = max(c_bw, bw_min_l, bw_min_r)
      bw.adj <- 1
    }
    
    ### Step 1: d_bw
    C_d_l = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=q+1, nu=q+1, o_B=q+2, h_V=c_bw, h_B=range_l, 0, vce, nnmatch, kernel, dups_l, dupsid_l, covs_drop_coll, ginv.tol)
    C_d_r = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=q+1, nu=q+1, o_B=q+2, h_V=c_bw, h_B=range_r, 0, vce, nnmatch, kernel, dups_r, dupsid_r, covs_drop_coll, ginv.tol)
    
    #### TWO
    if  (bwselect=="msetwo" |  bwselect=="certwo" | bwselect=="msecomb2" | bwselect=="cercomb2" )  {		
      d_bw_l = c((  C_d_l$V              /   C_d_l$B^2             )^C_d_l$rate)
      d_bw_r = c((  C_d_r$V              /   C_d_r$B^2             )^C_d_l$rate)
      if (isTRUE(bwrestrict)) {
        d_bw_l <- min(d_bw_l, bw_max_l)
        d_bw_r <- min(d_bw_r, bw_max_r)
      }
      
      if (!is.null(bwcheck)) {
        d_bw_l  <- max(d_bw_l, bw_min_l)
        d_bw_r  <- max(d_bw_r, bw_min_r)
      }
      C_b_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_l, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l, covs_drop_coll, ginv.tol)
      b_bw_l = c((  C_b_l$V              /   (C_b_l$B^2 + scaleregul*C_b_l$R)        )^C_b_l$rate)
      C_b_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_r, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r, covs_drop_coll, ginv.tol)
      b_bw_r = c((  C_b_r$V              /   (C_b_r$B^2 + scaleregul*C_b_r$R)        )^C_b_l$rate)
      
      if (isTRUE(bwrestrict)) {
        b_bw_l <- min(b_bw_l, bw_max_l)
        b_bw_r <- min(b_bw_r, bw_max_r)
      }
      
      C_h_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_l, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l, covs_drop_coll, ginv.tol)
      h_bw_l = c((  C_h_l$V              /   (C_h_l$B^2 + scaleregul*C_h_l$R)         )^C_h_l$rate)
      C_h_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_r, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r, covs_drop_coll, ginv.tol)
      h_bw_r = c((  C_h_r$V              /   (C_h_r$B^2 + scaleregul*C_h_r$R)         )^C_h_l$rate)
      
      if (isTRUE(bwrestrict)) {
        h_bw_l <- min(h_bw_l, bw_max_l)
        h_bw_r <- min(h_bw_r, bw_max_r) 
      }
      
    }
    
    ### SUM
    if  (bwselect=="msesum" | bwselect=="cersum" |  bwselect=="msecomb1" | bwselect=="msecomb2" |  bwselect=="cercomb1" | bwselect=="cercomb2")  {
      
      d_bw_s = c(( (C_d_l$V + C_d_r$V)  /  (C_d_r$B + C_d_l$B)^2 )^C_d_l$rate)
      
      if (isTRUE(bwrestrict)) {
        d_bw_s <- min(d_bw_s, bw_max)
      }
      
      if (!is.null(bwcheck)) d_bw_s  <-  max(d_bw_s, bw_min_l, bw_min_r)
      
      C_b_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_s, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l, covs_drop_coll, ginv.tol)
      C_b_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_s, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r, covs_drop_coll, ginv.tol)
      b_bw_s = c(( (C_b_l$V + C_b_r$V)  /  ((C_b_r$B + C_b_l$B)^2 + scaleregul*(C_b_r$R+C_b_l$R)) )^C_b_l$rate)
      
      if (isTRUE(bwrestrict)) {
        b_bw_s <- min(b_bw_s, bw_max)
      }
      
      C_h_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_s, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l, covs_drop_coll, ginv.tol)
      C_h_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_s, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r, covs_drop_coll, ginv.tol)
      h_bw_s = c(( (C_h_l$V + C_h_r$V)  /  ((C_h_r$B + C_h_l$B)^2 + scaleregul*(C_h_r$R + C_h_l$R)) )^C_h_l$rate)
      
      if (isTRUE(bwrestrict)) {
        h_bw_s <- min(h_bw_s, bw_max)
      }
    }
    
    ### RD
    if  (bwselect=="mserd" | bwselect=="cerrd" | bwselect=="msecomb1" | bwselect=="msecomb2" | bwselect=="cercomb1" | bwselect=="cercomb2" | bwselect=="" ) {
      d_bw_d = c(( (C_d_l$V + C_d_r$V)  /  (C_d_r$B - C_d_l$B)^2 )^C_d_l$rate)
      
      if (isTRUE(bwrestrict)) {
        d_bw_d <- min(d_bw_d, bw_max)
      }
      
      if (!is.null(bwcheck)) d_bw_d  <- max(d_bw_d, bw_min_l, bw_min_r)
      
      C_b_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_d, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l, covs_drop_coll, ginv.tol)
      C_b_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=q, nu=p+1, o_B=q+1, h_V=c_bw, h_B=d_bw_d, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r, covs_drop_coll, ginv.tol)
      b_bw_d = c(( (C_b_l$V + C_b_r$V)  /  ((C_b_r$B - C_b_l$B)^2 + scaleregul*(C_b_r$R + C_b_l$R)) )^C_b_l$rate)
      
      if (isTRUE(bwrestrict)) {
        b_bw_d <- min(b_bw_d, bw_max)
      }
      
      C_h_l  = rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_d, scaleregul, vce, nnmatch, kernel, dups_l, dupsid_l, covs_drop_coll, ginv.tol)
      C_h_r  = rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c=c, o=p, nu=deriv, o_B=q, h_V=c_bw, h_B=b_bw_d, scaleregul, vce, nnmatch, kernel, dups_r, dupsid_r, covs_drop_coll, ginv.tol)
      h_bw_d = c(( (C_h_l$V + C_h_r$V)  /  ((C_h_r$B - C_h_l$B)^2 + scaleregul*(C_h_r$R + C_h_l$R)) )^C_h_l$rate)
      
      if (isTRUE(bwrestrict)) {
        h_bw_d <- min(h_bw_d, bw_max)
      }
      
    }	
    
    if  (bwselect=="mserd" | bwselect=="cerrd" | bwselect=="msecomb1" | bwselect=="msecomb2" | bwselect=="cercomb1" | bwselect=="cercomb2" | bwselect=="" ) {
      h_mserd = x_sd*h_bw_d
      b_mserd = x_sd*b_bw_d
    }	
    if  (bwselect=="msesum" | bwselect=="cersum" |  bwselect=="msecomb1" | bwselect=="msecomb2" |  bwselect=="cercomb1" | bwselect=="cercomb2" )  {
      h_msesum = x_sd*h_bw_s
      b_msesum = x_sd*b_bw_s
    }
    if  (bwselect=="msetwo" |  bwselect=="certwo" | bwselect=="msecomb2" | bwselect=="cercomb2")  {		
      h_msetwo_l = x_sd*h_bw_l
      h_msetwo_r = x_sd*h_bw_r
      b_msetwo_l = x_sd*b_bw_l
      b_msetwo_r = x_sd*b_bw_r
    }
    if  (bwselect=="msecomb1" | bwselect=="cercomb1" ) {
      h_msecomb1 = min(c(h_mserd,h_msesum))
      b_msecomb1 = min(c(b_mserd,b_msesum))
    }
    if  (bwselect=="msecomb2" | bwselect=="cercomb2") {
      h_msecomb2_l = median(c(h_mserd,h_msesum,h_msetwo_l))
      h_msecomb2_r = median(c(h_mserd,h_msesum,h_msetwo_r))
      b_msecomb2_l = median(c(b_mserd,b_msesum,b_msetwo_l))
      b_msecomb2_r = median(c(b_mserd,b_msesum,b_msetwo_r))
    }
    
    
    cer_h = N^(-(p/((3+p)*(3+2*p))))
    
    if (!is.null(cluster)) {
      cer_h = (g_l+g_r)^(-(p/((3+p)*(3+2*p))))
    }
    
    cer_b = 1
    if  (bwselect=="cerrd"){
      h_cerrd = h_mserd*cer_h
      b_cerrd = b_mserd*cer_b
    }
    if  (bwselect=="cersum"){
      h_cersum = h_msesum*cer_h
      b_cersum=  b_msesum*cer_b
    }
    if  (bwselect=="certwo"){
      h_certwo_l   = h_msetwo_l*cer_h
      h_certwo_r   = h_msetwo_r*cer_h
      b_certwo_l   = b_msetwo_l*cer_b
      b_certwo_r   = b_msetwo_r*cer_b
    }
    if  (bwselect=="cercomb1"){
      h_cercomb1 = h_msecomb1*cer_h
      b_cercomb1 = b_msecomb1*cer_b
    }
    if  (bwselect=="cercomb2"){
      h_cercomb2_l = h_msecomb2_l*cer_h
      h_cercomb2_r = h_msecomb2_r*cer_h
      b_cercomb2_l = b_msecomb2_l*cer_b
      b_cercomb2_r = b_msecomb2_r*cer_b
    }
    
    bw_list = bwselect
    bws = matrix(NA,1,4)
    colnames(bws)=c("h (left)","h (right)","b (left)","b (right)")
    rownames(bws)=bwselect
    if  (bwselect=="mserd" | bwselect=="") bws[1,] = c(h_mserd,      h_mserd,      b_mserd,      b_mserd)
    if  (bwselect=="msetwo")               bws[1,] = c(h_msetwo_l,   h_msetwo_r,   b_msetwo_l,   b_msetwo_r)
    if  (bwselect=="msesum")               bws[1,] = c(h_msesum,     h_msesum,     b_msesum,     b_msesum)
    if  (bwselect=="msecomb1")             bws[1,] = c(h_msecomb1,   h_msecomb1,   b_msecomb1,   b_msecomb1)
    if  (bwselect=="msecomb2")             bws[1,] = c(h_msecomb2_l, h_msecomb2_r, b_msecomb2_l, b_msecomb2_r) 
    if  (bwselect=="cerrd")                bws[1,] = c(h_cerrd,      h_cerrd,      b_cerrd,      b_cerrd)
    if  (bwselect=="certwo")               bws[1,] = c(h_certwo_l,   h_certwo_r,   b_certwo_l,   b_certwo_r)
    if  (bwselect=="cersum")               bws[1,] = c(h_cersum,     h_cersum,     b_cersum,     b_cersum)
    if  (bwselect=="cercomb1")             bws[1,] = c(h_cercomb1,   h_cercomb1,   b_cercomb1,   b_cercomb1)
    if  (bwselect=="cercomb2")             bws[1,] = c(h_cercomb2_l, h_cercomb2_r, b_cercomb2_l, b_cercomb2_r)
    
    h_l = c(bws[1]); b_l = c(bws[3])
    h_r = c(bws[2]); b_r = c(bws[4])
    
    if (!is.null(rho)) {
      b_l = h_l/rho
      b_r = h_r/rho
    }
    
  } else{
    if (length(h)==1) h_l = h_r = h
    if (length(h)==2) {
      h_l = h[1]
      h_r = h[2]
    }
    if (is.null(b)) {
      b_l = h_l
      b_r = h_r
    } else {
      if (length(b)==1) b_l = b_r = b
      if (length(b)==2) {
        b_l = b[1]
        b_r = b[2]
      }  
    }  
  }
  
  if (isTRUE(stdvars)) { 
    c = c*x_sd
    X_l = X_l*x_sd;	X_r = X_r*x_sd
    Y_l = Y_l*y_sd;	Y_r = Y_r*y_sd
  }
  
  ### end BW selection
  
  
  ### Estimation
  w_h_l <- rdrobust_kweight(X_l,c,h_l,kernel);	w_h_r <- rdrobust_kweight(X_r,c,h_r,kernel)
  w_b_l <- rdrobust_kweight(X_l,c,b_l,kernel);	w_b_r <- rdrobust_kweight(X_r,c,b_r,kernel)
  
  if (!is.null(weights)) {
    w_h_l <- fw_l*w_h_l;	w_h_r <- fw_r*w_h_r
    w_b_l <- fw_l*w_b_l;	w_b_r <- fw_r*w_b_r			
  }
  
  ind_h_l <- w_h_l> 0;		ind_h_r <- w_h_r> 0
  ind_b_l <- w_b_l> 0;		ind_b_r <- w_b_r> 0
  N_h_l <- sum(ind_h_l); N_b_l <- sum(ind_b_l)
  N_h_r <- sum(ind_h_r); N_b_r <- sum(ind_b_r)
  
  ind_l = ind_b_l; ind_r = ind_b_r
  if (h_l>b_l) ind_l = ind_h_l   
  if (h_r>b_r) ind_r = ind_h_r   
  
  eN_l = sum(ind_l); eN_r = sum(ind_r)
  eY_l  = Y_l[ind_l,,drop=FALSE];	eY_r  = Y_r[ind_r,,drop=FALSE]
  eX_l  = X_l[ind_l,,drop=FALSE];	eX_r  = X_r[ind_r,,drop=FALSE]
  W_h_l = w_h_l[ind_l];	W_h_r = w_h_r[ind_r]
  W_b_l = w_b_l[ind_l];	W_b_r = w_b_r[ind_r]
  
  edups_l = edupsid_l = edups_r = edupsid_r = 0
  
  if (vce=="nn") {
    edups_l   = dups_l[ind_l] 
    edups_r   = dups_r[ind_r]
    edupsid_l = dupsid_l[ind_l]
    edupsid_r = dupsid_r[ind_r]
  }          
  
  u_l <- (eX_l-c)/h_l;	u_r <-(eX_r-c)/h_r
  R_q_l = matrix(NA,eN_l,(q+1)); R_q_r = matrix(NA,eN_r,(q+1))
  for (j in 1:(q+1))  {
    R_q_l[,j] = (eX_l-c)^(j-1)
    R_q_r[,j] = (eX_r-c)^(j-1)
  }
  R_p_l = R_q_l[,1:(p+1)]; R_p_r = R_q_r[,1:(p+1)]
  
  
  #print(Sys.time()-start_time)
  #print("Computing RD estimates.")
  #start_time <- Sys.time()
  
  L_l = crossprod(R_p_l*W_h_l,u_l^(p+1)); L_r = crossprod(R_p_r*W_h_r,u_r^(p+1)) 
  invG_q_l  = qrXXinv((sqrt(W_b_l)*R_q_l));	invG_q_r  = qrXXinv((sqrt(W_b_r)*R_q_r))
  invG_p_l  = qrXXinv((sqrt(W_h_l)*R_p_l));	invG_p_r  = qrXXinv((sqrt(W_h_r)*R_p_r))
  e_p1 = matrix(0,(q+1),1); e_p1[p+2]=1
  e_v  = matrix(0,(p+1),1); e_v[deriv+1]=1
  Q_q_l = t(t(R_p_l*W_h_l) - h_l^(p+1)*(L_l%*%t(e_p1))%*%t(t(invG_q_l%*%t(R_q_l))*W_b_l))
  Q_q_r = t(t(R_p_r*W_h_r) - h_r^(p+1)*(L_r%*%t(e_p1))%*%t(t(invG_q_r%*%t(R_q_r))*W_b_r))
  D_l = eY_l; D_r = eY_r
  
  eT_l = eT_r = NULL
  if (!is.null(fuzzy)) {
    
    if (perf_comp==TRUE | sharpbw==TRUE) {
      dT = 1
      T_l  = fuzzy[x<c,,drop=FALSE];  T_r  = fuzzy[x>=c,,drop=FALSE]
    }
    
    eT_l = T_l[ind_l,,drop=FALSE]; D_l  = cbind(D_l,eT_l)
    eT_r = T_r[ind_r,,drop=FALSE]; D_r  = cbind(D_r,eT_r)
  }
  
  eZ_l = eZ_r = NULL
  if (!is.null(covs)) {
    eZ_l  = Z_l[ind_l,,drop=FALSE]; D_l   = cbind(D_l,eZ_l)
    eZ_r  = Z_r[ind_r,,drop=FALSE]; D_r   = cbind(D_r,eZ_r)
    U_p_l = crossprod(R_p_l*W_h_l,D_l); U_p_r = crossprod(R_p_r*W_h_r,D_r)
  }
  
  eC_l = eC_r = NULL
  if (!is.null(cluster)) {
    eC_l  = C_l[ind_l]; eC_r  = C_r[ind_r]
  }
  
  beta_p_l  = invG_p_l%*%crossprod(R_p_l*W_h_l,D_l) 
  beta_p_r  = invG_p_r%*%crossprod(R_p_r*W_h_r,D_r) 
  beta_q_l  = invG_q_l%*%crossprod(R_q_l*W_b_l,D_l)
  beta_q_r  = invG_q_r%*%crossprod(R_q_r*W_b_r,D_r)
  beta_bc_l = invG_p_l%*%crossprod(Q_q_l,D_l) 
  beta_bc_r = invG_p_r%*%crossprod(Q_q_r,D_r)
  beta_p    = beta_p_r  - beta_p_l
  beta_q    = beta_q_r  - beta_q_l
  beta_bc   = beta_bc_r - beta_bc_l
  
  gamma_p = NULL
  if (is.null(covs)) {	
    tau_cl = tau_Y_cl = scalepar*factorial(deriv)*beta_p[(deriv+1),1]
    tau_bc = tau_Y_bc = scalepar*factorial(deriv)*beta_bc[(deriv+1),1]
    s_Y = 1
    
    tau_Y_cl_l = scalepar*factorial(deriv)*beta_p_l[(deriv+1),1]
    tau_Y_cl_r = scalepar*factorial(deriv)*beta_p_r[(deriv+1),1]
    tau_Y_bc_l = scalepar*factorial(deriv)*beta_bc_l[(deriv+1),1]
    tau_Y_bc_r = scalepar*factorial(deriv)*beta_bc_r[(deriv+1),1]
    bias_l = tau_Y_cl_l-tau_Y_bc_l
    bias_r = tau_Y_cl_r-tau_Y_bc_r 
    
    beta_Y_p_l = scalepar*factorial(deriv)*beta_p_l[,1]
    beta_Y_p_r = scalepar*factorial(deriv)*beta_p_r[,1]
    
    
    if (!is.null(fuzzy)) {
      tau_T_cl = factorial(deriv)*beta_p[(deriv+1),2]
      tau_T_bc = factorial(deriv)*beta_bc[(deriv+1),2]
      tau_cl   = tau_Y_cl/tau_T_cl
      s_Y      = c(1/tau_T_cl , -(tau_Y_cl/tau_T_cl^2))
      B_F      = c(tau_Y_cl-tau_Y_bc , tau_T_cl-tau_T_bc)
      tau_bc   = tau_cl - t(s_Y)%*%B_F
      sV_T     = c(0 , 1)
      
      tau_T_cl_l = factorial(deriv)*beta_p_l[(deriv+1),2]
      tau_T_cl_r = factorial(deriv)*beta_p_r[(deriv+1),2]
      tau_T_bc_l = factorial(deriv)*beta_bc_l[(deriv+1),2]
      tau_T_bc_r = factorial(deriv)*beta_bc_r[(deriv+1),2]
      B_F_l = c(tau_Y_cl_l-tau_Y_bc_l, tau_T_cl_l-tau_T_bc_l)
      B_F_r = c(tau_Y_cl_r-tau_Y_bc_r, tau_T_cl_r-tau_T_bc_r)
      bias_l = t(s_Y)%*%B_F_l
      bias_r = t(s_Y)%*%B_F_r
      
      beta_T_p_l = scalepar*factorial(deriv)*beta_p_l[,2]
      beta_T_p_r = scalepar*factorial(deriv)*beta_p_r[,2]
      
    }	
  } else {	
    ZWD_p_l  = crossprod(eZ_l*W_h_l,D_l)
    ZWD_p_r  = crossprod(eZ_r*W_h_r,D_r)
    colsZ = (2+dT):max(c(2+dT+dZ-1,(2+dT)))
    UiGU_p_l =  crossprod(matrix(U_p_l[,colsZ],nrow=p+1),invG_p_l%*%U_p_l) 
    UiGU_p_r =  crossprod(matrix(U_p_r[,colsZ],nrow=p+1),invG_p_r%*%U_p_r) 
    ZWZ_p_l = ZWD_p_l[,colsZ] - UiGU_p_l[,colsZ] 
    ZWZ_p_r = ZWD_p_r[,colsZ] - UiGU_p_r[,colsZ]     
    ZWY_p_l = ZWD_p_l[,1:(1+dT)] - UiGU_p_l[,1:(1+dT)] 
    ZWY_p_r = ZWD_p_r[,1:(1+dT)] - UiGU_p_r[,1:(1+dT)]     
    ZWZ_p = ZWZ_p_r + ZWZ_p_l
    ZWY_p = ZWY_p_r + ZWY_p_l
    if (covs_drop_coll == 0) gamma_p = chol2inv(chol(ZWZ_p))%*%ZWY_p
    if (covs_drop_coll == 1) gamma_p = ginv(ZWZ_p, tol = ginv.tol)%*%ZWY_p
    s_Y = c(1 ,  -gamma_p[,1])
    
    if (is.null(fuzzy)) {
      tau_cl = scalepar*t(s_Y)%*%beta_p[(deriv+1),]
      tau_bc = scalepar*t(s_Y)%*%beta_bc[(deriv+1),]
      
      tau_Y_cl_l = scalepar*t(s_Y)%*%beta_p_l[(deriv+1),]
      tau_Y_cl_r = scalepar*t(s_Y)%*%beta_p_r[(deriv+1),]
      tau_Y_bc_l = scalepar*t(s_Y)%*%beta_bc_l[(deriv+1),]
      tau_Y_bc_r = scalepar*t(s_Y)%*%beta_bc_r[(deriv+1),]
      bias_l = tau_Y_cl_l-tau_Y_bc_l
      bias_r = tau_Y_cl_r-tau_Y_bc_r 
      
      beta_Y_p_l = scalepar*tcrossprod(s_Y,beta_p_l)
      beta_Y_p_r = scalepar*tcrossprod(s_Y,beta_p_r)
      
    } else {
      s_T  = c(1,    -gamma_p[,2])
      sV_T = c(0, 1, -gamma_p[,2])
      tau_Y_cl = c(scalepar*factorial(deriv)*t(s_Y)%*%c(beta_p[(deriv+1),1], beta_p[(deriv+1),colsZ]))
      tau_Y_bc = c(scalepar*factorial(deriv)*t(s_Y)%*%c(beta_bc[(deriv+1),1],beta_bc[(deriv+1),colsZ]))
      tau_T_cl = c(         factorial(deriv)*t(s_T)%*%c(beta_p[(deriv+1),2], beta_p[(deriv+1),colsZ]))
      tau_T_bc = c(         factorial(deriv)*t(s_T)%*%c(beta_bc[(deriv+1),2],beta_bc[(deriv+1),colsZ]))
      
      tau_Y_cl_l = c(scalepar*factorial(deriv)*t(s_Y)%*%c(beta_p_l[(deriv+1),1], beta_p_l[(deriv+1),colsZ]))
      tau_Y_cl_r = c(scalepar*factorial(deriv)*t(s_Y)%*%c(beta_p_r[(deriv+1),2], beta_p_r[(deriv+1),colsZ]))
      tau_Y_bc_l = c(scalepar*factorial(deriv)*t(s_Y)%*%c(beta_bc_l[(deriv+1),1],beta_bc_l[(deriv+1),colsZ]))
      tau_Y_bc_r = c(scalepar*factorial(deriv)*t(s_Y)%*%c(beta_bc_r[(deriv+1),2],beta_bc_r[(deriv+1),colsZ]))
      
      tau_T_cl_l = c(factorial(deriv)*t(s_T)%*%c(beta_p_l[(deriv+1),1], beta_p_l[(deriv+1), colsZ]))
      tau_T_cl_r = c(factorial(deriv)*t(s_T)%*%c(beta_p_r[(deriv+1),2], beta_p_r[(deriv+1), colsZ]))
      tau_T_bc_l = c(factorial(deriv)*t(s_T)%*%c(beta_bc_l[(deriv+1),1],beta_bc_l[(deriv+1),colsZ]))
      tau_T_bc_r = c(factorial(deriv)*t(s_T)%*%c(beta_bc_r[(deriv+1),2],beta_bc_r[(deriv+1),colsZ]))
      
      beta_Y_p_l = scalepar*factorial(deriv)*t(s_Y)%*%t(cbind(beta_p_l[,1], beta_p_l[,colsZ]))
      beta_Y_p_r = scalepar*factorial(deriv)*t(s_Y)%*%t(cbind(beta_p_r[,1], beta_p_r[,colsZ]))
      beta_T_p_l =          factorial(deriv)*t(s_T)%*%t(cbind(beta_p_l[,2], beta_p_l[,colsZ]))
      beta_T_p_r =          factorial(deriv)*t(s_T)%*%t(cbind(beta_p_r[,2], beta_p_r[,colsZ]))
      
      tau_cl = tau_Y_cl/tau_T_cl
      B_F   = c(tau_Y_cl-tau_Y_bc,     tau_T_cl-tau_T_bc)
      B_F_l = c(tau_Y_cl_l-tau_Y_bc_l, tau_T_cl_l-tau_T_bc_l)
      B_F_r = c(tau_Y_cl_r-tau_Y_bc_r, tau_T_cl_r-tau_T_bc_r)
      
      s_Y = c(1/tau_T_cl , -(tau_Y_cl/tau_T_cl^2))
      tau_bc = tau_cl - t(s_Y)%*%B_F
      
      bias_l = t(s_Y)%*%B_F_l
      bias_r = t(s_Y)%*%B_F_r
      
      s_Y = c(1/tau_T_cl , -(tau_Y_cl/tau_T_cl^2) , -(1/tau_T_cl)*gamma_p[,1] + (tau_Y_cl/tau_T_cl^2)*gamma_p[,2])
    }
  }
  
  #print(Sys.time()-start_time)
  #print("Computing variance-covariance matrix.")
  #start_time <- Sys.time()
  
  hii_l = hii_r = predicts_p_l = predicts_p_r = predicts_q_l = predicts_q_r = 0
  if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
    predicts_p_l=R_p_l%*%beta_p_l
    predicts_p_r=R_p_r%*%beta_p_r
    predicts_q_l=R_q_l%*%beta_q_l
    predicts_q_r=R_q_r%*%beta_q_r
    if (vce=="hc2" | vce=="hc3") {
      hii_l = rowSums((R_p_l%*%invG_p_l)*(R_p_l*W_h_l))
      hii_r = rowSums((R_p_r%*%invG_p_r)*(R_p_r*W_h_r))
    }
  }
  
  res_h_l = rdrobust_res(eX_l, eY_l, eT_l, eZ_l, predicts_p_l, hii_l, vce, nnmatch, edups_l, edupsid_l, p+1)
  res_h_r = rdrobust_res(eX_r, eY_r, eT_r, eZ_r, predicts_p_r, hii_r, vce, nnmatch, edups_r, edupsid_r, p+1)
  
  if (vce=="nn") {
    res_b_l = res_h_l;	res_b_r = res_h_r
  } 	else {
    res_b_l = rdrobust_res(eX_l, eY_l, eT_l, eZ_l, predicts_q_l, hii_l, vce, nnmatch, edups_l, edupsid_l, q+1)
    res_b_r = rdrobust_res(eX_r, eY_r, eT_r, eZ_r, predicts_q_r, hii_r, vce, nnmatch, edups_r, edupsid_r, q+1)
  }
  
  V_Y_cl_l = invG_p_l%*%rdrobust_vce(dT+dZ, s_Y, as.matrix(R_p_l*W_h_l), res_h_l, eC_l)%*%invG_p_l
  V_Y_cl_r = invG_p_r%*%rdrobust_vce(dT+dZ, s_Y, as.matrix(R_p_r*W_h_r), res_h_r, eC_r)%*%invG_p_r
  V_Y_rb_l = invG_p_l%*%rdrobust_vce(dT+dZ, s_Y, as.matrix(Q_q_l),       res_b_l, eC_l)%*%invG_p_l
  V_Y_rb_r = invG_p_r%*%rdrobust_vce(dT+dZ, s_Y, as.matrix(Q_q_r),       res_b_r, eC_r)%*%invG_p_r
  V_tau_cl = scalepar^2*factorial(deriv)^2*(V_Y_cl_l+V_Y_cl_r)[deriv+1,deriv+1]
  V_tau_rb = scalepar^2*factorial(deriv)^2*(V_Y_rb_l+V_Y_rb_r)[deriv+1,deriv+1]
  se_tau_cl = sqrt(V_tau_cl);	se_tau_rb = sqrt(V_tau_rb)
  
  if (!is.null(fuzzy)) {
    V_T_cl_l = invG_p_l%*%rdrobust_vce(dT+dZ, sV_T, as.matrix(R_p_l*W_h_l), res_h_l, eC_l)%*%invG_p_l
    V_T_cl_r = invG_p_r%*%rdrobust_vce(dT+dZ, sV_T, as.matrix(R_p_r*W_h_r), res_h_r, eC_r)%*%invG_p_r
    V_T_rb_l = invG_p_l%*%rdrobust_vce(dT+dZ, sV_T, as.matrix(Q_q_l), res_b_l, eC_l)%*%invG_p_l
    V_T_rb_r = invG_p_r%*%rdrobust_vce(dT+dZ, sV_T, as.matrix(Q_q_r), res_b_r, eC_r)%*%invG_p_r
    V_T_cl = factorial(deriv)^2*(V_T_cl_l+V_T_cl_r)[deriv+1,deriv+1]
    V_T_rb = factorial(deriv)^2*(V_T_rb_l+V_T_rb_r)[deriv+1,deriv+1]
    se_tau_T_cl = sqrt(V_T_cl);	se_tau_T_rb = sqrt(V_T_rb)
  }
  
  #print(Sys.time()-start_time)
  
  if (is.null(fuzzy)) {
    if (is.null(covs)) {
      if      (deriv==0) rdmodel = "Sharp RD estimates using local polynomial regression." 
      else if (deriv==1) rdmodel = "Sharp Kink RD estimates using local polynomial regression."	
      else               rdmodel = "Sharp RD estimates using local polynomial regression. Derivative of order d" 
    }
    else {
      if      (deriv==0) rdmodel = "Covariate-adjusted Sharp RD estimates using local polynomial regression." 
      else if (deriv==1) rdmodel = "Covariate-adjusted Sharp Kink RD estimates using local polynomial regression."	
      else               rdmodel = paste("Covariate-adjusted Sharp RD estimates using local polynomial regression. Derivative of order ", deriv, ".")	
    }
  } else {
    if (is.null(covs)) {
      if      (deriv==0) rdmodel = "Fuzzy RD estimates using local polynomial regression." 
      else if (deriv==1) rdmodel = "Fuzzy Kink RD estimates using local polynomial regression."	
      else               rdmodel = paste("Fuzzy RD estimates using local polynomial regression. Derivative of order ", deriv, ".")	
    }
    else {
      if      (deriv==0) rdmodel = "Covariate-adjusted Fuzzy RD estimates using local polynomial regression." 
      else if (deriv==1) rdmodel = "Covariate-adjusted Fuzzy Kink RD estimates using local polynomial regression."	
      else               rdmodel = paste("Covariate-adjusted Fuzzy RD estimates using local polynomial regression. Derivative of order ", deriv, ".")			
    }
  }
  
  tau = c(tau_cl, tau_bc, tau_bc)
  se  = c(se_tau_cl,se_tau_cl,se_tau_rb)
  t   =  tau/se
  pv  = 2*pnorm(-abs(t))
  ci  = matrix(NA,nrow=3,ncol=2)
  rownames(ci) = c("Conventional","Bias-Corrected","Robust")
  colnames(ci) = c("Lower","Upper")
  ci[1,] = c(tau[1] - quant*se[1], tau[1] + quant*se[1])
  ci[2,] = c(tau[2] - quant*se[2], tau[2] + quant*se[2])
  ci[3,] = c(tau[3] - quant*se[3], tau[3] + quant*se[3])
  
  if (!is.null(fuzzy)) {  
    tau_T = c(tau_T_cl, tau_T_bc, tau_T_bc)
    se_T  = c(se_tau_T_cl, se_tau_T_cl, se_tau_T_rb)
    t_T   = tau_T/se_T
    pv_T  = 2*pnorm(-abs(t_T))
    ci_T  = matrix(NA,nrow=3,ncol=2)
    ci_T[1,] = c(tau_T[1] - quant*se_T[1], tau_T[1] + quant*se_T[1])
    ci_T[2,] = c(tau_T[2] - quant*se_T[2], tau_T[2] + quant*se_T[2])
    ci_T[3,] = c(tau_T[3] - quant*se_T[3], tau_T[3] + quant*se_T[3])
  }
  
  coef = matrix(tau,3,1)
  se   = matrix(se, 3,1)
  z    = matrix(t,  3,1)
  pv   = matrix(pv, 3,1)
  ci   = ci
  
  bws=matrix(c(h_l,b_l,h_r,b_r),2,2)
  colnames(bws)=c("left","right")
  rownames(bws)=c("h","b")
  
  rownames(coef)=rownames(se)=rownames(se)=rownames(z)=rownames(pv)=c("Conventional","Bias-Corrected","Robust")
  colnames(coef)="Coeff"
  colnames(se)="Std. Err."
  colnames(z)="z"
  colnames(pv)="P>|z|"
  colnames(bws)=c("left","right")
  rownames(ci)=c("Conventional","Bias-Corrected","Robust")
  colnames(ci)=c("CI Lower","CI Upper")
  
  Estimate=matrix(NA,1,4)
  colnames(Estimate)=c("tau.us","tau.bc","se.us","se.rb")
  Estimate[1,] <- c(tau_cl,tau_bc, se_tau_cl, se_tau_rb) 
  
  if (is.null(fuzzy)) { 
    out=list(Estimate=Estimate, bws=bws, coef=coef, se=se, z=z, pv=pv, ci=ci,
             beta_Y_p_l = beta_Y_p_l, beta_Y_p_r= beta_Y_p_r, Y = y, X = x,
             V_cl_l=V_Y_cl_l, V_cl_r=V_Y_cl_r, V_rb_l=V_Y_rb_l, V_rb_r=V_Y_rb_r,
             N=c(N_l,N_r), N_h=c(N_h_l,N_h_r), N_b=c(N_b_l,N_b_r), M=c(M_l,M_r),
             tau_cl=c(tau_Y_cl_l,tau_Y_cl_r), tau_bc=c(tau_Y_bc_l,tau_Y_bc_r),
             c=c, p=p, q=q, bias=c(bias_l,bias_r), kernel=kernel_type, all=all,
             vce=vce_type, bwselect=bwselect, level=level, masspoints=masspoints,
             rdmodel=rdmodel)
  } else {  
    out=list(Estimate=Estimate, bws=bws, coef=coef, se=se, z=z, pv=pv, ci=ci,
             beta_Y_p_l = beta_Y_p_l, beta_Y_p_r= beta_Y_p_r, Y = y, X = x,
             beta_T_p_l = beta_T_p_l, beta_T_p_r= beta_T_p_r,
             tau_T = tau_T, se_T  = se_T, t_T   = t_T, pv_T  = pv_T, ci_T  = ci_T,
             V_cl_l=V_Y_cl_l, V_cl_r=V_Y_cl_r, V_rb_l=V_Y_rb_l, V_rb_r=V_Y_rb_r,
             N=c(N_l,N_r), N_h=c(N_h_l,N_h_r), N_b=c(N_b_l,N_b_r), M=c(M_l,M_r),
             tau_cl=c(tau_Y_cl_l,tau_Y_cl_r), tau_bc=c(tau_Y_bc_l,tau_Y_bc_r),
             c=c, p=p, q=q, bias=c(bias_l,bias_r), kernel=kernel_type, all=all,
             vce=vce_type, bwselect=bwselect, level=level, masspoints=masspoints,
             rdmodel=rdmodel, beta_covs=gamma_p)
  }
  
  out$call <- match.call()
  class(out) <- "rdrobust"
  return(out)
  
  
}


qrXXinv = function(x, ...) {
  #tcrossprod(solve(qr.R(qr(x, tol = 1e-10)), tol = 1e-10))
  #tcrossprod(solve(qr.R(qr(x))))
  chol2inv(chol(crossprod(x))) 
}

qrreg = function(x,y,w,s2=0,var.comp=TRUE, ...) {
  M.X = sqrt(w)*x
  X.M.X_inv = qrXXinv(M.X) 
  X.M.Y = crossprod(M.X,sqrt(w)*y)
  beta.hat = X.M.X_inv%*%X.M.Y
  Psi.hat=Sigma.hat=0
  if (var.comp==TRUE) {
    Psi.hat = crossprod((w*s2*w)*x,x)
    Sigma.hat = crossprod(Psi.hat%*%X.M.X_inv,X.M.X_inv)
  }
  output = list(X.M.X_inv=X.M.X_inv, X.M.Y=X.M.Y, beta.hat=beta.hat, Psi.hat=Psi.hat, Sigma.hat=Sigma.hat)
  return(output)
}

rdrobust_kweight = function(X, c,  h,  kernel){
  u = (X-c)/h
  if (kernel=="epanechnikov" | kernel=="epa") {
    w = (0.75*(1-u^2)*(abs(u)<=1))/h
  }
  else if (kernel=="uniform" | kernel=="uni") {
    w = (0.5*(abs(u)<=1))/h
  }
  else {
    w = ((1-abs(u))*(abs(u)<=1))/h
  }
  return(w)	
}

rdrobust_res = function(X, y, T, Z, m, hii, vce, matches, dups, dupsid, d) {
  n = length(y)
  dT=dZ=0
  if (!is.null(T)) dT = 1
  if (!is.null(Z)) dZ = ncol(Z)
  res = matrix(NA,n,1+dT+dZ)  	
  
  if (vce=="nn") {
    for (pos in 1:n) {
      rpos = dups[pos] - dupsid[pos]
      lpos = dupsid[pos] - 1
      while (lpos+rpos < min(c(matches,n-1))) {
        if (pos-lpos-1 <= 0) rpos = rpos + dups[pos+rpos+1]
        else if (pos+rpos+1>n) lpos = lpos + dups[pos-lpos-1]
        else if ((X[pos]-X[pos-lpos-1]) > (X[pos+rpos+1]-X[pos])) rpos = rpos + dups[pos+rpos+1]
        else if ((X[pos]-X[pos-lpos-1]) < (X[pos+rpos+1]-X[pos])) lpos = lpos + dups[pos-lpos-1]
        else {
          rpos = rpos + dups[pos+rpos+1]
          lpos = lpos + dups[pos-lpos-1]
        }
      }
      ind_J = max(c(0,(pos-lpos))):min(c(n,(pos+rpos)))
      y_J   = sum(y[ind_J])-y[pos]
      Ji = length(ind_J)-1
      res[pos,1] = sqrt(Ji/(Ji+1))*(y[pos] - y_J/Ji)
      if (!is.null(T)) {
        T_J = sum(T[ind_J])-T[pos]
        res[pos,2] = sqrt(Ji/(Ji+1))*(T[pos] - T_J/Ji)
      }
      if (!is.null(Z)) {
        for (i in 1:dZ) {
          Z_J = sum(Z[ind_J,i])-Z[pos,i]
          res[pos,1+dT+i] = sqrt(Ji/(Ji+1))*(Z[pos,i] - Z_J/Ji)
        }
      }
    }		
  }
  else {
    if (vce=="hc0") w = 1
    else if (vce=="hc1") w = sqrt(n/(n-d))
    else if (vce=="hc2") w = sqrt(1/(1-hii))
    else                 w =      1/(1-hii)
    res[,1] = w*(y-m[,1])
    if (dT==1) res[,2] = w*(T-m[,2])
    if (dZ>0) {
      for (i in 1:dZ) {
        res[,1+dT+i] = w*(Z[,i]-m[,1+dT+i])
      }
    }
  }
  return(res)
}


rdrobust_bw = function(Y, X, T, Z, C, W, c, o, nu, o_B, h_V, h_B, scale, vce, nnmatch, kernel, dups, dupsid, covs_drop_coll, ginv.tol){
  dT = dZ = dC = 0
  w = rdrobust_kweight(X, c, h_V, kernel)
  dW = length(W)
  if (dW>1) {
    w = W*w
  }
  
  ind_V = w> 0; eY = Y[ind_V];eX = X[ind_V];eW = w[ind_V]
  n_V = sum(ind_V)
  D_V = eY
  R_V = matrix(NA,n_V,o+1)
  for (j in 1:(o+1)) R_V[,j] = (eX-c)^(j-1)
  invG_V = qrXXinv(R_V*sqrt(eW))
  e_v = matrix(0,(o+1),1); e_v[nu+1]=1
  s = 1
  eT=eC=eZ=NULL
  if (!is.null(T)) {
    dT = 1
    eT = T[ind_V]
    D_V = cbind(D_V,eT)
  }
  if (!is.null(Z)) {
    dZ = ncol(Z)
    eZ = Z[ind_V,,drop=FALSE]
    D_V = cbind(D_V,eZ)
    U = crossprod(R_V*eW,D_V)
    ZWD  = crossprod(eZ*eW,D_V)
    colsZ = (2+dT):max(c(2+dT+dZ-1,(2+dT)))
    UiGU =  crossprod(matrix(U[,colsZ],nrow=o+1),invG_V%*%U) 
    ZWZ = ZWD[,colsZ] - UiGU[,colsZ] 
    ZWY = ZWD[,1:(1+dT)] - UiGU[,1:(1+dT)] 
    #if (covs_drop_coll==0) gamma = chol2inv(chol(ZWZ))%*%ZWY
    if (covs_drop_coll==1) {
      gamma = ginv(ZWZ, tol=ginv.tol)%*%ZWY
    }
    else {
      gamma = chol2inv(chol(ZWZ))%*%ZWY
    }
    s = c(1 , -gamma[,1])
  }
  if (!is.null(C)) {
    dC = 1
    eC =  C[ind_V] 
  }
  beta_V = invG_V%*%crossprod(R_V*eW,D_V)	
  if (is.null(Z) & !is.null(T)) {	
    tau_Y = c(factorial(nu)*beta_V[nu+1,1])
    tau_T = c(factorial(nu)*beta_V[nu+1,2])
    s = c(1/tau_T , -(tau_Y/tau_T^2))
  }
  if (!is.null(Z) & !is.null(T)) {	
    s_T = c(1 , -gamma[,2])
    tau_Y = c(factorial(nu)*t(s)%*%  c(beta_V[nu+1,1],beta_V[nu+1,colsZ]))
    tau_T = c(factorial(nu)*t(s_T)%*%c(beta_V[nu+1,2],beta_V[nu+1,colsZ]))
    s = c(1/tau_T , -(tau_Y/tau_T^2) , -(1/tau_T)*gamma[,1] + (tau_Y/tau_T^2)*gamma[,2])
  }	
  dups_V=dupsid_V=predicts_V=0
  
  if (vce=="nn") {
    dups_V   = dups[ind_V]
    dupsid_V = dupsid[ind_V]
  }
  
  if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
    predicts_V=R_V%*%beta_V
    if (vce=="hc2" | vce=="hc3") {
      hii = rowSums((R_V%*%invG_V)*(R_V*eW))
    }
  }	
  res_V = rdrobust_res(eX, eY, eT, eZ, predicts_V, hii, vce, nnmatch, dups_V, dupsid_V, o+1)
  V_V = (invG_V%*%rdrobust_vce(dT+dZ, s, R_V*eW, res_V, eC)%*%invG_V)[nu+1,nu+1]
  v = crossprod(R_V*eW,((eX-c)/h_V)^(o+1))
  Hp = 0
  for (j in 1:(o+1)) Hp[j] = h_V^((j-1))
  BConst = (Hp*(invG_V%*%v))[nu+1]
  
  w = rdrobust_kweight(X, c, h_B, kernel)
  if (dW>1) {
    w = W*w
  }
  ind = w> 0 
  n_B = sum(ind)
  eY = Y[ind];eX = X[ind];eW = w[ind]
  D_B = eY
  R_B = matrix(NA,n_B,o_B+1)
  for (j in 1:(o_B+1)) R_B[,j] = (eX-c)^(j-1)
  invG_B = qrXXinv(R_B*sqrt(eW))
  eT=eC=eZ=NULL
  if (!is.null(T)) {
    eT = T[ind]
    D_B = cbind(D_B,eT)
  }
  if (!is.null(Z)) {
    eZ = Z[ind,,drop=FALSE]
    D_B = cbind(D_B,eZ)
  }
  if (!is.null(C)) {
    eC=C[ind]
  }	
  beta_B = invG_B%*%crossprod(R_B*eW,D_B)	
  BWreg=0
  if (scale>0) {
    e_B = matrix(0,(o_B+1),1); e_B[o+2]=1
    dups_B=dupsid_B=hii=predicts_B=0
    if (vce=="nn") {
      dups_B   = dups[ind]
      dupsid_B = dupsid[ind]
    }
    if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
      predicts_B=R_B%*%beta_B
      if (vce=="hc2" | vce=="hc3") {
        hii = rowSums((R_B%*%invG_B)*(R_B*eW))
      }
    }	
    res_B = rdrobust_res(eX, eY, eT, eZ, predicts_B, hii, vce, nnmatch, dups_B, dupsid_B,o_B+1)
    V_B = (invG_B%*%rdrobust_vce(dT+dZ, s, R_B*eW, res_B, eC)%*%invG_B)[o+2,o+2]
    BWreg = 3*BConst^2*V_B
  }
  B =  sqrt(2*(o+1-nu))*BConst%*%(t(s)%*%(beta_B[o+2,]))
  V = (2*nu+1)*h_V^(2*nu+1)*V_V
  R = scale*(2*(o+1-nu))*BWreg
  rate = 1/(2*o+3)
  output = list(V=V,B=B,R=R,rate=rate)
  return(output)
}

rdrobust_vce = function(d, s, RX, res, C) {	
  k = ncol(as.matrix(RX))
  M = matrix(0,k,k)
  n  = length(C)
  if (is.null(C)) {
    w = 1
    if (d==0){
      M  = crossprod(c(res)*RX)
    }
    else {
      for (i in 1:(1+d)) {
        SS = res[,i]*res
        for (j in 1:(1+d)) {
          M = M + crossprod(RX*(s[i]*s[j])*SS[,j],RX)
        }
      }
    }
  }
  else {	
    clusters = unique(C)
    g     = length(clusters)
    w=((n-1)/(n-k))*(g/(g-1))
    if (d==0){
      for (i in 1:g) {
        ind=C==clusters[i]
        Xi = RX[ind,,drop=FALSE]
        ri = res[ind,,drop=FALSE]
        Xr = t(crossprod(Xi,ri))
        M = M + crossprod(Xr,Xr)
      }
    }
    else {
      for (i in 1:g) {
        ind=C==clusters[i]
        Xi = RX[ind,,drop=FALSE]
        ri = res[ind,,drop=FALSE]
        MHolder = matrix(0,1+d,k)
        for (l in 1:(1+d)) {	
          MHolder[l,] = t(crossprod(Xi,s[l]*ri[,l]))
        }	
        summedvalues = t(colSums(MHolder))
        M = M + crossprod(summedvalues,summedvalues)
      }
    }
  }
  return(w*M)		
}

J.fun = function(B,V,n) {ceiling((((2*B)/V)*n)^(1/3))}








bwconst = function(p,v,kernel){
  if (kernel=="epanechnikov" | kernel=="epa" | kernel==3) {
    K.fun = function(u) {(0.75*(1-u^2)*(abs(u)<=1))}
  }
  else if (kernel=="uniform" | kernel=="uni" | kernel==2) {
    K.fun = function(u) {(0.5*(abs(u)<=1))}
  }
  else  {
    K.fun = function(u) {((1-abs(u))*(abs(u)<=1))}
  }
  p1 = p+1  
  Gamma_p = Phi_p = matrix(NA,p1,p1)
  Omega_pq = matrix(NA,p1,1)
  for (i in 1:p1) {
    Omega.fun = function(u) {K.fun(u)*(u^(p1))*(u^(i-1))}
    Omega_pq[i] = integrate(Omega.fun,lower=0,upper=1)$value
    for (j in 1:p1) {
      Gamma.fun = function(u) {K.fun(u)*(u^(i-1))*(u^(j-1))}
      Phi.fun   = function(u) {(K.fun(u)^2)*(u^(i-1))*(u^(j-1))}
      Gamma_p[i,j] = integrate(Gamma.fun,lower=0,upper=1)$value
      Phi_p[i,j] = integrate(Phi.fun,lower=0,upper=1)$value
    }
  }
  B_const = solve(Gamma_p)%*%Omega_pq
  V_const = solve(Gamma_p)%*%Phi_p%*%solve(Gamma_p)
  C1 = B_const[v+1,1]
  C2 = V_const[v+1,v+1]
  return(c(C1,C2))
}

rdvce= function(X,y,frd=NULL,p,h,matches,vce,kernel){
  m = matches+1
  n = length(X)
  p1 = p+1
  sigma = matrix(0,n,1)
  if (vce=="resid") {
    for (k in 1:n) {
      cutoff = matrix(X[k],n,1)
      cutoff1 = X[k]
      W = rdrobust_kweight(X,cutoff1,h,"kernel")
      ind=W>0
      if (sum(ind)>5) {
        w.p=W[ind]; X.p=X[ind]; y.p=y[ind]
        XX.p = matrix(c((X.p-cutoff1)^0, poly(X.p-cutoff1,degree=p,raw=T)),length(X.p),p+1)
        mu0_phat_y = qr.coef(qr(XX.p*sqrt(w.p), tol = 1e-10), sqrt(w.p)*y.p)[1]
        if (is.null(frd)) {
          sigma[k] = (y[k] - mu0_phat_y)^2
        }
        else if (!is.null(frd)) {
          z.p=frd[ind]
          out=qrreg(XX.p, z.p, w.p, var.comp=FALSE) 
          mu0_phat_z = out$beta.hat[1]
          sigma[k] = (y[k] - mu0_phat_y)*(frd[k] - mu0_phat_z)
        }
      }
    }
  }
  else  {
    #y_match_avg = z_match_avg = matrix(NA,n,1)
    for (k in 1:n) {
      diffx = abs(X - X[k])
      m.group = sort(unique(diffx))[2:m]
      ind = which(diffx %in% m.group)
      y_match_avg = z_match_avg = mean(y[ind])
      Ji = length(ind)
      if (is.null(frd)) {
        sigma[k] = (Ji/(Ji+1))*(y[k] - y_match_avg)^2
      } 
      else if (!is.null(frd)) {
        z_match_avg = mean(frd[ind])
        sigma[k] = (Ji/(Ji+1))*(y[k] - y_match_avg)*(frd[k] - z_match_avg)
      }
    }
  }
  return(sigma)
}

regconst = function(d,h){
  d2 = 2*d+1
  d1 = d+1
  mu = matrix(0,d2, 1)
  mu[1] = 1
  XX = matrix(0,d1,d1)
  for (j in 2:d2) {
    i = j-1
    if (j%%2==1) {
      mu[j] = (1/(i+1))*(h/2)^i
    }
  }
  for (j in 1:d1) {
    XX[j,] = t(mu[j:(j+d)])
  }
  invXX =solve(XX)
  return(invXX)
}

covs_drop_fun <- function(z) {
  z          <- as.matrix(z)
  ncovs      <- ncol(z)
  df         <- data.frame(z = z)
  constant   <- rep(1,nrow(df))
  tmp        <- lm(constant ~ ., data=df)
  to_keep    <- tmp$coefficients[!is.na(tmp$coefficients)]
  ncovs_keep <- length(to_keep)
  to_keep    <- names(to_keep[2:ncovs_keep])
  ncovs_keep <- ncovs_keep-1
  covs <- as.matrix(df[to_keep])
  #qr.X <- qr(x,  LAPACK = FALSE, tol = 1e-2) 
  #(rnkX <- qr.X$rank)
  #(keep <- qr.X$pivot[seq_len(rnkX)])
  #xx <- as.matrix(x[,keep])
  #output = list(xx=xx,rank_covs=rnkX)
  output = list(covs=covs, ncovs=ncovs_keep)
  return(output)
}