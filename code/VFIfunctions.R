### Value function iteration - Deaton 
bellman_operator <- function(grid, w){
  Valfunc = approxfun(grid, w, rule = 2)
  optimizer <- function(x){
    if (x==B) {return(u(cmin) + beta*mean(Valfunc(B)))} else{
      objective <- function(c) {
        xtplus1 = R*(x - c) + y
        xtplus1 = if_else(xtplus1<B, B, xtplus1)
        r = u(c) + beta*mean(Valfunc(xtplus1))
        return(r)
      }
      l <- optimize(objective, interval = c(1, x - B), maximum = TRUE)
      return(l$objective)
    }
  }
  Tw = rep(NA, length(grid))
  Tw = mclapply(grid, optimizer, mc.cores = 6)
  return(unlist(Tw))
}

### Policy function Deaton
policyfunc <- function(x, Vfx){
  if (x==B) {return(cmin)} else{
    objective <- function(c) {
      xtplus1 = R*(x - c) + y
      xtplus1 = if_else(xtplus1<B, B, xtplus1)
      r = u(c) + beta*mean(Vfx(xtplus1))
      return(r)
    }
    l <- optimize(objective, interval = c(1, x - B), maximum = TRUE)
    return(l$maximum)
  }
}

### Optimization setup - Banerjee Duflo
ROI_setup <- function(Vfx, x){
  objective <- function(ck) {
    c <- ck[1]
    k <- ck[2]
    xtplus1 = R*(x - c - k) + f(k)
    r = u(c) + beta*mean(Vfx(xtplus1))
    return(r)
  }
  ### has to hold with certainty, so min of production
  constraint <- function(ck){
    c <- ck[1]
    k <- ck[2]
    return(B - R*(x - c - k) - min(f(k)))
  }
  nlprob <- OP(F_objective(objective, 2),
               F_constraint(constraint, dir = "<=", rhs = 0),
               maximum = TRUE,
               bounds = V_bound(li = c(1,2), lb = c(1,1)))
  r <- ROI_solve(nlprob, solver = "nloptr.cobyla", start = c(1.1, 1.1))
  return(objective(solution(r)))
}

### Value function iteration
bellman_operator_BD <- function(grid, w){
  Vfx = approxfun(grid, w, rule = 2)
  Tw = rep(NA, length(grid))
  Tw = mclapply(xgrid, ROI_setup, Vfx = Vfx, mc.cores = 6, mc.preschedule = FALSE)
  return(unlist(Tw))
}

####

VFI <- function(grid, vinit, tol = 1e-13, maxiter = 300, bd = FALSE){
  w = matrix(0, length(grid), 1)
  w[,1] = vinit
  d = 1
  i = 2
  while (d > tol & i < maxiter){
    w = cbind(w, rep(0, length(grid)))
    if (bd=="bd") {w[,i] = bellman_operator_BD(grid, w[,i-1])} else {
      w[,i] = bellman_operator(grid, w[,i-1])}
    d = sqrt(sum((w[,i] - w[,i-1])^2))
    i = i+1
  }
  return(w)
}

dfplot <- function(V){
  V <- as.data.frame(V)
  V$grid <- xgrid
  colnames(V) <- c(1:(dim(V)[2] - 1), "grid")
  V <- pivot_longer(V, cols = -grid, names_to = "iteration", values_to = "fx")
  V$iteration <- as.integer(V$iteration)
  return(V)
}

bufferstock <- function(hhid, Vlist){
  ### pull right parameters
  c = df$food_consumption[df$hhid==hhid]
  mu = df$mu[df$hhid==hhid]
  i = which(unlist(lapply(Vlist, function(x) x$mu))==mu)
  cfx = Vlist[[i]]$cfx
  xgrid = Vlist[[i]]$xgrid
  
  ### find root of policy function
  cfxrt = function(x){cfx(x) - c}
  if (c>cfx(max(xgrid))) {r = max(xgrid)} else {
    r = uniroot(cfxrt, interval = c(1, max(xgrid)), extendInt = "upX")$root
  }
  ### pre and post aid buffer stock - aid = 300000?
  bsx_pre = r - df$quake_aid[df$hhid==hhid]
  bsx_post = r + 50000
  bsx_actual = r
  
  ### Calculate counterfactual consumption
  cdist_pre = cfx(bsx_pre)
  cdist_post = cfx(bsx_post)
  return(c(cdist_pre, cdist_post, bsx_pre, bsx_post, bsx_actual))
}
