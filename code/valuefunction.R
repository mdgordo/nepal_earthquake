library(tidyverse)
library(parallel)

##read in data, set parameters
df.hh <- read_csv("~/project/WBHRVS/full_panel.csv", guess_max = 7500)
gamma <- 3
beta <- .95

value <- function(hhid){
  ### find 3 year average of log of household income
  income = mean(df.hh$avg_income_q[df.hh$hhid==hhid], na.rm = TRUE)
  ### vector of consumption values for each year
  consumption = df.hh$consumption[df.hh$hhid==hhid]
  ### vector of amounts of earth quake aid received each year
  quake_aid = df.hh$quake_aid[df.hh$hhid==hhid]
  
  ### skip household if income or consumption is NA
  if (income <= 0 | sum(is.na(income) + is.na(consumption)>0)) {
    res <- list("x" = NA, "vx" = NA, "vaiddelta25" = NA, "vaiddelta100" = NA, "caiddelta25" = NA, "caiddelta100" = NA,
                "maxiterE" = NA, "maxiterV" = NA)
  } else {
    ### generate 200 random draws from lognormal distribution with mu and sigma based on household specific data
    y <- rlnorm(200, meanlog = income, sdlog = income/44)
    ### grid for buffer stock
    xgrid <- seq(1, max(y)*10, len = 3000)
    ### initial guesses (I found using xgrid works better and faster than a vector of zeros)
    v0 <- xgrid
    
    ### marginal utility function for CRRA
    mu <- function(x) {1/x^gamma}
    
    ### function takes grid, guess, and vector of y values and returns new guess
    euler_operator <- function(grid, w, y){
      interp <- approxfun(grid, w, rule = 2)
      Tfx <- function(x) {
        r <- max(c(beta*mean(mu(interp(x - interp(x) + y))), mu(x)))
        r <- r^(-1/gamma)
        return(r)
      }
      Tw <- sapply(grid, Tfx)
      return(Tw)
    }
    
    ### Euler function iterator
    EFI <- function(grid, vinit, y, tol = 1e-5, maxiter = 1000){
      w <- matrix(0, length(grid), 1)
      w[,1] <- vinit
      d <- 1
      i <- 2
      while (d > tol & i < maxiter){
        w <- cbind(w, rep(0, length(grid)))
        w[,i] <- euler_operator(grid, w[,i-1], y)
        d <- sqrt(sum((w[,i] - w[,i-1])^2))
        i <- i+1
      }
      return(w)
    }
    
    V <- EFI(xgrid, v0, y)
    
    ### return final function - consumption as a function of income + buffer stock
    Vfin <- V[,dim(V)[2]]
    cfx <- approxfun(xgrid, Vfin, rule = 2)
    
    ### flag to indicate if HH hit the max iterations
    if(dim(V)[2]>998) maxiterE <- 1 else maxiterE <- 0
    
    ### Use final function and actual consumption data to solve for size of buffer stock in each period
    x <- rep(0, length(consumption))
    maxflag <- rep(0, length(consumption))
    
    for (i in 1:length(consumption)) {
      c <- consumption[i]
      if(c > max(Vfin)){
        x[i] <- max(xgrid)
        maxflag[i] <- 1
      } else {
        cfxrt <- function(c){cfx(c) - consumption[i]}
        x[i] <- uniroot(cfxrt, interval = c(1, max(xgrid)), extendInt = "upX")$root
      }
    }
    
    ### Subtract any earthquake aid the household received in each period to get the value of the 'pre-aid' buffer stock
    x <- x - quake_aid
    
    ### initial guess for value function
    v0b <- rep(0, length(xgrid))
    
    ### CRRA utility
    u <- function(x){
      u <- x^(1-gamma)/(1-gamma)
      return(u)
    }
    
    ### Solve for value function
    bellman_operator <- function(grid, w, y){
      interp <- approxfun(grid,w, rule = 2)
      Tw <- rep(0, length(grid))
      for (i in 1:length(grid)) {
        xi <- grid[i]
        objective <- function(c) {
          r <- u(c) + beta*mean(interp(xi - c + y))
          return(r)
        }
        res <- optimize(objective, interval = c(1e-16, xi), maximum = TRUE)
        Tw[i] <- res$objective
      }
      return(Tw)
    }
    
    VFI <- function(grid, vinit, y, tol = 1e-13, maxiter = 300){
      w <- matrix(0, length(grid), 1)
      w[,1] <- vinit
      d <- 1
      i <- 2
      while (d > tol & i < maxiter){
        w <- cbind(w, rep(0, length(grid)))
        w[,i] <- bellman_operator(grid, w[,i-1], y)
        d <- sqrt(sum((w[,i] - w[,i-1])^2))
        i <- i+1
      }
      return(w)
    }
    
    Vb <- VFI(xgrid, v0b, y)
    
    vfx <- approxfun(xgrid, Vb[,dim(Vb)[2]], rule = 2)
    
    ### flag for max iterations
    if(dim(Vb)[2]>298) maxiterV <- 1 else maxiterV <- 0
    
    ### return x (buffer stock in each period); vx (value in utility of that buffer stock); a25 (value of 25000 rs in aid); a100 (value of 100000 rs in aid);
    ### c25 (instantaneous utility of extra consumption induced by 25000 rs in aid); c100 (same for 100000 rs aid)
    vx <- vfx(x)
    a25 = vfx(x + 25000) - vfx(x)
    a100 = vfx(x + 100000) - vfx(x)
    c25 = u(cfx(x + 25000)) - u(cfx(x))
    c100 = u(cfx(x + 100000)) - u(cfx(x))
    
    res <- list("x" = x, "vx" = vx, "vaiddelta25" = a25, "vaiddelta100" = a100, "caiddelta25" = c25, "caiddelta100" = c100, "maxflag" = maxflag,
                "maxiterE" = maxiterE, "maxiterV" = maxiterV)
  }
  return(res)
}

vallist <- mclapply(unique(df.hh$hhid), value, mc.preschedule = FALSE, mc.cores = 46)

saveRDS(vallist, "vallist.rds")


