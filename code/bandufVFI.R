library(tidyverse)
library(parallel)
library(ROI)
library(ROI.plugin.alabama)
library(ROI.plugin.nloptr)
library(ROI.plugin.deoptim)

### parameters
gamma <- 3
beta <- .9
alpha <- .3
R <- 1.03

### utility and production functional forms
u <- function(x){x^(1-gamma)/(1-gamma)}
f <- function(k) {epsilon*k^alpha}

### Need to calibrate these
epsilon <- rlnorm(1000, 10, .35)
B <- -min(f(1))/(R-1)

### Initial grid and guess
xgrid <- seq(B + 3, 5000000, len = 5000)
v0 <- rep(0, length(xgrid))

### Optimization setup
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
bellman_operator <- function(grid, w){
  Vfx = approxfun(grid, w, rule = 2)
  Tw = rep(NA, length(grid))
  Tw = mclapply(xgrid, ROI_setup, Vfx = Vfx, mc.cores = 6, mc.preschedule = FALSE)
  return(unlist(Tw))
}

#grid = xgrid; vinit = v0; tol = 1e-13; maxiter = 300

VFI <- function(grid, vinit, tol = 1e-13, maxiter = 300){
  w = matrix(0, length(grid), 1)
  w[,1] = vinit
  d = 1
  i = 2
  while (d > tol & i < maxiter){
    w = cbind(w, rep(0, length(grid)))
    w[,i] = bellman_operator(grid, w[,i-1])
    d = sqrt(sum((w[,i] - w[,i-1])^2))
    i = i+1
  }
  return(w)
}

V <- VFI(xgrid, v0)

### Plot value function
dfplot <- function(V){
  V <- as.data.frame(V)
  V$grid <- xgrid
  colnames(V) <- c(1:(dim(V)[2] - 1), "grid")
  V <- pivot_longer(V, cols = -grid, names_to = "iteration", values_to = "fx")
  return(V)
}

Vplot <- dfplot(V)

ggplot(filter(Vplot, grid!=1 & iteration!=1)) +
  geom_line(aes(x = grid, y = fx, color = as.integer(iteration), group = as.integer(iteration)))

ggplot(filter(Vplot, grid!=1 & as.integer(iteration)==232)) +
  geom_line(aes(x = grid, y = fx))

### Policy Functions
vfxfin <- approxfun(xgrid, V[,dim(V)[2]], rule = 2)

optimlist <- mclapply(xgrid, function(x) constrOptim(c(.9*x, .1*x), objective, method = "Nelder-Mead", ui = matrix(c(-1, 1, 0, -1, 0, 1), nrow = 3, ncol = 2),
                                       ci = c(B - x, 0 , 0), Vfx = vfxfin, x = x, control = list(fnscale = -1, maxit = 300)), mc.cores = 6)

consumpfx <- unlist(lapply(optimlist, function(x) x$par[1]))
kapfx <- unlist(lapply(optimlist, function(x) x$par[2]))

ggplot() +
  geom_line(aes(x = xgrid, y = consumpfx))

ggplot() +
  geom_line(aes(x = xgrid, y = kapfx))


