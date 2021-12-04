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

source(paste(getwd(), "/code/VFIfunctions.r", sep = ""))
#grid = xgrid; vinit = v0; tol = 1e-13; maxiter = 300

V <- VFI(xgrid, v0, "bd")

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

ggplot(filter(Vplot, grid!=1 & as.integer(iteration)==dim(V)[2])) +
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


