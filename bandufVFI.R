library(tidyverse)
library(parallel)

### parameters
gamma <- 3
beta <- .9
alpha <- .3
R <- 1.03
epsilon <- rlnorm(10000, 0, 1)
B <- -.2*300000

### Initial grid and guess
xgrid <- seq(1, 5000000, len = 3000)
v0 <- rep(0, length(xgrid))

### utility and production functional forms
u <- function(x){x^(1-gamma)/(1-gamma)}
uprime <- function(x) {1/x^gamma}
f <- function(k) {epsilon*k^alpha}
fprime <- function(k) {epsilon*alpha*k^(1-alpha)}

### Value function iteration
objective <- function(ck, Vfx, x) {
  c <- ck[1]
  k <- ck[2]
  xtplus1 = R*(x - c - k) + f(k)
  r = u(c) + beta*mean(Vfx(xtplus1))
  return(r)
}

bellman_operator <- function(grid, w){
  Valfunc = approxfun(grid, w, rule = 2)
  Tw = rep(NA, length(grid))
  Tw = mclapply(grid, function(x) constrOptim(c(.9*x, .1*x), objective, method = "Nelder-Mead", ui = matrix(c(-1, 1, 0, -1, 0, 1), nrow = 3, ncol = 2),
                                          ci = c(B - x, 0 , 0), Vfx = Valfunc, x = x, control = list(fnscale = -1, maxit = 300)), mc.cores = 6)
  Tw = unlist(lapply(Tw, function(x) x$value))
  return(Tw)
}

### Other approach
optimobj <- function(xi, k, b, Valfunc){
  Valfunc = Valfunc[[1]]
  o = optimize(objective, interval = c(0, b), x = xi, k = k, Vfx = Valfunc, maximum = TRUE)
  return(o$objective)
}

bellman_operator <- function(grid, w){
  Valfunc = approxfun(grid, w, rule = 2)
  Valfunc = list(Valfunc)
  Tw = rep(NA, length(grid))
  for (i in 1:length(grid)) {
    xi = grid[i]
    kgrid = seq(0, xi, len = 3000)
    bounds = B + kgrid - xi
    Twk = mcmapply(optimobj, k = kgrid, b = bounds, x = xi, Valfunc = Valfunc, mc.cores = 6)
    Tw[i] = max(unlist(Twk))
  }
  return(Tw)
}
#####

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


