library(tidyverse)
library(parallel)

### some reasonable parameters
gamma <- 3
beta <- .95
alpha <- .8

u <- function(x){
  u <- x^(1-gamma)/(1-gamma)
  return(u)
}

y <- rlnorm(1000, meanlog = 11, sdlog = .35)
mu = exp(11 + .35^2/2)
## mean = 63656
## sd = sqrt((exp(.35^2) - 1)*exp(2*11 + .35^2)) = 22979

### initial grid and guess
xgrid <- seq(1, 5000000, len = 3000)
v0 <- rep(0, length(xgrid))

### bellman and value function iteration functions
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

### this might take a while
V <- VFI(xgrid, v0, y)

### final value function 
vfx <- approxfun(xgrid, V[,dim(V)[2]])

### Assume post earthquake we observe hh buffer stock = 500000, and want to find "WTP" for 100000 in aid
vfx(500000)

### now we have to calculate value functions for a grid of potential means for the income process (hold the log sd constant for now)
means <- seq(10.5, 11, by = .01)
ylist <- lapply(means, rlnorm, n = 1000, sdlog = .35)

### This will definitely take a long time and might break your computer
Vlist <- mclapply(ylist, VFI, grid = xgrid, vinit = v0, mc.cores  = detectCores() - 4)

### extract the value functions
vfxtract <- function(V){
  vfx <- approxfun(xgrid, V[,dim(V)[2]])
  return(vfx)
}

vfxlist <- lapply(Vlist, vfxtract)

### calculate vfx(600000) for each and find the one that is closest to vfx(500000) with the original income
val4aid <- lapply(vfxlist, function(f) f(600000))

vi <- which.min((unlist(val4aid) - vfx(500000))^2)
muprime <- exp(means[vi] + .35^2/2)

### calculate infinite discounted sum of income difference
(mu - muprime)/(1 - beta)


ggplot() +
  geom_line(aes(x = xgrid[50:1000], y = vfx(xgrid[50:1000])), color = "red") +
  geom_line(aes(x = xgrid[50:1000], y = vfxlist[[vi]](xgrid[50:1000])), color = "blue") +
  geom_line(aes(x = xgrid[50:1000], y = vfxlist[[vi + 2]](xgrid[50:1000])), color = "grey") +
  geom_line(aes(x = xgrid[50:1000], y = vfxlist[[vi + 1]](xgrid[50:1000])), color = "grey") +
  geom_line(aes(x = xgrid[50:1000], y = vfxlist[[vi - 1]](xgrid[50:1000])), color = "grey") +
  theme_bw() + labs(x = "Assets + Income", y = "V(x)")
  

wtp100k <- function(x){
  val4aid <- lapply(vfxlist, function(f) f(x+100000))
  vi <- which.min((unlist(val4aid) - vfx(x))^2)
  muprime <- exp(means[vi] + .35^2/2)
  return((mu - muprime)/(1 - beta))
}

xvals <- seq(10000, 300000, 10000)
yvals <- sapply(xvals, wtp100k)

ggplot() +
  geom_line(aes(x = xvals, y = yvals)) +
  theme_bw() + labs(x = "Assets + Income", y = "WTP")



##### Production

alpha <- .5
f <- 1
p <- .1
e <- rlnorm(100000, 11, .35)

yces <- function(k){
  y <- f*(alpha*e^p + (1-alpha)*k^p)^(1/p)
}

d1 <- yces(50)
d2 <- yces(100)

ggplot() + 
  geom_histogram(aes(x = d1), bins=100, fill = "red", position = "identity", alpha = .3) + 
  geom_histogram(aes(x = d2), bins=100, fill = "blue", position = "identity", alpha = .3)

