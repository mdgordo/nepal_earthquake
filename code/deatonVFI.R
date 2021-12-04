library(tidyverse)
library(parallel)

### parameters
gamma <- 3
beta <- .95
R <- 1.03
y <- rlnorm(1000, meanlog = 11, sd = .35)
# mean = exp(11 + .35^2/2) = 63656
# sd = sqrt((exp(.35^2) - 1)*exp(2*11 + .35^2)) = 22979
#ylist <- rep(NA, 1000)
#for (i in c(1:1000)) {
#  ylist[i] <- min(rlnorm(10000, meanlog = 11, sd = .35))
#}  
# Expectation of minimum of 100 draws is ~26150, 1000 draws ~19400, 10000 draws is 15600
#w = readRDS("w.rds")

### natural borrowing limit - not real b/c depends on # of draws from y
B <- -.2*mean(y)
cmin <- 10000

### Initial grid and guess
xgrid <- seq(B, 5000000, len = 5000)
v0 <- rep(0, length(xgrid))

### utility and production functional forms
u <- function(x){x^(1-gamma)/(1-gamma)}

### VFI functions
source(paste(getwd(), "/code/VFIfunctions.r", sep = ""))

V <- VFI(xgrid, v0)
#V <- readRDS("V.rds")

### Plots
Vp <- dfplot(V)

ggplot(filter(Vp, grid>B & grid < 10000 & iteration!=1)) +
  geom_line(aes(x = grid, y = fx, color = iteration, group = iteration))

ggplot(filter(Vp, grid>0 & iteration!=1)) +
  geom_line(aes(x = grid, y = fx, color = iteration, group = iteration))

### Get policy functions
Vfx <- approxfun(xgrid, V[,dim(V)[2]], rule = 2)

cfx <- mclapply(xgrid, policyfunc, Vfx = Vfx, mc.cores = 6)
cfx <- approxfun(xgrid, cfx, rule = 2)

ggplot() +
  geom_line(aes(x = xgrid, y = cfx(xgrid))) +
  xlim(B, 500000) + ylim(0, 100000)

### Effect of aid on consumption

a <- cfx(xgrid + 100000) - cfx(xgrid)

ggplot() +
  geom_line(aes(x = xgrid, y = a)) +
  xlim(B, 500000) + ylim(0, 100000)

### Benefits of Aid

b <- Vfx(xgrid + 100000) - Vfx(xgrid)

ggplot() +
  geom_line(aes(x = xgrid, y = b)) +
  xlim(-10000, 500000) + ylim(0, 1e-8)

### Benefits of resilience
uc <- u(cfx(xgrid + 100000)) - u(cfx(xgrid))

ggplot() +
  geom_line(aes(x = xgrid, y = b), color = "blue") +
  geom_line(aes(x = xgrid, y = uc), color = "red") +
  xlim(-10000, 500000) + ylim(0, 1e-8)
