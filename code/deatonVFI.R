library(tidyverse)

gamma <- 3
beta <- .95
alpha <- .8
y <- rlnorm(1000, meanlog = 0, sd = 1)
#y <- rlnorm(1000, meanlog = 11, sd = .25)
# mean = exp(11 + .25^2/2) = 61774
# sd = sqrt((exp(.25^2) - 1)*exp(2*11 + .25^2)) = 15688
y1 <- rlnorm(1000, meanlog = 11, sd = .35)
# mean = exp(11 + .35^2/2) = 63656
# sd = sqrt((exp(.35^2) - 1)*exp(2*11 + .35^2)) = 22979
y2 <- rlnorm(1000, meanlog = 12, sd = .35)
# mean = exp(12 + .35^2/2) = 173035
# sd = sqrt((exp(.35^2) - 1)*exp(2*12 + .35^2)) = 62465
xgrid <- seq(1, 5000000, len = 3000)
v0 <- xgrid

u <- function(x) {1/x^gamma}

### y*x^alpha makes it so consumption rises w buffer stock but not proportionally
### could also do sectoral choice between risky and not risky sector

euler_operator <- function(grid, w, y){
  interp <- approxfun(grid, w, rule = 2)
  Tfx <- function(x) {
    r <- max(c(beta*mean(u(interp(x - interp(x) + y*(x - interp(x))^alpha))), u(x)))
    #r <- max(c(beta*mean(u(interp(x - interp(x) + y))), u(x)))
    r <- r^(-1/gamma)
    return(r)
  }
  Tw <- sapply(grid, Tfx)
  return(Tw)
}

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
V1 <- EFI(xgrid, v0, y1)
V2 <- EFI(xgrid, v0, y2)

dfplot <- function(V){
  V <- as.data.frame(V)
  V$grid <- xgrid
  colnames(V) <- c(1:(dim(V)[2] - 1), "grid")
  V <- pivot_longer(V, cols = -grid, names_to = "iteration", values_to = "fx")
  return(V)
}

V <- dfplot(V)
V1 <- dfplot(V1)
V2 <- dfplot(V2)

ggplot(filter(V, grid!=1 & iteration!=1)) +
  geom_line(aes(x = grid, y = fx, color = as.integer(iteration), group = as.integer(iteration))) +
  xlim(0,100000) + ylim(0,100000)

Vfinall <- rbind(filter(V, iteration==max(as.numeric(V$iteration))) %>% mutate(meansd = "62k, 16k"),
                 filter(V1, iteration==max(as.numeric(V1$iteration))) %>% mutate(meansd = "64k, 22k"),
                 filter(V2, iteration==max(as.numeric(V2$iteration))) %>% mutate(meansd = "172k, 62k"))

ggplot(Vfinall) +
  geom_line(aes(x = grid, y = fx, color = meansd, group = meansd)) +
  geom_abline(aes(slope = 1, intercept = 0), color = "grey", alpha = .2) +
  theme_bw() + labs(x = "Income + Assets", y = "Consumption", color = "mean, sd") +
  xlim(0, 500000) + theme(text = element_text(size=20)) +
  geom_vline(aes(xintercept = exp(11 + .25^2/2)), color = "grey68") + 
  geom_vline(aes(xintercept = exp(11 + .35^2/2)), color = "grey68") + 
  geom_vline(aes(xintercept = exp(12 + .35^2/2)), color = "grey68")
ggsave("/Users/mdgordo/Dropbox (Yale_FES)/Nepal/presentation_materials/figures/fxcons.png")

### Effect of aid on consumption

aidfx <- approxfun(xgrid, Vfinall$fx[Vfinall$meansd=="62k, 16k"])
aidfx1 <- approxfun(xgrid, Vfinall$fx[Vfinall$meansd=="64k, 22k"])
aidfx2 <- approxfun(xgrid, Vfinall$fx[Vfinall$meansd=="172k, 62k"])

a <- aidfx(xgrid + 100000) - aidfx(xgrid)
a1 <- aidfx1(xgrid + 100000) - aidfx1(xgrid)
a2 <- aidfx2(xgrid + 100000) - aidfx2(xgrid)

ggplot() +
  geom_line(aes(x = xgrid, y = a), color = "green") + 
  geom_line(aes(x = xgrid, y = a1), color = "blue") + 
  geom_line(aes(x = xgrid, y = a2), color = "red") +
  theme_bw() + labs(x = "Income + assets", y = "Increased consumption", title = "Predicted effect of 100k npr transfer") +
  theme(text = element_text(size=20)) + xlim(0, 1000000)
ggsave("/Users/mdgordo/Dropbox (Yale_FES)/Nepal/presentation_materials/figures/fxaid.png")

### Effect of aid on utility - do for different Ys
v0 <- rep(0, length(xgrid))

u <- function(x){
  u <- x^(1-gamma)/(1-gamma)
  return(u)
}

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

V <- VFI(xgrid, v0, y)
V1 <- VFI(xgrid, v0, y1)
V2 <- VFI(xgrid, v0, y2)

V <- dfplot(V)
V1 <- dfplot(V1)
V2 <- dfplot(V2)

ggplot(filter(V, grid > 10000)) +
  geom_line(aes(x = grid, y = fx, color = as.integer(iteration), group = as.integer(iteration)))

Vfinall <- rbind(filter(V, iteration==max(as.numeric(V$iteration))) %>% mutate(meansd = "62k, 16k"),
                 filter(V1, iteration==max(as.numeric(V1$iteration))) %>% mutate(meansd = "64k, 22k"),
                 filter(V2, iteration==max(as.numeric(V2$iteration))) %>% mutate(meansd = "172k, 62k"))

ggplot(filter(Vfinall, grid > 10000)) +
  geom_line(aes(x = grid, y = fx, color = meansd, group = meansd)) +
  theme_bw() + labs(x = "Income + Assets", y = "V(x)", color = "Mean, Std Dev") +
  xlim(0, 1000000)

### Benefits of Aid

vfxaid <- approxfun(xgrid, Vfinall$fx[Vfinall$meansd=="62k, 16k"])
vfxaid1 <- approxfun(xgrid, Vfinall$fx[Vfinall$meansd=="64k, 22k"])
vfxaid2 <- approxfun(xgrid, Vfinall$fx[Vfinall$meansd=="172k, 62k"])

r <- vfxaid(xgrid + 100000) - vfxaid(xgrid)
r1 <- vfxaid1(xgrid + 100000) - vfxaid1(xgrid)
r2 <- vfxaid2(xgrid + 100000) - vfxaid2(xgrid)

ggplot() +
  geom_line(aes(x = xgrid[xgrid>100000], y = r[xgrid>100000]), color = "red") +
  geom_line(aes(x = xgrid[xgrid>100000], y = r1[xgrid>100000]), color = "green") +
  geom_line(aes(x = xgrid[xgrid>100000], y = r2[xgrid>100000]), color = "blue") +
  theme_bw() + labs(x = "Income + Assets", y = "Benefits of 100,000 in Aid", color = "Mean, Std Dev") + 
  xlim(0, 1000000)

### Benefits of change in consumption compared to overall benefits (not sure this is right)
uc <- u(aidfx(xgrid + 100000)) - u(aidfx(xgrid))

df.bens <- data.frame("x" = rep(xgrid, 2),
                      "y" = c(r, uc),
                      "b" = c(rep("total", length(xgrid)), rep("consumption", length(xgrid))))

ggplot(filter(df.bens, x>25000)) +
  geom_line(aes(x = x, y = y, group = b, color = b)) + 
  theme_bw() + labs(x = "Income + Assets", y = "V(x), U(c)", color = "", title = "Benefits of 100k npr transfer") + xlim(0, 1000000) + 
  theme(text = element_text(size = 24), axis.text.y = element_blank())
ggsave("/Users/mdgordo/Dropbox (Yale_FES)/Nepal/presentation_materials/figures/resilience_benefits.png", height = 12, width = 14)

uc1 <- u(aidfx1(xgrid + 100000)) - u(aidfx1(xgrid))

ggplot() +
  geom_line(aes(x = xgrid[xgrid>25000], y = r1[xgrid>25000]), color = "red") + 
  geom_line(aes(x = xgrid[xgrid>25000], y = uc1[xgrid>25000]), color = "blue") +
  theme_bw() + labs(x = "Income + Assets", y = "V(x), U(c)", title = "Full Benefits of Aid compared to Benefits of Current Period Consumption") +
  xlim(0, 1000000)

uc2 <- u(aidfx2(xgrid + 100000)) - u(aidfx2(xgrid))

ggplot() +
  geom_line(aes(x = xgrid[xgrid>25000], y = r2[xgrid>25000]), color = "red") + 
  geom_line(aes(x = xgrid[xgrid>25000], y = uc2[xgrid>25000]), color = "blue") +
  theme_bw() + labs(x = "Income + Assets", y = "V(x), U(c)", title = "Full Benefits of Aid compared to Benefits of Current Period Consumption") +
  xlim(0, 1000000)

