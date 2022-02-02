### packages and data
library(tidyverse)
library(parallel)
library(nloptr)
library(flatlandr)
library(rdrobust)
library(gmm)
options('nloptr.show.inequality.warning'=FALSE)
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))
source(paste(getwd(), "/code/rdhelpers.r", sep = ""))

df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500) 
#df.hh <- read_csv("/home/mdg59/project/WBHRVS/full_panel.csv", guess_max = 7500)

b <- optbw("consumption", b = "dist_2_seg13", df = df.hh, fuzzy = TRUE, dist.exclude = "none")[1,3]

df.hh <- df.hh %>%
  filter(abs(dist_2_seg13) < b & designation !="none") %>%
  select(hhid, wave, wt_hh, food_consumption, total_income, var_inc, avg_inc, home_value, home_investment, imputed_bufferstock, quake_aid)  %>%
  group_by(hhid) %>%
  mutate(lag_h = lag(home_value, order_by = wave)) 

### Shocks - Gaussian Quadrature points
gqpts = c(-2.651961, -1.673552, -.8162879, 0, .8162879, 1.673552, 2.651961)
gqwts = c(.0009717812, .0545155828, .4256072526, .8102646176, .4256072526, .0545155828, .0009717812)/sqrt(pi)

### Grid size
xn = hn = 30; yn = 15
ygrid = seq(8, 13, length.out = yn)

### Final Value Function
g = readRDS("g.rds")
theta = theta0
gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]

### statespace
Blist = lapply(ygrid, function(y) -lambda*exp(y))
cminlist = lapply(ygrid, function(y) cbar*exp(y))
hminlist = lapply(ygrid, function(y) hbar*exp(y))
statespace = lapply(c(1:length(ygrid)), create.statespace, Blist = Blist, cminlist = cminlist, 
                    hminlist = hminlist)
statespace = do.call(rbind, statespace)

### Initial guess 
v0 = guesser(statespace, theta)
V = VFI(v0, statespace, theta = theta)
V = compact(readRDS("/users/mdgordo/Desktop/w.rds"))

## Checks and plots
tolchecker(V)
do.call(rbind, lapply(ygrid, satisficer, w = V, sigma = theta[7]))

iterplotter(V, hidx = statespace$h[1], yidx = statespace$y[1])
iterplotter(V, hidx = statespace$h[10], yidx = statespace$y[10])
iterplotter(V, xidx = statespace$x[1], yidx = statespace$y[1])
iterplotter(V, xidx = statespace$x[31], yidx = statespace$y[31])

library(rayshader)
threedplotter(V, 1, "Tw", d3 = FALSE)
threedplotter(V, 1, "cfx", d3 = FALSE)
threedplotter(V, 1, "ifx", d3 = FALSE)

### Policy Induced Surface
aidamt = 250000

threedplotter(V, 1, "aidcfx", d3=FALSE, aidamt = aidamt)
threedplotter(V, 1, "aidifx", d3=FALSE, aidamt = aidamt)

### Quantile Reg moments
Vfx = V[[length(V)]]

df.hh <- mutate(df.hh, x_noaid = imputed_bufferstock - quake_aid,
                x_aid = x_noaid + aidamt,
                lag_h = if_else(is.na(lag_h), home_value/delta - home_investment, lag_h)) %>%
  drop_na(x_noaid, x_aid, avg_inc, lag_h) %>%
  mutate(lag_h = if_else(lag_h<0, 0, lag_h))

cfacts <- mclapply(c(1:nrow(df.hh)), counterfactualizer, vfx = Vfx, statespace = statespace, df.hh = df.hh, 
                   mc.cores = detectCores())
cfacts <- do.call(rbind, cfacts)
colnames(cfacts) <- c("Consumption_noAid", "Consumption_Aid", "Investment_noAid", "Investment_Aid",
                      "Borrowing_noAid", "Borrowing_Aid")

pdfmaker <- function(var, ctpts){
  xvar = as.vector(cfacts[,var])
  pvec = sapply(ctpts, function(x) weighted.mean(xvar<x, df.hh$wt_hh, na.rm = TRUE))
  return(pvec)
}

ctpts = seq(0, 200000, 20000)
ggplot() +
  geom_line(aes(x = ctpts, y = pdfmaker("Consumption_noAid", ctpts)), color = "red") +
  geom_line(aes(x = ctpts, y = pdfmaker("Consumption_Aid", ctpts)), color = "blue") + 
  theme_bw() + labs(x = "food consumption", y = "P(X<x)")

ggplot() +
  geom_line(aes(x = ctpts, y = pdfmaker("Investment_noAid", ctpts)), color = "red") +
  geom_line(aes(x = ctpts, y = pdfmaker("Investment_Aid", ctpts)), color = "blue") + 
  theme_bw() + labs(x = "home investment", y = "P(X<x)")

ctpts = seq(-200000, 200000, 20000)
ggplot() +
  geom_line(aes(x = ctpts, y = pdfmaker("Borrowing_noAid", ctpts)), color = "red") +
  geom_line(aes(x = ctpts, y = pdfmaker("Borrowing_Aid", ctpts)), color = "blue") + 
  theme_bw() + labs(x = "borrowing", y = "P(X<x)")

### WTP for aid
df.hh$wtp <- mapply(wtpeliciter, x = df.hh$x_noaid, y = log(df.hh$avg_inc), h = df.hh$lag_h, 
                    aidamt = aidamt, vfx = list(Vfx), statespace = list(statespace))

ggplot() +
  geom_histogram(aes(x = df.hh$wtp))

### Plot WTP surface
yspace = filter(statespace, y==ygrid[10])
yspace$wtp = mclapply(wtpeliciter, x = yspace$x, y = yspace$y, h = yspace$h, 
                        aidamt = aidamt, vfx = list(Vfx), statespace = list(statespace),
                      mc.cores = detectCores())

wtpfunk = interpolater.creater(ygrid[10], yspace, yspace$wtp)
df.test = crossing(x = seq(Blist[[10]], 20*exp(ygrid[10]), length.out = 200), h = seq(0, 20*exp(ygrid[10]), length.out = 200))
df.test$projwtp = mclapply(wtpfunk, x = df.test$x, h = df.test$h, mc.cores = detectCores())
p = ggplot(df.test) + geom_tile(aes(x = x, y = h, fill = projwtp)) + scale_fill_viridis_c()
plot_gg(p, multicore=TRUE)

