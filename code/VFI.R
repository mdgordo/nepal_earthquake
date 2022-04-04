### packages and data
library(tidyverse)
library(parallel)
library(nloptr)
library(flatlandr)
library(Rearrangement)
library(fastGHQuad)
library(rdrobust)
options('nloptr.show.inequality.warning'=FALSE)
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))
source(paste(getwd(), "/code/rdhelpers.R", sep = ""))

### Read in final Value function + statespace
aidamt = 300000 
theta = readRDS(paste(getwd(), "/data/model_output/theta.rds", sep = ""))
V = readRDS(paste(getwd(), "/data/model_output/V.rds", sep = ""))
finalV = V[[length(V)]]
#finalVr = lapply(ygrid, monotonizer, w = finalV)
#finalVr = do.call(rbind, finalVr)
#V[[length(V)+1]] <- finalVr

gamma = theta[1]; beta = theta[2]; R = theta[3]; cbar = theta[4]; hbar = theta[5]
lambda = theta[6]; sigma = theta[7]; alpha = theta[8]; delta = theta[9]
ygrid = unique(finalV$y)
gqpts = gaussHermiteData(31)

## Checks and plots
tolchecker(V)
do.call(rbind, lapply(ygrid, satisficer, w = V, sigma = theta[7]))

iterplotter(V, hidx = finalV$h[1], yidx = finalV$y[1], xfilt = finalV$x[1], var = "cfx", ifilt = c(3:4))
iterplotter(V, xidx = finalV$x[5651], yidx = finalV$y[5851]) + xlim(1e5, 6e9) + ylim(-3e-19,0)

#library(rayshader)
threedplotter(V, 10, "Tw", d3 = FALSE, lbm = sort(unique(finalV$B))[10])
pcfx = threedplotter(V, 10, "cfx", d3 = FALSE, lbm = sort(unique(finalV$B))[10])  + theme(text = element_text(size = 17)) 
pifx = threedplotter(V, 10, "ifx", d3 = FALSE, lbm = sort(unique(finalV$B))[10]) + theme(text = element_text(size = 17)) 
cowplot::plot_grid(pcfx, pifx)
ggsave(paste(getwd(), "/presentation_materials/policyfx.png", sep = ""), width = 16, height = 7)

yplotter(V, 0, d3 = FALSE, yvals = c(log(50000), log(100000)), xvals = c(-5000, 500000)) + 
  geom_contour(aes(x = x, y = y, z = -1*log(-proj)), color = "grey30") +
  labs(y = "Avg Income", fill = "V(x,mu)") + scale_fill_viridis_c(labels = NULL) + theme(text = element_text(size = 20)) 
ggsave(paste(getwd(), "/presentation_materials/wtpcurves.png", sep = ""), width = 10, height = 7)


### Policy Induced Surface

p = threedplotter(V, 5, "aidcfx", d3=FALSE, aidamt = aidamt) 
p1 = p +
  labs(fill = "Change in Consumption") + theme_bw() + theme(text = element_text(size = 17)) 
p = threedplotter(V, 10, "aidifx", d3=FALSE, aidamt = aidamt)
p2 = p + 
  labs(fill = "Change in Housing Investment") + theme_bw() + theme(text = element_text(size = 17)) 
cowplot::plot_grid(p1, p2)
ggsave(paste(getwd(), "/presentation_materials/policysurf.png", sep = ""), width = 16, height = 7)

### Read in Data

df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500) 

b <- optbw("consumption", b = "dist_2_seg13", df = df.hh, fuzzy = TRUE, dist.exclude = "none")[1,3]

df.hh <- df.hh %>% 
  filter(abs(dist_2_seg13) < b & designation !="none") %>%
  select(hhid, wave, wt_hh, food_consumption, total_income, var_inc, avg_inc, home_value, home_investment, 
         imputed_bufferstock, quake_aid, quake_losses, NGO_transfers)  %>%
  group_by(hhid) %>%
  mutate(lag_h = lag(home_value, order_by = wave)) %>% ungroup() %>%
  mutate(x_noaid = imputed_bufferstock - quake_aid,
         x_aid = x_noaid + aidamt,
         lag_h = if_else(is.na(lag_h), home_value/delta - home_investment, lag_h)) %>%  ### this gives post earthquake housing/could estimate pre quake w/losses
  mutate(lag_h = if_else(lag_h<0, 0, lag_h)) %>%
  filter(avg_inc>0 & food_consumption>0 & lag_h + home_investment>0) 

df.eg = filter(df.hh, hhid %in% c(735, 2927)) %>%
  mutate(wave = as.integer(wave), hhid = as.factor(hhid))

eg1 <- ggplot(df.eg) + geom_point(aes(x = wave, y = food_consumption, color = as.factor(hhid)), size = 3) + 
  theme_bw() + theme(legend.position = "none", text = element_text(size = 20)) + 
  scale_x_continuous(breaks = c(1,2,3)) + scale_color_manual(values = c("red", "blue"))
eg2 <- ggplot(df.eg) + geom_point(aes(x = wave, y = total_income, color = as.factor(hhid)), size = 3) + 
  theme_bw() + theme(legend.position = "none", text = element_text(size = 20)) + 
  scale_x_continuous(breaks = c(1,2,3)) + scale_color_manual(values = c("red", "blue"))
cowplot::plot_grid(eg1, eg2)
ggsave(paste(getwd(), "/presentation_materials/eghhs.png", sep = ""), width = 10, height = 7)

df.hh <- df.hh %>%
  drop_na(x_noaid, x_aid, avg_inc, lag_h)
#df.hh <- filter(df.hh, quake_losses>0)

ggplot(filter(df.hh, wave==1)) +
  geom_point(aes(x = imputed_bufferstock, y = quake_losses))

ggplot(df.hh) +
  geom_point(aes(x = imputed_bufferstock, y = home_value))

ggplot(df.hh) +
  geom_point(aes(x = imputed_bufferstock, y = avg_inc))

### Quantile Reg moments
Vfx = V[[length(V)]]

cfacts <- mclapply(c(1:nrow(df.hh)), counterfactualizer, vfx = Vfx, df.hh = df.hh, 
                   mc.cores = detectCores())
cfacts <- do.call(rbind, cfacts)
colnames(cfacts) <- c("Consumption_noAid", "Consumption_Aid", "Investment_noAid", "Investment_Aid",
                      "Borrowing_noAid", "Borrowing_Aid")

pdfmaker <- function(var, ctpts){
  xvar = as.vector(cfacts[,var])
  pvec = sapply(ctpts, function(x) weighted.mean(xvar<x, df.hh$wt_hh, na.rm = TRUE))
  return(pvec)
}

rearranger <- function(q, y) {
  r <- Rearrangement::rearrangement(x = as.data.frame(q), y)
  return(as.vector(r))
}

q = c(0, 55000, 69000, 82000, 93000, 104000, 117000, 133000, 154000, 187000)
simdat = data.frame(q = q,
                    true_aid = rearranger(q, c(0, 0, 0, 0, 0, 0, .2, .45, .57, .76)),
                    true_noaid = rearranger(q, c(0, .15, .27, .36, .52, .61, .6, .85, .8, .92)),
                    sim_noaid = pdfmaker("Consumption_noAid", q),
                    sim_aid = pdfmaker("Consumption_Aid", q)) %>%
  pivot_longer(-q, names_to = "series", values_to = "pdf") %>%
  mutate(aid = if_else(series %in% c("sim_aid", "true_aid"), "aid", "none"),
         data = if_else(series %in% c("sim_aid", "sim_noaid"), "simulated", "actual"))

psim1 = ggplot(simdat) +
  geom_line(aes(x = q, y = pdf, color = aid, linetype = data)) +
  theme_bw() + theme(text = element_text(size = 20)) + labs(x = "food consumption", y = "P(X<x)", fill = "data")

q = c(0, 15000, 135000, 500000)
simdat = data.frame(q = q,
                    true_aid = rearranger(q, c(.68, .74, .77, .9)),
                    true_noaid = rearranger(q, c(.91, 1, 1, 1)),
                    sim_noaid = pdfmaker("Investment_noAid", q),
                    sim_aid = pdfmaker("Investment_Aid", q)) %>%
  pivot_longer(-q, names_to = "series", values_to = "pdf") %>%
  mutate(aid = if_else(series %in% c("sim_aid", "true_aid"), "aid", "none"),
         data = if_else(series %in% c("sim_aid", "sim_noaid"), "simulated", "actual"))

psim2 = ggplot(simdat) +
  geom_line(aes(x = q, y = pdf, color = aid, linetype = data)) +
  theme_bw()  + theme(text = element_text(size = 20)) + labs(x = "home investment", y = "P(X<x)", fill = "data") + ylim(.5,1.1)

cowplot::plot_grid(psim1, psim2)
ggsave(paste(getwd(), "/presentation_materials/simulations.png", sep = ""), width = 12, height = 7)

q = c(-182000, -82000, -37000, -10000, 0, 5000, 20000, 50000, 121000)
simdat = data.frame(q = q,
                    true_aid = rearranger(q, c(0, 0, .12, .19, .34, .47, .6, .72, .83)),
                    true_noaid = rearranger(q, c(.09, .33, .36, .40, .45, .59, .74, .78, 1)),
                    sim_noaid = pdfmaker("Borrowing_noAid", q),
                    sim_aid = pdfmaker("Borrowing_Aid", q)) %>%
  pivot_longer(-q, names_to = "series", values_to = "pdf") %>%
  mutate(aid = if_else(series %in% c("sim_aid", "true_aid"), "aid", "none"),
         data = if_else(series %in% c("sim_aid", "sim_noaid"), "simulated", "actual"))

ggplot(simdat) +
  geom_line(aes(x = q, y = pdf, color = aid, linetype = data)) +
  theme_bw()  + theme(text = element_text(size = 20)) + labs(x = "borrowing", y = "P(X<x)", fill = "data")

### WTP for aid
df.hh$wtp <- mcmapply(wtpeliciter, x = df.hh$x_noaid, y = df.hh$avg_inc, h = df.hh$lag_h, 
                    aidamt = aidamt, vfx = list(Vfx), statespace = list(statespace), 
                    mc.cores = detectCores()-2, mc.preschedule = FALSE)
df.hh$tau <- df.hh$wtp * (1-beta)/exp(df.hh$avg_inc + (sigma*df.hh$avg_inc)^2/2)
#df.hh$wtp <- as.numeric(df.hh$wtp)

#df.hh <- filter(df.hh, wtp < quantile(df.hh$wtp, .99, na.rm = TRUE))

w1 = ggplot() +
  geom_histogram(aes(x = df.hh$wtp)) + xlim(0, 2e6) + theme_bw() + theme(text = element_text(size = 20)) + labs(x = "HH Surplus")
w2 = ggplot() +
  geom_histogram(aes(x = df.hh$tau)) + xlim(0,1) + theme_bw() + theme(text = element_text(size = 20)) + labs(x = "tau")
cowplot::plot_grid(w1, w2)
ggsave(paste(getwd(), "/presentation_materials/wtphist.png", sep = ""), width = 12, height = 7)

ggplot(df.hh) + 
  geom_point(aes(x = food_consumption, y = wtp, color = avg_inc)) +
  theme_bw() + scale_color_viridis_c()
ggplot(df.hh) + 
  geom_point(aes(x = home_value, y = wtp, color = avg_inc)) +
  theme_bw() + scale_color_viridis_c()
ggplot(df.hh) + 
  geom_point(aes(x = imputed_bufferstock, y = wtp, color = avg_inc)) +
  theme_bw() + scale_color_viridis_c()
ggplot(df.hh) + 
  geom_point(aes(x = quake_losses, y = wtp, color = avg_inc)) +
  theme_bw() + scale_color_viridis_c()

cor(df.hh$quake_losses, df.hh$wtp)
cor(df.hh$food_consumption, df.hh$wtp)

a = weighted.mean(df.hh$wtp, df.hh$wt_hh, na.rm = TRUE)
ct_act = weighted.mean(df.hh$wtp[df.hh$quake_aid>0], w = df.hh$wt_hh[df.hh$quake_aid>0], na.rm = TRUE)/a
ct_NGO = weighted.mean(df.hh$wtp[df.hh$NGO_transfers>0], w = df.hh$wt_hh[df.hh$NGO_transfers>0], na.rm = TRUE)/a

qd <- quantile(df.hh$quake_losses, 1-mean(df.hh$quake_aid>0))
ct_losses = weighted.mean(df.hh$wtp[df.hh$quake_losses>qd], w = df.hh$wt_hh[df.hh$quake_losses>qd], na.rm = TRUE)/a
qc <- quantile(df.hh$food_consumption, mean(df.hh$quake_aid>0))
weighted.mean(df.hh$wtp[df.hh$food_consumption<qc], w = df.hh$wt_hh[df.hh$food_consumption<qc], na.rm = TRUE)/a
qb <- quantile(df.hh$imputed_bufferstock, mean(df.hh$quake_aid>0))
weighted.mean(df.hh$wtp[df.hh$imputed_bufferstock<qb], w = df.hh$wt_hh[df.hh$imputed_bufferstock<qb], na.rm = TRUE)/a
qi <- quantile(df.hh$avg_inc, mean(df.hh$quake_aid>0))
weighted.mean(df.hh$wtp[df.hh$avg_inc<qi], w = df.hh$wt_hh[df.hh$avg_inc<qi], na.rm = TRUE)/a

wtpmod <- lm(wtp ~ poly(food_consumption,2) + poly(quake_losses,2) + poly(home_value,2), 
             df.hh, weights = wt_hh)

df.hh$wtppred <- wtpmod$fitted.values
qp <- quantile(df.hh$wtppred, mean(df.hh$quake_aid>0))
weighted.mean(df.hh$wtp[df.hh$wtppred<qp], w = df.hh$wt_hh[df.hh$wtppred<qp], na.rm = TRUE)/a

### Plot WTP surface
yspace = filter(statespace, y==ygrid[7])
yspace$wtp = mcmapply(wtpeliciter, x = yspace$x, y = yspace$y, h = yspace$h, 
                        aidamt = aidamt, vfx = list(Vfx), statespace = list(statespace),
                      mc.cores = detectCores()-2)

wtpfunk = interpolater.creater(ygrid[7], yspace, yspace$wtp)
df.testm = crossing(x = seq(-lambda*ygrid[7], 5e5, length.out = 200), h = seq(0, 1e5, length.out = 200))
df.testm$projwtp = mcmapply(wtpfunk, x = df.testm$x, h = df.testm$h, mc.cores = detectCores())
df.testm$y = ygrid[7]

yspace = filter(statespace, y==ygrid[5])
yspace$wtp = mcmapply(wtpeliciter, x = yspace$x, y = yspace$y, h = yspace$h, 
                      aidamt = aidamt, vfx = list(Vfx), statespace = list(statespace),
                      mc.cores = detectCores()-2)

wtpfunk = interpolater.creater(ygrid[5], yspace, yspace$wtp)
df.testl = crossing(x = seq(-lambda*ygrid[5], 5e5, length.out = 200), h = seq(0, 1e5, length.out = 200))
df.testl$projwtp = mcmapply(wtpfunk, x = df.testl$x, h = df.testl$h, mc.cores = detectCores())
df.testl$y = ygrid[5]

yspace = filter(statespace, y==ygrid[9])
yspace$wtp = mcmapply(wtpeliciter, x = yspace$x, y = yspace$y, h = yspace$h, 
                      aidamt = aidamt, vfx = list(Vfx), statespace = list(statespace),
                      mc.cores = detectCores()-2)

wtpfunk = interpolater.creater(ygrid[9], yspace, yspace$wtp)
df.testh = crossing(x = seq(-lambda*ygrid[9], 5e5, length.out = 200), h = seq(0, 1e5, length.out = 200))
df.testh$projwtp = mcmapply(wtpfunk, x = df.testh$x, h = df.testh$h, mc.cores = detectCores())
df.testh$y = ygrid[9]

plow = ggplot(df.testl) + geom_tile(aes(x = x, y = h, fill = log(projwtp))) + 
  labs(fill = "log(surplus)") + 
  scale_fill_viridis_c(limits = c(11,14)) + theme_bw()  + theme(text = element_text(size = 20)) 
pmed = ggplot(df.testm) + geom_tile(aes(x = x, y = h, fill = log(projwtp))) + 
  labs(fill = "log(surplus)") + scale_fill_viridis_c() + 
  #scale_fill_viridis_c(limits = c(11,14)) + 
  theme_bw() + theme(text = element_text(size = 20)) 
phigh = ggplot(df.testh) + geom_tile(aes(x = x, y = h, fill = log(projwtp))) + 
  labs(fill = "log(surplus)") + 
  scale_fill_viridis_c(limits = c(11,14)) + theme_bw()  + theme(text = element_text(size = 20)) 
cowplot::plot_grid(plow, pmed, phigh, nrow =1)
pmed
ggsave(paste(getwd(), "/presentation_materials/wtpsurf.png", sep = ""), width = 12, height = 7)
