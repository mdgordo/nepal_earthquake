### packages and data
library(tidyverse)
library(parallel)
library(nloptr)
library(fastGHQuad)
library(flatlandr)
library(gmm)
library(Rearrangement)
library(mgcv)
library(neuralnet)
options('nloptr.show.inequality.warning'=FALSE)
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))

df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500) 

df.hh <- df.hh %>%
  filter(designation !="none") %>%
  select(hhid, wave, wt_hh, food_consumption, total_income, var_inc, avg_inc, 
         home_value, home_investment, imputed_bufferstock, quake_aid, M_avg) %>% 
  mutate(home_value = home_value/M_avg,
         imputed_bufferstock = imputed_bufferstock/M_avg,
         food_consumption = food_consumption/M_avg,
         home_investment = home_investment/M_avg,
         quake_aid = quake_aid/M_avg,
         total_income = total_income/M_avg) %>%
  group_by(hhid) %>%
  mutate(lag_h = lag(home_value, order_by = wave),
         lag_x = lag(imputed_bufferstock, order_by = wave),
         lag_y = lag(total_income, order_by = wave),
         lag_c = lag(food_consumption, order_by = wave),
         lag_i = lag(home_investment, order_by = wave)) %>%
  filter(avg_inc>0 & food_consumption>0 & lag_h + home_investment>0)
df.hh <- df.hh[complete.cases(df.hh),]

### Grid size
xn = 40; hn = 40

### generate starting points - run 35 versions of this on cluster using dSQ
### dsq --job-file ~/project/HPC_WBHRVS_DSQ.txt -c 20 --mem-per-cpu 2g -t 24:00:00 --mail-type ALL

slurmseed <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(slurmseed); seeds = sample(c(1:1e6),9)
set.seed(seeds[1]); gamma <- runif(100, 1.01, 10)
set.seed(seeds[2]); beta <- runif(100, .6, .99)
set.seed(seeds[3]); R <- runif(100, 1, 1.4)
set.seed(seeds[4]); cbar <- runif(100, .01, .9)
set.seed(seeds[5]); hbar <- runif(100, .01, .9)
set.seed(seeds[6]); lambda <- runif(100, .01, 10)
set.seed(seeds[7]); sigma <- runif(100, .01, 1)
set.seed(seeds[8]); alpha <- runif(100, .1, .99)
set.seed(seeds[9]); delta <- runif(100, .5, 1)

iterations <- data.frame(gamma, beta, R, cbar, hbar,
                         lambda, sigma, alpha, delta)
iterations = filter(iterations, R*beta<1)

globalwrap <- function(ridx){
  theta = as.vector(t(iterations[ridx, ]))
  momentmat = gmmmomentmatcher(theta, df.hh)
  g = sum(colMeans(momentmat^2)); print(g)
  return(g)
}

f <- sapply(c(1:nrow(iterations)), globalwrap)
iterations$f <- f

saveRDS(iterations, paste(getwd(), "/data/model_output/iterations", slurmseed, ".rds", sep = ""))

### after all jobs finish, run
# setwd("data/model_output")
# ifiles <- list.files(pattern = "^iter")
# ilist <- lapply(ifiles, readRDS)
# iterations <- do.call(rbind, ilist)
# saveRDS(iterations, iterations.rds")
