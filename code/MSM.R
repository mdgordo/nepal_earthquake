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
library(AER)
library(fixest)
library(lhs)
options('nloptr.show.inequality.warning'=FALSE)
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))

df.adj <- readRDS(paste(getwd(), "/data/model_output/df.adj.rds", sep = ""))
df.adj <- df.adj[complete.cases(df.adj),]

### Grid size
xn = 50; hn = 50

iterationsr <- readRDS(paste(getwd(), "/data/model_output/iterations_raw.rds", sep = ""))
iframes <- split(iterationsr, rep(c(1:30), each = ceiling(nrow(iterationsr)/30), length.out = nrow(iterationsr)))

### dsq --job-file ~/project/HPC_WBHRVS_DSQ.txt -c 20 --mem-per-cpu 2g -t 24:00:00 --mail-type ALL
slurmseed <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
iterations <- iframes[[slurmseed+1]]

#gpenalty = 27946515 ## max iterations value

globalwrap <- function(ridx){ #chg back to t
  #if (t[2]*t[3] > 1 | t[4] > t[8] | t[5] > 1 - t[8]){
  #  return(gpenalty)
  #} else {
    t = as.vector(t(iterations[ridx, ]))
    momentmat = gmmmomentmatcher(t, df.adj)
    g = colSums(momentmat); print(g)
    return(g) ## change this to sum(g^2)
  #}
}

#sol <- crs2lm(x0 = c(1.2, .926, 1.054, .33, .19, .08, .42, .63, .93), 
#              fn = globalwrap, 
#              lower = c(1.01, .6, 1, .01, .01, .01, .01, .1, .5), 
#              upper = c(10, .99, 1.2, .9, .9, 10, 1.5, .99, 1),
#              xtol_rel = 1e-2)
#sol

f <- lapply(c(1:nrow(iterations)), globalwrap)
f <- do.call(rbind, f)
colnames(f) <- paste("m", seq(1,11,1), sep = "_")
iterations <- cbind(iterations, f)

saveRDS(iterations, paste(getwd(), "/data/model_output/iterations", slurmseed, ".rds", sep = ""))

### after all jobs finish, run
# setwd("data/model_output")
# ifiles <- list.files(pattern = "^iter")
# ifiles <- ifiles[ifiles!="iterations_raw.rds"]
# ilist <- lapply(ifiles, readRDS)
# iterations <- do.call(rbind, ilist)
# saveRDS(iterations, "iterations.rds")
