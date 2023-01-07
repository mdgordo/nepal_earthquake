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

globalwrap <- function(ridx){ 
    t = as.vector(t(iterations[ridx, ]))
    momentmat = gmmmomentmatcher(t, df.adj)
    g = colSums(momentmat); print(g)
    return(g) 
}

f <- lapply(c(1:nrow(iterations)), globalwrap)
f <- do.call(rbind, f)
colnames(f) <- paste("m", seq(1,12,1), sep = "_")
iterations <- cbind(iterations, f)

saveRDS(iterations, paste(getwd(), "/data/model_output/iterations", slurmseed, ".rds", sep = ""))

### after all jobs finish, run
# setwd("data/model_output")
# ifiles <- list.files(pattern = "^iter")
# ifiles <- ifiles[ifiles!="iterations_raw.rds"]
# ilist <- lapply(ifiles, readRDS)
# iterations <- do.call(rbind, ilist)
# saveRDS(iterations, "iterations.rds")

#iterations <- readRDS(paste(getwd(), "/data/model_output/iterations.rds", sep = ""))
#iterations$obj <- apply(iterations[,c(10:21)], 1, function(x) sum(x^2))

#ggplot(iterations) + geom_point(aes(x = alpha, y = gamma, color = (m_1^2))) + scale_color_viridis_c()
