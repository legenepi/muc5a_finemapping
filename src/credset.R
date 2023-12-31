#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
df <- fread(args[1])
head(df)
dim(df)
df$cumsum_pip <- cumsum(df$PIP)
credset <- df %>% filter(cumsum_pip <= 0.95)
dim(credset)
fwrite(credset,args[2],row.names=F,quote=F)

credset_bysusie <- df %>% filter(CREDIBLE_SET == 1)
fwrite(credset_bysusie,args[3],row.names=F,quote=F)