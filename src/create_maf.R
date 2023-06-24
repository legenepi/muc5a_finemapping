#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
df <- fread("/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/modsev_sumstat")
df <- df %>% mutate(MAF=ifelse(df$A1_freq <= 0.5, df$A1_freq, 1-df$A1_freq))
fwrite(df,"/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/modsev_sumstat_maf",sep=" ",quote=F,row.names=F)