#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
df <- fread("/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/gbmi_eur_sumstat")
df <- df %>% mutate(MAF=ifelse(df$A1_freq <= 0.5, df$A1_freq, 1-df$A1_freq))
#Rename N-sample size col:
colnames(df)[10] <- "N"
df <- df %>% filter(Chromosome == 11)
write.table(df,"/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/b38_chr11_gbmi_eur_sumstat_maf", na="NA", sep=" ",quote=F,row.names=F)