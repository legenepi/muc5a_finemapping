#usage: Rscript susie.R <genotype matrix> <z scores> <plot name>
#remotes::install_github("stephenslab/susieR")

args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(data.table)
library(susieR)

print("Reading data...")

data <- as.matrix(fread(args[1], header=FALSE))

print("Data read.")

print("Creating correlation matrix...")

data.cor = cor(data)

print("Correlation matrix complete.")

sumstat <- read.table(args[2])
colnames(sumstat) <- c("")

print("Running Susie...")

#coverage default is 0.95:
fitted_rss <- susie_rss(bhat = sumstat $V1, shat = sumstat$V2, R = data.cor, L = 10, n = as.numeric(args[3]), coverage=0.90)

write.table(summary(fitted_rss),args[4],row.names=FALSE,quote=FALSE,col.names=T,sep="\t")

jpeg(args[5])

susie_plot(fitted_rss, y="PIP")

dev.off()