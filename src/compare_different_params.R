#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(readxl)
library(VennDiagram)
library(gplots)

args = commandArgs(trailingOnly=TRUE)

df1 <- fread(args[1])
df1_type <- as.factor(args[2])

df2 <- fread(args[3])
df2_type <- as.factor(args[4])

df3 <- fread(args[5])
df3_type <- as.factor(args[6])

df4 <- fread(args[7])
df4_type <- as.factor(args[8])

venn.diagram(
   x = list(
     df1 %>% select(BP) %>% distinct() %>% unlist(),
     df2 %>%  select(BP) %>% distinct() %>% unlist(),
     df3 %>%  select(BP) %>% distinct() %>% unlist(),
     df4 %>%  select(BP) %>% distinct() %>% unlist()
    ),
   category.names = c(df1_type, df2_type, df3_type, df4_type),
   filename = args[9],
   output = TRUE ,
           imagetype="png" ,
           height = 1200 ,
           width = 1500 ,
           resolution = 400,
           compression = "lzw",
           lwd = 1,
           col=c("#440154ff", '#21908dff', '#fde725ff', 'Dark Olive Green 4'),
           fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3), alpha('Dark Olive Green 4',0.3)),
           cex = 0.5,
           fontfamily = "sans",
           cat.cex = 0.3,
           cat.default.pos = "outer")

#use venn() in package gplots to obtain the list of participants for each intersection:
listInput <- list(df1 %>% select(BP) %>% distinct() %>% unlist(),
     df2 %>%  select(BP) %>% distinct() %>% unlist(),
     df3 %>%  select(BP) %>% distinct() %>% unlist(),
     df4 %>%  select(BP) %>% distinct() %>% unlist())
ItemsList <- venn(listInput, show.plot = FALSE, category.names = c(df1_type, df2_type, df3_type, df4_type))
#write into an output file the intersection:
file.create(args[8])
for (i in seq(1,length(attributes(ItemsList)$intersections))) {
write.table(attributes(ItemsList)$intersections[i],args[10],append=T,row.names = FALSE, quote=FALSE)}
