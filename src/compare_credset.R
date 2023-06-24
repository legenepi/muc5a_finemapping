#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(readxl)
library(VennDiagram)
library(gplots)

chsp <- read_xlsx("/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/chronic_sputum_suppl.xlsx", sheet="S10-Credible sets", skip = 1, col_names=T)
modsev <- fread("/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/output/polyfun_susie_modsev_credset")
gbmi <- fread("/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/output/polyfun_susie_gbmi_eur_credset_Nmedian")
colnames(chsp)[2] <- "to_separate"
chsp <- chsp %>% separate(to_separate, c("CHR", "BP"), sep=":")
colnames(chsp)[5] <- "to_separate2"
chsp <- chsp %>% separate(to_separate2, c("A1", "A2"), sep="/")
chsp <- chsp %>% filter(Locus == "MUC2")
chsp <- chsp %>% rename(SNP = RSID)
chsp <- chsp %>% rename(PIP = 'Posterior probability')
#Select cols of interest:
chsp <- chsp %>% select("CHR","BP","SNP","A1","A2","P","PIP")
colnames(chsp) <- paste("chsp", colnames(chsp), sep = "_")
modsev <- modsev %>% select("CHR","BP","SNP","A1","A2","P","PIP")
colnames(modsev) <- paste("modsev", colnames(modsev), sep = "_")
gbmi <- gbmi %>% select("CHR","BP","SNP","A1","A2","P","PIP")
colnames(gbmi) <- paste("gbmi", colnames(gbmi), sep = "_")
#join the dataset to have a unique file with the results:
chsp <- chsp %>% rename(BP = chsp_BP)
chsp$BP <- as.numeric(chsp$BP)
modsev <- modsev %>% rename(BP = modsev_BP)
modsev$BP <- as.numeric(modsev$BP)
gbmi <- gbmi %>% rename(BP = gbmi_BP)
gbmi$BP <- as.numeric(gbmi$BP)

chsp <- chsp %>% rename(SNP = chsp_SNP)
modsev <- modsev %>% rename(SNP = modsev_SNP)
gbmi <- gbmi %>% rename(SNP = gbmi_SNP)

chsp <- chsp %>% rename(CHR = chsp_CHR)
chsp$CHR <- as.numeric(chsp$CHR)
modsev <- modsev %>% rename(CHR = modsev_CHR)
modsev$CHR <- as.numeric(modsev$CHR)
gbmi <- gbmi %>% rename(CHR = gbmi_CHR)
gbmi$CHR <- as.numeric(gbmi$CHR)

chsp_modsev <- full_join(chsp,modsev,by=c("BP","CHR","SNP"))
chsp_modsev_gbmi <- full_join(chsp_modsev,gbmi,by=c("BP","CHR","SNP"))
chsp_modsev_gbmi <- chsp_modsev_gbmi %>% arrange(BP)
write.table(chsp_modsev_gbmi,"/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/output/credset_3_gwas.txt", na="NA", sep=" ",quote=F,row.names=F)

venn.diagram(
   x = list(
     chsp %>% select(BP) %>% distinct() %>% unlist(),
     modsev %>%  select(BP) %>% distinct() %>% unlist(),
     gbmi %>%  select(BP) %>% distinct() %>% unlist()
    ),
   category.names = c("chronic sputum", "moderate-to-severe", "GBMI-all comer asthma"),
   filename = "/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/output/Venn_asthma_credset.png",
   output = TRUE ,
           imagetype="png" ,
           height = 900 ,
           width = 1000 ,
           resolution = 400,
           compression = "lzw",
           lwd = 1,
           col=c("#440154ff", '#21908dff', '#fde725ff'),
           fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
           cex = 0.5,
           fontfamily = "sans",
           cat.cex = 0.3,
           cat.default.pos = "outer")
#use venn() in package gplots to obtain the list of participants for each intersection:
listInput <- list(chsp %>% select(BP) %>% distinct() %>% unlist(),
     modsev %>%  select(BP) %>% distinct() %>% unlist(),
     gbmi %>%  select(BP) %>% distinct() %>% unlist())
ItemsList <- venn(listInput, show.plot = FALSE, category.names = c("chronic sputum", "moderate-to-severe", "GBMI-all comer asthma"))
#write into an output file the intersection:
file.create("${path_dir}/output/Intersection_asthma_credset_BP_list")
for (i in seq(1,length(attributes(ItemsList)$intersections))) {
write.table(attributes(ItemsList)$intersections[i],"/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/output/Intersection_asthma_credset_BP_list",append=T,row.names = FALSE, quote=FALSE)}
