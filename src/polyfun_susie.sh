#!/bin/bash


module load python/gcc/3.9.10
conda activate polyfun
module load gcc/9.3

#install parquet-tools to read parquet file in bash: python3 -m pip install parquet-tools
#Step1. Create a munged summary statistics file in a PolyFun-friendly parquet format.

##MODERATE-TO-SEVERE:
#SNP	Chromosome	Position_b37	Coded Non_coded	Coded_freq	INFO	beta	SE_GC P_GC
zcat /rfs/TobinGroup/kaf19/Public_data/ukb_mod_sev_asthma/Shrine_30552067_moderate-severe_asthma.txt.gz | \
    sed 's/Position_b37/BP/g' | \
    sed 's/SE_GC/SE/g' | \
    sed 's/P_GC/P/g' | \
    sed 's/#SNP/SNP/g' | awk -F "\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13}' > /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/modsev_sumstat

sed -i 's/Coded/A1/g' /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/modsev_sumstat
sed -i 's/Non_coded/A2/g' /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/modsev_sumstat
sed -i 's/Coded_freq/A1_freq/g' /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/modsev_sumstat

#cretae MAF in R:
module unload R/4.2.1
module load R/4.1.0

#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
df <- fread("input/modsev_sumstat")
df <- df %>% mutate(MAF=ifelse(df$A1_freq <= 0.5, df$A1_freq, 1-df$A1_freq))
fwrite(df,"input/modsev_sumstat_maf",sep=" ",quote=F,row.names=F)


#Run actual polyfun:
##Create parquet file:
python /home/n/nnp5/software/polyfun/munge_polyfun_sumstats.py \
  --sumstats /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/modsev_sumstat_maf \
  --n 30810 \
  --out /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/modsev_sumstats_munged.parquet \
  --min-info 0.85 \
  --min-maf 0.01

#Approach 1 (UK Biobank pre-computed prior causal probabilities); no information for HLA region !
##A.Extract pre-computed prior causal probabilities in UK Biobank - White British:
python /home/n/nnp5/software/polyfun/extract_snpvar.py \
    --sumstats /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/modsev_sumstats_munged.parquet \
    --allow-missing \
    --out /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/modsev_SNPs_PriCauPro

#Sentinel variant: rs11603634
#BP37: 1136478
#+/- 500000
start=$((1136478-500000))
end=$((1136478+500000))

#Download pre-computed LD matrix from British ancestry UK Biobank individuals:
wget https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/chr11_1_3000001.npz
wget https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/chr11_1_3000001.gz

#Run the actual FIFO:
python /home/n/nnp5/software/polyfun/finemapper.py \
    --ld /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/chr11_1_3000001 \
    --sumstats /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/modsev_SNPs_PriCauPro \
    --n 30810 \
    --chr 11 \
    --start 636478 \
    --end 1636478 \
    --method susie \
    --max-num-causal 10 \
    --allow-missing \
    --out /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/output/polyfun_susie_modsev

#Extract credset 95% in R
library(tidyverse)
library(data.table)
df <- fread("output/polyfun_susie_modsev")
head(df)
dim(df)
df$cumsum_pip <- cumsum(df$PIP)
credset <- df %>% filter(cumsum_pip <= 0.95)
dim(credset)
fwrite(credset,"output/polyfun_susie_modsev_credset",row.names=F,quote=F)





#GBMI ALL-COMER ASTHMA IN EUROPEAN:
#/rfs/TobinGroup/kaf19/Public_data/gbmi_asthma/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz

##ALT is the coded allele:
#CHR	POS	REF	ALT	rsid	all_meta_AF	inv_var_meta_beta	inv_var_meta_sebeta	inv_var_meta_p	inv_var_het_p	direction	N_case	N_ctrl	n_dataset	n_bbk	is_strand_flip	is_diff_AF_gnomAD
zcat /rfs/TobinGroup/kaf19/Public_data/gbmi_asthma/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz | \
    awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12+$13}' | sed 's/#CHR/Chromosome/g' | \
    sed 's/POS/BP/g' | sed 's/REF/A2/g' | sed 's/ALT/A1/g' | sed 's/rsid/SNP/g' | \
    sed 's/all_meta_AF/A1_freq/g' | sed 's/inv_var_meta_beta/beta/g' | \
    sed 's/inv_var_meta_sebeta/SE/g' | sed 's/inv_var_meta_p/P/g' > /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/gbmi_eur_sumstat

#cretae MAF in R:
module unload R/4.2.1
module load R/4.1.0

#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
df <- fread("/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/gbmi_eur_sumstat")
df <- df %>% mutate(MAF=ifelse(df$A1_freq <= 0.5, df$A1_freq, 1-df$A1_freq))
colnames(df)[10] <- "N" 
fwrite(df,"/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/gbmi_eur_sumstat_maf",sep=" ",quote=F,row.names=F)


#Run actual polyfun:
##Create parquet file:
python /home/n/nnp5/software/polyfun/munge_polyfun_sumstats.py \
  --sumstats /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/gbmi_eur_sumstat_maf \
  --out /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/gbmi_eur_sumstats_munged.parquet \
  --min-info 0.85 \
  --min-maf 0.01

#Approach 1 (UK Biobank pre-computed prior causal probabilities); no information for HLA region !
##A.Extract pre-computed prior causal probabilities in UK Biobank - White British:
python /home/n/nnp5/software/polyfun/extract_snpvar.py \
    --sumstats /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/gbmi_eur_sumstats_munged.parquet \
    --allow-missing \
    --out /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/gbmi_eur_SNPs_PriCauPro

#rs11245964
#1113278
start=$((1113278 - 500000))
end=$((1113278 + 500000))

#Download pre-computed LD matrix from British ancestry UK Biobank individuals:
wget https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/chr11_1_3000001.npz
wget https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/chr11_1_3000001.gz


#I got different N sample size for each SNPs, so I use the R package implementation that allows me to use different N.
#Not able to do it with the command line tool.
#R
install.packages("remotes")
remotes::install_github("RajLabMSSM/echofinemap")
library(echofinemap)
##!! Not able to download this package ! problem with github-memory acces.

locus_dir <- "/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/output"
dat <- "/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/gbmi_eur_SNPs_PriCauPro"
LD_matrix <- "/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/chr11_1_3000001"
dat2 <- echofinemap::POLYFUN(locus_dir=locus_dir,
                             dat=dat,
                             LD_matrix = LD_matrix,
                             method="SUSIE",
                             max_causal = 10,
                             conda_env = "polyfun",
                             verbose = T)
#Run the actual FIFO:
#median N-sample size: 1376100
#max N-sample size: 1376100
#min N-sample size: 26133
#USING MEDIAN VALUE FOR n-SAMPLE SIZE:
python /home/n/nnp5/software/polyfun/finemapper.py \
    --ld /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/chr11_1_3000001 \
    --sumstats /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/gbmi_eur_SNPs_PriCauPro \
    --n 1376100 \
    --chr 11 \
    --start $start \
    --end $end \
    --method susie \
    --max-num-causal 10 \
    --allow-missing \
    --out /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/output/polyfun_susie_gbmi_eur

python /home/n/nnp5/software/polyfun/finemapper.py \
    --ld /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/chr11_1_3000001 \
    --sumstats /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/input/gbmi_eur_SNPs_PriCauPro \
    --n 26133 \
    --chr 11 \
    --start $start \
    --end $end \
    --method susie \
    --max-num-causal 10 \
    --allow-missing \
    --out /home/n/nnp5/PhD/PhD_project/muc5a_finemapping/output/polyfun_susie_gbmi_eur_Nmin

#Extract credset 95% in R:
library(tidyverse)
library(data.table)
df <- fread("/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/output/polyfun_susie_gbmi_eur")
head(df)
dim(df)
df$cumsum_pip <- cumsum(df$PIP)
credset <- df %>% filter(cumsum_pip <= 0.95)
dim(credset)
fwrite(credset,"/home/n/nnp5/PhD/PhD_project/muc5a_finemapping/output/polyfun_susie_gbmi_eur_credset_Nmedian",row.names=F,quote=F)


#For now, rs11600289 looks to be the one found by all three fine-mapping.