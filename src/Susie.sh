#!/bin/bash

#PBS -N SuSie_muc5ac
#PBS -j oe
#PBS -o SuSie_muc5ac
#PBS -l walltime=32:0:0
#PBS -l vmem=50gb
#PBS -l nodes=1:ppn=1
#PBS -d .
#PBS -W umask=022

path_dir="/home/n/nnp5/PhD/PhD_project/muc5a_finemapping"
tmp_data="/scratch/gen1/nnp5/muc5a_finemapping/tmp_data"

cd ${path_dir}

module load gcc/9.3
module unload R/4.2.1
module load R/4.1.0
module load plink2

#input data:
#bgen:
#data/gen1/AIRPROM/imputation/chr11/chr11.hap.impute2.chunk01.bgen
#.sample
#awk '{print $1, $2, $3}' /data/gen1/AIRPROM/assoc/severe_asthma/airprom_pheno_array.sample \
#    > ${tmp_data}/modsev_all.sample
#awk 'NR == 1 || NR == 2 || $24 == 0 || $24 == 1 {print $1, $2, $3}' /data/gen1/AIRPROM/assoc/severe_asthma/airprom_pheno_array.sample \
#    > ${tmp_data}/modsev_casecontrol.sample

#z_score file and snp_list file in python from input/modsev_sumstats_munged.parquet:
#start: 636478
#end: 1636478
#zcat /rfs/TobinGroup/kaf19/Public_data/ukb_mod_sev_asthma/Shrine_30552067_moderate-severe_asthma.txt.gz | \
#    awk -F "\t" '$2 == 11 && $3 >= 636478 && $3 <= 1636478 {print $1, $8, $9}' > ${tmp_data}/modsev_snp_beta_se.txt

#Exclude multi-allelic variants and find the common SNP IDs for the genotyped matrix and the zscore input files:
#use the file for each regions created by FINEMAP.sh:
#grep -v -w -F -f /data/gen1/UKBiobank_500K/imputed/multiallelic.snps \
#    ${tmp_data}/modsev_snp_beta_se.txt | awk '{print $2, $3}' \
#    > ${tmp_data}/modsev_noma_beta_se.txt

#grep -v -w -F -f /data/gen1/UKBiobank_500K/imputed/multiallelic.snps \
#    ${tmp_data}/modsev_snp_beta_se.txt  | awk '{print $1}' \
#    > ${tmp_data}/modsev_noma_snp.txt

#Format region data for input to R
#plink2 \
#    --bgen /data/gen1/AIRPROM/imputation/chr11/chr11.hap.impute2.chunk01.bgen ref-first \
#    --sample ${tmp_data}/modsev_all.sample \
#    --keep ${tmp_data}/modsev_casecontrol.sample \
#    --export A \
#    --extract ${tmp_data}/modsev_noma_snp.txt \
#    --out ${tmp_data}/ukb_modsev_snp

#cut -f7- ${tmp_data}/ukb_modsev_snp.raw \
#    > ${tmp_data}/ukb_modsev_snp.cols.raw

#awk 'NR>1 {print}' ${tmp_data}/ukb_modsev_snp.cols.raw \
#    > ${tmp_data}/ukb_modsev_snp.cols_noheader.raw

Rscript src/susie.R \
    ${tmp_data}/ukb_modsev_snp.cols_noheader.raw \
    ${tmp_data}/modsev_noma_beta_se.txt \
    ${path_dir}/output/susie_modsev_qsub.txt \
    ${path_dir}/output/susie_modsev_qsub.jpeg

















#if there is some NAs in the genotype matrix (keep it because it is a nice script)
#awk '{ lines[NR] = $0; for (i = 1; i <= NF; i++) if ($i == "NA") skip[i] = 1;} END { for (i = 1; i <= NR; i++) {
#    nf = split(lines[i], fields);
#    for (j = 1; j <= nf; j++) if (!(j in skip)) printf("%s ", fields[j]);
#    printf("\n");
#    }
#    }' /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header.raw  > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header_noNA.raw