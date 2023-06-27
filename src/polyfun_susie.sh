#!/bin/bash

#set up local variable, environment to use:
path_dir="/home/n/nnp5/PhD/PhD_project/muc5a_finemapping"
module load python/gcc/3.9.10
conda activate polyfun
module load gcc/9.3
cd ${path_dir}/

#required install:
#install parquet-tools to read parquet file in bash: python3 -m pip install parquet-tools

##MODERATE-TO-SEVERE:
#Step1. Create a munged summary statistics file in a PolyFun-friendly parquet format.
#SNP	Chromosome	Position_b37	Coded Non_coded	Coded_freq	INFO	beta	SE_GC P_GC
zcat /rfs/TobinGroup/kaf19/Public_data/ukb_mod_sev_asthma/Shrine_30552067_moderate-severe_asthma.txt.gz | \
    sed 's/Position_b37/BP/g' | \
    sed 's/SE_GC/SE/g' | \
    sed 's/P_GC/P/g' | \
    sed 's/#SNP/SNP/g' | awk -F "\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $13}' > ${path_dir}/input/modsev_sumstat

sed -i 's/Coded/A1/g' ${path_dir}/input/modsev_sumstat
sed -i 's/Non_coded/A2/g' ${path_dir}/input/modsev_sumstat
sed -i 's/Coded_freq/A1_freq/g' ${path_dir}/input/modsev_sumstat

#cretae MAF in R:
module unload R/4.2.1
module load R/4.1.0
dos2unix ${path_dir}/src/create_maf.R
chmod o+x ${path_dir}/src/create_maf.R
Rscript ${path_dir}/src/create_maf.R

#Run Polyfun:
##Create parquet file:
python /home/n/nnp5/software/polyfun/munge_polyfun_sumstats.py \
  --sumstats ${path_dir}/input/modsev_sumstat_maf \
  --n 30810 \
  --out ${path_dir}/input/modsev_sumstats_munged.parquet \
  --min-info 0.85 \
  --min-maf 0.01

#Approach 1 (UK Biobank pre-computed prior causal probabilities); no information for HLA region !
##A.Extract pre-computed prior causal probabilities in UK Biobank - White British:
python /home/n/nnp5/software/polyfun/extract_snpvar.py \
    --sumstats ${path_dir}/input/modsev_sumstats_munged.parquet \
    --allow-missing \
    --out ${path_dir}/input/modsev_SNPs_PriCauPro

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
    --ld ${path_dir}/input/chr11_1_3000001 \
    --sumstats ${path_dir}/input/modsev_SNPs_PriCauPro \
    --n 30810 \
    --chr 11 \
    --start 636478 \
    --end 1636478 \
    --method susie \
    --max-num-causal 10 \
    --allow-missing \
    --out ${path_dir}/output/polyfun_susie_modsev

#Extract credset 95% in R
dos2unix ${path_dir}/src/credset.R
chmod o+x ${path_dir}/src/credset.R
Rscript ${path_dir}/src/credset.R \
    ${path_dir}/output/polyfun_susie_modsev \
    ${path_dir}/output/polyfun_susie_modsev_credset \
    ${path_dir}/output/polyfun_susie_modsev_credset_credsetbysusie


#####################################
#GBMI ALL-COMER ASTHMA IN EUROPEAN:
#/rfs/TobinGroup/kaf19/Public_data/gbmi_asthma/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz

##ALT is the coded allele:
#CHR	POS	REF	ALT	rsid	all_meta_AF	inv_var_meta_beta	inv_var_meta_sebeta	inv_var_meta_p	inv_var_het_p	direction	N_case	N_ctrl	n_dataset	n_bbk	is_strand_flip	is_diff_AF_gnomAD
zcat /rfs/TobinGroup/kaf19/Public_data/gbmi_asthma/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz | \
    awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $12+$13}' | sed 's/#CHR/Chromosome/g' | \
    sed 's/POS/BP/g' | sed 's/REF/A2/g' | sed 's/ALT/A1/g' | sed 's/rsid/SNP/g' | \
    sed 's/all_meta_AF/A1_freq/g' | sed 's/inv_var_meta_beta/beta/g' | \
    sed 's/inv_var_meta_sebeta/SE/g' | sed 's/inv_var_meta_p/P/g' > ${path_dir}/input/gbmi_eur_sumstat

##GRCh38: use liftover to convert in GRCh37

#Extract variants of interest and create MAF in R:
module unload R/4.2.1
module load R/4.1.0
dos2unix ${path_dir}/src/create_maf_gbmi.R
chmod o+x ${path_dir}/src/create_maf_gbmi.R
Rscript ${path_dir}/src/create_maf_gbmi.R



#convert sentinel SNPs interactively with liftOver (GRCh37 to GRCh38): https://genome.ucsc.edu/cgi-bin/hgLiftOver
#and then download the file for broad pheno sentinel vars.
awk '{print "chr"$1, $2, $2}' ${path_dir}/input/b38_chr11_gbmi_eur_sumstat_maf | \
    tail -n +2 > input/b38_chr11_gbmi.bed

#Successfully converted 1341801 records
#Conversion failed on 805 records.
##Explain failure messages:
#Deleted in new:
#    Sequence intersects no chains
#Partially deleted in new:
#    Sequence insufficiently intersects one chain
#Split in new:
#    Sequence insufficiently intersects multiple chains
#Duplicated in new:
#    Sequence sufficiently intersects multiple chains
#Boundary problem:
#    Missing start or end base in an exon

#Update the BP with GRCh37:
awk -F "\t|:|-" '{print $6}' ${path_dir}/input/b37_chr11_gbmi.bed | \
    grep -F -w -f - ${path_dir}/input/b38_chr11_gbmi_eur_sumstat_maf \
    > ${path_dir}/input/b38_with_b37_chr11_gbmi_sumstat

#in R:
dos2unix ${path_dir}/src/gbmi_b37.R
chmod o+x ${path_dir}/src/gbmi_b37.R
Rscript ${path_dir}/src/gbmi_b37.R

#Run Polyfun:
##Create parquet file:
python /home/n/nnp5/software/polyfun/munge_polyfun_sumstats.py \
  --sumstats ${path_dir}/input/b37_chr11_gbmi_eur_sumstat_maf \
  --out ${path_dir}/input/gbmi_eur_sumstats_munged.parquet \
  --min-info 0.85 \
  --min-maf 0.01

#Approach 1 (UK Biobank pre-computed prior causal probabilities); no information for HLA region !
##A.Extract pre-computed prior causal probabilities in UK Biobank - White British:
python /home/n/nnp5/software/polyfun/extract_snpvar.py \
    --sumstats ${path_dir}/input/gbmi_eur_sumstats_munged.parquet \
    --allow-missing \
    --out ${path_dir}/input/gbmi_eur_SNPs_PriCauPro


#Download pre-computed LD matrix from British ancestry UK Biobank individuals:
#rs11245964 is te sentinel variant from GBMI eur all-comer asthma.
#1119370 GRCh38
#1113278 GRCh37
start=$((1113278 - 500000))
end=$((1113278 + 500000))

#The same that I downloaded for moderate-to-severe (do not run it):
#wget https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/chr11_1_3000001.npz
#wget https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/chr11_1_3000001.gz

#Run the actual FIFO:
#median N-sample size: 1376100
#max N-sample size: 1376100
#min N-sample size: 26133
#USING MEDIAN VALUE FOR n-SAMPLE SIZE:
python /home/n/nnp5/software/polyfun/finemapper.py \
    --ld ${path_dir}/input/chr11_1_3000001 \
    --sumstats ${path_dir}/input/gbmi_eur_SNPs_PriCauPro \
    --n 1376100 \
    --chr 11 \
    --start $start \
    --end $end \
    --method susie \
    --max-num-causal 10 \
    --allow-missing \
    --out ${path_dir}/output/polyfun_susie_gbmi_eur_Nmedian

python /home/n/nnp5/software/polyfun/finemapper.py \
    --ld ${path_dir}/input/chr11_1_3000001 \
    --sumstats ${path_dir}/input/gbmi_eur_SNPs_PriCauPro \
    --n 26133 \
    --chr 11 \
    --start $start \
    --end $end \
    --method susie \
    --max-num-causal 10 \
    --allow-missing \
    --out ${path_dir}/output/polyfun_susie_gbmi_eur_Nmin

conda deactivate polyfun

#Extract credset 95% in R:
Rscript ${path_dir}/src/credset.R \
    ${path_dir}/output/polyfun_susie_gbmi_eur_Nmedian \
    ${path_dir}/output/polyfun_susie_gbmi_eur_credset_Nmedian \
    ${path_dir}/output/polyfun_susie_gbmi_eur_Nmedian_credsetbysusie

Rscript ${path_dir}/src/credset.R \
    ${path_dir}/output/polyfun_susie_gbmi_eur_Nmin \
    ${path_dir}/output/polyfun_susie_gbmi_eur_credset_Nmin \
    ${path_dir}/output/polyfun_susie_gbmi_eur_Nmin_credsetbysusie

#In R, make Venn diagram comparing the three fine-mapping credible sets:
dos2unix ${path_dir}/src/compare_credset.R
chmod o+x ${path_dir}/src/compare_credset.R
Rscript ${path_dir}/src/compare_credset.R


#Using the credible set as found by susie:
#In order for variables to be in a rho credible set, their PIP have to sum to greater than or equal to rho,
#and their "purity", defined by minimum absolute pairwise correlation, should be greater than or equal to r.
# Currently, rho is set to 0.95 and r is set to 0.5. You can lower either these numbers, or both,
# to be more lenient with credible sets, depending on your analysis context (current default is
# what we believe reasonable for genetic fine-mapping applications)

#QC:
#MODSEV: 441 variants not in fine-mapping
zcat input/modsev_SNPs_PriCauPro.miss.gz | awk -F "\t" '$2 == 11 && $3 >= 636478 && $3 <= 1636478 {print}' | wc -l
#2 variants not in fine-mapping for modsev and found in the credset of chronic sputum:
zcat input/modsev_SNPs_PriCauPro.miss.gz | awk -F "\t" '$2 == 11 && $3 >= 636478 && $3 <= 1636478 {print $3}' | \
    grep -F -f - input/chronic_sputum_credset_muc2
#GBMI: 389 variants not in fine-mapping
zcat input/gbmi_eur_SNPs_PriCauPro.miss.gz | awk -F "\t" '$2 >= 613278 && $2 <= 1613278 {print}' | wc -l
#2 variants not in fine-mapping for gbmi and found in the credset of chronic sputum (same x modsev):
zcat input/gbmi_eur_SNPs_PriCauPro.miss.gz | awk -F "\t" '$2 >= 613278 && $2 <= 1613278 {print $2}' | \
    grep -F -f - input/chronic_sputum_credset_muc2

#Look at sentinel variants in other study:
#moderate-to-severe sentinel variant: rs11603634
grep -w "rs11603634" /data/gen1/Phlegm/Results_files/GWAS/combined_results_INFO_MAC20
grep "1136478" ${path_dir}/input/chronic_sputum_credset_muc2
zgrep -w "rs11603634" /rfs/TobinGroup/kaf19/Public_data/gbmi_asthma/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz
grep -w "rs11603634" ${path_dir}/output/polyfun_susie_gbmi_eur_Nmedian
#chronic sputum sentinel variant: rs779167905, BP:1116931 (rs35606069 in moderate-to-severe)
zgrep "1116931" /rfs/TobinGroup/kaf19/Public_data/ukb_mod_sev_asthma/Shrine_30552067_moderate-severe_asthma.txt.gz
grep "1116931" ${path_dir}/output/polyfun_susie_modsev
zgrep -w "1116931" ${path_dir}/modsev_SNPs_PriCauPro.miss.gz
zgrep -w "rs779167905" /rfs/TobinGroup/kaf19/Public_data/gbmi_asthma/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz
grep -w "1116931" ${path_dir}/output/polyfun_susie_gbmi_eur_Nmedian
zgrep -w "1116931" /rfs/TobinGroup/kaf19/Public_data/gbmi_asthma/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz

#GBMI sentinel variant: rs11245964
zgrep "rs11245964" /rfs/TobinGroup/kaf19/Public_data/ukb_mod_sev_asthma/Shrine_30552067_moderate-severe_asthma.txt.gz
grep "rs11245964" ${path_dir}/output/polyfun_susie_modsev
grep -w "rs11245964" /data/gen1/Phlegm/Results_files/GWAS/combined_results_INFO_MAC20