#!/bin/bash

#Rationale: test different info thresholds and number of maximum causal variants for polyfun+susie analysis

#set up local variable, environment to use:
path_dir="/home/n/nnp5/PhD/PhD_project/muc5a_finemapping"
module load python/gcc/3.9.10
conda activate polyfun
module load gcc/9.3
cd ${path_dir}/

##MODERATE-TO-SEVERE:
#${path_dir}/input/modsev_sumstat_maf already created in polyfun_susie.sh
#Sentinel variant: rs11603634
#BP37: 1136478
#+/- 500000: 636478, 1636478
start=$((1136478-500000))
end=$((1136478+500000))

## INFO 085, --max-num-causal 1:
python /home/n/nnp5/software/polyfun/finemapper.py \
    --sumstats ${path_dir}/input/modsev_SNPs_PriCauPro \
    --n 30810 \
    --chr 11 \
    --start ${start} \
    --end  ${end} \
    --method susie \
    --max-num-causal 1 \
    --allow-missing \
    --out ${path_dir}/output/N1_polyfun_susie_modsev
#Extract credset 95% in R
Rscript ${path_dir}/src/credset.R ${path_dir}/output/N1_polyfun_susie_modsev ${path_dir}/output/N1_polyfun_susie_modsev_credset

##INFO 09: I need to create parquet file with --min-info 0.9 and extract snpvar again:
#Run Polyfun:
##Create parquet file:
python /home/n/nnp5/software/polyfun/munge_polyfun_sumstats.py \
  --sumstats ${path_dir}/input/modsev_sumstat_maf \
  --n 30810 \
  --out ${path_dir}/input/info09_modsev_sumstats_munged.parquet \
  --min-info 0.9 \
  --min-maf 0.01

#Approach 1 (UK Biobank pre-computed prior causal probabilities); no information for HLA region !
##A.Extract pre-computed prior causal probabilities in UK Biobank - White British:
python /home/n/nnp5/software/polyfun/extract_snpvar.py \
    --sumstats ${path_dir}/input/info09_modsev_sumstats_munged.parquet \
    --allow-missing \
    --out ${path_dir}/input/info09_modsev_SNPs_PriCauPro

#Run FIFO:
## INFO 09, --max-num-causal 10:
python /home/n/nnp5/software/polyfun/finemapper.py \
    --ld ${path_dir}/input/chr11_1_3000001 \
    --sumstats ${path_dir}/input/info09_modsev_SNPs_PriCauPro \
    --n 30810 \
    --chr 11 \
    --start ${start} \
    --end  ${end} \
    --method susie \
    --max-num-causal 10 \
    --allow-missing \
    --out ${path_dir}/output/info09_polyfun_susie_modsev
#Extract credset 95% in R
Rscript ${path_dir}/src/credset.R ${path_dir}/output/info09_polyfun_susie_modsev ${path_dir}/output/info09_polyfun_susie_modsev_credset

## INFO 09, --max-num-causal 1:
python /home/n/nnp5/software/polyfun/finemapper.py \
    --sumstats ${path_dir}/input/info09_modsev_SNPs_PriCauPro \
    --n 30810 \
    --chr 11 \
    --start 636478 \
    --end 1636478 \
    --method susie \
    --max-num-causal 1 \
    --allow-missing \
    --out ${path_dir}/output/N1_info09_polyfun_susie_modsev
#Extract credset 95% in R
Rscript ${path_dir}/src/credset.R ${path_dir}/output/N1_info09_polyfun_susie_modsev ${path_dir}/output/N1_info09_polyfun_susie_modsev_credset


#COMPARE DIFFERENT INFO AND MAXIMUM NUMBER OF CAUSAL VARIANTS FOR MODERATE-TO-SEVERE:
dos2unix ${path_dir}/src/compare_different_params.R
chmod o+x ${path_dir}/src/compare_different_params.R
Rscript ${path_dir}/src/compare_different_params.R \
    ${path_dir}/output/polyfun_susie_modsev_credset "info_085_N10" \
    ${path_dir}/output/N1_polyfun_susie_modsev_credset "info_085_N1" \
    ${path_dir}/output/info09_polyfun_susie_modsev_credset "info_090_N10" \
    ${path_dir}/output/N1_info09_polyfun_susie_modsev_credset "info_090_N1" \
    ${path_dir}/output/Venn_asthma_credset_modsev_params.png \
    ${path_dir}/output/Intersection_asthma_credset_modsev_params_BP_list


#####################################
#GBMI ALL-COMER ASTHMA IN EUROPEAN:
#rs11245964 is te sentinel variant from GBMI eur all-comer asthma.
#1119370 GRCh38
#1113278 GRCh37
start=$((1113278 - 500000))
end=$((1113278 + 500000))

#Run the actual FIFO:
#median N-sample size: 1376100
#max N-sample size: 1376100
#min N-sample size: 26133
#INFO 085, N-MEDIAN, --max-num-causal 1:
python /home/n/nnp5/software/polyfun/finemapper.py \
    --sumstats ${path_dir}/input/gbmi_eur_SNPs_PriCauPro \
    --n 1376100 \
    --chr 11 \
    --start $start \
    --end $end \
    --method susie \
    --max-num-causal 1 \
    --allow-missing \
    --out ${path_dir}/output/polyfun_susie_gbmi_eur_Nmedian_N1
#Extract credset 95% in R:
Rscript ${path_dir}/src/credset.R ${path_dir}/output/polyfun_susie_gbmi_eur_Nmedian_N1 ${path_dir}/output/polyfun_susie_gbmi_eur_credset_Nmedian_N1

#INFO 085, N-MIN, --max-num-causal 1:
python /home/n/nnp5/software/polyfun/finemapper.py \
    --sumstats ${path_dir}/input/gbmi_eur_SNPs_PriCauPro \
    --n 26133 \
    --chr 11 \
    --start $start \
    --end $end \
    --method susie \
    --max-num-causal 1 \
    --allow-missing \
    --out ${path_dir}/output/polyfun_susie_gbmi_eur_Nmin_N1
#Extract credset 95% in R:
Rscript ${path_dir}/src/credset.R ${path_dir}/output/polyfun_susie_gbmi_eur_Nmin_N1 ${path_dir}/output/polyfun_susie_gbmi_eur_credset_Nmin_N1

##INFO 09: I need to create parquet file with --min-info 0.9 and extract snpvar again:#Run Polyfun:
##Create parquet file:
python /home/n/nnp5/software/polyfun/munge_polyfun_sumstats.py \
  --sumstats ${path_dir}/input/b37_chr11_gbmi_eur_sumstat_maf \
  --out ${path_dir}/input/info09_gbmi_eur_sumstats_munged.parquet \
  --min-info 0.9\
  --min-maf 0.01

#Approach 1 (UK Biobank pre-computed prior causal probabilities); no information for HLA region !
##A.Extract pre-computed prior causal probabilities in UK Biobank - White British:
python /home/n/nnp5/software/polyfun/extract_snpvar.py \
    --sumstats ${path_dir}/input/info09_gbmi_eur_sumstats_munged.parquet \
    --allow-missing \
    --out ${path_dir}/input/info09_gbmi_eur_SNPs_PriCauPro

#INFO 09, N-MEDIAN, --max-num-causal 10:
python /home/n/nnp5/software/polyfun/finemapper.py \
    --ld ${path_dir}/input/chr11_1_3000001 \
    --sumstats ${path_dir}/input/info09_gbmi_eur_SNPs_PriCauPro \
    --n 1376100 \
    --chr 11 \
    --start $start \
    --end $end \
    --method susie \
    --max-num-causal 10 \
    --allow-missing \
    --out ${path_dir}/output/info09_polyfun_susie_gbmi_eur_Nmedian
#Extract credset 95% in R:
Rscript ${path_dir}/src/credset.R ${path_dir}/output/info09_polyfun_susie_gbmi_eur_Nmedian ${path_dir}/output/info09_polyfun_susie_gbmi_eur_credset_Nmedian

#INFO 09, N-MEDIAN, --max-num-causal 1:
python /home/n/nnp5/software/polyfun/finemapper.py \
    --sumstats ${path_dir}/input/info09_gbmi_eur_SNPs_PriCauPro \
    --n 1376100 \
    --chr 11 \
    --start $start \
    --end $end \
    --method susie \
    --max-num-causal 1 \
    --allow-missing \
    --out ${path_dir}/output/info09_polyfun_susie_gbmi_eur_Nmedian_N1
#Extract credset 95% in R:
Rscript ${path_dir}/src/credset.R ${path_dir}/output/info09_polyfun_susie_gbmi_eur_Nmedian_N1 ${path_dir}/output/info09_polyfun_susie_gbmi_eur_credset_Nmedian_N1

#INFO 09, N-MIN, --max-num-causal 10:
python /home/n/nnp5/software/polyfun/finemapper.py \
    --ld ${path_dir}/input/chr11_1_3000001 \
    --sumstats ${path_dir}/input/info09_gbmi_eur_SNPs_PriCauPro \
    --n 26133 \
    --chr 11 \
    --start $start \
    --end $end \
    --method susie \
    --max-num-causal 10 \
    --allow-missing \
    --out ${path_dir}/output/info09_polyfun_susie_gbmi_eur_Nmin
Rscript ${path_dir}/src/credset.R ${path_dir}/output/info09_polyfun_susie_gbmi_eur_Nmin ${path_dir}/output/info09_polyfun_susie_gbmi_eur_credset_Nmin

#INFO 09, N-MIN, --max-num-causal 1:
python /home/n/nnp5/software/polyfun/finemapper.py \
    --sumstats ${path_dir}/input/info09_gbmi_eur_SNPs_PriCauPro \
    --n 26133 \
    --chr 11 \
    --start $start \
    --end $end \
    --method susie \
    --max-num-causal 1 \
    --allow-missing \
    --out ${path_dir}/output/info09_polyfun_susie_gbmi_eur_Nmin_N1
Rscript ${path_dir}/src/credset.R ${path_dir}/output/info09_polyfun_susie_gbmi_eur_Nmin_N1 ${path_dir}/output/info09_polyfun_susie_gbmi_eur_credset_Nmin_N1

#COMPARE DIFFERENT INFO AND MAXIMUM NUMBER OF CAUSAL VARIANTS FOR GBMI:
#N_median:
Rscript ${path_dir}/src/compare_different_params.R \
    ${path_dir}/output/polyfun_susie_gbmi_eur_credset_Nmedian "info_085_Nmedian_N10" \
    ${path_dir}/output/polyfun_susie_gbmi_eur_credset_Nmedian_N1 "info_085_Nmedian_N1" \
    ${path_dir}/output/info09_polyfun_susie_gbmi_eur_credset_Nmedian "info_090_Nmedian_N10" \
    ${path_dir}/output/info09_polyfun_susie_gbmi_eur_credset_Nmedian_N1 "info_090_Nmedian_N1" \
    ${path_dir}/output/Venn_asthma_credset_gbmi_Nmedian_params.png \
    ${path_dir}/output/Intersection_asthma_credset_gbmi_Nmedian_params_BP_list

#N-min:
Rscript ${path_dir}/src/compare_different_params.R \
    ${path_dir}/output/polyfun_susie_gbmi_eur_credset_Nmin "info_085_Nmin_N10" \
    ${path_dir}/output/polyfun_susie_gbmi_eur_credset_Nmin_N1 "info_085_Nmin_N1" \
    ${path_dir}/output/info09_polyfun_susie_gbmi_eur_credset_Nmin "info_090_Nmin_N10" \
    ${path_dir}/output/info09_polyfun_susie_gbmi_eur_credset_Nmin_N1 "info_090_Nmin_N1" \
    ${path_dir}/output/Venn_asthma_credset_gbmi_Nmin_params.png \
    ${path_dir}/output/Intersection_asthma_credset_gbmi_Nmin_params_BP_list

#SUSiE did not fins any credible set in the moderate-to-severe with coverage set at 0.95;
# what about if I set the coverage at .90?
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
    --out ${path_dir}/output/polyfun_susie_modsev_coverage90

##but if rs35606069 is excluded because it is not found in the pre-calculated priors
##could you repeat with just Susie (not with PolyFun)?
#Run moderate-to-severe without prior so that the sumstats has the chronic sputum sentinel variant:
#sumstat: gbmi_eur_sumstats_munged.parquet
#--non-funct parameter for non-functionally informed fine-mapping
python /home/n/nnp5/software/polyfun/finemapper.py \
    --ld ${path_dir}/input/chr11_1_3000001 \
    --sumstats ${path_dir}/input/gbmi_eur_sumstats_munged.parquet \
    --n 30810 \
    --chr 11 \
    --start 636478 \
    --end 1636478 \
    --method susie \
    --non-funct \
    --max-num-causal 10 \
    --allow-missing \
    --out ${path_dir}/output/polyfun_susie_modsev_noprior
Rscript ${path_dir}/src/credset.R \
    ${path_dir}/output/polyfun_susie_modsev_noprior \
    ${path_dir}/output/polyfun_susie_modsev_noprior_credset \
    ${path_dir}/output/polyfun_susie_modsev_noprior_credsetbysusie
