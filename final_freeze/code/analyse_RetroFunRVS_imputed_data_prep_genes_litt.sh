module load StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16 r/4.2
path_genes_sz_bp=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/genes_litt/2024-06-25_selected_genes_ID_BP_SZ.txt
path_genes_mullins=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/2024-06-25_selected_genes_ID_Mullins2021.txt
path_rvs=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS
path_exons_litt=${path_rvs}/genes_litt
path_anno=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/annotations
#Should we only keep the variants with the consequences we are interested in ?
path_effects_to_keep=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/2024_09_30_annotations_effects_to_keep.txt

mkdir ${path_exons_litt}; mkdir ${path_exons_litt}/exons_only; mkdir ${path_exons_litt}/exons_only_strict

#Keep every position in our sequences found in these genes.
sed "s/\r//g" ${path_genes_sz_bp} | sed "s/^/|/g" | sed "s/$/|/g" > ${path_genes_sz_bp}_tmp
for chr in {1..22} 
do
  bcftools query -f 'ID=%ID, annotation=%ANN\n' ${path_anno}/seq_FINAL_chr${chr}_ccds_annoted_dbNSFP4_snpeff_with_mask.vcf.gz \
        | grep -f ${path_genes_sz_bp}_tmp >> ${path_exons_litt}/bp_sz_genes_effects_to_keep_seq_FINAL.txt
done
grep -f ${path_effects_to_keep} ${path_exons_litt}/bp_sz_genes_effects_to_keep_seq_FINAL.txt >> ${path_exons_litt}/exons_only/bp_sz_genes_effects_to_keep_seq_FINAL.txt
grep -f <(sed "/synonymous_variant/d;/^$/d" ${path_effects_to_keep}) ${path_exons_litt}/bp_sz_genes_effects_to_keep_seq_FINAL.txt > ${path_exons_litt}/exons_only_strict/bp_sz_genes_effects_to_keep_seq_FINAL.txt

#Create the .ped files and the .bed. -- All variants
for chr in {1..22} 
do 
  grep "chr${chr}:" ${path_exons_litt}/bp_sz_genes_effects_to_keep_seq_FINAL.txt > ${path_exons_litt}/var_to_extract_chr${chr}.txt
  while IFS= read -r gene
  do
    grep ${gene} ${path_exons_litt}/var_to_extract_chr${chr}.txt | cut -d " " -f 1 | sed "s/ID=//g" | sed "s/,//g"  > ${path_exons_litt}/var_to_extract.txt_tmp
    if [ -s "${path_exons_litt}/var_to_extract.txt_tmp" ]; then
      gene_out=$(echo ${gene} | sed "s/|//g")
      plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr} --extract ${path_exons_litt}/var_to_extract.txt_tmp --keep-allele-order --recode --out ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_gene_${gene_out}
      plink --file ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_gene_${gene_out} --keep-allele-order --nonfounders --freq --out ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_gene_${gene_out}
    fi
    rm ${path_exons_litt}/var_to_extract.txt_tmp
  done < ${path_genes_sz_bp}_tmp
done

#Create the .ped files and the .bed. -- Exons only
for chr in {1..22} 
do 
  grep "chr${chr}:" ${path_exons_litt}/exons_only/bp_sz_genes_effects_to_keep_seq_FINAL.txt > ${path_exons_litt}/exons_only/var_to_extract_chr${chr}.txt
  while IFS= read -r gene
  do
    grep ${gene} ${path_exons_litt}/exons_only/var_to_extract_chr${chr}.txt | cut -d " " -f 1 | sed "s/ID=//g" | sed "s/,//g"  > ${path_exons_litt}/exons_only/var_to_extract.txt_tmp
    if [ -s "${path_exons_litt}/exons_only/var_to_extract.txt_tmp" ]; then
      gene_out=$(echo ${gene} | sed "s/|//g")
      plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr} --extract ${path_exons_litt}/exons_only/var_to_extract.txt_tmp --keep-allele-order --recode --out ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_gene_${gene_out}
      plink --file ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_gene_${gene_out} --keep-allele-order --nonfounders --freq --out ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_gene_${gene_out}
    fi
    rm ${path_exons_litt}/exons_only/var_to_extract.txt_tmp
  done < ${path_genes_sz_bp}_tmp
done
rm ${path_genes_sz_bp}_tmp

#It is possible that some .ped/.map are empty. In fact, the annotated data also present non-rare variants
#So if the variants found in the genes are all non-rare, then the file is empty because we only keep rare variants.




#Prepare a dataset with every genes.
sed "s/\r//g" ${path_genes_sz_bp} | sed "s/^/|/g" | sed "s/$/|/g" > ${path_genes_sz_bp}_tmp

#Create the .ped file and the .bed. -- All variants
for chr in {1..22} 
do 
  grep "chr${chr}:" ${path_exons_litt}/bp_sz_genes_effects_to_keep_seq_FINAL.txt > ${path_exons_litt}/all_genes_var_to_extract_chr${chr}.txt
  cut -d " " -f 1 ${path_exons_litt}/all_genes_var_to_extract_chr${chr}.txt  | sed "s/ID=//g" | sed "s/,//g"  > ${path_exons_litt}/var_to_extract.txt_tmp_chr${chr}
  if [ -s "${path_exons_litt}/var_to_extract.txt_tmp_chr${chr}" ]; then
    plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr} --extract ${path_exons_litt}/var_to_extract.txt_tmp_chr${chr} --keep-allele-order --recode --out ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes
    plink --file ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes --keep-allele-order --nonfounders --freq --out ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes
  fi
  rm ${path_exons_litt}/var_to_extract.txt_tmp_chr${chr}
  if [[ -f "${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.ped" ]]; then
    echo ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.ped ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.map >> ${path_exons_litt}/all_genes_merge_list.txt
  fi
done
plink --merge-list ${path_exons_litt}/all_genes_merge_list.txt --keep-allele-order --recode --out ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene
rm ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.bim ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.bed ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.fam
rm ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_chr*_all_genes.ped
plink --file ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene --keep-allele-order --nonfounders --freq --out ${path_exons_litt}/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene
cat ${path_exons_litt}/all_genes_var_to_extract_chr*.txt > ${path_exons_litt}/all_genes_var_to_extract.txt; rm ${path_exons_litt}/all_genes_var_to_extract_chr*.txt

#Create the .ped file and the .bed. -- Exons only
for chr in {1..22} 
do 
 {
  grep "chr${chr}:" ${path_exons_litt}/exons_only/bp_sz_genes_effects_to_keep_seq_FINAL.txt > ${path_exons_litt}/exons_only/all_genes_var_to_extract_chr${chr}.txt
  cut -d " " -f 1 ${path_exons_litt}/exons_only/all_genes_var_to_extract_chr${chr}.txt | sed "s/ID=//g" | sed "s/,//g"  > ${path_exons_litt}/exons_only/var_to_extract.txt_tmp_chr${chr}
  if [ -s "${path_exons_litt}/exons_only/var_to_extract.txt_tmp_chr${chr}" ]; then
    plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr} --extract ${path_exons_litt}/exons_only/var_to_extract.txt_tmp_chr${chr} --keep-allele-order --recode --out ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes
    plink --file ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes --keep-allele-order --nonfounders --freq --out ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes
  fi
  rm ${path_exons_litt}/exons_only/var_to_extract.txt_tmp_chr${chr}
  if [[ -f "${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.ped" ]]; then
    echo ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.ped ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.map >> ${path_exons_litt}/exons_only/all_genes_merge_list.txt
  fi
 } &
done
wait
plink --merge-list ${path_exons_litt}/exons_only/all_genes_merge_list.txt --keep-allele-order --recode --out ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene
rm ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.bim ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.bed ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.fam
rm ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr*_all_genes.ped
plink --file ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene --keep-allele-order --nonfounders --freq --out ${path_exons_litt}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene
cat ${path_exons_litt}/exons_only/all_genes_var_to_extract_chr*.txt > ${path_exons_litt}/exons_only/all_genes_var_to_extract.txt; rm ${path_exons_litt}/exons_only/all_genes_var_to_extract_chr*.txt


#Create the .ped file and the .bed. -- Exons only without synonymous_variant
for chr in {1..22} 
do 
 {
  grep "chr${chr}:" ${path_exons_litt}/exons_only_strict/bp_sz_genes_effects_to_keep_seq_FINAL.txt > ${path_exons_litt}/exons_only_strict/all_genes_var_to_extract_chr${chr}.txt
  cut -d " " -f 1 ${path_exons_litt}/exons_only_strict/all_genes_var_to_extract_chr${chr}.txt | sed "s/ID=//g" | sed "s/,//g"  > ${path_exons_litt}/exons_only_strict/var_to_extract.txt_tmp_chr${chr}
  if [ -s "${path_exons_litt}/exons_only_strict/var_to_extract.txt_tmp_chr${chr}" ]; then
    plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr} --extract ${path_exons_litt}/exons_only_strict/var_to_extract.txt_tmp_chr${chr} --keep-allele-order --recode --out ${path_exons_litt}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes
    plink --file ${path_exons_litt}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes --keep-allele-order --nonfounders --freq --out ${path_exons_litt}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes
  fi
  rm ${path_exons_litt}/exons_only_strict/var_to_extract.txt_tmp_chr${chr}
  if [[ -f "${path_exons_litt}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.ped" ]]; then
    echo ${path_exons_litt}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.ped ${path_exons_litt}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.map >> ${path_exons_litt}/exons_only_strict/all_genes_merge_list.txt
  fi
 } &
done
wait
plink --merge-list ${path_exons_litt}/exons_only_strict/all_genes_merge_list.txt --keep-allele-order --recode --out ${path_exons_litt}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene
rm ${path_exons_litt}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.bim ${path_exons_litt}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.bed ${path_exons_litt}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.fam
rm ${path_exons_litt}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_chr*_all_genes.ped
plink --file ${path_exons_litt}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene --keep-allele-order --nonfounders --freq --out ${path_exons_litt}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene
cat ${path_exons_litt}/exons_only_strict/all_genes_var_to_extract_chr*.txt > ${path_exons_litt}/exons_only_strict/all_genes_var_to_extract.txt; rm ${path_exons_litt}/exons_only_strict/all_genes_var_to_extract_chr*.txt
rm ${path_genes_sz_bp}_tmp



