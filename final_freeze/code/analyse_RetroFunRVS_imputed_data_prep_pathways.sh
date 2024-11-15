module load StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16 r/4.2
path_rvs=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS
path_pathways=${path_rvs}/pathways
mkdir ${path_pathways}; mkdir ${path_pathways}/exons_only; mkdir ${path_pathways}/exons_only_strict
path_anno=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/annotations

#Should we only keep the variants with the consequences we are interested in ?
path_effects_to_keep=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/2024_09_30_annotations_effects_to_keep.txt

#To get the correspondance between Ensembl genes codes and chr:pos_begin:pos_end, go to https://www.ensembl.org/biomart/martview/eef451126d44dd44b71b8705c916e5dd
#Download the file and bring it to Beluga
#At first, the file is named "mart_export.txt". I renamed it biomart_syngo_genes_pos.txt
path_ensembl=${path_pathways}/biomart_syngo_genes_pos.txt

#The syngo ontology file at first is in .xlsx format. We don't want that. Use a converter to go from .xlsx to txt .format.
#Create a whole list with every genes hgnc_symbol
cut -f 7 ${path_pathways}/syngo_ontologies.txt | sed "s/, /\n/g" | sort | uniq > ${path_pathways}/syngo_ontologies_all_genes_symbols.txt

#Keep every position in our sequences found in these genes.
sed "s/\r//g" ${path_pathways}/syngo_ontologies_all_genes_symbols.txt | sed "s/^/|/g" | sed "s/$/|/g" > ${path_pathways}/syngo_ontologies_all_genes_symbols.txt_tmp
for chr in {1..22} 
do
  bcftools query -f 'ID=%ID, annotation=%ANN\n' ${path_anno}/seq_FINAL_chr${chr}_ccds_annoted_dbNSFP4_snpeff_with_mask.vcf.gz \
        | grep -f ${path_pathways}/syngo_ontologies_all_genes_symbols.txt_tmp >> ${path_pathways}/pathways_genes_effects_to_keep_seq_FINAL.txt
done
grep -f ${path_effects_to_keep} ${path_pathways}/pathways_genes_effects_to_keep_seq_FINAL.txt > ${path_pathways}/exons_only/pathways_genes_effects_to_keep_seq_FINAL.txt
grep -f <(sed "/synonymous_variant/d;/^$/d" ${path_effects_to_keep}) ${path_pathways}/pathways_genes_effects_to_keep_seq_FINAL.txt > ${path_pathways}/exons_only_strict/pathways_genes_effects_to_keep_seq_FINAL.txt

#How many rare variants for every line of the ontology.
cut -f2 ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr*.map > ${path_pathways}/rare_ID.txt
cut -f7 "${path_pathways}/syngo_ontologies.txt" | while IFS= read -r line
do
  echo ${line} | sed "s/, /\n/g" | sort | uniq | sed "s/^/|/g" | sed "s/$/|/g" > ${path_pathways}/genes_to_keep.tmp
  grep -f ${path_pathways}/genes_to_keep.tmp ${path_pathways}/pathways_genes_effects_to_keep_seq_FINAL.txt | grep -w -f ${path_pathways}/rare_ID.txt | wc -l >> ${path_pathways}/syngo_ontologies_rare_n.txt
  grep -f ${path_pathways}/genes_to_keep.tmp ${path_pathways}/exons_only/pathways_genes_effects_to_keep_seq_FINAL.txt | grep -w -f ${path_pathways}/rare_ID.txt | wc -l >> ${path_pathways}/exons_only/syngo_ontologies_rare_n.txt
  grep -f ${path_pathways}/genes_to_keep.tmp ${path_pathways}/exons_only_strict/pathways_genes_effects_to_keep_seq_FINAL.txt | grep -w -f ${path_pathways}/rare_ID.txt | wc -l >> ${path_pathways}/exons_only_strict/syngo_ontologies_rare_n.txt
done

#Create the .ped file and the .bed. -- All variants
for chr in {1..22} 
do 
 {
  grep "chr${chr}:" ${path_pathways}/pathways_genes_effects_to_keep_seq_FINAL.txt > ${path_pathways}/all_genes_var_to_extract_chr${chr}.txt
  cut -d " " -f 1 ${path_pathways}/all_genes_var_to_extract_chr${chr}.txt  | sed "s/ID=//g" | sed "s/,//g"  > ${path_pathways}/var_to_extract.txt_tmp_chr${chr}
  if [ -s "${path_pathways}/var_to_extract.txt_tmp_chr${chr}" ]; then
    plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr} --extract ${path_pathways}/var_to_extract.txt_tmp_chr${chr} --keep-allele-order --recode --out ${path_pathways}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes
    plink --file ${path_pathways}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes --keep-allele-order --nonfounders --freq --out ${path_pathways}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes
  fi
  rm ${path_pathways}/var_to_extract.txt_tmp_chr${chr}
  if [[ -f "${path_pathways}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.ped" ]]; then
    echo ${path_pathways}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.ped ${path_pathways}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.map >> ${path_pathways}/all_genes_merge_list.txt
  fi
 } &
done
wait
plink --merge-list ${path_pathways}/all_genes_merge_list.txt --keep-allele-order --recode --out ${path_pathways}/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene
rm ${path_pathways}/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.bim ${path_pathways}/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.bed ${path_pathways}/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.fam
rm ${path_pathways}/impute5_gigi2_combined_seq_RV_FINAL_chr*_all_genes.ped
plink --file ${path_pathways}/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene --keep-allele-order --nonfounders --freq --out ${path_pathways}/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene
cat ${path_pathways}/all_genes_var_to_extract_chr*.txt > ${path_pathways}/all_genes_var_to_extract.txt; rm ${path_pathways}/all_genes_var_to_extract_chr*.txt

#Create the .ped file and the .bed. -- Exons only
for chr in {1..22} 
do 
 {
  grep "chr${chr}:" ${path_pathways}/exons_only/pathways_genes_effects_to_keep_seq_FINAL.txt > ${path_pathways}/exons_only/all_genes_var_to_extract_chr${chr}.txt
  cut -d " " -f 1 ${path_pathways}/exons_only/all_genes_var_to_extract_chr${chr}.txt | sed "s/ID=//g" | sed "s/,//g"  > ${path_pathways}/exons_only/var_to_extract.txt_tmp_chr${chr}
  if [ -s "${path_pathways}/exons_only/var_to_extract.txt_tmp_chr${chr}" ]; then
    plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr} --extract ${path_pathways}/exons_only/var_to_extract.txt_tmp_chr${chr} --keep-allele-order --recode --out ${path_pathways}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes
    plink --file ${path_pathways}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes --keep-allele-order --nonfounders --freq --out ${path_pathways}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes
  fi
  rm ${path_pathways}/exons_only/var_to_extract.txt_tmp_chr${chr}
  if [[ -f "${path_pathways}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.ped" ]]; then
    echo ${path_pathways}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.ped ${path_pathways}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.map >> ${path_pathways}/exons_only/all_genes_merge_list.txt
  fi
 } &
done
wait
plink --merge-list ${path_pathways}/exons_only/all_genes_merge_list.txt --keep-allele-order --recode --out ${path_pathways}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene
rm ${path_pathways}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.bim ${path_pathways}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.bed ${path_pathways}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.fam
rm ${path_pathways}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_chr*_all_genes.ped
plink --file ${path_pathways}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene --keep-allele-order --nonfounders --freq --out ${path_pathways}/exons_only/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene
cat ${path_pathways}/exons_only/all_genes_var_to_extract_chr*.txt > ${path_pathways}/exons_only/all_genes_var_to_extract.txt; rm ${path_pathways}/exons_only/all_genes_var_to_extract_chr*.txt

#Create the .ped file and the .bed. -- Exons only without synonymous_variant
for chr in {1..22} 
do 
 {
  grep "chr${chr}:" ${path_pathways}/exons_only_strict/pathways_genes_effects_to_keep_seq_FINAL.txt > ${path_pathways}/exons_only_strict/all_genes_var_to_extract_chr${chr}.txt
  cut -d " " -f 1 ${path_pathways}/exons_only_strict/all_genes_var_to_extract_chr${chr}.txt | sed "s/ID=//g" | sed "s/,//g"  > ${path_pathways}/exons_only_strict/var_to_extract.txt_tmp_chr${chr}
  if [ -s "${path_pathways}/exons_only_strict/var_to_extract.txt_tmp_chr${chr}" ]; then
    plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr} --extract ${path_pathways}/exons_only_strict/var_to_extract.txt_tmp_chr${chr} --keep-allele-order --recode --out ${path_pathways}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes
    plink --file ${path_pathways}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes --keep-allele-order --nonfounders --freq --out ${path_pathways}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes
  fi
  rm ${path_pathways}/exons_only_strict/var_to_extract.txt_tmp_chr${chr}
  if [[ -f "${path_pathways}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.ped" ]]; then
    echo ${path_pathways}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.ped ${path_pathways}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_all_genes.map >> ${path_pathways}/exons_only_strict/all_genes_merge_list.txt
  fi
 } &
done
wait
plink --merge-list ${path_pathways}/exons_only_strict/all_genes_merge_list.txt --keep-allele-order --recode --out ${path_pathways}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene
rm ${path_pathways}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.bim ${path_pathways}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.bed ${path_pathways}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.fam
rm ${path_pathways}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_chr*_all_genes.ped
plink --file ${path_pathways}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene --keep-allele-order --nonfounders --freq --out ${path_pathways}/exons_only_strict/impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene
cat ${path_pathways}/exons_only_strict/all_genes_var_to_extract_chr*.txt > ${path_pathways}/exons_only_strict/all_genes_var_to_extract.txt; rm ${path_pathways}/exons_only_strict/all_genes_var_to_extract_chr*.txt
rm ${path_pathways}/syngo_ontologies_all_genes_symbols.txt_tmp
