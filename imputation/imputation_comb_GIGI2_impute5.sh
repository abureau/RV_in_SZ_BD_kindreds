#merge pop (Impute5) and fam (GIGI2) imputation
path_pop=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/imputation_pop
path_fam=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/imputation_fam
path_comb=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/imputation_comb
soft_comb=/lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples/ped_pop
mask=with_mask
mkdir ${path_comb}; mkdir ${path_comb}/logs

#Parameters file to use Mohammed's software.
echo "MAF 0.01" > ${path_comb}/param.txt
echo "Method PED_POP" >> ${path_comb}/param.txt

#Only keep the imputed subjects in the GIGI2 files. #REVOIR LE MATCHING POUR LES SUJETS SANS FAMID POUR IMPUTE5.
chr=22
bcftools query -l ${path_pop}/impute5_imputed_chr${chr}.vcf.gz > ${path_comb}/impute5_samples.txt
bcftools query -l ${path_fam}/chr${chr}/geno_impute_chr${chr}_ATCG.vcf.gz > ${path_comb}/gigi2_samples.txt
comm -1 -2 <(cut -d "_" -f 2 ${path_comb}/impute5_samples.txt | sort) <(cut -d "_" -f 2 ${path_comb}/gigi2_samples.txt | sort) > ${path_comb}/common_ID_samples.txt
egrep -f <(sed "s/^/_/g" ${path_comb}/common_ID_samples.txt | sed "s/$/\$/g") ${path_comb}/impute5_samples.txt > ${path_comb}/impute5_common_samples.txt
egrep -v -f <(sed "s/^/_/g" ${path_comb}/common_ID_samples.txt | sed "s/$/\$/g") ${path_comb}/impute5_samples.txt > ${path_comb}/impute5_only_samples.txt
egrep -f <(sed "s/^/_/g" ${path_comb}/common_ID_samples.txt | sed "s/$/\$/g") ${path_comb}/gigi2_samples.txt > ${path_comb}/gigi2_common_samples.txt
egrep -v -f <(sed "s/^/_/g" ${path_comb}/common_ID_samples.txt | sed "s/$/\$/g") ${path_comb}/gigi2_samples.txt > ${path_comb}/gigi2_only_samples.txt

#Sort the common lists to make sure that the subjects in the following file are found in the same order.
#We sort based on gigi2's list, where the family ID are complete. Then, we sort impute5's list according to gigi2's list.
sort -o ${path_comb}/gigi2_common_samples.txt ${path_comb}/gigi2_common_samples.txt
unsorted=${path_comb}/impute5_common_samples_unsorted.txt
sorted=${path_comb}/impute5_common_samples.txt
mv ${sorted} ${unsorted}
cut -d "_" -f 2 ${path_comb}/gigi2_common_samples.txt | sed "s/^/_/g" | sed "s/$/\$/g" | while read line; do
  egrep "$line" "$unsorted" >> "$sorted"
done
#Verification
diff <(cut -d "_" -f 2 ${path_comb}/gigi2_common_samples.txt) <(cut -d "_" -f 2 ${path_comb}/impute5_common_samples.txt)

#Convert dosage to genotypes, we need to create a fam file in the sorted order
sort -k 2 <(awk 'NR==FNR {samples[$1]; next} $2 in samples' ${path_comb}/gigi2_common_samples.txt ${path_fam}/chr${chr}/geno_impute_chr${chr}.tfam) > ${path_comb}/impute5_gigi2_dose.fam

for chr in {1..22}
do
  #Using -S on gigi2_common_samples.txt, a sorted list, on GIGI2 and Impute5 vcf assure that the samples are in the same order in both files.
  bcftools view -S ${path_comb}/gigi2_common_samples.txt -O z -o ${path_comb}/gigi2_imputed_chr${chr}.vcf.gz ${path_fam}/chr${chr}/geno_impute_chr${chr}_ATCG.vcf.gz

  #For Impute5, we wish to only keep the variants imputed in GIGI2 as well.
  bcftools query -f '%ID\n' ${path_fam}/chr${chr}/geno_impute_chr${chr}_ATCG.vcf.gz > ${path_comb}/gigi2_variants_chr${chr}.txt
  bcftools view -S ${path_comb}/impute5_common_samples.txt -i "ID=@${path_comb}/gigi2_variants_chr${chr}.txt" -O z -o ${path_comb}/impute5_imputed_chr${chr}.vcf.gz ${path_pop}/impute5_imputed_chr${chr}.vcf.gz

  #Prepare the prob and INFO files for Impute5. As we are using ${path_comb}/impute5_imputed_chr${chr}.vcf.gz, the prob and INFO are sorted by ID.
  bcftools query -f '[%GP ]\n' ${path_comb}/impute5_imputed_chr${chr}.vcf.gz | sed 's/,/ /g' > ${path_comb}/impute5_prob_chr${chr}.txt
  paste -d " " <(bcftools query -f '%ID\n' ${path_comb}/impute5_imputed_chr${chr}.vcf.gz) <(bcftools query -f '%INFO\n' ${path_comb}/impute5_imputed_chr${chr}.vcf.gz | grep -oP 'INFO=\K[^;]+') > ${path_comb}/impute5_INFO_chr${chr}.txt

  #GIGI2 prob files are not sorted. To get them in the right order, we match the GIGI2 sample ID with the gigi2_common_samples.txt list and
  #we print the line number of every match. Then, we sort by sample ID and we get the right prob using the line number.
  grep -nFxf "${path_comb}/gigi2_common_samples.txt" "${path_comb}/gigi2_samples.txt" | sed "s/:/ /g" | sort -k 2 | cut -d " " -f1 > ${path_comb}/matching_IDs.txt
  awk '
  NR==FNR {
      indices[NR] = $1; num_indices = NR; next;
  }
  {
      output = "";
      for (i = 1; i <= num_indices; i++) {
          id = indices[i]; output = output (output ? " " : "") $(3*id-1) " " $(3*id) " " $(3*id+1);
      }
      print output;
  }' ${path_comb}/matching_IDs.txt ${path_fam}/chr${chr}/GIGI2_imputed_chr${chr}_by_fam.prob > ${path_comb}/gigi2_prob_chr${chr}.txt
  rm ${path_comb}/matching_IDs.txt

  #We running the software for the first time, I understood that my prob files needs to be transposed. 
  #I wrote an executable R scripts to do so.
  Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples/call/transpose_3prob.R -prob_type var_subject -prob ${path_comb}/gigi2_prob_chr${chr}.txt -out ${path_comb}/gigi2_prob_transposed_chr${chr}.txt
  Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples/call/transpose_3prob.R -prob_type var_subject -prob ${path_comb}/impute5_prob_chr${chr}.txt -out ${path_comb}/impute5_prob_transposed_chr${chr}.txt
  rm ${path_comb}/gigi2_prob_chr${chr}.txt ${path_comb}/impute5_prob_chr${chr}.txt

  #Use Mohammed's software. Tranpose the output.
  ${soft_comb}/PED_POP ${path_comb}/gigi2_prob_transposed_chr${chr}.txt ${path_comb}/impute5_prob_transposed_chr${chr}.txt ${path_comb}/impute5_INFO_chr${chr}.txt ${path_comb}/param.txt ${path_comb}/impute5_gigi2_prob_transposed_chr${chr}.txt ${path_comb}/impute5_gigi2_dose_transposed_chr${chr}.txt
  Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples/call/transpose_dosage.R -dosage ${path_comb}/impute5_gigi2_dose_transposed_chr${chr}.txt -out ${path_comb}/impute5_gigi2_dose_chr${chr}_step1.txt
  gzip ${path_comb}/impute5_gigi2_dose_transposed_chr${chr}.txt; gzip ${path_comb}/impute5_gigi2_prob_transposed_chr${chr}.txt
  gzip ${path_comb}/gigi2_prob_transposed_chr${chr}.txt; gzip ${path_comb}/impute5_prob_transposed_chr${chr}.txt

  #To convert dosage to genotypes, we need to add SNPID, A1, A2 for each. We put the allele in order %ALT %REF because --import-dosage interpret it that way.
  paste -d " " <(bcftools query -f '%ID %ALT %REF\n' ${path_comb}/impute5_imputed_chr${chr}.vcf.gz) ${path_comb}/impute5_gigi2_dose_chr${chr}_step1.txt > ${path_comb}/impute5_gigi2_dose_chr${chr}.txt
  rm ${path_comb}/impute5_gigi2_dose_chr${chr}_step1.txt

  #Convert to genotypes. A cutoff of 0.1 on the dosage is used to call the genotypes as Impute2 uses a cutoff of 0.9 on the probabilities when calling the genotypes.
  plink2 --import-dosage ${path_comb}/impute5_gigi2_dose_chr${chr}.txt 'noheader' --fam ${path_comb}/impute5_gigi2_dose_chr${chr}.fam --hard-call-threshold 0.1 --make-bed --out ${path_comb}/impute5_gigi2_combined_chr${chr}
  mv ${path_comb}/impute5_gigi2_combined_chr${chr}.bim ${path_comb}/impute5_gigi2_combined_chr${chr}.bim~
  awk -v OFS="\t" '{print $1, $3, $4, $2, $5, $6}' <(paste <(cut -f 2 ${path_comb}/impute5_gigi2_combined_chr${chr}.bim~ | cut -d ":" -f 1,2 | sed "s/:/\t/g") <(cut -f 2,3,5,6 ${path_comb}/impute5_gigi2_combined_chr${chr}.bim~)) > ${path_comb}/impute5_gigi2_combined_chr${chr}.bim
  rm ${path_comb}/impute5_gigi2_combined_chr${chr}.bim~
done

#Compare the output to the Impute5 and GIGI2 files to verify if the output makes sense...
mkdir ${path_comb}/comp/
for chr in {1..22}
do
  plink2 --vcf ${path_comb}/gigi2_imputed_chr${chr}.vcf.gz --make-bed --out ${path_comb}/comp/gigi2_chr${chr}
  plink2 --vcf ${path_comb}/impute5_imputed_chr${chr}.vcf.gz --make-bed --out ${path_comb}/comp/impute5_chr${chr}
  rm ${path_comb}/comp/impute5_chr${chr}.fam; cp ${path_comb}/comp/gigi2_chr${chr}.fam ${path_comb}/comp/impute5_chr${chr}.fam #Because samples aren't all correcly labeled in Impute5 data. 
  rm ${path_comb}/comp/gigi2_chr${chr}.bim; cp ${path_comb}/comp/impute5_chr${chr}.bim ${path_comb}/comp/gigi2_chr${chr}.bim #Because INDELS alleles are not labelled the same between Impute5 and Gigi2 (we used Impute5 as the header reference when combining). 
  plink --bfile ${path_comb}/impute5_gigi2_combined_chr${chr} --bmerge ${path_comb}/comp/gigi2_chr${chr} --merge-mode 7 --keep-allele-order --out ${path_comb}/comp/gigi2_chr${chr}_comp
  plink --bfile ${path_comb}/impute5_gigi2_combined_chr${chr} --bmerge ${path_comb}/comp/impute5_chr${chr} --merge-mode 7 --keep-allele-order --out ${path_comb}/comp/impute5_chr${chr}_comp
  rm ${path_comb}/comp/gigi2_chr${chr}.*; rm ${path_comb}/comp/impute5_chr${chr}.*; gzip ${path_comb}/comp/impute5_chr${chr}_comp.diff; gzip ${path_comb}/comp/gigi2_chr${chr}_comp.diff
  grep "rate" ${path_comb}/comp/impute5_chr${chr}_comp.log | cut -d " " -f 8 | sed "s/.$//" >> ${path_comb}/comp/impute5_comp.txt
  grep "rate" ${path_comb}/comp/gigi2_chr${chr}_comp.log | cut -d " " -f 8 | sed "s/.$//" >> ${path_comb}/comp/gigi2_comp.txt

  #Rare variants (based on GIGI2 to compare the same set of variants)
  plink2 --vcf ${path_comb}/gigi2_imputed_chr${chr}.vcf.gz --max-maf 0.01 --make-bed --out ${path_comb}/comp/gigi2_chr${chr}
  cut -f2 ${path_comb}/comp/gigi2_chr${chr}.bim > ${path_comb}/comp/gigi2_chr${chr}_rare.snplist
  plink2 --vcf ${path_comb}/impute5_imputed_chr${chr}.vcf.gz --extract ${path_comb}/comp/gigi2_chr${chr}_rare.snplist --make-bed --out ${path_comb}/comp/impute5_chr${chr}
  rm ${path_comb}/comp/impute5_chr${chr}.fam; cp ${path_comb}/comp/gigi2_chr${chr}.fam ${path_comb}/comp/impute5_chr${chr}.fam #Because samples aren't all correcly labeled in Impute5 data. 
  rm ${path_comb}/comp/gigi2_chr${chr}.bim; cp ${path_comb}/comp/impute5_chr${chr}.bim ${path_comb}/comp/gigi2_chr${chr}.bim #Because INDELS alleles are not labelled the same between Impute5 and Gigi2 (we used Impute5 as the header reference when combining). 
  plink --bfile ${path_comb}/impute5_gigi2_combined_chr${chr} --bmerge ${path_comb}/comp/gigi2_chr${chr} --merge-mode 7 --keep-allele-order --out ${path_comb}/comp/gigi2_chr${chr}_comp_rare
  plink --bfile ${path_comb}/impute5_gigi2_combined_chr${chr} --bmerge ${path_comb}/comp/impute5_chr${chr} --merge-mode 7 --keep-allele-order --out ${path_comb}/comp/impute5_chr${chr}_comp_rare
  rm ${path_comb}/comp/gigi2_chr${chr}.*; rm ${path_comb}/comp/impute5_chr${chr}.*; gzip ${path_comb}/comp/impute5_chr${chr}_comp_rare.diff; gzip ${path_comb}/comp/gigi2_chr${chr}_comp_rare.diff
  grep "rate" ${path_comb}/comp/impute5_chr${chr}_comp_rare.log | cut -d " " -f 8 | sed "s/.$//" >> ${path_comb}/comp/impute5_comp_rare.txt
  grep "rate" ${path_comb}/comp/gigi2_chr${chr}_comp_rare.log | cut -d " " -f 8 | sed "s/.$//" >> ${path_comb}/comp/gigi2_comp_rare.txt

  #Common variants
  plink2 --vcf ${path_comb}/impute5_imputed_chr${chr}.vcf.gz --exclude ${path_comb}/comp/gigi2_chr${chr}_rare.snplist --make-bed --out ${path_comb}/comp/impute5_chr${chr}
  plink2 --vcf ${path_comb}/gigi2_imputed_chr${chr}.vcf.gz --exclude ${path_comb}/comp/gigi2_chr${chr}_rare.snplist --make-bed --out ${path_comb}/comp/gigi2_chr${chr}
  rm ${path_comb}/comp/impute5_chr${chr}.fam; cp ${path_comb}/comp/gigi2_chr${chr}.fam ${path_comb}/comp/impute5_chr${chr}.fam #Because samples aren't all correcly labeled in Impute5 data. 
  rm ${path_comb}/comp/gigi2_chr${chr}.bim; cp ${path_comb}/comp/impute5_chr${chr}.bim ${path_comb}/comp/gigi2_chr${chr}.bim #Because INDELS alleles are not labelled the same between Impute5 and Gigi2 (we used Impute5 as the header reference when combining). 
  plink --bfile ${path_comb}/impute5_gigi2_combined_chr${chr} --bmerge ${path_comb}/comp/gigi2_chr${chr} --merge-mode 7 --keep-allele-order --out ${path_comb}/comp/gigi2_chr${chr}_comp_common
  plink --bfile ${path_comb}/impute5_gigi2_combined_chr${chr} --bmerge ${path_comb}/comp/impute5_chr${chr} --merge-mode 7 --keep-allele-order --out ${path_comb}/comp/impute5_chr${chr}_comp_common
  rm ${path_comb}/comp/gigi2_chr${chr}.*; rm ${path_comb}/comp/impute5_chr${chr}.*; gzip ${path_comb}/comp/impute5_chr${chr}_comp_common.diff; gzip ${path_comb}/comp/gigi2_chr${chr}_comp_common.diff
  grep "rate" ${path_comb}/comp/impute5_chr${chr}_comp_common.log | cut -d " " -f 8 | sed "s/.$//" >> ${path_comb}/comp/impute5_comp_common.txt
  grep "rate" ${path_comb}/comp/gigi2_chr${chr}_comp_common.log | cut -d " " -f 8 | sed "s/.$//" >> ${path_comb}/comp/gigi2_comp_common.txt
done
#NOTE: in .diff files, NEW = genotypes from GIGI2 or Impute5, OLD = genotypes from the combining method.

for chr in {1..22}; do wc -l ${path_comb}/comp/gigi2_chr${chr}_rare.snplist | awk '{print $1}' >> ${path_comb}/comp/n_rare.txt; done
awk '{ sum += $1 } END { print sum }' ${path_comb}/comp/n_rare.txt
for chr in {1..22}; do zcat ${path_comb}/comp/gigi2_chr${chr}_comp.diff.gz | wc -l | awk '{print $1}' >> ${path_comb}/comp/n_gigi2_diff.txt; done
awk '{ sum += $1 } END { print sum }' ${path_comb}/comp/n_gigi2_diff.txt
for chr in {1..22}; do zcat ${path_comb}/comp/impute5_chr${chr}_comp.diff.gz | wc -l  | awk '{print $1}' >> ${path_comb}/comp/n_impute5_diff.txt; done
awk '{ sum += $1 } END { print sum }' ${path_comb}/comp/n_impute5_diff.txt
for chr in {1..22}; do zcat ${path_comb}/comp/gigi2_chr${chr}_comp_rare.diff.gz | wc -l | awk '{print $1}' >> ${path_comb}/comp/n_gigi2_rare_diff.txt; done
awk '{ sum += $1 } END { print sum }' ${path_comb}/comp/n_gigi2_rare_diff.txt
for chr in {1..22}; do zcat ${path_comb}/comp/impute5_chr${chr}_comp_rare.diff.gz | wc -l  | awk '{print $1}' >> ${path_comb}/comp/n_impute5_rare_diff.txt; done
awk '{ sum += $1 } END { print sum }' ${path_comb}/comp/n_impute5_rare_diff.txt

awk '{ sum += $1 } END { print sum/NR }' ${path_comb}/comp/impute5_comp.txt
awk '{ sum += $1 } END { print sum/NR }' ${path_comb}/comp/gigi2_comp.txt
awk '{ sum += $1 } END { print sum/NR }' ${path_comb}/comp/impute5_comp_rare.txt
awk '{ sum += $1 } END { print sum/NR }' ${path_comb}/comp/gigi2_comp_rare.txt
awk '{ sum += $1 } END { print sum/NR }' ${path_comb}/comp/impute5_comp_common.txt
awk '{ sum += $1 } END { print sum/NR }' ${path_comb}/comp/gigi2_comp_common.txt


#We wish to obtain a complete imputed file of our RELATED GENOTYPES OR SEQUENCED subjects.
#The combination is used on the intercection of GIGI2 (related subjects) and Impute5 (genotypes subjects) imputed subjects (related genotyped subjects).
#Thus, to the combined imputation file, will need to add the sequenced subjects.
#At the end, we should have every subject in this list /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/iid_geno_plus_seq_only_related.txt
mkdir ${path_comb}/merged_with_seq
path_phasing=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/phasing
for chr in {22..1}
do
  #Make a vcf out of the bfiles 
  plink --bfile ${path_comb}/impute5_gigi2_combined_chr${chr} --keep-allele-order --recode vcf-iid --out ${path_comb}/merged_with_seq/impute5_gigi2_combined_chr${chr}
  bgzip ${path_comb}/merged_with_seq/impute5_gigi2_combined_chr${chr}.vcf
  bcftools index -t -f ${path_comb}/merged_with_seq/impute5_gigi2_combined_chr${chr}.vcf.gz
  
  #First, we need to remove the imputed genotypes of the sequenced subjects
  comm -1 -2 <(bcftools query -l ${path_phasing}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf | grep -v "Sample" | cut -d "_" -f 2 | sort) <(bcftools query -l ${path_comb}/merged_with_seq/impute5_gigi2_combined_chr${chr}.vcf.gz | cut -d "_" -f 2 | sort) > ${path_comb}/merged_with_seq/seq_ID_to_exclude_chr${chr}.txt
  egrep -f <(sed "s/^/_/g" ${path_comb}/merged_with_seq/seq_ID_to_exclude_chr${chr}.txt | sed "s/$/\$/g") <(bcftools query -l ${path_comb}/merged_with_seq/impute5_gigi2_combined_chr${chr}.vcf.gz) > ${path_comb}/merged_with_seq/seq_samples_to_exclude_chr${chr}.txt
  bcftools view -S ^${path_comb}/merged_with_seq/seq_samples_to_exclude_chr${chr}.txt -O z -o ${path_comb}/merged_with_seq/impute5_gigi2_combined_without_seq_chr${chr}.vcf.gz ${path_comb}/merged_with_seq/impute5_gigi2_combined_chr${chr}.vcf.gz
  bcftools index -t -f ${path_comb}/merged_with_seq/impute5_gigi2_combined_without_seq_chr${chr}.vcf.gz
  
  #Second, from the sequencing data, we remove the CaG samples and the Intercepts
  bcftools query -l ${path_phasing}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf | egrep "Sample|_19" > ${path_comb}/merged_with_seq/CaG_Intercept_samples_to_exclude_chr${chr}.txt
  bcftools view -S ^${path_comb}/merged_with_seq/CaG_Intercept_samples_to_exclude_chr${chr}.txt -i "ID=@${path_comb}/gigi2_variants_chr${chr}.txt" -O z -o  ${path_comb}/merged_with_seq/seq_FINAL_only_related_PHASED_${mask}_in_genmap_chr${chr}.vcf.gz  ${path_phasing}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf
  bcftools index -t -f ${path_comb}/merged_with_seq/seq_FINAL_only_related_PHASED_${mask}_in_genmap_chr${chr}.vcf.gz

  #Third, we merge this imputed dataset with the whole sequencing dataset
  bcftools merge -O z -o ${path_comb}/merged_with_seq/impute5_gigi2_combined_seq_chr${chr}.vcf.gz ${path_comb}/merged_with_seq/impute5_gigi2_combined_without_seq_chr${chr}.vcf.gz ${path_comb}/merged_with_seq/seq_FINAL_only_related_PHASED_${mask}_in_genmap_chr${chr}.vcf.gz
  bcftools index -t -f ${path_comb}/merged_with_seq/impute5_gigi2_combined_seq_chr${chr}.vcf.gz

  #Finally, we remove the intermediate dataset 
  rm ${path_comb}/merged_with_seq/impute5_gigi2_combined_without_seq_chr${chr}.vcf.gz ${path_comb}/merged_with_seq/impute5_gigi2_combined_without_seq_chr${chr}.vcf.gz.tbi
  rm ${path_comb}/merged_with_seq/seq_FINAL_only_related_PHASED_${mask}_in_genmap_chr${chr}.vcf.gz ${path_comb}/merged_with_seq/seq_FINAL_only_related_PHASED_${mask}_in_genmap_chr${chr}.vcf.gz.tbi
  rm ${path_comb}/merged_with_seq/*_to_exclude_chr${chr}.txt
  rm ${path_comb}/merged_with_seq/impute5_gigi2_combined_chr${chr}.*
done
