path_phasing=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/phasing
path_impu=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/imputation_pop
path_data=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/genotyping
mask=with_mask
mkdir ${path_impu}/logs

#Compute the recombination rate and generate the genetic map in the correct format
gen_map=/lustre03/project/6033529/IBD_epilepsy/data/genetic_map_grch38
mkdir=${path_impu}/gen_map
if [ ! -s "${path_impu}/gen_map" ]
then
  Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/gen_map_phaseit2_format.R -ref_dir ${gen_map} -out ${path_impu}/gen_map
fi

#Genotyping data
puce=omni
for chr in {22..1}
do
  impute5 --h ${path_phasing}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf \
	  --g ${path_data}/${puce}_PHASED_grch38_in_genmap_chr${chr}.bcf \
	  --m ${path_impu}/gen_map/plink.chr${chr}.GRCh38.map \
	  --l ${path_impu}/logs/log_${puce}_chr${chr} \
	  --r ${chr} \
	  --o ${path_impu}/${puce}_chr${chr}_IMPUTED.vcf.gz \
          --out-gp-field
  bcftools index -t -f ${path_impu}/${puce}_chr${chr}_IMPUTED.vcf.gz
done

puce=gsa
for chr in {22..1}
do
  impute5 --h ${path_phasing}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf \
	  --g ${path_data}/${puce}_PHASED_grch38_in_genmap_chr${chr}.bcf \
	  --m ${path_impu}/gen_map/plink.chr${chr}.GRCh38.map \
	  --l ${path_impu}/logs/log_${puce}_chr${chr} \
	  --r ${chr} \
	  --o ${path_impu}/${puce}_chr${chr}_IMPUTED.vcf.gz \
          --out-gp-field
  bcftools index -t -f ${path_impu}/${puce}_chr${chr}_IMPUTED.vcf.gz
done

#Combine the imputations for the 2 chips.
for chr in {22..1}
do
  bcftools merge ${path_impu}/omni_chr${chr}_IMPUTED.vcf.gz ${path_impu}/gsa_chr${chr}_IMPUTED.vcf.gz -O z -o ${path_impu}/impute5_imputed_chr${chr}.vcf.gz
  bcftools index -t -f ${path_impu}/impute5_imputed_chr${chr}.vcf.gz
done

#Removed the indep. chips files
rm ${path_impu}/gsa_chr*_IMPUTED.* ${path_impu}/omni_chr*_IMPUTED.*

#Get the INFO scores of every variants
for chr in {22..1}
do
  final_data=${path_impu}/impute5_imputed_chr${chr}.vcf.gz
  paste <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${final_data}) <(bcftools query -f '%INFO\n' ${final_data} | tr ";" "\n" | grep INFO) > ${path_impu}/impute5_imputed_chr${chr}_INFO.txt
done

#PAS ROULÉ ENCORE.
#We wish to create datasets where the imputed data are merged with the sequenced data of our cohort
mkdir ${path_impu}/merged
for chr in {22..1}
do
  #First, we need to remove the imputed genotypes of the sequenced subjects
  comm -1 -2 <(bcftools query -l ${path_phasing}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf | grep -v "Sample" | cut -d "_" -f 2 | sort) <(bcftools query -l ${path_impu}/impute5_imputed_chr${chr}.vcf.gz | cut -d "_" -f 2 | sort) > ${path_impu}/merged/seq_ID_to_exclude_chr${chr}.txt
  egrep -f <(sed "s/^/_/g" ${path_impu}/merged/seq_ID_to_exclude_chr${chr}.txt | sed "s/$/\$/g") <(bcftools query -l ${path_impu}/impute5_imputed_chr${chr}.vcf.gz) > ${path_impu}/merged/seq_samples_to_exclude_chr${chr}.txt
  bcftools view -S ^${path_impu}/merged/seq_samples_to_exclude_chr${chr}.txt -O z -o ${path_impu}/merged/impute5_imputed_without_seq_chr${chr}.vcf.gz ${path_impu}/impute5_imputed_chr${chr}.vcf.gz
  bcftools index -t -f ${path_impu}/merged/impute5_imputed_without_seq_chr${chr}.vcf.gz
  
  #Second, from the sequencing data, we remove the CaG samples
  bcftools query -l ${path_phasing}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf | grep "Sample" > ${path_impu}/merged/CaG_samples_to_exclude_chr${chr}.txt
  bcftools view -S ^${path_impu}/merged/CaG_samples_to_exclude_chr${chr}.txt -O z -o  ${path_impu}/merged/seq_FINAL_without_CaG_PHASED_${mask}_in_genmap_chr${chr}.vcf.gz  ${path_phasing}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf
  bcftools index -t -f ${path_impu}/merged/seq_FINAL_without_CaG_PHASED_${mask}_in_genmap_chr${chr}.vcf.gz

  #Third, we merge this imputed dataset with the whole sequencing dataset
  bcftools merge -O z -o ${path_impu}/merged/impute5_imputed_seq_chr${chr}.vcf.gz ${path_impu}/merged/impute5_imputed_without_seq_chr${chr}.vcf.gz ${path_impu}/merged/seq_FINAL_without_CaG_PHASED_${mask}_in_genmap_chr${chr}.vcf.gz
  bcftools index -t -f ${path_impu}/merged/impute5_imputed_seq_chr${chr}.vcf.gz

  #Finally, we remove the intermediate dataset 
  rm ${path_impu}/merged/impute5_imputed_without_seq_chr${chr}.vcf.gz ${path_impu}/merged/impute5_imputed_without_seq_chr${chr}.vcf.gz.tbi
  rm ${path_impu}/merged/seq_FINAL_without_CaG_PHASED_${mask}_in_genmap_chr${chr}.vcf.gz ${path_impu}/merged/seq_FINAL_without_CaG_PHASED_${mask}_in_genmap_chr${chr}.vcf.gz.tbi
  rm ${path_impu}/merged/*_to_exclude_chr${chr}.txt
done