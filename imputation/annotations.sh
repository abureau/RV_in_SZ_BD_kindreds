export JAVA_TOOL_OPTIONS="-Xmx8g"
export _JAVA_OPTIONS="-Xmx8g"

#Path to data and objects
path_data=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag
path=${path_data}/annotations; mkdir ${path}
ccds=/lustre03/project/6033529/GENOMES/ccds_grch38.bed
fasta=/lustre03/project/6033529/ctb-pacoss/GENOME_REF/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna
snpeff=/lustre03/project/6033529/pregene_2023/results/snpEff.config

#Do we wish to annotate all variants? if false, only the ones not found in CaG are annotated. (true or false)
anno_all=true

for chr in {22..1}
do

  echtvar_file=/lustre03/project/6033529/GENOMES/echtvar/gnomad.genomes.v3.1.2.sites.chr${chr}.vcf.echtvar.zip

  if [ "$anno_all" = "true" ]; then
    #If we wish to annotate all variants:
    out_name=seq_FINAL
    plink2 --vcf ${path_data}/QC/seq_FINAL_chr${chr}_with_mask.vcf.gz --make-bed --out ${path}/${out_name}_chr${chr}_with_mask
  else
    out_name=${out_name}
    #From the final sequenced dataset, we only keep the variants that aren't in CARTaGENE because we already have the annotation for them
    awk '{print $1, $2, $2, NR}' ${path_data}/RV/cag/seq_pos_NOT_in_CaG_chr${chr}_with_mask.txt > ${path}/seq_pos_NOT_in_CaG_chr${chr}_with_mask.extract
    plink2 --vcf ${path_data}/QC/seq_FINAL_chr${chr}_with_mask.vcf.gz --extract range ${path}/seq_pos_NOT_in_CaG_chr${chr}_with_mask.extract --make-bed --out ${path}/${out_name}_chr${chr}_with_mask
  fi

  #Check
  wc -l ${path}/${out_name}_chr${chr}_with_mask.bim
  wc -l ${path}/seq_pos_NOT_in_CaG_chr${chr}_with_mask.extract

  #step 1: take only the coding regions
  sed -i 's/^/chr/' ${path}/${out_name}_chr${chr}_with_mask.bim
  cut -f-3 ${ccds} | awk '{print $1, $2, $3, NR}' > ${path}/ccds_grch38.extract
  plink2 --bfile ${path}/${out_name}_chr${chr}_with_mask --extract range ${path}/ccds_grch38.extract --recode vcf bgz --out ${path}/${out_name}_chr${chr}_ccds_with_mask
  bcftools index -t -f ${path}/${out_name}_chr${chr}_ccds_with_mask.vcf.gz

  #Transform the vcf in bcf to make it work for echtvar
  bcftools norm -m - ${path}/${out_name}_chr${chr}_ccds_with_mask.vcf.gz -w 10000 -O b -o ${path}/${out_name}_chr${chr}_ccds_with_mask.bcf

  #Use echtvar
  /lustre03/project/6033529/SOFT/echtvar anno -e ${echtvar_file} ${path}/${out_name}_chr${chr}_ccds_with_mask.bcf ${path}/${out_name}_chr${chr}_ccds_annoted_with_mask.vcf.gz

  #Take the file created in by etchvar and add dbNSFP4.3a annotation on top of echtvar and CCDS
  #chr by chr.
  dbNSFP=/lustre03/project/6033529/pregene_2023/data/dbNSFP/dbNSFP4.3a/dbNSFP4.3a_variant.chr${chr}.gz
  plink2 --vcf ${path}/${out_name}_chr${chr}_ccds_annoted_with_mask.vcf.gz --chr ${chr} --recode vcf bgz --out ${path}/${out_name}_chr${chr}_ccds_annoted_chr_with_mask
  java -Xmx8g -jar  $EBROOTSNPEFF/SnpSift.jar dbnsfp -v -db ${dbNSFP} ${path}/${out_name}_chr${chr}_ccds_annoted_with_mask.vcf.gz > ${path}/${out_name}_chr${chr}_ccds_annoted_dbNSFP4_with_mask.vcf
  bgzip ${path}/${out_name}_chr${chr}_ccds_annoted_dbNSFP4_with_mask.vcf
  bcftools index -t -f ${path}/${out_name}_chr${chr}_ccds_annoted_dbNSFP4_with_mask.vcf.gz

  ## Extra step
  java -Xmx8g -jar $EBROOTSNPEFF/snpEff.jar -c ${snpeff} -v GRCh38.99 ${path}/${out_name}_chr${chr}_ccds_annoted_dbNSFP4_with_mask.vcf.gz  > ${path}/${out_name}_chr${chr}_ccds_annoted_dbNSFP4_snpeff_with_mask.vcf
  bcftools view -Oz ${path}/${out_name}_chr${chr}_ccds_annoted_dbNSFP4_snpeff_with_mask.vcf > ${path}/${out_name}_chr${chr}_ccds_annoted_dbNSFP4_snpeff_with_mask.vcf.gz
  bcftools index -t -f ${path}/${out_name}_chr${chr}_ccds_annoted_dbNSFP4_snpeff_with_mask.vcf.gz
  
  # Remove unnecessary files
  rm ${path}/${out_name}_chr${chr}_with_mask.bim ${path}/${out_name}_chr${chr}_with_mask.bed ${path}/${out_name}_chr${chr}_with_mask.fam ${path}/${out_name}_chr${chr}_with_mask.log
  rm ${path}/ccds_grch38.extract
  rm ${path}/${out_name}_chr${chr}_ccds_annoted_dbNSFP4_snpeff_with_mask.vcf
done

#Not to overuse the disk.
rm ${path}/${out_name}_chr*_ccds_annoted_dbNSFP4_with_mask.vcf.gz
rm ${path}/${out_name}_chr*_ccds_annoted_dbNSFP4_with_mask.vcf.gz.tbi 
rm ${path}/${out_name}_chr*_ccds_annoted_chr_with_mask.vcf.gz
rm ${path}/${out_name}_chr*_ccds_annoted_chr_with_mask.log
rm ${path}/${out_name}_chr*_ccds_annoted_with_mask.vcf.gz
rm ${path}/${out_name}_chr*_ccds_with_mask.vcf.gz
rm ${path}/${out_name}_chr*_ccds_with_mask.vcf.gz.tbi
rm ${path}/${out_name}_chr*_ccds_with_mask.bcf
rm ${path}/${out_name}_chr*_ccds_with_mask.log