#Find the rare variants in CARTaGENE among our sequenced variants using EUR subjects only.

#Load modules
module load r/4.2 StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16

#Results will be exported here:
pathData=/lustre03/project/6033529/quebec_10x/data/freeze/RV

#Our QCed sequencing data are found here:
#A loop will be used to complete the real dataset name: seq_FINAL_without_mask.vcf.gz and seq_FINAL_with_mask.vcf.gz
seqQcData=/lustre03/project/6033529/quebec_10x/data/freeze/QC/seq_FINAL

#CARTaGENE sequencing data are found here:
pathCaG=/lustre03/project/6033529/cartagene/wgs_qc_freeze1/vcfs_freeze1_20221206

#We need to keep the 1756 EUR subjects only ine the CARTaGENE data. 
cut -f 1 /lustre03/project/6033529/cartagene/wgs_qc_freeze1/summary_statistics/list_sample_FINALVCF_FrenchCanada.ID > ${pathData}/EUR_CaG_sample.txt

for mask in without_mask with_mask
do
  seqQcData_i=${seqQcData}_${mask}.vcf.gz

  #First, we need to extract our SNPs only in the CaG data.
  #Note that in CaG data, the 23th chromosome is written as "X", and there are "chr" prefixes. We need to adjust that in the extraction list.
  zcat ${seqQcData_i} | cut -f-2 | grep -v "^#" | sed 's/^23/X/' | sed 's/^/chr/' > ${pathData}/pos_to_extract_${mask}.txt

  #Create the desired subsets of the original CaG data.
  #Print a list of every variants in our data that our found in CaG
  for chr in {1..22} X
  do
    bcftools view -S ${pathData}/EUR_CaG_sample.txt -T ${pathData}/pos_to_extract_${mask}.txt -O z -o ${pathData}/chr${chr}_CaG_FINAL_${mask}.vcf.gz ${pathCaG}/chr${chr}_FINAL.vcf.gz
    zcat ${pathData}/chr${chr}_CaG_FINAL_${mask}.vcf.gz | cut -f-2 | grep -v "^#" >> ${pathData}/seq_pos_in_CaG_${mask}.txt
  done
  awk -F'\t' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${pathData}/seq_pos_in_CaG_${mask}.txt ${pathData}/pos_to_extract_${mask}.txt > ${pathData}/seq_pos_NOT_in_CaG_${mask}.txt

  #Isolate the rare variants.
  for chr in {1..22} X
  do
    plink -vcf ${pathData}/chr${chr}_CaG_FINAL_${mask}.vcf.gz --maf 0.01 --write-snplist --double-id --out ${pathData}/CV_CaG_chr${chr}_${mask}
    plink -vcf ${pathData}/chr${chr}_CaG_FINAL_${mask}.vcf.gz --max-maf 0.01 --write-snplist --double-id --out ${pathData}/RV_CaG_chr${chr}_${mask}
  done
done

mv ${pathData}/*.log ${pathData}/log
rm ${pathData}/*.nosex               