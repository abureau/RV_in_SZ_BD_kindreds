#Find the rare variants in CARTaGENE among our sequenced variants using the EUR subjects only.

#Load modules
module load r/4.2 StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16

#Results will be exported here:
path_data=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/RV
mkdir ${path_data}

#Our QCed sequencing data are found here:
#A loop will be used to complete the real dataset name: seq_FINAL_without_mask.vcf.gz and seq_FINAL_with_mask.vcf.gz
seq_data_qc=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/QC/seq_FINAL

#CARTaGENE sequencing data are found here:
path_cag=/lustre03/project/6033529/cartagene/wgs_qc_freeze1/vcfs_freeze1_20221206

#We need to keep the EUR subjects only in the CaG data. 
#In Claudia's CaG freeze, we find 2173 subjects while CaG published a total of 2185. 10 EUR are missing.
grep -f /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/cag_samples/2023-11-30_1889_canadians_ID.txt <(bcftools query -l ${path_cag}/chr${chr}_FINAL.vcf.gz) > ${path_data}/EUR_cag_sample.txt

#On 2024-01-04, We do not generate the list for the "without_mask" version of the WGS. It needs to be done if ever we wish to use it.
#for mask in without_mask with_mask
for mask in with_mask
do
  seq_data_qc_i=${seq_data_qc}_${mask}.vcf.gz

  #First, we need to extract our SNPs only in the CaG data.
  #Note that in CaG data, the 23th chromosome is written as "X", and there are "chr" prefixes. We need to adjust that in the extraction list.
  zcat ${seq_data_qc_i} | cut -f-2 | grep -v "^#" | sed 's/^23/X/' | sed 's/^/chr/' > ${path_data}/pos_to_extract_${mask}.txt

  #Create the desired subsets of the original CaG data.
  #Print a list of every variants in our data that is found in CaG
  for chr in {1..22} X
  do
    sbatch --time=03:00:00 --mem=30G --export=path_data=${path_data},mask=${mask},path_cag=${path_cag},chr=${chr} /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/cag_data_seq_pos.sh
  done
  for chr in {1..22} X; do zcat ${path_data}/chr${chr}_CaG_FINAL_${mask}.vcf.gz | cut -f-2 | grep -v "^#" >> ${path_data}/seq_pos_in_CaG_${mask}.txt; done
  awk -F'\t' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${path_data}/seq_pos_in_CaG_${mask}.txt ${path_data}/pos_to_extract_${mask}.txt > ${path_data}/seq_pos_NOT_in_CaG_${mask}.txt

  #Isolate the rare variants.
  for chr in {1..22} X
  do
    plink -vcf ${path_data}/chr${chr}_CaG_FINAL_${mask}.vcf.gz --maf 0.01 --write-snplist --double-id --out ${path_data}/CV_CaG_chr${chr}_${mask}
    plink -vcf ${path_data}/chr${chr}_CaG_FINAL_${mask}.vcf.gz --max-maf 0.01 --write-snplist --double-id --out ${path_data}/RV_CaG_chr${chr}_${mask}
  done
done

mkdir ${path_data}/log
mv ${path_data}/*.log ${path_data}/log
rm ${path_data}/*.nosex