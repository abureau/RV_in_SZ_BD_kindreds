#We wish to verify the MAF of the SNPs not found in CaG with the gnomAD data.

#Load modules
module load r/4.2 StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16

#Results will be exported here:
path_data=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/RV
mkdir ${path_data}

#We first work on seq_FINAL_with_mask.vcf.gz
mask=with_mask
seq_data_qc=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/QC/seq_FINAL
seq_data_qc_i=${seq_data_qc}_${mask}.vcf.gz

#Path to gnomAD data
path_gnomad=/lustre03/project/6033529/GENOMES

#We need to extract our SNPs only in the CaG data.
#Note that the chromosomes are specified in the same way as in the CaG data. We can reuse our last extraction list.
pos_to_extract=${path_data}/pos_to_extract_${mask}.txt

#We create a new dataset from the original gnomAD data containing our variants only.
#As the MAF are part of the vcf file and it doesn't need to be computed, we drop the individual genotypes (-G).
#Then, we create a list of rare variants in gnomAD.
for chr in {1..22} X
do
  sbatch --time=5:00:00 --mem=50G --export=path_data=${path_data},pos_to_extract=${pos_to_extract},mask=${mask},path_gnomad=${path_gnomad},chr=${chr} /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/gnomad_data_seq_pos.sh
done
for chr in {1..22} X; do cat ${path_data}/chr${chr}_gnomAD_FINAL_${mask}.vcf | cut -f-2 | grep -v "^#" >> ${path_data}/seq_pos_in_gnomAD_${mask}.txt; done
awk -F'\t' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${path_data}/seq_pos_in_gnomAD_${mask}.txt ${pos_to_extract} > ${path_data}/seq_pos_NOT_in_gnomAD_${mask}.txt

#We want to know what happens to the variants not in CaG...
for chr in {1..22} X
do
  bcftools view -G -T /lustre03/project/6033529/quebec_10x/data/freeze/RV/seq_pos_NOT_in_CaG_${mask}.txt -O v -o ${path_data}/chr${chr}_gnomAD_NOT_in_CaG_${mask}.vcf ${path_data}/chr${chr}_gnomAD_FINAL_${mask}.vcf
  bcftools view -i 'AF_nfe < 0.01' ${path_data}/chr${chr}_gnomAD_NOT_in_CaG_${mask}.vcf | cut -f 1,2,4,5 | grep -v "^#" > ${path_data}/RV_gnomAD_NOT_in_CaG_chr${chr}_${mask}.snplist
  bcftools view -i 'AF_nfe >= 0.01' ${path_data}/chr${chr}_gnomAD_NOT_in_CaG_${mask}.vcf | cut -f 1,2,4,5 | grep -v "^#" > ${path_data}/CV_gnomAD_NOT_in_CaG_chr${chr}_${mask}.snplist
  cat ${path_data}/chr${chr}_gnomAD_NOT_in_CaG_${mask}.vcf | cut -f-2 | grep -v "^#" >> ${path_data}/seq_pos_in_gnomAD_NOT_in_CaG_${mask}.txt
done
awk -F'\t' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${path_data}/seq_pos_in_gnomAD_NOT_in_CaG_${mask}.txt ${path_data}/seq_pos_NOT_in_CaG_${mask}.txt > ${path_data}/seq_pos_NOT_in_gnomAD_NOT_in_CaG_${mask}.txt

#Use R to compare the lists easily.
Rscript gnomad.R