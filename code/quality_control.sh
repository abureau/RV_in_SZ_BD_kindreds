#WGS Quality Control
#JR, 2023-11-17
#I follow Guillaume Lettre's steps.
#In this code, I modify the data.
#I understand that a first QC is done on the individuals using PLINK, after that, a small QC is done on the SNPs.

#Load modules
module load r/4.2 StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16

#Define paths
#The freeze can be found here:
seq_data=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/freeze_24112023_schizo.vcf.gz
#The QCed data will be generated here:
path_data=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/QC
mkdir ${path_data}
mkdir ${path_data}/logs

#Modify the sample label problems in the initial dataset
#M-2452 becomes M-2338; M-2773 becomes M-2733; M-2733 becomes M-2773.
sample_label_change=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/freeze_24112023_sample_to_change.txt
echo "M-2452 M-2338" > ${sample_label_change}
echo "M-2773 M-2733" >> ${sample_label_change}
echo "M-2733 M-2773" >> ${sample_label_change}
bcftools reheader -o /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/freeze_24112023_schizo_changed.vcf.gz --sample ${sample_label_change} ${seq_data}

#To be sure that the change was made correctly.
bcftools query -l /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/freeze_24112023_schizo_changed.vcf.gz >> /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/freeze_24112023_sample_list_changed.txt
diff /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/freeze_24112023_sample_list.txt /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/freeze_24112023_sample_list_changed.txt
#IF THE MODIFICATIONS WERE MADE CORRECTLY, replace the old seq file by the new one.
mv ${seq_data} /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/freeze_24112023_schizo_old.vcf.gz
mv /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/freeze_24112023_schizo_changed.vcf.gz ${seq_data}
rm /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/freeze_24112023_schizo_old.vcf.gz

#Pedigree for the sequenced samples is here:
seq_ped=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/manip/seq_pedigree.txt
#Lets first add the pedigree information for the new subject M-2338.
echo -e "2338\tM-2338\t121\t2454\t2460\t2" >> ${seq_ped}

#Create multi-allelic PLINK files using the VCFs created by Illumina DRAGEN, filtering again on variant genotyping percent (at least 90%).
#These multi-allelic variants are not used for quality control. 
#We modify the pedigree right away and set as missing the Mendel errors.
sbatch --export=seq_data=${seq_data},path_data=${path_data},split=multi /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/seq_allelic_split.sh
plink --vcf ${path_data}/seq_multiallelic.vcf --geno 0.1 --set-missing-var-ids @:#:\$1:\$2 --make-bed --keep-allele-order --chr 1-22, X --allow-extra-chr --out ${path_data}/seq_multiallelic; mv ${path_data}/seq_multiallelic.log ${path_data}/logs/seq_multiallelic.log
Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/update_seq_fam.R -fam ${path_data}/seq_multiallelic.fam -ref ${seq_ped}
plink --bfile ${path_data}/seq_multiallelic --set-me-missing --mendel-duos --mendel-multigen --make-bed --out ${path_data}/seq_multiallelic; mv ${path_data}/seq_multiallelic.log ${path_data}/logs/seq_multiallelic_mendel_errors.log; rm ${path_data}/seq_multiallelic.*~
#We create a .vcf that will be used later to create the final QCed dataset.
plink --bfile ${path_data}/seq_multiallelic --recode vcf --out ${path_data}/seq_multiallelic
rm ${path_data}/seq_multiallelic.log ${path_data}/seq_multiallelic.bim ${path_data}/seq_multiallelic.bed ${path_data}/seq_multiallelic.fam ${path_data}/seq_multiallelic.nosex

#Filter VCFs created by Illumina DRAGEN on variant genotyping percent (at least 90%) and keep variants with 2 alleles.
#Set as missing the Mendel errors. We need the pedigree to be updated. 
sbatch --export=seq_data=${seq_data},path_data=${path_data},split=bi /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/seq_allelic_split.sh
plink --vcf ${path_data}/seq.vcf --geno 0.1 --set-missing-var-ids @:#:\$1:\$2 --make-bed --keep-allele-order --chr 1-22, X --allow-extra-chr --out ${path_data}/seq; mv ${path_data}/seq.log ${path_data}/logs/seq.log
Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/update_seq_fam.R -ref ${seq_ped} -fam ${path_data}/seq.fam
plink --bfile ${path_data}/seq --set-me-missing --mendel-duos --mendel-multigen --make-bed --out ${path_data}/seq; mv ${path_data}/seq.log ${path_data}/logs/seq_mendel_errors.log; rm ${path_data}/seq.*~
#We create a .vcf that will be used later to create the final QCed dataset.
plink --bfile ${path_data}/seq --recode vcf --out ${path_data}/seq

#The following steps are informative only.
#Remove variants with "NON_REF" as alternate allele. (11 variants)
grep NON_REF ${path_data}/seq.bim | cut -f 2 >> ${path_data}/var_to_remove_non_ref.txt
plink --bfile ${path_data}/seq --exclude ${path_data}/var_to_remove_non_ref.txt --make-bed --keep-allele-order --out ${path_data}/seq; mv ${path_data}/seq.log ${path_data}/logs/seq_non_ref.log; rm ${path_data}/seq.*~
grep NON_REF ${path_data}/seq_multiallelic.vcf | cut -f 2 | wc -l

#Compute variant missingness and individual missingness.
plink --bfile ${path_data}/seq --missing --out ${path_data}/seq; mv ${path_data}/seq.log ${path_data}/logs/seq_missingness.log
awk '$6 > 0.05' ${path_data}/seq.imiss 
#None, final freeze

#Compute Hardy-Weinberg equilibrium.  Not used for filtering because of mixed populations in the dataset.
plink --bfile ${path_data}/seq --hardy midp --out ${path_data}/seq; mv ${path_data}/seq.log ${path_data}/logs/seq_HW.log
awk '$9 < 0.05' ${path_data}/seq.hwe | wc -l 

#Create a pruned version of the dataset.
plink --bfile ${path_data}/seq --indep-pairwise 500 100 0.3 --out ${path_data}/seq_indep; mv ${path_data}/seq_indep.log ${path_data}/logs/seq_indep.log
plink --bfile ${path_data}/seq --extract ${path_data}/seq_indep.prune.in --make-bed --out ${path_data}/seq_pruned; mv ${path_data}/seq_pruned.log ${path_data}/logs/seq_pruned.log

#Compute heterozygozity on the pruned version.
plink --bfile ${path_data}/seq_pruned --het --out ${path_data}/seq_pruned
awk '$6 > 0.4 || $6 < -0.4' ${path_data}/seq_pruned.het 
#None, final freeze

#Create QC'ed version.  Remove monomorphic variants.
#add to "list_sample_to_remove.txt" the above problematic ID
#final freeze, we wish to keep every subject.
plink --bfile ${path_data}/seq --geno 0.05 --maf 0.0000001 --make-bed --out ${path_data}/seq; mv ${path_data}/seq.log ${path_data}/logs/seq_mono.log; rm ${path_data}/seq.*~

#Check for sex concordance
plink --bfile ${path_data}/seq --check-sex --out ${path_data}/seq_sex_check

#From here, we create the final vcf.
#First, we remove the problematic samples identified in the PLINK QC above, variants with missingness > 0.05, variants with NON_REF as alt allele, and monomorphic variants.
#We also set the variant IDs, so that it's easy to tell which variants were split from multi-allelic ones to create bi-allelic.
plink --vcf ${path_data}/seq.vcf --geno 0.05 --recode vcf --out ${path_data}/seq_bial_step1; mv ${path_data}/seq_bial_step1.log ${path_data}/logs/seq_bial_step1.log
bcftools norm -m -both -O v -o ${path_data}/seq_bial_step2.vcf ${path_data}/seq_bial_step1.vcf
bcftools view -i 'ALT[0] != "NON_REF"' -c 1 -O v -o ${path_data}/seq_bial_step3.vcf ${path_data}/seq_bial_step2.vcf
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -O v -o ${path_data}/seq_bial_step4.vcf ${path_data}/seq_bial_step3.vcf
bcftools sort -O v -o ${path_data}/seq_bial_QCed.vcf ${path_data}/seq_bial_step4.vcf
bgzip ${path_data}/seq_bial_QCed.vcf
tabix -p vcf ${path_data}/seq_bial_QCed.vcf.gz

#If everything went right, to save space, remove step1 to step3 files.
rm ${path_data}/seq_bial_step*.vcf

#Now we deal with the multi-allelic file
plink --vcf ${path_data}/seq_multiallelic.vcf --geno 0.05 --recode vcf --out ${path_data}/seq_multial_step1; mv ${path_data}/seq_multial_step1.log ${path_data}/logs/seq_multial_step1.log
bcftools norm -m -both -O v -o ${path_data}/seq_multial_step2.vcf ${path_data}/seq_multial_step1.vcf
bcftools view -i 'ALT[0] != "NON_REF"' -c 1 -O v -o ${path_data}/seq_multial_step3.vcf ${path_data}/seq_multial_step2.vcf
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT\_m' -O v -o ${path_data}/seq_multial_step4.vcf ${path_data}/seq_multial_step3.vcf
bcftools sort -O v -o ${path_data}/seq_multial_QCed.vcf ${path_data}/seq_multial_step4.vcf
bgzip ${path_data}/seq_multial_QCed.vcf
tabix -p vcf ${path_data}/seq_multial_QCed.vcf.gz

#Again, if everything went right, to save space, remove step1 to step3 files.
rm ${path_data}/seq_multial_step*.vcf
rm ${path_data}/*.nosex

#Merge the bi-allelic part and the multi-allelic split part.
bcftools concat -a ${path_data}/seq_bial_QCed.vcf.gz ${path_data}/seq_multial_QCed.vcf.gz -O v -o ${path_data}/seq_merged.vcf 
bcftools sort -T /scratch/jasmric/tmp_bcftools_sort -O z -o ${path_data}/seq_merged.vcf.gz -m 9G ${path_data}/seq_merged.vcf 
rm ${path_data}/seq_merged.vcf 

#Verify if variants with NON_REF as an alternate allele are still remaining (maybe a bug in bcftools).
zcat ${path_data}/seq_merged.vcf.gz | grep -v -w '<NON_REF>' > ${path_data}/seq_merged_nonref.vcf
bgzip ${path_data}/seq_merged_nonref.vcf
tabix -p vcf ${path_data}/seq_merged_nonref.vcf.gz

#First, we remove variants located in the centromeres or in the ENCODE blacklist.
centro=/lustre03/project/6033529/GENOMES/centromeres_hg38.tsv
blacklist=/lustre03/project/6033529/GENOMES/ENCODE_blacklist_hg38.tsv
sed '1d' ${centro} | cut -f 2-5 | sed 's/^chr//' > ${path_data}/centro_blacklist_regions.txt
sed '1d' ${blacklist} | cut -f 1-4 | sed 's/^chr//' >> ${path_data}/centro_blacklist_regions.txt
plink --vcf ${path_data}/seq_merged_nonref.vcf.gz --exclude range ${path_data}/centro_blacklist_regions.txt --recode vcf --out ${path_data}/seq_FINAL_without_mask; mv ${path_data}/seq_FINAL_without_mask.log ${path_data}/logs/seq_FINAL_without_mask.log
bgzip ${path_data}/seq_FINAL_without_mask.vcf
tabix -p vcf ${path_data}/seq_FINAL_without_mask.vcf.gz
zcat ${path_data}/seq_FINAL_without_mask.vcf.gz | cut -f 1,2,4,5 | grep -v "^#" >> ${path_data}/seq_FINAL_without_mask.snplist

#Keep variants located in the following mask only.
mask=/lustre03/project/6033529/GENOMES/20160622.allChr.mask.bed
cut -f 1-3 ${mask} | sed 's/^chr//' | sed 's/^X/23/' > ${path_data}/mask_to_keep.txt
bcftools view -R ${path_data}/mask_to_keep.txt -O v -o ${path_data}/seq_FINAL_with_mask.vcf ${path_data}/seq_FINAL_without_mask.vcf.gz
bgzip ${path_data}/seq_FINAL_with_mask.vcf
tabix -p vcf ${path_data}/seq_FINAL_with_mask.vcf.gz
zcat ${path_data}/seq_FINAL_with_mask.vcf.gz | cut -f 1,2,4,5 | grep -v "^#" >> ${path_data}/seq_FINAL_with_mask.snplist

#The last dataset is most likely the one we will use. Let's generate the bfiles format
#Be sure that the pedigree is complete, it is lost when PLINK is used to go from bfiles to vcf.
plink --vcf ${path_data}/seq_FINAL_with_mask.vcf.gz --make-bed --out ${path_data}/seq_FINAL_with_mask
Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/update_seq_fam.R -fam ${path_data}/seq_FINAL_with_mask.fam -ref ${seq_ped} -match_var "iidGeno"
