#WGS Quality Control
#I follow Guillaume Lettre's steps.
#I understand that a first QC is done on the individuals using PLINK, after that, a small QC is done on the SNPs.

#Load modules
module load r/4.2 StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16

#Define paths
#The freeze can be found here:
seqData=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/freezes/freeze_16052023/freeze_16052023.vcf.gz
pathData=/lustre03/project/6033529/quebec_10x/data/freeze/QC

#Filter VCFs created by Illumina DRAGEN on variant genotyping percent (at least 90%) and keep variants with 2 alleles.
#Create a PLINK version
bcftools view -i 'N_ALT < 3' -O v -o ${pathData}/seq.vcf ${seqData};
plink --vcf ${pathData}/seq.vcf --geno 0.1 --set-missing-var-ids @:#:\$1:\$2 --make-bed --keep-allele-order --chr 1-22, X, Y, XY --allow-extra-chr --out ${pathData}/seq;
plink --bfile ${pathData}/seq --recode vcf --out ${pathData}/seq

#Create multi-allelic PLINK files, filtering again on variant genotyping percent (at least 90%).  These multi-allelic variants are not used for quality control.
bcftools view -i 'N_ALT > 2' -O v -o ${pathData}/seqMultiallelic.vcf ${seqData};
plink --vcf ${pathData}/seqMultiallelic.vcf --geno 0.1 --set-missing-var-ids @:#:\$1:\$2 --recode vcf --keep-allele-order --chr 1-22, X, Y, XY --allow-extra-chr --out ${pathData}/seqMultiallelic;

#Remove variants with "NON_REF" as alternate allele.
grep NON_REF ${pathData}/seq.bim 
#None, freeze 3

#Compute variant missingness and individual missingness.
plink --bfile ${pathData}/seq --missing --out ${pathData}/seq
awk '$6 > 0.05' ${pathData}/seq.imiss 
#None, freeze 3

#Compute Hardy-Weinberg equilibrium.  Not used for filtering because of mixed populations in the dataset.
plink --bfile ${pathData}/seq --hardy midp --out ${pathData}/seq

#Create a pruned version of the dataset.
plink --bfile ${pathData}/seq --indep-pairwise 500 100 0.3 --out ${pathData}/seq_indep
plink --bfile ${pathData}/seq --extract ${pathData}/seq_indep.prune.in --make-bed --out ${pathData}/seq_pruned

#Compute heterozygozity on the pruned version.
plink --bfile ${pathData}/seq_pruned --het --out ${pathData}/seq_pruned
awk '$6 > 0.4 || $6 < -0.4' ${pathData}/seq_pruned.het 
#None, freeze 3

#Create QC'ed version.  Remove monomorphic variants.
#add to "list_sample_to_remove.txt" the above problematic ID
#Freeze 3, we remove samples suffering from low coverage, as mentionned by Claudia
rm ${pathData}/list_sample_to_remove.txt
echo 'M-8305 M-8305' > ${pathData}/list_sample_to_remove.txt
echo 'M-0541 M-0541' >> ${pathData}/list_sample_to_remove.txt
echo 'M-4812 M-4812' >> ${pathData}/list_sample_to_remove.txt
echo 'M-4796 M-4796' >> ${pathData}/list_sample_to_remove.txt
echo 'M-4802 M-4802' >> ${pathData}/list_sample_to_remove.txt
echo 'M-4830 M-4830' >> ${pathData}/list_sample_to_remove.txt
echo 'M-3266 M-3266' >> ${pathData}/list_sample_to_remove.txt
plink --bfile ${pathData}/seq --remove ${pathData}/list_sample_to_remove.txt --geno 0.05 --maf 0.0000001 --make-bed --out ${pathData}/seq

#Check for sex concordance
plink --bfile ${pathData}/seq --check-sex --out ${pathData}/seq_sexCheck

#Create final version of VCFs.
#First, we remove the problematic samples identified in the PLINK QC above, variants with missingness > 0.05, variants with NON_REF as alt allele, and monomorphic variants.
#We also set the variant IDs, so that it's easy to tell which variants were split from multi-allelic ones to create bi-allelic.
plink --vcf ${pathData}/seq.vcf --remove ${pathData}/list_sample_to_remove.txt --geno 0.05 --recode vcf --out ${pathData}/seq_bial_step1
bcftools norm -m -both -O v -o ${pathData}/seq_bial_step2.vcf ${pathData}/seq_bial_step1.vcf
bcftools view -i 'ALT[0] != "NON_REF"' -c 1 -O v -o ${pathData}/seq_bial_step3.vcf ${pathData}/seq_bial_step2.vcf
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -O v -o ${pathData}/seq_bial_step4.vcf ${pathData}/seq_bial_step3.vcf
bcftools sort -O v -o ${pathData}/seq_bial_QCed.vcf ${pathData}/seq_bial_step4.vcf
bgzip ${pathData}/seq_bial_QCed.vcf
tabix -p vcf ${pathData}/seq_bial_QCed.vcf.gz

#If everything went right, to save space, remove step1 to step3 files.
rm ${pathData}/seq_bial_step*.vcf

#Now we deal with the multi-allelic file
plink --vcf ${pathData}/seqMultiallelic.vcf --remove ${pathData}/list_sample_to_remove.txt --geno 0.05 --recode vcf --out ${pathData}/seq_multial_step1
bcftools norm -m -both -O v -o ${pathData}/seq_multial_step2.vcf ${pathData}/seq_multial_step1.vcf
bcftools view -i 'ALT[0] != "NON_REF"' -c 1 -O v -o ${pathData}/seq_multial_step3.vcf ${pathData}/seq_multial_step2.vcf
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT\_m' -O v -o ${pathData}/seq_multial_step4.vcf ${pathData}/seq_multial_step3.vcf
bcftools sort -O v -o ${pathData}/seq_multial_QCed.vcf ${pathData}/seq_multial_step4.vcf
bgzip ${pathData}/seq_multial_QCed.vcf
tabix -p vcf ${pathData}/seq_multial_QCed.vcf.gz

#Again, if everything went right, to save space, remove step1 to step3 files.
rm ${pathData}/seq_multial_step*.vcf

#Merge the bi-allelic part and the multi-allelic split part.
bcftools concat -a ${pathData}/seq_bial_QCed.vcf.gz ${pathData}/seq_multial_QCed.vcf.gz -O v -o ${pathData}/seq_merged.vcf 
bcftools sort -T /scratch/jasmric/tmp_bcftools_sort -O z -o ${pathData}/seq_merged.vcf.gz -m 9G ${pathData}/seq_merged.vcf 

#Verify if variants with NON_REF as an alternate allele are still remaining (maybe a bug in bcftools).
zcat ${pathData}/seq_merged.vcf.gz | grep -v -w '<NON_REF>' > ${pathData}/seq_FINAL.vcf

#Compress and index the final VCF.
bgzip ${pathData}/seq_FINAL.vcf
tabix -p vcf ${pathData}/seq_FINAL.vcf.gz

