#Run this code using Arcturus.
#Here, we prepare the input for lg_auto.
#This time, we want a larger framework. I do this by:
#1- removing the critera on MAF greater than 0.48 for the clumping.
#2- doing a clumping BASED on MAF, as proposed by Privé (https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html).
#To get more variants and less inconsistencies as the authors of GIGI2 propose, I set the 1-MAF threshold to 0.65 (MAF greater than 0.35) (values selected by trial and error in order to get a variant every 0.5 cM)

pathData=/home/jricard/lg_auto_data/larger_framework
variant_number=0
for chr in {1..22}
do
  /mnt-biostats/Plink/plink --bfile ${pathData}/genotypes-chr${chr}-Commun --freq --out ${pathData}/clumping_freq
  awk 'BEGIN {OFS="\t"} NR>1 {$5 = 1 - $5} 5' ${pathData}/clumping_freq.frq > ${pathData}/clumping_freq_formated.frq
  /mnt-biostats/Plink/plink --bfile ${pathData}/genotypes-chr${chr}-Commun --remove ${pathData}/non_related_to_remove.txt --clump ${pathData}/clumping_freq_formated.frq --clump-p1 0.65 --clump-p2 0.65 --clump-r2 0.1 --clump-kb 250 --clump-snp-field SNP --clump-field MAF --out ${pathData}/GSA_genotypes_clumped_chr${chr} 
  head -n 1 ${pathData}/genotypes-chr${chr}-Commun.bim | cut -f 2 > ${pathData}/clumped_plus_extremite_to_keep.snplist
  awk 'NR>1 {print $3}' ${pathData}/GSA_genotypes_clumped_chr${chr}.clumped | cat >> ${pathData}/clumped_plus_extremite_to_keep.snplist
  tail -n 1 ${pathData}/genotypes-chr${chr}-Commun.bim | cut -f 2 >> ${pathData}/clumped_plus_extremite_to_keep.snplist
  /mnt-biostats/Plink/plink --bfile ${pathData}/genotypes-chr${chr}-Commun --remove ${pathData}/non_related_to_remove.txt --extract ${pathData}/clumped_plus_extremite_to_keep.snplist --make-bed --recode --out ${pathData}/GSA_genotypes_clumped_chr${chr}
  /mnt-biostats/Plink/plink --bfile ${pathData}/GSA_genotypes_clumped_chr${chr} --freq --out ${pathData}/GSA_genotypes_clumped_chr${chr}
  variant_number_chr=$(wc -l ${pathData}/GSA_genotypes_clumped_chr${chr}.bim | awk '{print $1}')
  variant_number=$(expr ${variant_number} + ${variant_number_chr}) 
  #Check for mendelien errors
  /mnt-biostats/Plink/plink --bfile ${pathData}/GSA_genotypes_clumped_chr${chr} --set-me-missing --mendel-duos --make-bed --out ${pathData}/GSA_genotypes_clumped_mendel_chr${chr}
  /mnt-biostats/Plink/plink --bfile ${pathData}/GSA_genotypes_clumped_mendel_chr${chr} --recode --out ${pathData}/GSA_genotypes_clumped_mendel_chr${chr}
done
rm ${pathData}/clumping_freq*
rm ${pathData}/*.nosex
echo $variant_number
#21,422 SNPs remaining

Rscript fichier_morgan_SNPs.R

for chr in {1..22}
do
cd /home/jricard/lg_auto_data/larger_framework/chr${chr}
nohup /home/jricard/2dossier-biostats/sboies/MORGAN/MORGAN_V34_Release/Autozyg/gl_auto glauto.par > glauto_console_output.txt &
done
