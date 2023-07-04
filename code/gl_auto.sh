#Prepare the bfiles for gl_auto.
#Needs to be executed on Arcturus.
module load r/4.2 StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16

pathData=/home/jricard/gl_auto_data
variant_number=0
for chr in {1..22}
do
  #Clump the data using MAF to create a bfile of quasi-independant variants.
  #Then, the first and the last genotyped variants are added to the bfile.
  /mnt-biostats/Plink/plink --bfile ${pathData}/genotypes-chr${chr}-Commun --freq --out ${pathData}/clumping_freq
  /mnt-biostats/Plink/plink --bfile ${pathData}/genotypes-chr${chr}-Commun --remove ${pathData}/non_related_to_remove.txt --write-snplist --maf 0.48 --clump ${pathData}/clumping_freq.frq --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --clump-snp-field SNP --clump-field MAF --out ${pathData}/GSA_genotypes_clumped_chr${chr} 
  head -n 1 ${pathData}/genotypes-chr${chr}-Commun.bim | cut -f 2 > ${pathData}/clumped_plus_extremite_to_keep.snplist
  cat ${pathData}/GSA_genotypes_clumped_chr${chr}.snplist >> ${pathData}/clumped_plus_extremite_to_keep.snplist
  tail -n 1 ${pathData}/genotypes-chr${chr}-Commun.bim | cut -f 2 >> ${pathData}/clumped_plus_extremite_to_keep.snplist
  /mnt-biostats/Plink/plink --bfile ${pathData}/genotypes-chr${chr}-Commun --remove ${pathData}/non_related_to_remove.txt --extract ${pathData}/clumped_plus_extremite_to_keep.snplist --make-bed --recode --out ${pathData}/GSA_genotypes_clumped_chr${chr}
  /mnt-biostats/Plink/plink --bfile ${pathData}/GSA_genotypes_clumped_chr${chr} --freq --out ${pathData}/GSA_genotypes_clumped_chr${chr}

  #Output the number of variants remaining.
  variant_number_chr=$(wc -l ${pathData}/GSA_genotypes_clumped_chr${chr}.bim | awk '{print $1}')
  variant_number=$(expr ${variant_number} + ${variant_number_chr}) 

  #Check for mendelien errors.
  #Create a bfile setting missing values to mendelien errors.
  /mnt-biostats/Plink/plink --bfile ${pathData}/GSA_genotypes_clumped_chr${chr} --set-me-missing --mendel-duos --make-bed --out ${pathData}/GSA_genotypes_clumped_mendel_chr${chr}
  /mnt-biostats/Plink/plink --bfile ${pathData}/GSA_genotypes_clumped_mendel_chr${chr} --recode --out ${pathData}/GSA_genotypes_clumped_mendel_chr${chr}
done
echo $variant_number
rm ${pathData}/clumping_freq*
rm ${pathData}/*.nosex

Rscript fichier_morgan_SNPs.R

#Then, we run gl_auto.
for chr in {1..22}
do
cd /home/jricard/gl_auto_data/chr${chr}
nohup /home/jricard/2dossier-biostats/sboies/MORGAN/MORGAN_V34_Release/Autozyg/gl_auto glauto.par > glauto_console_output.txt &
done
