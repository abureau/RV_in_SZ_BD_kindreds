#Run this code using Arcturus.
#Here, we prepare the input for lg_auto.
#This time, we want a larger framework. I do this by:
#1- removing the criteria on MAF greater than 0.48 for the clumping.
#2- doing a clumping BASED on MAF, as proposed by Privé (https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html).
#To get more variants and less inconsistencies as the authors of GIGI2 propose, I set the 1-MAF threshold to 0.65 (MAF greater than 0.35) (values selected by trial and error in order to get a variant every 0.5 cM)

path=/home/jricard/gl_auto_data
path_data=${path}/corrected_pedigree
glauto_data=/home/jricard/gl_auto_data/corrected_pedigree

#Genotyping data on variants present in both GSA and OmniExpress (lifted over to GRCh38) are found here.
#We need to remove the subjects without genotypes found in the OmniExpress data set (lets use a --mind threshold of 0.9).
#We also need to modify the .fam following the detection of different errors by Jasmin and Mylene in Decembre 2023.
#We wish to keep one subject 4705 only, the one presenting less missing values (4705 4705 presents more missing values than 151 4705)
#2347 -> 2760, 2760 -> 2347, 2413 -> 2461, 650 -> 560, 2848 removed
geno=/mnt-biostats/genotypes_omni_gsa/common_var_GRCh38/omni_gsa_common_var_grch38
echo "4705 4705" > ${path_data}/remove_4705.txt
/mnt-biostats/Plink/plink --bfile ${geno} --mind 0.9 --remove ${path_data}/remove_4705.txt --allow-extra-chr --chr 1-22 --make-bed --out ${path_data}/omni_gsa_common_var
R
  fam <- data.table::fread("/home/jricard/gl_auto_data/corrected_pedigree/omni_gsa_common_var.fam")
  data.table::fwrite(fam, "/home/jricard/gl_auto_data/corrected_pedigree/omni_gsa_common_var.fam~", col.names = FALSE, row.names = FALSE, sep = "\t")
  ped_ref <- data.table::fread("/home/jricard/gl_auto_data/GCbroad.pre")
  fam_check <- merge(x = fam, y = ped_ref[,paste0("V", 1:4)], by.x = "V2", by.y = "V2", all.x = TRUE, sort = FALSE)
  fam_check[fam_check$V1.x != fam_check$V1.y, c("V1.x","V3.x","V4.x")] <- fam_check[fam_check$V1.x != fam_check$V1.y, c("V1.y","V3.y","V4.y")]
  fam <- setNames(fam_check[,c("V1.x", "V2", "V3.x","V4.x", "V5", "V6")], paste0("V", 1:6))
  old_2347 <- which(fam$V2 == "2347"); old_2760 <- which(fam$V2 == "2760"); new_2760 <- fam[old_2347,]; new_2347 <- fam[old_2760,]
  fam[old_2347,] <- new_2347; fam[old_2760,] <- new_2760
  fam[which(fam$V2 == "1002413"), "V2"] <- 2461
  fam[which(fam$V2 == "650"), "V1"] <- 250; fam[which(fam$V2 == "650"), "V3"] <- 540; fam[which(fam$V2 == "650"), "V4"] <- 541; fam[which(fam$V2 == "650"), "V5"] <- 2; fam[which(fam$V2 == "650"), "V2"] <- 560
  data.table::fwrite(fam, "/home/jricard/gl_auto_data/corrected_pedigree/omni_gsa_common_var.fam", col.names = FALSE, row.names = FALSE, sep = "\t")
quit()
n

#Add the missing subject in the existing list of related to keep in the dataset
#egrep " 2347| 2760| 2413| 2461| 650| 560| 2848" ${path}/related_to_keep.txt
egrep -v " 650| 2848" ${path}/related_to_keep.txt > ${path}/2023-12-19_related_to_keep.txt
echo -e "119 2347\n122 2461\n250 560" >> ${path}/2023-12-19_related_to_keep.txt
comm -1 -3 <(sed 's/ /\t/g' ${path}/2023-12-19_related_to_keep.txt | cut -f-2 | uniq | sort) <(cut -f-2 ${path_data}/omni_gsa_common_var.fam | uniq | sort) > ${path}/2023-12-19_non_related_to_remove.txt
/mnt-biostats/Plink/plink --bfile ${path_data}/omni_gsa_common_var --remove ${path}/2023-12-19_non_related_to_remove.txt --make-bed --out ${path_data}/omni_gsa_common_var_related

#Set mendelien errors to missing.
/mnt-biostats/Plink/plink --bfile ${path_data}/omni_gsa_common_var_related --set-me-missing --mendel-duos --mendel-multigen --make-bed --out ${path_data}/omni_gsa_common_var_mendel

#Break the dataset by chromosome.
for chr in {1..22}; do /mnt-biostats/Plink/plink --bfile ${path_data}/omni_gsa_common_var_mendel --chr ${chr} --make-bed --out ${path_data}/omni_gsa_common_var_mendel_chr${chr}; done
rm ${path_data}/*.nosex

variant_number=0
for chr in {1..22}
do
  /mnt-biostats/Plink/plink --bfile ${path_data}/omni_gsa_common_var_mendel_chr${chr} --freq --nonfounders --out ${path_data}/clumping_freq
  awk 'BEGIN {OFS="\t"} NR>1 {$5 = 1 - $5} 5' ${path_data}/clumping_freq.frq > ${path_data}/clumping_freq_formated.frq
  /mnt-biostats/Plink/plink --bfile ${path_data}/omni_gsa_common_var_mendel_chr${chr} --clump ${path_data}/clumping_freq_formated.frq --clump-p1 0.65 --clump-p2 0.65 --clump-r2 0.1 --clump-kb 300 --clump-snp-field SNP --clump-field MAF --out ${path_data}/clumps_chr${chr} 
  head -n 1 ${path_data}/omni_gsa_common_var_mendel_chr${chr}.bim | cut -f 2 > ${path_data}/clumps_plus_extremite_to_keep.snplist
  awk 'NR>1 {print $3}' ${path_data}/clumps_chr${chr}.clumped | cat >> ${path_data}/clumps_plus_extremite_to_keep.snplist
  tail -n 1 ${path_data}/omni_gsa_common_var_mendel_chr${chr}.bim | cut -f 2 >> ${path_data}/clumps_plus_extremite_to_keep.snplist
  /mnt-biostats/Plink/plink --bfile ${path_data}/omni_gsa_common_var_mendel_chr${chr} --extract ${path_data}/clumps_plus_extremite_to_keep.snplist --freq --make-bed --recode --out ${path_data}/omni_gsa_common_var_mendel_clumped_chr${chr}
  variant_number_chr=$(wc -l ${path_data}/omni_gsa_common_var_mendel_clumped_chr${chr}.bim | awk '{print $1}'); variant_number=$(expr ${variant_number} + ${variant_number_chr}) 
done
rm ${path_data}/clumping_freq*; rm ${path_data}/*.nosex
echo $variant_number
#21,843 SNPs remaining

Rscript fichier_morgan_SNPs.R

for chr in {1..22}
do
  cd ${glauto_data}/chr${chr}
  nohup /home/jricard/2dossier-biostats/sboies/MORGAN/MORGAN_V34_Release/Autozyg/gl_auto glauto.par > glauto_console_output.txt &
done

#It is possible that some mendel errors remain in the datasets because PLINK struggles to set to missing mendel errors with non-genotyped subjects.
#gl_auto, with the help of the whole pedigree used in its input, detect them and stops
#Hence, gl_auto may stop on certain chromosome.

#The only way I found do correct this easily is by hand... Look at the log file and set to missing the error.
for chr in {1..22}; do echo ${chr}; grep "Mendelian inconsistency" ${glauto_data}/chr${chr}/glauto_console_output.txt; done

#chr5 - 220_4724 - marker 878
grep "220_4724" ${glauto_data}/chr5/framework_gcbroad.map_af_geno | cut -d " " -f 1756,1757
mv ${glauto_data}/chr5/framework_gcbroad.map_af_geno ${glauto_data}/chr5/framework_gcbroad.map_af_geno_error
awk '/220_4724/ {$1756=0; $1757=0} 1' ${glauto_data}/chr5/framework_gcbroad.map_af_geno_error > ${glauto_data}/chr5/framework_gcbroad.map_af_geno

#chr12 - 235_4730, 235_4729, 235_4728 - marker 443
grep "235_4730" ${glauto_data}/chr12/framework_gcbroad.map_af_geno | cut -d " " -f 886,887
grep "235_4729" ${glauto_data}/chr12/framework_gcbroad.map_af_geno | cut -d " " -f 886,887
grep "235_4728" ${glauto_data}/chr12/framework_gcbroad.map_af_geno | cut -d " " -f 886,887
mv ${glauto_data}/chr12/framework_gcbroad.map_af_geno ${glauto_data}/chr12/framework_gcbroad.map_af_geno_error
awk '/235_4730/ {$886=0; $887=0} 1' ${glauto_data}/chr12/framework_gcbroad.map_af_geno_error > ${glauto_data}/chr12/framework_gcbroad.map_af_geno_tmp1
awk '/235_4729/ {$886=0; $887=0} 1' ${glauto_data}/chr12/framework_gcbroad.map_af_geno_tmp1 > ${glauto_data}/chr12/framework_gcbroad.map_af_geno_tmp2
awk '/235_4728/ {$886=0; $887=0} 1' ${glauto_data}/chr12/framework_gcbroad.map_af_geno_tmp2 > ${glauto_data}/chr12/framework_gcbroad.map_af_geno
rm ${glauto_data}/chr12/framework_gcbroad.map_af_geno_tmp*

#chr14 - 220_4724 - marker 396
grep "220_4724" ${glauto_data}/chr14/framework_gcbroad.map_af_geno | cut -d " " -f 792,793
mv ${glauto_data}/chr14/framework_gcbroad.map_af_geno ${glauto_data}/chr14/framework_gcbroad.map_af_geno_error
awk '/220_4724/ {$792=0; $793=0} 1' ${glauto_data}/chr14/framework_gcbroad.map_af_geno_error > ${glauto_data}/chr14/framework_gcbroad.map_af_geno

#chr18 - 228_779 - marker 251
grep "228_4715" ${glauto_data}/chr18/framework_gcbroad.map_af_geno | cut -d " " -f 502,503
grep "228_4716" ${glauto_data}/chr18/framework_gcbroad.map_af_geno | cut -d " " -f 502,503
#chr18 - 235_4730, 235_4729, 235_4728 - marker 251
grep "235_4730" ${glauto_data}/chr18/framework_gcbroad.map_af_geno | cut -d " " -f 502,503
grep "235_4729" ${glauto_data}/chr18/framework_gcbroad.map_af_geno | cut -d " " -f 502,503
grep "235_4728" ${glauto_data}/chr18/framework_gcbroad.map_af_geno | cut -d " " -f 502,503
mv ${glauto_data}/chr18/framework_gcbroad.map_af_geno ${glauto_data}/chr18/framework_gcbroad.map_af_geno_error
awk '/228_4715/ {$502=0; $503=0} 1' ${glauto_data}/chr18/framework_gcbroad.map_af_geno_error > ${glauto_data}/chr18/framework_gcbroad.map_af_geno_tmp1
awk '/228_4716/ {$502=0; $503=0} 1' ${glauto_data}/chr18/framework_gcbroad.map_af_geno_tmp1 > ${glauto_data}/chr18/framework_gcbroad.map_af_geno_tmp2
awk '/235_4730/ {$502=0; $503=0} 1' ${glauto_data}/chr18/framework_gcbroad.map_af_geno_tmp2 > ${glauto_data}/chr18/framework_gcbroad.map_af_geno_tmp3
awk '/235_4729/ {$502=0; $503=0} 1' ${glauto_data}/chr18/framework_gcbroad.map_af_geno_tmp3 > ${glauto_data}/chr18/framework_gcbroad.map_af_geno_tmp4
awk '/235_4728/ {$502=0; $503=0} 1' ${glauto_data}/chr18/framework_gcbroad.map_af_geno_tmp4 > ${glauto_data}/chr18/framework_gcbroad.map_af_geno
rm ${glauto_data}/chr18/framework_gcbroad.map_af_geno_tmp*

#Run it again
for chr in 5 12 14 18
do
  cd ${glauto_data}/chr${chr}
  nohup /home/jricard/2dossier-biostats/sboies/MORGAN/MORGAN_V34_Release/Autozyg/gl_auto glauto.par > glauto_console_output.txt &
done