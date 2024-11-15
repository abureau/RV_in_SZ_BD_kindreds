#To run GIGI2 by family
for chr in {1..22}
do
  path_impu=${path}/imputation_fam/chr${chr}

  #gl_auto produced a check.oped where the pedigree starts on the 7th line. Pedigree must start on the 5th line if not, results are awkward.
  #We remove the 3rd ("input pedigree record father mother") and the 4th lines ("input pedigree record gender present") which aren't present in the GIGI2 example.
  egrep -v "input pedigree record father mother|input pedigree record gender present" ${path_impu}/check.oped > ${path_impu}/check.oped.cor

  #Create a temporary file where GIGI2 will output its results for every family.
  mkdir ${path_impu}/tmp_fam
  path_impu_tmp=${path_impu}/tmp_fam

  #We iterate on a list of every unique family ID found in the data file.
  #To create a merged file, we set an index in our loop
  first_fam=1
  head -n 1 ${path_impu}/dense_geno_long.txt | tr ' ' '\n' | cut -d "_" -f1 | tail -n +2 | uniq | sort | while read fam
  do
    echo "******************Imputation of fam ${fam}******************"

    #We input the current family information only to GIGI2.
    head -n 4 ${path_impu}/check.oped.cor > ${path_impu_tmp}/check.oped.cor_fam${fam} #Get the header of the check.oped.cor
    fam_n=$(grep "${fam}_" ${path_impu}/check.oped.cor | wc -l)
    grep "${fam}_" ${path_impu}/check.oped.cor >> ${path_impu_tmp}/check.oped.cor_fam${fam} #Paste the pedigree of the current family in the check.oped.cor
    sed -i 's/input pedigree size 1557/input pedigree size '"$fam_n"'/' ${path_impu_tmp}/check.oped.cor_fam${fam} #In the header, correct the number of subject in the family
    grep "${fam}_" ${path_impu}/framework.IVs > ${path_impu_tmp}/framework.IVs_fam${fam} #Keep our current family IVs only
    fam_idx=$(head -n 1 ${path_impu}/dense_geno_long.txt | tr ' ' '\n' | egrep -n "${fam}_|id" | cut -d: -f1) #Get the index of the columns where we find the current family's genotypes.
    cut -d " " -f $(echo $fam_idx | tr ' ' ',') ${path_impu}/dense_geno_long.txt > ${path_impu_tmp}/dense_geno_long.txt_fam${fam} #Keep our current family's genotypes only.

    #Create the parameters file
    echo "--ped ${path_impu_tmp}/check.oped.cor_fam${fam}" > ${path_impu_tmp}/param_gigi2.txt
    echo "--meiosis ${path_impu_tmp}/framework.IVs_fam${fam}" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--iter 1000" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--genocall 2 0.8 0.9" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--smap ${path_impu}/framework_pos.txt" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--dmap ${path_impu}/dense_pos_vf.txt" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--geno ${path_impu_tmp}/dense_geno_long.txt_fam${fam}" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--afreq ${path_impu}/dense_freq.txt" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--out ${path_impu_tmp}/GIGI2_imputed_chr${chr}_fam${fam}" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--seed 1234" >> ${path_impu_tmp}/param_gigi2.txt

    #Imputation using GIGI2
    /lustre03/project/6033529/SOFT/GIGI2/src/GIGI2 ${path_impu_tmp}/param_gigi2.txt

    #Get every subject ID from the pedigree, duplicate the values, add a colname for the variants ID and add this whole header to the imputed genotypes file.
    mv ${path_impu_tmp}/GIGI2_imputed_chr${chr}_fam${fam}.geno ${path_impu_tmp}/GIGI2_imputed_chr${chr}_fam${fam}_no_ID.geno
    cat <(sed 's/^[ \t]*//' ${path_impu_tmp}/check.oped.cor_fam${fam} | tail -n +5 | cut -d " " -f1 | sed 'p' | paste -s -d' ' | sed 's/^/id /') ${path_impu_tmp}/GIGI2_imputed_chr${chr}_fam${fam}_no_ID.geno > ${path_impu_tmp}/GIGI2_imputed_chr${chr}_fam${fam}.geno
    rm ${path_impu_tmp}/GIGI2_imputed_chr${chr}_fam${fam}_no_ID.geno

    #To merge the result...
    if [ $first_fam -eq 1 ]; then
        cp ${path_impu_tmp}/GIGI2_imputed_chr${chr}_fam${fam}.geno ${path_impu}/GIGI2_imputed_chr${chr}_by_fam.geno
        cp ${path_impu_tmp}/GIGI2_imputed_chr${chr}_fam${fam}.prob ${path_impu}/GIGI2_imputed_chr${chr}_by_fam.prob
        ((first_fam++))
    else
        mv ${path_impu}/GIGI2_imputed_chr${chr}_by_fam.geno ${path_impu}/GIGI2_imputed_chr${chr}_by_fam_tmp.geno
        paste -d ' ' ${path_impu}/GIGI2_imputed_chr${chr}_by_fam_tmp.geno <(cut -f2- -d' ' ${path_impu_tmp}/GIGI2_imputed_chr${chr}_fam${fam}.geno) > ${path_impu}/GIGI2_imputed_chr${chr}_by_fam.geno
        rm ${path_impu}/GIGI2_imputed_chr${chr}_by_fam_tmp.geno
        mv ${path_impu}/GIGI2_imputed_chr${chr}_by_fam.prob ${path_impu}/GIGI2_imputed_chr${chr}_by_fam_tmp.prob
        paste -d ' ' ${path_impu}/GIGI2_imputed_chr${chr}_by_fam_tmp.prob <(cut -f2- -d' ' ${path_impu_tmp}/GIGI2_imputed_chr${chr}_fam${fam}.prob) > ${path_impu}/GIGI2_imputed_chr${chr}_by_fam.prob
        rm ${path_impu}/GIGI2_imputed_chr${chr}_by_fam_tmp.prob
    fi
  done
done

for chr in {1..22}
do
  #The pedigree is turned into a .tfam file for PLINK
  awk '{print "0", $0}' <(sed '1,6d' ${path_impu}/check.oped) > ${path_impu}/geno_impute_chr${chr}.tfam

  #The imputed genotypes are turned into a .tped file for PLINK
  sed 's/ / 0 0 /2' <(awk '{print '"${chr}"', $0}' <(paste --delimiters=' ' ${path_impu}/liste_dense.txt <(cut -d' ' -f2- <(tail +2 ${path_impu}/GIGI2_imputed_chr${chr}_by_fam.geno)))) > ${path_impu}/geno_impute_chr${chr}_step1.tped
  awk 'BEGIN {FS=OFS=" "} {split($2, parts, ":"); $4 = parts[2]; print}' ${path_impu}/geno_impute_chr${chr}_step1.tped > ${path_impu}/geno_impute_chr${chr}.tped
  rm ${path_impu}/geno_impute_chr${chr}_step1.tped; gzip ${path_impu}/GIGI2_imputed_chr${chr}_by_fam.geno

  Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples/call/post_imputation_tped_to_full_missing.R $chr
  
  #Generating the .ped and .fam files for the imputed data
  plink --tped ${path_impu}/geno_impute_chr${chr}_full_missing.tped --tfam ${path_impu}/geno_impute_chr${chr}.tfam --recode --missing --out ${path_impu}/geno_impute_chr${chr}
 
  #Back to A1 and A2 
  #Create the file to update the alleles
  awk '{print $2, 1, 2, $5, $6}' ${path_impu}/seq_FINAL_PHASED_with_mask_in_genmap_chr${chr}.bim > ${path_impu}/update_allele.txt
 
  #Generating the .ped and .fam files for the imputed data
  plink --file ${path_impu}/geno_impute_chr${chr} --update-alleles ${path_impu}/update_allele.txt --recode vcf --out ${path_impu}/geno_impute_chr${chr}_ATCG
  bgzip ${path_impu}/geno_impute_chr${chr}_ATCG.vcf
 
  #The vcf.gz samples names were merged, so they all present a "0_", we remove it.
  paste <(bcftools query -l ${path_impu}/geno_impute_chr${chr}_ATCG.vcf.gz) <(bcftools query -l ${path_impu}/geno_impute_chr${chr}_ATCG.vcf.gz | sed "s/^0_//g") > ${path_impu}/header_change.txt
  bcftools reheader -o ${path_impu}/geno_impute_chr${chr}_ATCG_cor.vcf.gz --sample ${path_impu}/header_change.txt ${path_impu}/geno_impute_chr${chr}_ATCG.vcf.gz
  rm ${path_impu}/geno_impute_chr${chr}_ATCG.vcf.gz; mv ${path_impu}/geno_impute_chr${chr}_ATCG_cor.vcf.gz ${path_impu}/geno_impute_chr${chr}_ATCG.vcf.gz
  bcftools index -t -f ${path_impu}/geno_impute_chr${chr}_ATCG.vcf.gz
done


#We want to know which positions are imputed with IMPUTE5 and not with GIGI2.
for chr in {1..22};do bcftools query -f '%ID\n' ${path}/imputation_pop/impute5_imputed_chr${chr}.vcf.gz >> ${path}/imputation_pop/imputed_positions.txt ;done
for chr in {1..22};do bcftools query -f '%ID\n' ${path}/imputation_fam/chr${chr}/geno_impute_chr${chr}_ATCG.vcf.gz >> ${path}/imputation_fam/imputed_positions.txt ;done
diff ${path}/imputation_pop/imputed_positions.txt ${path}/imputation_fam/imputed_positions.txt | grep "<" | sed "s/< //g" > ${path}/imputation_pop/imputed_popositions_only_impute5.txt
mkdir ${path}/imputation_pop/maf_positions_only_impute5
for chr in {1..22}
do
  plink --vcf ${path}/imputation_pop/impute5_imputed_chr${chr}.vcf.gz --extract ${path}/imputation_pop/imputed_popositions_only_impute5.txt --keep <(sed "s/_/ /g" ${path}/imputation_comb/impute5_common_samples.txt) --freq --nonfounders --out ${path}/imputation_pop/maf_positions_only_impute5/impute5_imputed_chr${chr}
done
rm ${path}/imputation_pop/maf_positions_only_impute5/impute5_imputed_chr*.nosex

#GIGI2 results take a lot of place on the disk, we need to make some cleaning
for chr in {1..22}
do
  #geno_impute_chr${chr}_full_missing.tped are the same as /geno_impute_chr${chr}.tped
  rm ${path}/imputation_fam/chr${chr}/geno_impute_chr${chr}_full_missing.tped
done

#This task takes a long time. We recommand doing it in parallel.
for chr in {1..22}
do
  #Archive the temporary files
  tar -zcvf ${path}/imputation_fam/chr${chr}/tmp_fam.tar.gz ${path}/imputation_fam/chr${chr}/tmp_fam &
done
wait