#Do not forget to import gl_auto data from Arcturus to Beluga.
#They are found in their respective chromosome folder in /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/imputation

module load r/4.2 StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16
path=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples

#Run the R code `imputation_data.R` to prepare the data for GIGI2.
for chr in {1..6}; do cd ${path}/imputation/chr${chr}; sbatch --time=6:00:00 --mem=100G --export=chr=${chr} /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/imputation_data_split_call.sh; done
for chr in {7..12}; do cd ${path}/imputation/chr${chr}; sbatch --time=04:00:00 --mem=100G --export=chr=${chr} /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/imputation_data_call.sh; done
for chr in {13..22}; do cd ${path}/imputation/chr${chr}; sbatch --time=02:00:00 --mem=50G --export=chr=${chr} /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/imputation_data_call.sh; done

#To run the following code using sbatch, run this line
#sbatch --export=path=${path} /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/imputation_GIGI2.sh

for chr in {22..1}
do
  path_impu=${path}/imputation/chr${chr}

  #gl_auto produced a check.oped where the pedigree starts on the 7th line. Pedigree must start on the 5th line if not, results are awkward.
  #We remove the 3rd ("input pedigree record father mother") and the 4th lines ("input pedigree record gender present") which aren't present in the GIGI2 example.
  egrep -v "input pedigree record father mother|input pedigree record gender present" ${path_impu}/check.oped > ${path_impu}/check.oped.cor

  #Genotype file from wide to long
  if [ ! -e "${path_impu}/dense_geno_long.txt" ]; then
    /lustre03/project/6033529/SOFT/GIGI2/src/Wide2Long ${path_impu}/dense_geno.txt ${path_impu}/dense_geno_long.txt
  fi

  #Create the parameters file
  echo "--ped ${path_impu}/check.oped.cor" > ${path_impu}/param_gigi2.txt
  echo "--meiosis ${path_impu}/framework.IVs" >> ${path_impu}/param_gigi2.txt
  echo "--iter 1000" >> ${path_impu}/param_gigi2.txt
  echo "--genocall 2 0.8 0.9" >> ${path_impu}/param_gigi2.txt
  echo "--smap ${path_impu}/framework_pos.txt" >> ${path_impu}/param_gigi2.txt
  echo "--dmap ${path_impu}/dense_pos_vf.txt" >> ${path_impu}/param_gigi2.txt
  echo "--geno ${path_impu}/dense_geno_long.txt" >> ${path_impu}/param_gigi2.txt
  echo "--afreq ${path_impu}/dense_freq.txt" >> ${path_impu}/param_gigi2.txt
  echo "--out ${path_impu}/GIGI2_imputed_chr${chr}" >> ${path_impu}/param_gigi2.txt
  echo "--seed 1234" >> ${path_impu}/param_gigi2.txt

  #Imputation using GIGI2
  /lustre03/project/6033529/SOFT/GIGI2/src/GIGI2 ${path_impu}/param_gigi2.txt
done

#To run the following code using sbatch, run this line
#sbatch --export=path=${path} /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/imputation_GIGI2_by_fam.sh

#To run GIGI2 by family
for chr in {1..22}
do
  path_impu=${path}/imputation/chr${chr}

  #gl_auto produced a check.oped where the pedigree starts on the 7th line. Pedigree must start on the 5th line if not, results are awkward.
  #We remove the 3rd ("input pedigree record father mother") and the 4th lines ("input pedigree record gender present") which aren't present in the GIGI2 example.
  egrep -v "input pedigree record father mother|input pedigree record gender present" ${path_impu}/check.oped > ${path_impu}/check.oped.cor

  #Create a temporary file where GIGI2 will output its results for every family.
  mkdir ${path_impu}/tmp_fam
  path_impu_tmp=${path_impu}/tmp_fam

  #We iterate on a list of every unique family ID found in the data file.
  #To create a merged file, we set an index in our loop
  first_fam=1
  head -n 1 ${path_impu}/dense_geno_long.txt | tr ' ' '\n' | cut -d "_" -f1 | uniq | sort | while read fam
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
        ((first_fam++))
    else
        mv ${path_impu}/GIGI2_imputed_chr${chr}_by_fam.geno ${path_impu}/GIGI2_imputed_chr${chr}_by_fam_tmp.geno
        paste -d ' ' ${path_impu}/GIGI2_imputed_chr${chr}_by_fam_tmp.geno <(cut -f2- -d' ' ${path_impu_tmp}/GIGI2_imputed_chr${chr}_fam${fam}.geno) > ${path_impu}/GIGI2_imputed_chr${chr}_by_fam.geno
        rm ${path_impu}/GIGI2_imputed_chr${chr}_by_fam_tmp.geno
    fi
  done
done