#Do not forget to import gl_auto data from Arcturus to Beluga.
#They are found in their respective chromosome folder in /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/imputation

module load r/4.2 StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16

#Run the R code `imputation_data.R` to prepare the data for GIGI2.
for chr in {1..21}; do sbatch --time=04:00:00 --mem=50G --export=chr=${chr} /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/imputation_data_call.sh; done

path=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples
for chr in {1..22}
do
  path_impu=${path}/imputation/chr${chr}

  #gl_auto produced a check.oped where the pedigree starts on the 7th line. Pedigree must start on the 5th line if not, results are awkward.
  #We remove the 3rd ("input pedigree record father mother") and the 4th lines ("input pedigree record gender present") which aren't present in the GIGI2 example.
  egrep -v "input pedigree record father mother|input pedigree record gender present" ${path_impu}/check.oped > ${path_impu}/check.oped.cor

  #Genotype file from wide to long
  /lustre03/project/6033529/SOFT/GIGI2/src/Wide2Long ${path_impu}/dense_geno.txt ${path_impu}/dense_geno_long.txt

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

#To run GIGI2 by family
for chr in {1..22}
do
  path_impu=${path}/imputation/chr${chr}

  #gl_auto produced a check.oped where the pedigree starts on the 7th line. Pedigree must start on the 5th line if not, results are awkward.
  #We remove the 3rd ("input pedigree record father mother") and the 4th lines ("input pedigree record gender present") which aren't present in the GIGI2 example.
  egrep -v "input pedigree record father mother|input pedigree record gender present" ${path_impu}/check.oped > ${path_impu}/check.oped.cor

  mkdir ${path_impu}/tmp_fam
  tail -n +5 /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/imputation/chr22/check.oped.cor | awk -F'_' '{gsub(/[^0-9]/, "", $1); print $1}' | uniq | sort | while read fam
  do
    echo "******************Imputation of fam ${fam}******************"
    path_impu_tmp=${path_impu}/tmp_fam

    #Empty the temporary directory
    #rm ${path_impu_tmp}/*

    #We input the current family information only to GIGI2.
    head -n 4 ${path_impu}/check.oped.cor > ${path_impu_tmp}/check.oped.cor
    fam_n=$(grep "${fam}_" ${path_impu}/check.oped.cor | wc -l)
    grep "${fam}_" ${path_impu}/check.oped.cor >> ${path_impu_tmp}/check.oped.cor
    sed -i 's/input pedigree size 1557/input pedigree size '"$fam_n"'/' ${path_impu_tmp}/check.oped.cor > ${path_impu_tmp}/check.oped.cor
    grep "${fam}_" ${path_impu}/framework.IVs > ${path_impu_tmp}/framework.IVs
    grep "${fam}_" ${path_impu}/dense_geno.txt > ${path_impu_tmp}/dense_geno.txt

    #Genotype file from wide to long
    /lustre03/project/6033529/SOFT/GIGI2/src/Wide2Long ${path_impu_tmp}/dense_geno.txt ${path_impu_tmp}/dense_geno_long.txt
    rm ${path_impu_tmp}/dense_geno.txt

    #Create the parameters file
    echo "--ped ${path_impu_tmp}/check.oped.cor" > ${path_impu_tmp}/param_gigi2.txt
    echo "--meiosis ${path_impu_tmp}/framework.IVs" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--iter 1000" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--genocall 2 0.8 0.9" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--smap ${path_impu}/framework_pos.txt" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--dmap ${path_impu}/dense_pos_vf.txt" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--geno ${path_impu_tmp}/dense_geno_long.txt" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--afreq ${path_impu}/dense_freq.txt" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--out ${path_impu_tmp}/GIGI2_imputed_chr${chr}_fam${fam}" >> ${path_impu_tmp}/param_gigi2.txt
    echo "--seed 1234" >> ${path_impu_tmp}/param_gigi2.txt

    #Imputation using GIGI2
    /lustre03/project/6033529/SOFT/GIGI2/src/GIGI2 ${path_impu_tmp}/param_gigi2.txt
  done
done