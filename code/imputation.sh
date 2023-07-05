#Do not forget to import gl_auto data from Arcturus to Beluga.
#They are found in their respective chromosome folder in /lustre03/project/6033529/quebec_10x/data/freeze/imputation
module load r/4.2 StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16

#Split the QCed data by chromosome
pathData=/lustre03/project/6033529/quebec_10x/data/freeze/QC
for chr in {1..22}
do
plink --vcf ${pathData}/seq_FINAL_with_mask.vcf.gz --chr ${chr} --recode --out ${pathData}/seq_FINAL_${chr}_with_mask
done

Rscript imputation_data.R

for chr in {1..22}
do
pathImpu=/lustre03/project/6033529/quebec_10x/data/freeze/imputation/chr${chr}

#Genotype file from wide to long
/lustre03/project/6033529/SOFT/GIGI2/src/Wide2Long ${pathImpu}/dense_geno.txt ${pathImpu}/dense_geno_long.txt

#Create the parameters file
echo "--ped ${pathImpu}/GCbroad_v2.ped" > ${pathImpu}/param_gigi2.txt
echo "--meiosis ${pathImpu}/framework.IVs" >> ${pathImpu}/param_gigi2.txt
echo "--iter 1000" >> ${pathImpu}/param_gigi2.txt
echo "--genocall 2 0.8 0.9" >> ${pathImpu}/param_gigi2.txt
echo "--smap ${pathImpu}/framework_pos.txt" >> ${pathImpu}/param_gigi2.txt
echo "--dmap ${pathImpu}/dense_pos_vf.txt" >> ${pathImpu}/param_gigi2.txt
echo "--geno ${pathImpu}/dense_geno_long.txt" >> ${pathImpu}/param_gigi2.txt
echo "--afreq ${pathImpu}/dense_freq.txt" >> ${pathImpu}/param_gigi2.txt
echo "--out ${pathImpu}/GIGI2_imputed_chr${chr}" >> ${pathImpu}/param_gigi2.txt
echo "--seed 1234" >> ${pathImpu}/param_gigi2.txt

#Imputation using GIGI2
/lustre03/project/6033529/SOFT/GIGI2/src/GIGI2 ${pathImpu}/param_gigi2.txt
done