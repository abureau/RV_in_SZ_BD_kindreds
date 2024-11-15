module load StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16 r/4.2
path_impu=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/imputation_comb
path_rvs=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS
mkdir ${path_rvs}

#The first half of this script is coded in R
R
library(data.table); library(dplyr)

#From our combined imputation files, keep the GnomaAD/CaG rare variants only
#Remove monomorphic variants
path_impu <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/imputation_comb"
path_rvs <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS"
RV <- fread("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RV/cag_gnomad/RV_CaG_gnomad.snplist")
if (!file.exists(paste0(path_impu, "/merged_with_seq/impute5_gigi2_combined_seq.bim"))){
  for (chr in 1:22){
    system(paste0("bcftools query -f '%CHROM %POS %REF %ALT %ID\n' ", path_impu, "/merged_with_seq/impute5_gigi2_combined_seq_chr", chr, ".vcf.gz >> ", path_impu, "/merged_with_seq/impute5_gigi2_combined_seq.bim"))
  }
}
impu <- fread(paste0(path_impu, "/merged_with_seq/impute5_gigi2_combined_seq.bim"), header = FALSE) 
RV$V1 <- as.numeric(gsub("chr", "", RV$V1))
impu$V1 <- as.numeric(gsub("chr", "", impu$V1)); impu$V2 <- as.numeric(impu$V2)
match <- lassosum:::matchpos(tomatch = RV, ref.df = impu, auto.detect.tomatch = F, auto.detect.ref = F, chr = "V1", ref.chr = "V1", pos = "V2", ref.pos = "V2", ref = "V3", ref.ref = "V3", alt = "V4", ref.alt = "V4", exclude.ambiguous = T, silent = F, rm.duplicates = T)
RV_keep <- data.table(impu$V5[match$ref.extract])
fwrite(RV_keep, paste0(path_rvs, "/RV_CaG_gnomad_impu.snplist"), col.names = FALSE, row.names = FALSE, sep = "\t")
for (chr in 1:22){
  system(paste0("bcftools view -i 'ID=@", path_rvs, "/RV_CaG_gnomad_impu.snplist' -O z -o ", path_rvs, "/impute5_gigi2_combined_seq_RV_chr", chr, ".vcf.gz ", path_impu, "/merged_with_seq/impute5_gigi2_combined_seq_chr", chr, ".vcf.gz"))
}

#Add pedigree information to the ped. We use the chr22 to match the values but apply the result on every chrm.
ped <- data.table(system(paste0("bcftools query -l ", path_rvs, "/impute5_gigi2_combined_seq_RV_chr22.vcf.gz"), intern = TRUE))
ped_list <- strsplit(ped$V1, "_", fixed = TRUE)
ped$V1 <- sapply(ped_list, "[[", 1); ped$V2 <- as.character(sapply(ped_list, "[[", 2))
ref <- fread("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/Photo2021.pre", header = TRUE)[,1:5]
ref$`Individual ID` <- as.character(ref$`Individual ID`)
match <- merge(x = ped, y = ref, by.x = "V2", by.y = "Individual ID", all.x = TRUE, sort = FALSE)
match <- match[,c("Family ID", "V2", "Paternal ID", "Maternal ID", "Sex (1=male; 2=female; other=unknown)")]
match[match$V2 == "2620b", c("Family ID", "Paternal ID", "Maternal ID", "Sex (1=male; 2=female; other=unknown)")] <- match[match$V2 == "2620", c("Family ID", "Paternal ID", "Maternal ID", "Sex (1=male; 2=female; other=unknown)")]
fwrite(data.table(match[,c("Family ID", "V2", "Paternal ID", "Maternal ID")]), paste0(path_rvs, "/update_parent_IDs.txt"), col.names = FALSE, row.names = FALSE, sep = "\t")
fwrite(data.table(match[,c("Family ID", "V2", "Sex (1=male; 2=female; other=unknown)")]), paste0(path_rvs, "/update_sex.txt"), col.names = FALSE, row.names = FALSE, sep = "\t")
for (chr in 1:22){
  #PLINK will change the order of the alleles (A1 will become the minor allele and A2, the major). We keep the .bim file to get the correspondance when using the ped file in 1-2 format.
  system(paste0('plink --vcf ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.vcf.gz --update-parents ', path_rvs, '/update_parent_IDs.txt --update-sex ', path_rvs, '/update_sex.txt --nonfounders --maf 0.000000001 --make-bed --out ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr))
  system(paste0('plink --bfile ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, ' --keep-allele-order --freq --nonfounders --recode 12  --out ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr))
  system(paste0('rm ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.fam ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.bed ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.nosex'))
}
quit()
n

#Remove variants that are still frequent in our data and that appear in more than 20% of the family.
cd ${path_rvs}
for chr in {1..22}; do sbatch --job-name=CV_20prct_fam_${chr} --time=03:00:00 --mem=20G --export=chr=${chr} /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples/call/CV_20prct_fam.sh; done

#Produce the files by TADs and compute variants frequency
mkdir ${path_rvs}/TADs
cd ${path_rvs}/TADs
for chr in {1..22}
do
  awk -v value="chr${chr}" '$1 == value' /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/Subset_TADs/TADs_all_chrs.bed > ${path_rvs}/TADs/TADs_list_chr${chr}.bed
  n_TADs=$(wc -l < "${path_rvs}/TADs/TADs_list_chr${chr}.bed")
  for ((TAD=1; TAD<=n_TADs; TAD++))
  do
    awk '{print $0, "\t" "1"}' <(sed -n ${TAD}p ${path_rvs}/TADs/TADs_list_chr${chr}.bed | cut -f-3) > ${path_rvs}/TADs/TADs_tmp.bed
    plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr} --extract range ${path_rvs}/TADs/TADs_tmp.bed --keep-allele-order --recode --out ${path_rvs}/TADs/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_TAD_${TAD}
    plink --file ${path_rvs}/TADs/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_TAD_${TAD} --keep-allele-order --nonfounders --freq --out ${path_rvs}/TADs/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_TAD_${TAD}
  done
  rm ${path_rvs}/TADs/TADs_tmp.bed
done

#Files by CHRs in 0 TAD
mkdir ${path_rvs}/TADs/overlap_0
cd ${path_rvs}/TADs/overlap_0
for chr in {1..22}
do
  for file in /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/CRHs_overlap_0_TAD/chr${chr}/*
  do
    #Extract the filename without .txt
    filename=$(basename "$file" .bed)
    #get the values in the filename (chr and CRH)
    IFS='_' read -r crh neu ipsc value_crh value_chr <<< "$filename"
    value_chr="${value_chr:3}"
    cp $file ${path_rvs}/TADs/overlap_0/chr${value_chr}_CRH_${value_crh}.txt
    plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${value_chr} --extract range <(sed "s/\r//g" $file | awk 'NR > 1 {print $1, $2, $3, NR}' | sed "s/chr//") --keep-allele-order --recode --out ${path_rvs}/TADs/overlap_0/impute5_gigi2_combined_seq_RV_FINAL_chr${value_chr}_CRH_${value_crh}
  done
done

#Files by CHRs in 2 TADs
mkdir ${path_rvs}/TADs/overlap_2
cd ${path_rvs}/TADs/overlap_2
for chr in {1..22}
do
  for file in /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/CRHs_overlap_2_TADs/chr${chr}/*
  do
    #Extract the filename without .txt
    filename=$(basename "$file" .bed)
    #get the values in the filename (chr and CRH)
    IFS='_' read -r crh neu ipsc value_crh value_chr <<< "$filename"
    value_chr="${value_chr:3}"
    cp $file ${path_rvs}/TADs/overlap_2/chr${value_chr}_CRH_${value_crh}.txt
    plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${value_chr} --extract range <(sed "s/\r//g" $file | awk 'NR > 1 {print $1, $2, $3, NR}' | sed "s/chr//") --keep-allele-order --recode --out ${path_rvs}/TADs/overlap_2/impute5_gigi2_combined_seq_RV_FINAL_chr${value_chr}_CRH_${value_crh}
  done
done