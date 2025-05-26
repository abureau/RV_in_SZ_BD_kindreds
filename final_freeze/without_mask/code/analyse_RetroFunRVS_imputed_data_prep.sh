module load StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16 r/4.2
path_impu=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag_without_mask/imputation_comb
path_rvs=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag_without_mask/RetroFunRVS
mkdir ${path_rvs}

#The first half of this script is coded in R
R
library(data.table); library(dplyr)

#From our combined imputation files, keep the GnomaAD/CaG rare variants only
#Remove monomorphic variants
path_impu <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag_without_mask/imputation_comb"
path_rvs <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag_without_mask/RetroFunRVS"
RV_cag_gnomad <- fread("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag_without_mask/RV/RV_CaG_gnomad_matched_seq_without_mask.snplist")
RV_cag_only <- fread("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag_without_mask/RV/seq_not_in_gnomad_RV_in_CaG_without_mask.snplist")
if (!file.exists(paste0(path_impu, "/merged_with_seq/impute5_gigi2_combined_seq.bim"))){
  for (chr in 1:22){
    system(paste0("bcftools query -f '%CHROM %POS %REF %ALT %ID\n' ", path_impu, "/merged_with_seq/impute5_gigi2_combined_seq_chr", chr, ".vcf.gz >> ", path_impu, "/merged_with_seq/impute5_gigi2_combined_seq.bim"))
  }
}
all_seq <- fread("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag_without_mask/RV/seq_without_mask.bim")
impu <- fread(paste0(path_impu, "/merged_with_seq/impute5_gigi2_combined_seq.bim"), header = FALSE)
all_seq$chr_pos <- paste0(all_seq$V1, "_", all_seq$V2)
impu$chr_pos <- paste0(impu$V1, "_", impu$V2)
all_seq_equiv <- merge(x = impu[,c("V5", "chr_pos")], y = all_seq, by = "chr_pos", all.y = TRUE, sort = FALSE)
fwrite(data.table("ID_seq" = all_seq_equiv$V5.y, "ID_impu" = all_seq_equiv$V5.x),
       paste0(path_impu, "/imputation_ID_in_seq_equivalence.txt"), col.names = TRUE, row.names = FALSE, sep = "\t")
RV <- data.table(c(RV_cag_gnomad$ID_seq, RV_cag_only$ID_seq))
RV_list <- strsplit(RV$V1, ":", fixed = TRUE)
RV_bim <- data.frame("V1" = sapply(RV_list, "[[", 1), "V2" = as.numeric(sapply(RV_list, "[[", 2)), "V3" = sapply(RV_list, "[[", 3), "V4" = sapply(RV_list, "[[", 4))
match <- lassosum:::matchpos(tomatch = RV_bim, ref.df = all_seq_equiv, auto.detect.tomatch = F, auto.detect.ref = F, chr = "V1", ref.chr = "V6", pos = "V2", ref.pos = "V2", ref = "V3", ref.ref = "V3", alt = "V4", ref.alt = "V4",
                             exclude.ambiguous = F, silent = F, rm.duplicates = T)
RV_keep <- data.table(all_seq_equiv$`V5.x`[match$ref.extract])
fwrite(RV_keep, paste0(path_rvs, "/RV_CaG_gnomad_plus_seq_only_impu.snplist"), col.names = FALSE, row.names = FALSE, sep = "\t") #AVANT: RV_CaG_gnomad_impu.snplist
for (chr in 1:22){
  system(paste0("bcftools view -i 'ID=@", path_rvs, "/RV_CaG_gnomad_plus_seq_only_impu.snplist' -O z -o ", path_rvs, "/impute5_gigi2_combined_seq_RV_chr", chr, ".vcf.gz ", path_impu, "/merged_with_seq/impute5_gigi2_combined_seq_chr", chr, ".vcf.gz"))
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
system(paste0("mkdir ", path_rvs, "/log"))
for (chr in 1:22){
  #PLINK will change the order of the alleles (A1 will become the minor allele and A2, the major). We keep the .bim file to get the correspondance when using the ped file in 1-2 format.
  system(paste0('plink --vcf ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.vcf.gz --update-parents ', path_rvs, '/update_parent_IDs.txt --update-sex ', path_rvs, '/update_sex.txt --nonfounders --maf 0.000000001 --keep-allele-order --make-bed --out ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr))
  system(paste0('mv ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.log ', path_rvs, '/log'))
  system(paste0('plink --bfile ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, ' --keep-allele-order --freq --nonfounders --recode 12  --out ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr))
  system(paste0('rm ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.fam ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.bed ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.nosex ', path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.log'))
}
quit()
n
grep "variants removed due to minor allele" ${path_rvs}/log/impute5_gigi2_combined_seq_RV_chr*.log > ${path_rvs}/log/impute5_gigi2_combined_seq_RV_variants_removed.txt
grep "people pass filters and QC" ${path_rvs}/log/impute5_gigi2_combined_seq_RV_chr*.log > ${path_rvs}/log/impute5_gigi2_combined_seq_RV_variants_remaining.txt
cut -f2 -d ":" ${path_rvs}/log/impute5_gigi2_combined_seq_RV_variants_removed.txt | cut -f1 -d " " | awk '{ sum += $1 } END { print sum }'
cut -f2 -d ":" ${path_rvs}/log/impute5_gigi2_combined_seq_RV_variants_remaining.txt | cut -f1 -d " " | awk '{ sum += $1 } END { print sum }'

#Remove variants that are still frequent in our data and that appear in more than 20% of the family.
cd ${path_rvs}
R
library(data.table); library(dplyr)
path_rvs <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag_without_mask/RetroFunRVS"
max_fam <- ceiling(0.20*48) #We have 48 families.
args <- commandArgs(TRUE)
chr <- as.numeric(args[1])
for(chr in 4:22){
  to_exclude <- c()
  print(paste0("chr : ", chr))
  system(paste0("plink --file ", path_rvs, "/impute5_gigi2_combined_seq_RV_chr", chr, " --keep-allele-order --nonfounders --freq --family --out ", path_rvs, "/impute5_gigi2_combined_seq_RV_chr", chr, "_family"))
  system(paste0("rm ", path_rvs, "/impute5_gigi2_combined_seq_RV_chr", chr, "_family.log"))

  freq_fam <- fread(paste0(path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '_family.frq.strat'))
  map <- fread(paste0(path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.map'))
  freq <- fread(paste0(path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.frq'))
  CV <- which(freq$MAF>0.01)
  map <- map[CV]; freq <- freq[CV]
  freq_fam <- freq_fam[freq_fam$SNP %in% map$V2,]
  n_var <- length(CV)
  for(var in 1:n_var){
    var_i <- map$V2[var]
    freq_fam_i <- freq_fam[freq_fam$SNP == var_i,]
    n_fam_minor_present <- sum(freq_fam_i$MAC>0)
    if(n_fam_minor_present >= max_fam){to_exclude <- c(to_exclude, var_i)}
  }
  #Remove the concerned variants.
  fwrite(data.table(to_exclude), paste0(path_rvs, "/filtered_CV_to_exclude_chr_", chr, ".snplist"), row.names = FALSE, col.names = FALSE)
  system(paste0("plink --file ", path_rvs, "/impute5_gigi2_combined_seq_RV_chr", chr, " --exclude ", path_rvs, "/filtered_CV_to_exclude_chr_", chr, ".snplist --keep-allele-order --nonfounders --freq --recode --out ", path_rvs, "/impute5_gigi2_combined_seq_RV_FINAL_chr", chr))
  system(paste0("mv ", path_rvs, "/impute5_gigi2_combined_seq_RV_FINAL_chr", chr, ".log ", path_rvs, "/log/impute5_gigi2_combined_seq_RV_freq_fam_FINAL_chr", chr, ".log"))
  system(paste0("rm ", path_rvs, "/impute5_gigi2_combined_seq_RV_chr", chr, ".ped ", path_rvs, "/impute5_gigi2_combined_seq_RV_chr", chr, ".map"))
  system(paste0("rm ", path_rvs, "/impute5_gigi2_combined_seq_RV_chr", chr, "_family.frq.strat"))
}
quit
n
grep "variants remaining" ${path_rvs}/log/impute5_gigi2_combined_seq_RV_freq_fam_FINAL_chr*.log > ${path_rvs}/log/impute5_gigi2_combined_seq_RV_freq_fam_FINAL_variants_remaining.txt
cut -f2 -d " " ${path_rvs}/log/impute5_gigi2_combined_seq_RV_freq_fam_FINAL_variants_remaining.txt | awk '{ sum += $1 } END { print sum }'

#Run the following code to format the files.
#/lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples/call/TADs_file_liftOver_executable_Loic_format.sh

#We need to create a softlink for the directory /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped
mkdir ${path_rvs}/objets_ped
ln -s /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped/* ${path_rvs}/objets_ped

#Produce the files by TADs and compute variants frequency
mkdir ${path_rvs}/TADs
cd ${path_rvs}/TADs
for chr in {1..22}
do
  awk -v value="chr${chr}" '$1 == value' /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/Subset_TADs/TADs_all_chrs.bed > ${path_rvs}/TADs/TADs_list_chr${chr}.bed
  n_TADs=$(wc -l < "${path_rvs}/TADs/TADs_list_chr${chr}.bed")
  for ((TAD=1; TAD<=n_TADs; TAD++))
  do
    awk '{print $0, "\t" "1"}' <(sed -n ${TAD}p ${path_rvs}/TADs/TADs_list_chr${chr}.bed | cut -f-3) > ${path_rvs}/TADs/TADs_tmp_chr${chr}_TAD${TAD}.bed
    plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr} --extract range ${path_rvs}/TADs/TADs_tmp_chr${chr}_TAD${TAD}.bed --keep-allele-order --recode --out ${path_rvs}/TADs/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_TAD_${TAD}
    plink --file ${path_rvs}/TADs/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_TAD_${TAD} --keep-allele-order --nonfounders --freq --out ${path_rvs}/TADs/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_TAD_${TAD}
    rm ${path_rvs}/TADs/TADs_tmp_chr${chr}_TAD${TAD}.bed
  done
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
for chr in {11..22}
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