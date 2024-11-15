module load StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16 r/4.2
#Find the exons associated with the CRHs with the annotated dataset.
path_rvs=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS
path_exons=${path_rvs}/genes
path_anno=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/annotations
path_effects_to_keep=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/2024_09_30_annotations_effects_to_keep.txt

mkdir ${path_exons}
#This file is found on Loic's Github. We need to download it on $path_exons.
cut -f14 ${path_exons}/EnhancerPredictions_NEU_iPSC_liftoverhg38.txt | tail -n +2 | sort | uniq > ${path_exons}/CRHs_list.txt
cut -f8 ${path_exons}/EnhancerPredictions_NEU_iPSC_liftoverhg38.txt | sort | uniq > ${path_exons}/genes_list.txt

#In our vcf, the gene names are found between "|".
sed -i "s/^/|/g" ${path_exons}/genes_list.txt
sed -i "s/$/|/g" ${path_exons}/genes_list.txt

#In some cases, the suffixe/prefix are used. Most of them are "AS".
#It appears that they are useful. We want to keep it.
grep "\-" ${path_exons}/genes_list.txt

tr -d '\r' < ${path_exons}/CRHs_list.txt | while read crh
do
   #For the actual CRH, get all the gene codes and the chromosome(s) where they are found.
   awk -v crh="$crh" '$14 == crh' "${path_exons}/EnhancerPredictions_NEU_iPSC_liftoverhg38.txt" | awk '{print "|"$8"|"}' | sort | uniq > ${path_exons}/genes_list_CRH${crh}.txt
   chrs=$(awk -v crh="$crh" '$14 == crh' "${path_exons}/EnhancerPredictions_NEU_iPSC_liftoverhg38.txt" | cut -f3 | sort | uniq | sed "s/chr//")
   for chr in ${chrs}
   do
     #Only keep the annotations regarding the gene(s) in question + only keep the variants in the exons (function found in the list ${path_effects_to_keep}).
     zcat ${path_anno}/seq_FINAL_chr${chr}_ccds_annoted_dbNSFP4_snpeff_with_mask.vcf.gz | grep -v "#" | cut -f-8 | grep -f ${path_exons}/genes_list_CRH${crh}.txt | grep -f ${path_effects_to_keep} > ${path_exons}/seq_FINAL_annotations_chr${chr}_CRH${crh}_genes_list.txt;
     awk '{print $1, $2, $2, NR}' ${path_exons}/seq_FINAL_annotations_chr${chr}_CRH${crh}_genes_list.txt >> ${path_exons}/seq_FINAL_annotations_CRH${crh}_genes_list.bed
   done
done

#MEMBERSHIP NUMBER PROBLEM
#Verify the matching between position in the file previously used
for chr in {1..22}
do
  file_chr=$(grep chr${chr} "${path_exons}/EnhancerPredictions_NEU_iPSC_liftoverhg38.txt")
  #Overlap 1 TAD
  for file in /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/TADs_for_CRHs_overlap_1_TAD/chr${chr}/*
  do
    while IFS= read -r line
    do
      field2=$(echo "$line" | cut -f2,3)
      cat <(grep "$field2" <(echo "$file_chr")) <(echo "$line") >> /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped/check_CRH_membership.txt
    done < <(tail -n +2 $file)
  done
  #Overlap 0 TAD
  #There are little modifications in the code as the files were note the same.
  for file in /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/CRHs_overlap_0_TAD/chr${chr}/*
  do
    CRH=$(basename "$file" | cut -d "_" -f4)
    while IFS= read -r line
    do
      field2=$(echo "$line" | awk 'BEGIN { OFS = "\t" } { print $2, $3}')
      cat <(grep "$field2" <(echo "$file_chr")) <(paste <(echo "$line")) >> /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped/check_CRH_membership.txt
    done < <(tail -n +2 $file | sed "s/\r//g" | awk -v CRH="$CRH" '{ print $0 "\t" CRH}')
  done
  #Overlap 2 TADs
  for file in /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/CRHs_overlap_2_TADs/chr${chr}/*
  do
    CRH=$(basename "$file" | cut -d "_" -f4)
    while IFS= read -r line
    do
      field2=$(echo "$line" | awk 'BEGIN { OFS = "\t" } { print $2, $3}')
      cat <(grep "$field2" <(echo "$file_chr")) <(paste <(echo "$line")) >> /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped/check_CRH_membership.txt
    done < <(tail -n +2 $file | sed "s/\r//g" | awk -v CRH="$CRH" '{ print $0 "\t" CRH}')
  done
done
Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples/call/membership_check.R

#Add exons to the TAD files
#First, we need to list every CRH in each TAD
mkdir ${path_rvs}/TADs/with_exons
R
path_impu <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs"
membership_equi <- data.table::fread("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped/CRH_problem_membership_equivalence.txt")
CRHs_by_TADs <- data.frame()
path_0 <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/CRHs_overlap_0_TAD"
path_1 <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/TADs_for_CRHs_overlap_1_TAD"
path_2 <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/CRHs_overlap_2_TADs"
for(chr in 1:22){
  files_0 <- list.files(paste0(path_0, "/chr", chr))
  files_1 <- list.files(paste0(path_1, "/chr", chr))
  files_2 <- list.files(paste0(path_2, "/chr", chr))
  for(file_1 in files_1){
    i <- which(files_1 == file_1)
    CRHs_by_TAD <- read.table(paste0(path_1, "/chr", chr, "/", file_1), header=TRUE)
    CRHs_by_TADs <- rbind(CRHs_by_TADs, data.frame("chr" = chr,  "CRH_membership" = unique(CRHs_by_TAD$name), "overlap" = 1, "TAD" = i))
  }
  for(file_0 in files_0){
    i <- which(files_0 == file_0)
    CRHs_by_TAD <- read.table(paste0(path_0, "/chr", chr, "/", file_0), header=TRUE)
    CRHs_by_TADs <- rbind(CRHs_by_TADs, data.frame("chr" = chr, "CRH_membership" = unlist(strsplit(file_0, "_", fixed = TRUE))[[4]], "overlap" = 0, "TAD" = i))
  }
  for(file_2 in files_2){
    i <- which(files_2 == file_2)
    CRHs_by_TAD <- read.table(paste0(path_2, "/chr", chr, "/", file_2), header=TRUE)
    CRHs_by_TADs <- rbind(CRHs_by_TADs, data.frame("chr" = chr, "CRH_membership" = unlist(strsplit(file_2, "_", fixed = TRUE))[[4]], "overlap" = 2, "TAD" = i))
  }
}
CRHs_by_TADs <- merge(x = CRHs_by_TADs, y = membership_equi, by.x = "CRH_membership", by.y = "membership_EnhancerPredictions_NEU_iPSC_liftoverhg38_file", all.x = TRUE, sort = FALSE)
colnames(CRHs_by_TADs) <- c("CRH_membership", "chr", "overlap", "TAD", "CRH_membership_equivalence")
data.table::fwrite(CRHs_by_TADs, "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/with_exons/CRHs_by_TADs.txt", row.names = FALSE, col.names = TRUE, sep = "\t")
quit()
n

#Then, we need to check if some CRH are empty (genes mentionned in EnhancerPredictions_NEU_iPSC_liftoverhg38.txt not found in our annotation) 
#but also which TADs will be affected by the addition of the genes. In fact, for a given TAD, if the genes position are included in the TAD's border, it will not change anything
Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples/call/2024-09-17_CRH_empty_gene_border.R
#Affected TADs: /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/with_exons/CRH_affecting_TAD_by_gene_addition.txt
for chr in {1..22}
do
   awk '{print $4}' ${path_rvs}/TADs/TADs_list_chr${chr}.bed | while read TAD
  do
    awk -v chr="$chr" -v TAD="$TAD" '$2 == chr && $3 == 1 && $4 == TAD' ${path_rvs}/TADs/with_exons/CRHs_by_TADs.txt | cut -f5 | while read CRH
    do
      if [ $(grep -c -w "$CRH" ${path_rvs}/TADs/with_exons/CRH_affecting_TAD_by_gene_addition.txt) -eq 1 ]
      then
        cat ${path_rvs}/genes/seq_FINAL_annotations_CRH${CRH}_genes_list.bed | cut -d " " -f-3 >> ${path_rvs}/TADs/with_exons/chr${chr}_TAD${TAD}_CRHs_genes.tmp
      fi
    done
    if [ -e ${path_rvs}/TADs/with_exons/chr${chr}_TAD${TAD}_CRHs_genes.tmp ]
    then
      sort -n -k 2 ${path_rvs}/TADs/with_exons/chr${chr}_TAD${TAD}_CRHs_genes.tmp | uniq | awk '{print $0, NR}' > ${path_rvs}/TADs/with_exons/chr${chr}_TAD${TAD}_CRHs_genes.bed; rm ${path_rvs}/TADs/with_exons/chr${chr}_TAD${TAD}_CRHs_genes.tmp

      awk -v chr="$chr" -v TAD="$TAD" '$1 == "chr"chr && $4 == TAD' ${path_rvs}/TADs/TADs_list_chr${chr}.bed | sed "s/chr//"  >> ${path_rvs}/TADs/with_exons/chr${chr}_TAD${TAD}_CRHs_genes.bed
      awk -i inplace '{print $1, $2, $3, NR}' ${path_rvs}/TADs/with_exons/chr${chr}_TAD${TAD}_CRHs_genes.bed
    fi
  done
done


#---- CHRs - Overlap 1 TAD ----

#Correct the TADs
for chr in {1..22}
do
  for file in ${path_rvs}/TADs/with_exons/chr${chr}_*.bed
  do
    TAD=$(echo $file | cut -d "/" -f12| cut -d "_" -f2 | sed "s/TAD//")
    plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr} --extract range $file --keep-allele-order --recode --out ${path_rvs}/TADs/with_exons/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_TAD_${TAD}
    plink --file ${path_rvs}/TADs/with_exons/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_TAD_${TAD} --keep-allele-order --nonfounders --freq --out ${path_rvs}/TADs/with_exons/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_TAD_${TAD}
  done
done

#Adding the exons in the CRHs files /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/TADs_for_CRHs_overlap_1_TAD
cp -R /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/TADs_for_CRHs_overlap_1_TAD ${path_rvs}/TADs/with_exons
cd  ${path_rvs}/TADs/with_exons/TADs_for_CRHs_overlap_1_TAD
for chr in {1..22}
do
  for TAD in ${path_rvs}/TADs/with_exons/TADs_for_CRHs_overlap_1_TAD/chr${chr}/*.bed
  do
    cp ${TAD} ${TAD}_backup
    cut -f 4 ${TAD}_backup | tail -n +2 | uniq | while read CRH
    do
      CRH_equi=$(awk -v CRH="$CRH" '$2 == CRH' /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped/CRH_problem_membership_equivalence.txt | cut -f1)
      cat ${TAD} <(awk -v CRH="$CRH" 'BEGIN{OFS="\t"}{print "chr"$1, $2, $3, CRH}' ${path_exons}/seq_FINAL_annotations_CRH${CRH_equi}_genes_list.bed) > ${TAD}_${CRH}_tmp.txt
      mv ${TAD}_${CRH}_tmp.txt ${TAD}
    done
    sort -k 2 -n ${TAD} > ${TAD}_tmp.txt
    mv ${TAD}_tmp.txt ${TAD}
    rm ${TAD}_backup
  done
done

#There might be cases where some intervals in the .bed overlap.
R
for (chr in 1:22){
  CRH_files <- list.files(paste0("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/with_exons/TADs_for_CRHs_overlap_1_TAD/chr", chr))
  CRH_files <- CRH_files[!grepl("_backup", CRH_files, fixed = TRUE)]
  overlap <- data.frame()
  for (file in CRH_files){
    CRH_og <- data.table::fread(paste0("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/with_exons/TADs_for_CRHs_overlap_1_TAD/chr", chr, "/", file))
    CRH <- data.frame(CRH_og)
    overlap_CRH <- 0
    if(nrow(CRH)>1){
      for (i in 2:nrow(CRH)){
        last_upper <- as.numeric(CRH[i-1, 3])
        present_lower <- as.numeric(CRH[i, 2])
        prob <- ifelse(present_lower < last_upper, yes = TRUE, no = FALSE)
        if(!prob){
          next
        } else {
          CRH[i-1, 3] <- CRH[i, 3]
          CRH[i, 2] <- NA
          overlap_CRH <- overlap_CRH+1
        }
      }
    }
    CRH <- CRH[!is.na(CRH[[2]]),]
    overlap <- rbind(overlap, data.frame("chr" = chr, "CRH_file" = file, "overlap" = overlap_CRH))
    system(paste0("mv /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/with_exons/TADs_for_CRHs_overlap_1_TAD/chr", chr, "/", file, " /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/with_exons/TADs_for_CRHs_overlap_1_TAD/chr", chr, "/", file, "_overlap"))
    data.table::fwrite(data.table::as.data.table(CRH), paste0("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/with_exons/TADs_for_CRHs_overlap_1_TAD/chr", chr, "/", file), col.names = FALSE, row.names = FALSE, sep = "\t")
  }
  data.table::fwrite(data.table::as.data.table(overlap), paste0("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/with_exons/TADs_for_CRHs_overlap_1_TAD/chr", chr, "/overlaps.txt"), col.names = FALSE, row.names = FALSE, sep = "\t")
}
quit()
n

#---- CHRs - Overlap 0 TAD ----
#/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/CRHs_overlap_0_TAD
mkdir ${path_rvs}/TADs/overlap_0/with_exons

#Adding the exons in the CRHs files
cd  ${path_rvs}/TADs/overlap_0/with_exons
for CRH_file in ${path_rvs}/TADs/overlap_0/*.txt
do
  CRH_basename=$(basename ${CRH_file})
  CRH=$(echo ${CRH_file} | cut -d "/" -f 12 | cut -d "_" -f3 | sed "s/.txt//g")

  CRH_equi=$(awk -v CRH="$CRH" '$2 == CRH' /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped/CRH_problem_membership_equivalence.txt | cut -f1)
  cat <(sed "s/\r//g" ${CRH_file}) <(awk '{print "chr"$1, $2, $3}' ${path_exons}/seq_FINAL_annotations_CRH${CRH_equi}_genes_list.bed) | awk -v CRH="$CRH" 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, CRH}' > ${path_rvs}/TADs/overlap_0/with_exons/${CRH_basename}_tmp
  cat <(echo -e "name\tchrom\tchromStart\tchromEnd") <(cat ${path_rvs}/TADs/overlap_0/with_exons/${CRH_basename}_tmp) | sort -k 2 -n > ${path_rvs}/TADs/overlap_0/with_exons/${CRH_basename}
done
rm ${path_rvs}/TADs/overlap_0/with_exons/*_tmp

#There might be cases where some intervals in the .bed overlap.
R
CRH_files <- list.files("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/overlap_0/with_exons")
CRH_files <- CRH_files[grepl(".txt", CRH_files, fixed = TRUE)]
overlap <- data.frame()
for (file in CRH_files){
  CRH_og <- data.table::fread(paste0("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/overlap_0/with_exons/", file))
  CRH <- data.frame(CRH_og)
  overlap_CRH <- 0
  if(nrow(CRH)>1){
    for (i in 2:nrow(CRH)){
      last_upper <- as.numeric(CRH[i-1, 3])
      present_lower <- as.numeric(CRH[i, 2])
      prob <- ifelse(present_lower < last_upper, yes = TRUE, no = FALSE)
      if(!prob){
        next
      } else {
        CRH[i-1, 3] <- CRH[i, 3]
        CRH[i, 2] <- NA
        overlap_CRH <- overlap_CRH+1
      }
    }
  }
  CRH <- CRH[!is.na(CRH[[2]]),]
  overlap <- rbind(overlap, data.frame("CRH_file" = file, "overlap" = overlap_CRH))
  system(paste0("mv /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/overlap_0/with_exons/", file, " /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/overlap_0/with_exons/", file, "_overlap"))
  data.table::fwrite(data.table::as.data.table(CRH), paste0("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/overlap_0/with_exons/", file), col.names = FALSE, row.names = FALSE, sep = "\t")
}
data.table::fwrite(data.table::as.data.table(overlap), paste0("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/overlap_0/with_exons/overlaps.txt"), col.names = FALSE, row.names = FALSE, sep = "\t")
quit()
n

#Then, create the .ped files
for file in ${path_rvs}/TADs/overlap_0/with_exons/*.txt
do
  chr=$(cat $file | cut -f1 | sed "s/chr//g" | uniq)
  crh=$(echo $file | cut -d "/" -f13 | cut -d "_" -f3 | cut -d "." -f1)
  plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr} --extract range $file --keep-allele-order --recode --out ${path_rvs}/TADs/overlap_0/with_exons/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_CRH_${crh}
  plink --file ${path_rvs}/TADs/overlap_0/with_exons/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_CRH_${crh} --keep-allele-order --nonfounders --freq --out ${path_rvs}/TADs/overlap_0/with_exons/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_CRH_${crh}
done

#---- CHRs - Overlap 2 TADs ----
#/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/CRHs_overlap_2_TADs
mkdir ${path_rvs}/TADs/overlap_2/with_exons

#Adding the exons in the CRHs files
cd  ${path_rvs}/TADs/overlap_2/with_exons
for CRH_file in ${path_rvs}/TADs/overlap_2/*.txt
do
  CRH_basename=$(basename ${CRH_file})
  CRH=$(echo ${CRH_file} | cut -d "/" -f 12 | cut -d "_" -f3 | sed "s/.txt//g")

  CRH_equi=$(awk -v CRH="$CRH" '$2 == CRH' /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped/CRH_problem_membership_equivalence.txt | cut -f1)
  cat <(sed "s/\r//g" ${CRH_file}) <(awk '{print "chr"$1, $2, $3}' ${path_exons}/seq_FINAL_annotations_CRH${CRH_equi}_genes_list.bed) | awk -v CRH="$CRH" 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, CRH}' > ${path_rvs}/TADs/overlap_2/with_exons/${CRH_basename}_tmp
  cat <(echo -e "name\tchrom\tchromStart\tchromEnd") <(cat ${path_rvs}/TADs/overlap_2/with_exons/${CRH_basename}_tmp) | sort -k 2 -n > ${path_rvs}/TADs/overlap_2/with_exons/${CRH_basename}
done
rm ${path_rvs}/TADs/overlap_2/with_exons/*_tmp

#There might be cases where some intervals in the .bed overlap.
R
CRH_files <- list.files("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/overlap_2/with_exons")
CRH_files <- CRH_files[grepl(".txt", CRH_files, fixed = TRUE)]
overlap <- data.frame()
for (file in CRH_files){
  CRH_og <- data.table::fread(paste0("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/overlap_2/with_exons/", file))
  CRH <- data.frame(CRH_og)
  overlap_CRH <- 0
  if(nrow(CRH)>1){
    for (i in 2:nrow(CRH)){
      last_upper <- as.numeric(CRH[i-1, 3])
      present_lower <- as.numeric(CRH[i, 2])
      prob <- ifelse(present_lower < last_upper, yes = TRUE, no = FALSE)
      if(!prob){
        next
      } else {
        CRH[i-1, 3] <- CRH[i, 3]
        CRH[i, 2] <- NA
        overlap_CRH <- overlap_CRH+1
      }
    }
  }
  CRH <- CRH[!is.na(CRH[[2]]),]
  overlap <- rbind(overlap, data.frame("CRH_file" = file, "overlap" = overlap_CRH))
  system(paste0("mv /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/overlap_2/with_exons/", file, " /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/overlap_2/with_exons/", file, "_overlap"))
  data.table::fwrite(data.table::as.data.table(CRH), paste0("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/overlap_2/with_exons/", file), col.names = FALSE, row.names = FALSE, sep = "\t")
}
data.table::fwrite(data.table::as.data.table(overlap), paste0("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/TADs/overlap_2/with_exons/overlaps.txt"), col.names = FALSE, row.names = FALSE, sep = "\t")
quit()
n

#Then, create the .ped files
for file in ${path_rvs}/TADs/overlap_2/with_exons/*.txt
do
  chr=$(cat $file | cut -f1 | sed "s/chr//g" | uniq)
  crh=$(echo $file | cut -d "/" -f13 | cut -d "_" -f3 | cut -d "." -f1)
  plink --file ${path_rvs}/impute5_gigi2_combined_seq_RV_FINAL_chr${chr} --extract range $file --keep-allele-order --recode --out ${path_rvs}/TADs/overlap_2/with_exons/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_CRH_${crh}
  plink --file ${path_rvs}/TADs/overlap_2/with_exons/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_CRH_${crh} --keep-allele-order --nonfounders --freq --out ${path_rvs}/TADs/overlap_2/with_exons/impute5_gigi2_combined_seq_RV_FINAL_chr${chr}_CRH_${crh}
done

