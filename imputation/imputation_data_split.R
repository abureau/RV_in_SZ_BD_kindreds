#This code should be run using the shell code `imputation.sh`.
#When the .ped imported are too huge to be imported in R, we should import them in half and merge the final outputs.

library(data.table); library(dplyr)
path_data <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/phasing"
path_impu <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/imputation_fam"

#We use the WGS datasets generated for the phasing. It is already separated by chromosome, mendel errors are fixed to 0 and GRCh38 genetic positions are already interpolated.
#Note that variants that are beyond the first and the last physical position in the genetic map reference are removed (1,937 among the 16,430,520 variants) 
#Family information was already verified in the datasets.

#For parallel computation...
args <- commandArgs(TRUE)
chr <- as.numeric(args[1])
sep <- as.numeric(args[2])

#Path where we can find gl_auto results by chromosome and where GIGI2 will be used
system(paste0("mkdir ", path_impu, "/chr", chr))
path_chr <- paste0(path_impu, "/chr", chr)
setwd(path_chr)

#Number of variants
n_var <- system(paste0("bcftools query -f '%ID\n' ", path_data, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".bcf | wc -l"), intern = TRUE) %>%
  as.numeric()
n_var_half_1 <- floor(n_var/sep)

#Samples
samples <- system(paste0("bcftools query -l ", path_data, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".bcf"), intern = TRUE)

#Convert alleles to 0 and 1.
fwrite(data.table(samples[!grepl("^0", samples)]), paste0(path_chr, "/seq_FINAL.tokeep"), col.names = FALSE, row.names = FALSE)
fwrite(data.table(samples[grepl("^0_Sample", samples)]), paste0(path_chr, "/seq_FINAL_CaG.tokeep"), col.names = FALSE, row.names = FALSE)
system(paste0("bcftools view -S ", path_chr, "/seq_FINAL.tokeep -O v -o ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".vcf ", path_data, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".bcf"))
system(paste0("bcftools view -S ", path_chr, "/seq_FINAL_CaG.tokeep -O v -o ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, "_CaG.vcf ", path_data, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".bcf"))
system(paste0("sed -i 's/_/ /g' ", path_chr, "/seq_FINAL.tokeep"))
system(paste0("sed -i 's/0_/0 /g' ", path_chr, "/seq_FINAL_CaG.tokeep"))
system(paste0("plink --vcf ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".vcf --keep-allele-order --nonfounders --freq --make-bed --recode 12 --out ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ""))
system(paste0("plink --vcf ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, "_CaG.vcf --keep-allele-order --nonfounders --freq --double-id --out ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, "_CaG"))
system(paste0("rm ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".vcf ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, "_CaG.vcf"))

#Compute genetic positions
gen_map <- "/lustre03/project/6033529/IBD_epilepsy/data/genetic_map_grch38"
system(paste0("Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/gen_pos_interpolation_grch38.R -bim ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".map -ref_dir ", gen_map, " -method 'linear'"))

#We use a loop.
for(half in 1:sep){
  
  #Family, father and mother IDs need to be added to the sequencing data
  fam <- fread("check.oped", skip = 6, colClasses = 'character')
  fam$FID <- strsplit(fam$V1, "_", fixed = TRUE) %>% sapply(., head, n=1)
  fam$IID <- strsplit(fam$V1, "_", fixed = TRUE) %>% sapply(., tail, n=1)
  fam$FTR <- strsplit(fam$V2, "_", fixed = TRUE) %>% sapply(., tail, n=1)
  fam$MTR <- strsplit(fam$V3, "_", fixed = TRUE) %>% sapply(., tail, n=1)
  colnames(fam)[4] <- "SEX"
  fam <- fam[,c("FID", "IID", "FTR", "MTR", "SEX")]
  
  #WGS data data
  #We will create a new .ped/.map where only a half of the variants are found.
  if(half==1){idx_map <- 1:n_var_half_1}else if(half==sep){idx_map <- (1+max(idx_map)):n_var}else{idx_map <- (1+max(idx_map)):(max(idx_map)+length(idx_map))}
  wgs_map <- fread(paste0(path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".map"), colClasses = 'character', data.table = FALSE)
  wgs_map <- wgs_map[idx_map,]
  fwrite(data.table(wgs_map$V2), "tmp.snplist", row.names = FALSE, col.names = FALSE)
  system(paste0("plink --file ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, " --keep-allele-order --extract ", path_chr, "/tmp.snplist --recode --out ", path_chr, "/tmp_seq"))
  wgs_ped <- fread(paste0(path_chr, "/tmp_seq.ped"), colClasses = 'character', data.table = FALSE)
  wgs_map$V3 <- as.numeric(wgs_map$V3)
  
  #Import subject IDs.
  id.fam.suj <- fread("check.oped", skip = 6, header = FALSE, colClasses = 'character')[[1]]
  id.suj <- strsplit(id.fam.suj, "_", fixed = TRUE) %>% sapply(., tail, n=1)
  ids.mat <- data.frame(id.suj, id.fam.suj)
  
  #Clean the IDs
  pre <- merge(x = ids.mat, y = wgs_ped, by.x = "id.suj", by.y = "V2", all.x = FALSE, all.y = TRUE, sort = FALSE)
  pre <- pre[,-1]; pre[,2] <- pre[,1]
  
  #Remove NA values, subjects not found in the pedigree (501 to 464 subjects (-37) : Intercept subjects 19xxx and 2620b) or duplicated subjects (none).
  pre <- pre[!is.na(pre[,1]),]
  pre <- pre[!duplicated(pre[,1]),]
  
  #Import framework markers
  #We need to keep the variants that are included in the markers framework only
  frame_pos <- fread(paste0("framework_pos.txt"))
  which.keep <- which(wgs_map$V3 <= max(frame_pos[,1]) & wgs_map$V3 >= min(frame_pos[,1]))
  pos_dens <- wgs_map[which.keep, c("V2","V3")]
  wgs_map <- wgs_map[which.keep, ]
  
  #GIGI2 wants the variants to be sorted by genetic position, we ensure that every output respects this order
  out_order <- order(wgs_map$V3)
  fwrite(pos_dens[out_order,], paste0("dense.pos_", half, ".txt"), row.names=F, quote=F, col.names=F, sep = "\t")
   
  #Remove variants not in framework marker
  which.geno.keep <- c(2*which.keep-1, 2*which.keep)
  which.geno.keep <- which.geno.keep[order(which.geno.keep)]
  geno <- pre[, 7:ncol(pre)]
  geno <- geno[, which.geno.keep]
  geno_order <- (2*rep(out_order, each = 2))+rep(c(-1,0), length(out_order))
  geno_out <- geno[, geno_order]
  pre_new <- data.frame(pre[,1:6], geno_out)
  write(apply(pre_new[, c(1,7:dim(pre_new)[2])], 1, paste, collapse=" "), paste0("dense_geno_", half, ".txt"))
  
  #Object need to the MAF output
  geno <- pre_new[,-(1:6)]
  tmp <- matrix(as.numeric(as.vector(as.matrix(geno))),ncol=ncol(geno)/2)
  est.freqs.all2 <- apply(tmp,2,function(x) mean(x[x!=0]==2))
  saveRDS(est.freqs.all2, paste0("est.freqs.all2_", half, ".RDS"))
  
  #Remove the temporary files
  system(paste0("rm ", path_chr, "/tmp*"))

  #To merge the dense_geno txt file, we need to get them in the long format, remove the header of the second file and merge them.
  system(paste0("/lustre03/project/6033529/SOFT/GIGI2/src/Wide2Long ", path_chr, "/dense_geno_", half, ".txt ", path_chr, "/dense_geno_long_", half, ".txt"))
  if(half!=1){system(paste0("sed -i '1d' ", path_chr, "/dense_geno_long_", half, ".txt"))}
}

#Merge the dense.pos outputs
wording_pos <-  paste0(path_chr, "/dense.pos_", 1:half, ".txt", collapse = " ")
system(paste0("cat ", wording_pos, " > ", path_chr, "/dense.pos.txt"))

#To merge the dense_geno txt file, we need to get them in the long format, remove the header of the second file and merge them.
#After, we can set the variants ID to rs0, rs1, rs2,...
wording_geno <-  paste0(path_chr, "/dense_geno_long_", 1:half, ".txt", collapse = " ")
system(paste0("cat ", wording_geno, " | awk 'BEGIN {OFS=\" \"} NR==1 {print; next} { $1=\"rs\" (NR-2); print }'> ", path_chr, "/dense_geno_long.txt"))

#Import the objects
pos_dens <- fread(paste0(path_chr, "/dense.pos.txt")) %>% setNames(., c("V2", "V3"))
out_order <- order(pos_dens$V3)
est.freqs.all2 <- c()
for(half in 1:sep){est.freqs.all2 <- c(est.freqs.all2, readRDS(paste0("est.freqs.all2_", half, ".RDS")))}

freqfile <- fread(paste0("seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, "_CaG.frq"))
freqfile <- freqfile[freqfile$SNP %in% pos_dens$V2,]
freq_out <- data.frame("MAF_A1" = freqfile$MAF, "MAF_A2" = 1-freqfile$MAF)
change_A1 <- which(freq_out$MAF_A1 == 0)
change_A2 <- which(freq_out$MAF_A2 == 0)
freq_out$MAF_A1[change_A1] <- 0.001
freq_out$MAF_A2[change_A1] <- 0.999
freq_out$MAF_A1[change_A2] <- 0.999
freq_out$MAF_A2[change_A2] <- 0.001
write.table(freq_out[out_order,], "dense_freq.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

#Final marker list for GIGI2
write.table(pos_dens[out_order, "V2"], "liste_dense.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(pos_dens[out_order, "V3"], "dense_pos_vf.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#only keeping the .bim allele is 1 (A1) and 2 (A2)
system(paste0("rm ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".ped ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".map ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".map~"))
system(paste0("rm ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".bed ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".fam ", path_chr, "/seq_FINAL_PHASED_with_mask_in_genmap_chr", chr, ".log"))
