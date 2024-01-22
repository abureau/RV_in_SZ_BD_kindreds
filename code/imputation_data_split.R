#This code should be run using the shell code `imputation.sh`.
#When the .ped imported are too huge to be imported in R, we should importe them in half and merge the final outputs.

library(data.table); library(dplyr)
path_data <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/phasing"
path_impu <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/imputation"

#We use the WGS datasets generated for the phasing. It is already separated by chromosome, mendel errors are fixed to 0 and GRCh38 genetic positions are already interpolated.
#Note that variants that are beyond the first and the last physical position in the genetic map reference are removed (1,937 among the 16,430,520 variants) 
#Family information was already verified in the datasets.

args <- commandArgs(TRUE)
chr <- as.numeric(args[1])

#Path where we can find gl_auto results by chromosome and where GIGI2 will be used
system(paste0("mkdir ", path_impu, "/chr", chr))
path_chr <- paste0(path_impu, "/chr", chr)
setwd(path_chr)

#Number of variants
n_var <- system(paste0("wc -l ", path_data, "/seq_FINAL_with_mask_gen_map_var_chr", chr, ".map | awk '{print $1}'"), intern = TRUE) %>%
  as.numeric()
n_var_half_1 <- floor(n_var/2)
n_var_half_2 <- n_var-n_var_half_1

#We use a loop.
for(half in 1:2){
  
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
  if(half==1){idx_map <- 1:n_var_half_1}else{idx_map <- (1+n_var_half_1):n_var}
  wgs_map <- fread(paste0(path_data, "/seq_FINAL_with_mask_gen_map_var_chr", chr, ".map"), colClasses = 'character', data.table = FALSE)
  wgs_map <- wgs_map[idx_map,]
  fwrite(data.table(wgs_map$V2), "tmp.snplist", row.names = FALSE, col.names = FALSE)
  system(paste0("plink --file ", path_data, "/seq_FINAL_with_mask_gen_map_var_chr", chr, " --extract ", path_chr, "/tmp.snplist --recode --out ", path_chr, "/tmp_seq"))
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
  
  #Convert allele in numbers. Monomorphic SNPs are coded as 1
  pre_start <- pre[,1:6]
  geno <- pre[,-(1:6)]; geno[geno=="0"] <- "Z"
  tmp <- matrix(as.vector(as.matrix(geno)),ncol=ncol(geno)/2)
  geno.new <- matrix(as.vector(as.matrix(apply(tmp,2,function(x) as.numeric(as.factor(x))))),ncol=ncol(geno))
  geno.new[geno.new==3] <- 0
  pre <- data.frame(pre_start, geno.new)
  
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
}

#Merge the dense.pos outputs
system(paste0("cat ", path_chr, "/dense.pos_1.txt ", path_chr, "/dense.pos_2.txt > ", path_chr, "/dense.pos.txt"))

#To merge the dense_geno txt file, we need to get them in the long format, remove the header of the second file and merge them.
#After, we can set the variants ID to rs0, rs1, rs2,...
system(paste0("/lustre03/project/6033529/SOFT/GIGI2/src/Wide2Long ", path_chr, "/dense_geno_1.txt ", path_chr, "/dense_geno_long_1.txt"))
system(paste0("/lustre03/project/6033529/SOFT/GIGI2/src/Wide2Long ", path_chr, "/dense_geno_2.txt ", path_chr, "/dense_geno_long_2.txt"))
system(paste0("sed -i '1d' ", path_chr, "/dense_geno_long_2.txt"))
system(paste0("cat ", path_chr, "/dense_geno_long_1.txt ", path_chr, "/dense_geno_long_2.txt | awk 'BEGIN {OFS=\" \"} NR==1 {print; next} { $1=\"rs\" (NR-2); print }'> ", path_chr, "/dense_geno_long.txt"))

#Import the objects
pos_dens <- fread(paste0(path_chr, "/dense.pos.txt")) %>% setNames(., c("V2", "V3"))
out_order <- order(pos_dens$V3)
est.freqs.all2_1 <- readRDS(paste0("est.freqs.all2_1.RDS"))
est.freqs.all2_2 <- readRDS(paste0("est.freqs.all2_2.RDS"))
est.freqs.all2 <- c(est.freqs.all2_1, est.freqs.all2_2)

#MAF extraction computed in CARTaGENE
if(!file.exists("CaG.frq")){
  system(paste0("plink --vcf /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/RV/chr", chr, "_CaG_FINAL_with_mask.vcf.gz --freq --double-id --out ", path_chr, "/CaG"))
  system(paste0("rm ", path_chr, "/*.nosex"))
}
freqfile <- fread("CaG.frq", header = TRUE)
freqfile$SNP <- gsub("_", ":", freqfile$SNP)
freqfile$pos <- sapply(strsplit(freqfile$SNP, ":", fixed = TRUE), `[`, 2)

#We keep all SNPs, but after the merge, we set a freq of 0.001 to missing values or values of 0
pos_dens$CHR <- as.numeric(chr)
pos_dens$pos <- sapply(strsplit(pos_dens$V2, ":", fixed = TRUE), `[`, 2)
pos_dens$A1 <- sapply(strsplit(pos_dens$V2, ":", fixed = TRUE), `[`, 3)
pos_dens$A2 <- sapply(strsplit(pos_dens$V2, ":", fixed = TRUE), `[`, 4)
#freqfile <- freqfile[as.numeric(freqfile$pos) >= min(as.numeric(pos_dens$pos)) & as.numeric(freqfile$pos) <= max(as.numeric(pos_dens$pos))]
SNP <- as.character(freqfile$SNP)
freq <- as.numeric(freqfile$MAF)
match <- lassosum:::matchpos(tomatch = freqfile, ref.df = pos_dens,
                             chr = "CHR", ref.chr = "CHR", pos = "pos", ref.pos = "pos",
                             ref = "A1", ref.ref = "A1", alt = "A2", ref.alt = "A2",
                             exclude.ambiguous = F, silent = T, rm.duplicates = F)
freq_to_add <- rep(NA, nrow(pos_dens))
freq_to_add[match(match$order, 1:nrow(pos_dens))] <- match$order
freq_dense <- data.frame(pos_dens, freqfile[freq_to_add, "MAF"])

#Verifications
all(freq_dense[,1] == pos_dens$V2)
freq_dense[, "MAF"] <- as.numeric(freq_dense[, "MAF"])
sum(freq_dense[, "MAF"] == 0, na.rm=T)
sum(is.na(freq_dense[, "MAF"]))
freq_dense[is.na(freq_dense[, "MAF"]) | freq_dense[, "MAF"] == 0, "MAF"] <- 0.001
write.table(pos_dens[out_order, "V3"], "dense_pos_vf.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#If A2 has a freq <= 0.5, we set the freq of A1 to the major allele freq and the freq of A2 to the minor allele.
freq <- freq_dense$MAF
freq1 <- rep(NA,length(freq))
freq2 <- rep(NA,length(freq))
freq1[est.freqs.all2<=0.5] <- 1-freq[est.freqs.all2<=0.5]
freq2[est.freqs.all2<=0.5] <- freq[est.freqs.all2<=0.5]
#If A2 has a freq > 0.5, we set the freq of A1 to the minor allele freq and the freq of A2 to the major allele.
freq1[est.freqs.all2>0.5] <- freq[est.freqs.all2>0.5]
freq2[est.freqs.all2>0.5] <- 1-freq[est.freqs.all2>0.5]
write.table(data.frame(freq1, freq2)[out_order,], "dense_freq.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

#Final marker list for GIGI2
write.table(freq_dense[out_order, 1], "liste_dense.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
