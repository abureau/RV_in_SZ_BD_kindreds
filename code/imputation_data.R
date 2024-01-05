#This code should be run using the shell code `imputation.sh`.

library(data.table); library(dplyr)
path_data <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/phasing"
path_impu <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/imputation"

#We use the WGS datasets generated for the phasing. It is already separated by choromosome, mendel errors are fixed to 0 and GRCh38 genetic positions are already interpolated.
#Note that variants that are beyond the first and the last physical position in the genetic map reference are removed (1,937 among the 16,430,520 variants) 
#Family information was already verified in the datasets.

args <- commandArgs(TRUE)
chr <- as.numeric(args[1])

#Path where we can find gl_auto results by chromosome and where GIGI2 will be used
system(paste0("mkdir ", path_impu, "/chr", chr))
setwd(paste0(path_impu, "/chr", chr))

#Family, father and mother IDs need to be added to the sequencing data
fam <- fread("check.oped", skip = 6, colClasses = 'character')
fam$FID <- strsplit(fam$V1, "_", fixed = TRUE) %>% sapply(., head, n=1)
fam$IID <- strsplit(fam$V1, "_", fixed = TRUE) %>% sapply(., tail, n=1)
fam$FTR <- strsplit(fam$V2, "_", fixed = TRUE) %>% sapply(., tail, n=1)
fam$MTR <- strsplit(fam$V3, "_", fixed = TRUE) %>% sapply(., tail, n=1)
colnames(fam)[4] <- "SEX"
fam <- fam[,c("FID", "IID", "FTR", "MTR", "SEX")]

#WGS data data
wgs_ped <- fread(paste0(path_data, "/seq_FINAL_with_mask_gen_map_var_chr", chr, ".ped"), colClasses = 'character', data.table = FALSE)
wgs_map <- fread(paste0(path_data, "/seq_FINAL_with_mask_gen_map_var_chr", chr, ".map"), colClasses = 'character', data.table = FALSE)

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
fwrite(pos_dens[out_order,], "dense.pos.txt", row.names=F, quote=F, col.names=F, sep = "\t")

#Convert allele in numbers. Monomorphic SNPs are coded as 1
pre_start <- pre[,1:6]
geno <- pre; geno[geno=="0"] <- "Z"
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
write(apply(pre_new[, c(1,7:dim(pre_new)[2])], 1, paste, collapse=" "), "dense_geno.txt")

#MAF extraction computed in CARTaGENE
if(!file.exists("CaG.frq")){
  system(paste0("plink --vcf /lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/RV/chr", chr, "_CaG_FINAL_with_mask.vcf.gz --freq --double-id --out ", path_impu, "/chr", chr, "/CaG"))
  system(paste0("rm ", path_impu, "/chr", chr, "/*.nosex"))
}
freqfile <- fread("CaG.frq", header = TRUE)
freqfile$SNP <- gsub("_", ":", freqfile$SNP)
freqfile$pos <- sapply(strsplit(freqfile$SNP, ":", fixed = TRUE), `[`, 2)
SNP <- as.character(freqfile$SNP)
freq <- as.numeric(freqfile$MAF)

#We keep all SNPs, but after the merge, we set a freq of 0.001 to missing values or values of 0
pos_dens$CHR <- as.numeric(chr)
pos_dens$pos <- sapply(strsplit(pos_dens$V2, ":", fixed = TRUE), `[`, 2)
pos_dens$A1 <- sapply(strsplit(pos_dens$V2, ":", fixed = TRUE), `[`, 3)
pos_dens$A2 <- sapply(strsplit(pos_dens$V2, ":", fixed = TRUE), `[`, 4)
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

freq <- freq_dense$MAF
geno <- pre_new[,-(1:6)]
tmp <- matrix(as.numeric(as.vector(as.matrix(geno))),ncol=ncol(geno)/2)
est.freqs.all2 <- apply(tmp,2,function(x) mean(x[x!=0]==2))
#If A2 has a freq <= 0.5, we set the freq of A1 to the major allele freq and the freq of A2 to the minor allele.
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
