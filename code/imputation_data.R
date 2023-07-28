library(data.table)
library(dplyr)
pathData <- "/lustre03/project/6033529/quebec_10x/data/freeze/QC"
pathimpu <- "/lustre03/project/6033529/quebec_10x/data/freeze/imputation"

chr_loop <- as.character(1:22)

for(chr in chr_loop){
  print(paste0("chr : ", chr))
  #Path where we can find gl_auto results by chromosome and where GIGI2 will be used
  setwd(paste0(pathimpu, "/chr", chr))

  #Split the QCed data by chromosome. 
  #system(paste0("plink --vcf ", pathData, "/seq_FINAL_with_mask.vcf.gz --chr ", chr, " --recode --out ", pathData, "/seq_FINAL_", chr, "_with_mask"))

  #Family, father and mother IDs need to be added to the sequencing data
  fam <- fread("check.oped", skip = 6)
  fam$FID <- strsplit(fam$V1, "_", fixed = TRUE) %>% sapply(., head, n=1)
  fam$IID <- strsplit(fam$V1, "_", fixed = TRUE) %>% sapply(., tail, n=1)
  fam$FTR <- strsplit(fam$V2, "_", fixed = TRUE) %>% sapply(., tail, n=1)
  fam$MTR <- strsplit(fam$V3, "_", fixed = TRUE) %>% sapply(., tail, n=1)
  colnames(fam)[4] <- "SEX"
  fam <- fam[,c("FID", "IID", "FTR", "MTR", "SEX")]

  #ped data
  pre.scan <- fread(paste0(pathData, "/seq_FINAL_", chr, "_with_mask.ped"))
  pre1 <- as.data.frame(pre.scan)
  colnames(pre1)[1] <- "id.suj"

  #Clean the IDs
  pre1[,1] <- strsplit(x = pre1[,1], split = "-", fixed = TRUE) %>%
	 	sapply(., FUN = tail, n = 1) %>% 
		stringr::str_replace(., "^0+", "")
  pre1[,2] <- strsplit(x = pre1[,2], split = "-", fixed = TRUE) %>%
	 	sapply(., FUN = tail, n = 1) %>% 
		stringr::str_replace(., "^0+", "")
  pre <- merge(x = pre1, y = fam, by.y = "IID", by.x = "id.suj", all.x = TRUE, all.y = FALSE)
  pre[,1] <- pre$FID
  pre[,3] <- pre$FTR
  pre[,4] <- pre$MTR
  pre[,5] <- pre$SEX
  pre <- pre[,1:(ncol(pre)-4)]
  pre[pre[,2] == "2007:R0", c(1,3,4,5)] <- 0
  pre[pre[,2] == "2620b", c(1,3,4,5)] <- 0
  pre <- pre[match(pre1[,2], pre[,2]),]
  fwrite(pre, paste0(pathData, "/mendel_corrected/seq_FINAL_full_ID_", chr, "_with_mask.ped"), row.names = FALSE, col.names = FALSE, sep = "\t")  

  #Create the resulting bfiles to set mendel error to missing values. We also remove SNPs where no correspondance was found from GRCh38 to hg19 (GRCh37)
  system(paste0("plink --map ", pathData , "/seq_FINAL_", chr, "_with_mask.map --ped ", pathData , "/mendel_corrected/seq_FINAL_full_ID_", chr, "_with_mask.ped --set-me-missing --mendel-duos --mendel-multigen --make-bed --out ", pathData , "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask"))
  system(paste0("plink --bfile ", pathData , "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask --recode --out ", pathData , "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask"))
  system(paste0("rm ", pathData , "/mendel_corrected/*.nosex"))

  #We perform a liftOver from GRCh38 to hg19 (GRCh37) to get genetic positions that match with the ones used in gl_auto
  #For liftOver information, see https://genome.sph.umich.edu/wiki/LiftOver#Lift_PLINK_format 
  liftOver_chain <- "/lustre03/project/6033529/SOFT/hg38ToHg19.over.chain.gz"
  liftOver_exec <- "/lustre03/project/6033529/SOFT/liftOver"
  system(paste0("python /lustre03/project/6033529/SOFT/liftOverPlink/liftOverPlink.py -m ", pathData , "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask.map -p ", pathData , "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask.ped -o ", pathData , "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask_hg19 -c ", liftOver_chain, " -e ", liftOver_exec))
  system(paste0("grep -v '^#' ", pathData , "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask_hg19.bed.unlifted | cut -f 4 > ", pathData , "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask_hg19_unlifted_ID.txt"))
  
  #Remove the SNPs where no matching was found between the two builds in the original dataset corrected for mendel errors
  system(paste0("plink --bfile ", pathData , "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask --exclude ", pathData , "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask_hg19_unlifted_ID.txt --make-bed --recode --out ", pathData , "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask"))
  system(paste0("rm ", pathData , "/mendel_corrected/*~"))

  #Import subject IDs.
  id.fam.suj <- fread("check.oped", skip = 6, header = FALSE)[[1]]
  id.suj <- strsplit(id.fam.suj, "_", fixed = TRUE) %>% sapply(., tail, n=1)
  ids.mat <- data.frame(id.suj, id.fam.suj)

  #map data for GRCh38 to hg19 (GRCh37)
  seq.tout <- read.table(paste0(pathData, "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask.map"))
  seq.tout_hg19 <- read.table(paste0(pathData, "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask_hg19.map"))

  #ped data without mendel errors
  pre.scan <- fread(paste0(pathData, "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask.ped"))
  pre1 <- as.data.frame(pre.scan)
  colnames(pre1)[2] <- "id.suj"

  #Clean the IDs
  pre <- merge(ids.mat, pre1, by="id.suj", all.x=F, all.y=T)
  pre <- pre[,-1]
  pre[,2] <- pre[,1]

  #Remove duplicated subjects or NA values (2 values removed)
  pre <- pre[!is.na(pre[,1]),]
  pre <- pre[!duplicated(as.character(pre[,1])),]

  #Genetic positions
  pos_seq <- setNames(fread(paste0(pathimpu, "/maps_for_hg19/chrom_", chr, "_map.txt"))[,c(1,3)], c("V4", "V3"))
  pos_seq$V1 <- as.character(chr); pos_seq$V2 <- "."
  pos_seq <- pos_seq[,c("V1", "V2", "V3", "V4")]

  #Add dummy SNPs at the start and the at the end of the chromosome
  bidon <- data.frame("V1" = as.character(chr), "V2" = "SNP_bidon", "V3" = 0, "V4" = 0)
  bidon2 <- data.frame("V1" = as.character(chr), "V2" = "SNP_bidon2", "V3" = pos_seq$V3[nrow(pos_seq)]+0.00001, "V4" = max(seq.tout_hg19$V4)+1)
  pos_seq <- rbind(bidon, pos_seq[,1:4], bidon2)

  #Interpolation genetic positions using hg19 physical positions
  #The variants are now sorted by physical position but we need to get them back to their initial order after
  temp <- seq.tout_hg19
  temp$order <- 1:nrow(temp)
  temp <- temp[order(temp$V4),]
  interpole.pos.gen <- function(x){
    which.juste.avant <- max(which(pos_seq$V4<=x))
    pos.phys.juste.avant <- pos_seq$V4[which.juste.avant]
    pos.phys.juste.apres <- pos_seq$V4[which.juste.avant+1]
    pos.gen.juste.avant <- pos_seq$V3[which.juste.avant]
    pos.gen.juste.apres <- pos_seq$V3[which.juste.avant+1]
    out <- pos.gen.juste.avant+(pos.gen.juste.apres-pos.gen.juste.avant)/(pos.phys.juste.apres-pos.phys.juste.avant)*(x-pos.phys.juste.avant)
    return(out)
  }
  temp$pos_gen_approx <- sapply(temp$V4, interpole.pos.gen)

  #Deal with duplicated genetic positions and missing values
  pos.tmp <- temp$pos_gen_approx
  y <- table(pos.tmp)
  z <- c(y[-1], y[length(y)])
  mat <- cbind(y, z, as.numeric(names(y)), as.numeric(names(z)))
  tmp <- apply(mat, 1, function(u){if(u[1]==1) return(u[3]); if(u[1]>1) return(c(u[3],u[3]+((u[4]-u[3])/(u[1]-1)/100)*(1:(u[1]-1))))})
  temp$pos_gen_fin <- as.numeric(unlist(tmp))
  temp <- temp[order(temp$order),]; temp$order <- NULL

  #Import framework markers
  frame_pos <- fread(paste0("framework_pos.txt"))

  #We need to keep the variants that are included in the markers framework only
  which.keep <- which(temp$pos_gen_fin <= max(frame_pos[,1]) & temp$pos_gen_fin >= min(frame_pos[,1]))
  pos_dens <- temp[which.keep, c("V2","pos_gen_fin")]
  temp <- temp[which.keep, ]

  #GIGI2 wants the variants to be sorted by genetic position, we ensure that every output respects this order
  out_order <- order(temp$pos_gen_fin)
  fwrite(pos_dens[out_order,], "dense.pos.txt", row.names=F, quote=F, col.names=F, sep = "\t")

  #Convert allele in numbers. Monomorphic SNPs are coded as 1
  pre_start <- pre[,1:6]
  geno <- apply(pre[, -(1:6)], 2, as.character)
  geno[geno=="0"] <- "Z"
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
    system(paste0("plink -vcf /lustre03/project/6033529/quebec_10x/data/freeze/RV/chr", chr, "_CaG_FINAL_with_mask.vcf.gz --freq --double-id --out ", pathimpu, "/chr", chr, "/CaG"))
    #system(paste0("plink --file ", pathData, "/mendel_corrected/seq_FINAL_mendel_", chr, "_with_mask --freq --out ", pathimpu, "/chr", chr, "/seq"))
    system(paste0("rm ", pathimpu, "/chr", chr, "/*.nosex"))
  }
  freqfile <- fread("CaG.frq", header=T)
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
  write.table(pos_dens[out_order, "pos_gen_fin"], "dense_pos_vf.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

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
}
