library(data.table)
library(dplyr)
pathData <- "/lustre03/project/6033529/quebec_10x/data/freeze/QC"
pathimpu <- "/lustre03/project/6033529/quebec_10x/data/freeze/imputation"

chr_loop <- as.character(1:22)

for(chr in chr_loop){
  dir.create(paste0(pathimpu, "/chr", chr))
  setwd(paste0(pathimpu, "/chr", chr))

  #Import subject IDs from one of the chromosomes files
  id.fam.suj <- fread("check.oped", skip = 6, header = FALSE)[[1]]
  id.suj <- strsplit(id.fam.suj, "_", fixed = TRUE) %>% sapply(., tail, n=1)
  ids.mat <- data.frame(id.suj, id.fam.suj)

  #map data
  seq.tout <- read.table(paste0(pathData, "/seq_FINAL_", chr, "_with_mask.map"))

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
  pre <- merge(ids.mat, pre1, by="id.suj", all.x=F, all.y=T)
  pre <- pre[,-1]
  pre[,2] <- pre[,1]

  #Remove duplicated subjects or NA values (from 349 subjects to 347)
  pre <- pre[!is.na(pre[,1]),]
  pre <- pre[!duplicated(as.character(pre[,1])),]

  #Genetic positions
  pos_seq <- data.frame(fread(paste0("/lustre03/project/6033529/quebec_10x/results/carte_genetique/plink.chr", chr, ".GRCh38.map")))
  pos_seq[,1] <- as.character(pos_seq[,1])

  #Add dummy SNPs at the start and the at the end of the chromosome
  bidon <- data.frame("V1" = as.character(chr), "V2" = "SNP_bidon", "V3" = 0, "V4" = 0)
  bidon2 <- data.frame("V1" = as.character(chr), "V2" = "SNP_bidon2", "V3" = pos_seq$V3[nrow(pos_seq)]+0.00001, "V4" = max(seq.tout$V4)+1)
  pos_seq <- rbind(bidon, pos_seq[,1:4], bidon2)

  #Interpolation
  temp <- seq.tout
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

  #Deal with duplicated genetic positions
  pos.tmp <- temp$pos_gen_approx
  y <- table(pos.tmp)
  z <- c(y[-1], y[length(y)])
  mat <- cbind(y, z, as.numeric(names(y)), as.numeric(names(z)))
  tmp <- apply(mat, 1, function(u){if(u[1]==1) return(u[3]); if(u[1]>1) return(c(u[3],u[3]+((u[4]-u[3])/(u[1]-1)/100)*(1:(u[1]-1))))})
  temp$pos_gen_fin <- as.numeric(unlist(tmp))

  #Import framework markers
  frame_pos <- fread(paste0("framework_pos.txt"))

  #We need to keep the variants that are included in the marqueurs framework only
  which.keep <- which(temp$pos_gen_fin <= max(frame_pos[,1]) & temp$pos_gen_fin >= min(frame_pos[,1]))
  pos_dens <- temp[which.keep,c("V2","pos_gen_fin")]
  fwrite(pos_dens, "dense.pos.txt", row.names=F, quote=F, col.names=F, sep = "\t")

  #Convert allele in numbers. Monomorphic SNPs are coded as 1
  geno <- pre[, -(1:6)] %>% apply(., 2, as.character)
  geno[geno=="0"] <- "Z"
  tmp <- matrix(as.vector(as.matrix(geno)),ncol=ncol(geno)/2)
  geno.new <- matrix(as.vector(as.matrix(apply(tmp,2,function(x) as.numeric(as.factor(x))))),ncol=ncol(geno))
  geno.new[geno.new==3] <- 0
  pre <- data.frame(pre[,1:6],geno.new)

  #Remove variants not in framework marker
  which.geno.keep <- c(2*which.keep-1, 2*which.keep)
  which.geno.keep <- which.geno.keep[order(which.geno.keep)]
  geno <- pre[,7:ncol(pre)]
  geno <- geno[,which.geno.keep]
  pre_new <- data.frame(pre[,1:6],geno)
  write(apply(pre_new[,c(1,7:dim(pre_new)[2])], 1, paste, collapse=" "), "dense_geno.txt")

  #MAF extraction computed in CARTaGENE
  if(!file.exists("CaG.frq")){
    system(paste0("plink -vcf /lustre03/project/6033529/quebec_10x/data/freeze/RV/chr", chr, "_CaG_FINAL_with_mask.vcf.gz --freq --double-id --out ", pathimpu, "/chr", chr, "/CaG"))
    system(paste0("rm ", pathimpu, "/chr", chr, "/CaG.nosex"))
  }
  freqfile <- fread("CaG.frq", header=T)
  freqfile$SNP <- gsub("_", ":", freqfile$SNP)
  SNP <- as.character(freqfile$SNP)
  freq <- as.numeric(freqfile$MAF)

  #We keep all SNPs, but after the merge, we set a freq of 0.001 to missing values or values of 0
  id.snp <- 1:nrow(pos_dens)
  pos_dens <- data.frame(pos_dens,id.snp)
  freq_dense <- merge(pos_dens, data.frame(SNP, freq), by.x = 'V2', by.y = 'SNP', all.x = TRUE, all.y = FALSE)
  freq_dense <- freq_dense[order(freq_dense$id.snp),]
  freq_dense <- freq_dense[, -which(colnames(freq_dense)=="id.snp")]

  #Verifications
  all(as.character(freq_dense[,1]) == as.character(pos_dens$V2))
  freq_dense[, ncol(freq_dense)] <- as.numeric(freq_dense[, ncol(freq_dense)])
  sum(freq_dense[, ncol(freq_dense)] == 0, na.rm=T)
  sum(is.na(freq_dense[, ncol(freq_dense)]))
  freq_dense[is.na(freq_dense[, ncol(freq_dense)]) | freq_dense[, ncol(freq_dense)] == 0, ncol(freq_dense)] <- 0.001
  write.table(pos_dens[, 2], "dense_pos_vf.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

  freq <- freq_dense$freq
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
  write.table(data.frame(freq1, freq2), "dense_freq.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

  #Final marker list for GIGI2
  write.table(freq_dense[, 1], "liste_dense.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
}
