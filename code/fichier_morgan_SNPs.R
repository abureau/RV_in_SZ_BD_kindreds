#Prepare the inputs for gl_auto.
#Needs to be executed on Arcturus.

for(chr in 1:22){
pathData <- "/home/jricard/lg_auto_data"
setwd(pathData)
dir.create(paste0("chr", chr))
setwd(paste0("chr", chr))
system(paste0("cp ../glauto.par ./"))

#Create a file where only the remaining variants MAF are found.
freq <- fread(paste0(pathData, "/GSA_genotypes_clumped_chr", chr, ".frq"))
bim <- fread(paste0(pathData, "/GSA_genotypes_clumped_chr", chr, ".bim"))
bim_not_clumped <- fread(paste0(pathData, "/genotypes-chr", chr, "-Commun.bim"))
map <- fread(paste0(pathData, "/chrom_", chr, "_map.txt"))
to_keep_cl <- which(bim_not_clumped$V2 %in% bim$V2)
map_cl <- map[to_keep_cl,]
fwrite(map_cl, paste0(pathData, "/GSA_genotypes_genetic_distance_clumped_chr", chr, ".txt"), col.names = FALSE, row.names = FALSE, sep = "\t")

#Import the complete pedigree file and add the following missing informations.
pre <- read.table("/home/jricard/lg_auto_data/GCbroad.pre")
pre[pre[, 1] == "245" & pre[, 3] == 2162, 3] <- 0 
pre[pre[, 1] == "245" & pre[, 4] == 2163, 4] <- 0 
pre[pre[, 3] %in% c(0, 3186, 2314, 2381), 3] <- 0 
pre[pre[, 4] %in% c(0, 3187, 2313, 2380), 4] <- 0

#Import the actual .ped.
ped <- read.table(paste0("/home/jricard/lg_auto_data/GSA_genotypes_clumped_mendel_chr", chr, ".ped"))

#Combine subject ID with family ID.
pre[,2]=paste(pre[,1],pre[,2],sep="_")
pere=pre[,3]
mere=pre[,4]
pre[,3]=paste(pre[,1],pere,sep="_")
pre[,4]=paste(pre[,1],mere,sep="_")
pre[pere==0,3]="0"
pre[mere==0,4]="0"
pedigree <- pre[,2:5]

#Create a dummy trait.
pedigree$V6 <- 0

#Create outputs.
cat("input pedigree size ",nrow(pedigree), "\ninput pedigree record names 3 integers 2\ninput pedigree record trait 1 integer 2\n*****\n", file="GCbroad_v2.ped")
write.table(pedigree,"GCbroad_v2.ped", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

#Importe the actual .map and the genetic positions.
op <- options(scipen=999)
position = read.table(paste0("/home/jricard/lg_auto_data/GSA_genotypes_clumped_mendel_chr", chr, ".map"))
genpos <- read.table(paste0("/home/jricard/lg_auto_data/GSA_genotypes_genetic_distance_clumped_chr", chr, ".txt"))
position[, 3] <- genpos[, 3]
names(position)[1:3] <- c("chr", "rs", "positionGenetique")
which.keep <- which(!duplicated(position$positionGenetique))
which.geno.keep <- c(2*which.keep-1, 2*which.keep)
which.geno.keep <- which.geno.keep[order(which.geno.keep)]
position <- position[which.keep, ]

#We set every MAF to 0.5.
nb_freq <- rep("0.5 0.5", nrow(position))

#Add the text that needs to be written in the frequency file.
#Create outputs.
text2 <- paste("set markers",1:length(nb_freq),"allele freq",nb_freq)
write.table(text2,"framework_freq.txt", quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(position$positionGenetique, "framework_gcbroad.map_af_geno", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(position$positionGenetique, "framework_pos.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#Convert alleles to numeric values
geno <- ped[, -(1:6)]
geno <- geno[, which.geno.keep]
geno <- apply(geno, 2, as.character)
geno[geno=="0"] <- "Z"
tmp <- matrix(as.vector(as.matrix(geno)), ncol = ncol(geno)/2)
geno.new <- matrix(as.vector(as.matrix(apply(tmp, 2, function(x) as.numeric(as.factor(x))))), ncol = ncol(geno))
geno.new[geno.new == 3] <- 0

#Create outputs.
text2 <- paste("set markers", 1:length(nb_freq), "allele freq", nb_freq)
write.table(as.matrix(c("map markers position\n", position$positionGenetique, text2, paste0("\nset markers data ", length(position$positionGenetique)))),
            file = "framework_gcbroad.map_af_geno", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(cbind(paste0(ped[, 1], "_", ped[, 2]), geno.new), file = "framework_gcbroad.map_af_geno", quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)

}



