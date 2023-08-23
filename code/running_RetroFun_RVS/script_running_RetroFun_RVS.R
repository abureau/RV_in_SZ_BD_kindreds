#module load r/4.2 StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16
library(stringr)
library(GenomicRanges)
library(RetroFunRVS)
options(scipen=999)

#Scripts pour RetroFun-RVS
#1ere étape: Du fait que certaines colonnes ne sont pas formattées correctement pour l'application de 
#RetroFun-RVS, on fait la correspondance entre les fichiers .ped et les données familiales 
#2ème étape: On préprocesse les fichiers .ped en utilisant la fonction agg.genos.by.fam de RetroFun-RVS
#3ème étape: On crée les fichiers d'annotations pour chaque TAD (à automatiser pour concilier les fichiers crées par Jasmin et les fichiers présents dans CRHs_by_TAD)
#4ème étape: Execution de la fonction RetroFun-RVS pour chaque TAD

pathimpu <- "/lustre03/project/6033529/quebec_10x/data/freeze/QC/mendel_corrected"
setwd(pathimpu)

#JR 2023-08-07, Il est important de retirer les variants � MAF==0 et ceux qui rencontre un probleme d'inconsistence (MAF==NA)
genome_results <- data.frame()
missing_TADs <- c()
for(chr in 1:22){
print(paste0("Traitement du chr ", chr))

TADs <- data.table::fread(paste0(pathimpu, "/TADs_list_chr", chr, ".bed"), header = FALSE)

#Pour certains TADs, aucun variant n'�tait trouv�. On les retire.
#JR 2023-08-18 ce bout n'est plus utile puisque Lo�c s'est assur� que ce n'est pas le cas.
fileschr <- list.files(pathimpu)
fileschrTAD <- fileschr[grep(pattern = paste0("chr_", chr, "_TAD.*\\.frq"), x = fileschr)]
fileschrTAD <- str_replace(fileschrTAD, ".frq", "")
TADs_to_keep <- sort(as.numeric(str_split_fixed(fileschrTAD, "_", 8)[,6]))
missing <- setdiff(TADs$V4, TADs_to_keep)
if(length(missing)!=0){missing_TADs <- c(missing_TADs, paste0("chr_", chr, "_TAD_", setdiff(TADs$V4, TADs_to_keep)))}

results_dataframe <- data.frame()
for(TAD in TADs_to_keep){
print(paste0("Traitement du TAD ", TAD))

datafile = paste0("seq_FINAL_chr_", chr, "_TAD_", TAD, "_with_mask")
MAF = data.table::fread(paste0(pathimpu, "/", datafile, ".frq"))

#Load le fichier .ped pour un TAD donné (ici le TAD1 pour le chromosome 22)
#2 doit etre l'allele mineur. Pour le moment, PLINK code l'allele mineur par 1 et le majeur par 2.
pedfile = read.table(paste0(pathimpu, "/", datafile, ".ped"))
pedfile[,c(7:ncol(pedfile))] <- dplyr::case_when(
  pedfile[,c(7:ncol(pedfile))] == 0 ~ 0,
  pedfile[,c(7:ncol(pedfile))] == 1 ~ 2,
  pedfile[,c(7:ncol(pedfile))] == 2 ~ 1)

#Load le fichier de phénotype ou les boucles de consanguinité n'ont pas été retirées (fichier ped2021.ped)
#Conversion en data.frame en vue de faciliter les étapes de manipulations
df.ped.2021 = read.table("/lustre03/project/6033529/quebec_10x/data/freeze/GCbroad_seq_inbred.pre", header=FALSE, sep=" ")
colnames(df.ped.2021) = c("famid","id","fid","mid","sex","affected")
subset.fam = c("103","105","110","115","119","121","124","125","126","129","131","133","151","182","207","210","211","212","217","220","224","228","230","233","234","235","238","255")
df.ped.2021 = df.ped.2021[df.ped.2021$famid%in%subset.fam,]

pedfile.fam.infos = data.frame(pedfile[,1], pedfile[,2], pedfile[,2], pedfile[,2])
colnames(pedfile.fam.infos) = c("famid", "id", "findex", "mindex")
pedfile.fam.infos.full = merge(pedfile.fam.infos, df.ped.2021[,c("famid","id","sex","affected")], by=c("famid","id"), sort = FALSE)

pedfile <- pedfile[pedfile$V2 %in% pedfile.fam.infos.full$id, ]
all(pedfile$V2 == pedfile.fam.infos.full$id)
pedfile[,1:6] = pedfile.fam.infos.full

agg.genos.by.fam = agg.genos.by.fam(pedfile.path=NULL, pedfile=pedfile, correction="none")

#Creation des fichiers d'annotations
#Ici on part directement des fichiers .map

mapfile = read.table(paste0(pathimpu, "/", datafile, ".map"), header=FALSE)
variants = str_split_fixed(mapfile$V2, ":", 4)

GRanges.variants =  GRanges(seqnames=variants[,1], ranges=IRanges(start=as.numeric(variants[,2]),end=as.numeric(variants[,2])))

#Correspondance entre le TAD et les CRHs
#Un fichier par TAD avec les CRHs présents
#Ici les CRHs correspondant au TAD 1 pour le chromosome 22 (à automatiser)

CRHs.by.TAD = read.table(paste0("/lustre03/project/6033529/quebec_10x/data/CRHs_iPSC_neurons/TADs_for_CRHs_overlap_1_TAD/chr", chr, "/chr", chr, "_", TADs$V2[TAD]+1, "_", TADs$V3[TAD], ".bed"), header=TRUE)
GRanges.CRHs.by.TAD = GRanges(seqnames=CRHs.by.TAD$chrom, ranges=IRanges(start=CRHs.by.TAD$chromStart, end=CRHs.by.TAD$chromEnd),name=CRHs.by.TAD$name)

#On split par CRH
split.by.CRH = split(GRanges.CRHs.by.TAD,GRanges.CRHs.by.TAD$name)

#On créée une matrice d'annotation, la première colonne constitue le Burden original

annotation.matrix = matrix(0, ncol=length(split.by.CRH), nrow=length(GRanges.variants))
for(i in 1:length(split.by.CRH)){
  c = countOverlaps(GRanges.variants,split.by.CRH[[i]])
  c[c>1] =1
  annotation.matrix[,i] = c 
}

annotation.matrix = cbind(1, annotation.matrix)

df.annotation = data.frame(annotation.matrix)
colnames(df.annotation) = c("Burden_Original", paste0("CRH",names(split.by.CRH)))

null.with.consanguinity = read.table("/lustre03/project/6033529/quebec_10x/data/CRHs_iPSC_neurons/launch_RetroFunRVS/null_var_seqped2021_with_consanguinity.txt", header=TRUE, sep="\t")
attributes(null.with.consanguinity)$distinguishHomo = TRUE
results = RetroFun.RVS(null.with.consanguinity, agg.genos.by.fam, Z_annot = df.annotation, W = rep(1, nrow(df.annotation)), independence=FALSE)

#Creer la sortie pour chaque TAD
n_CRH <- length(results)-2
names_CRH <- str_split_fixed(names(results)[1:n_CRH], "_", 2)[,2]
scores <- unname(unlist(results[1:n_CRH]))
out <- data.frame("chr" = rep(chr, n_CRH), "TAD" = rep(TAD, n_CRH), "TAD_name" = names_CRH, "score" = scores, "ACAT" = rep(results[[n_CRH+1]], n_CRH), "Fisher" = rep(results[[n_CRH+2]], n_CRH))

results_dataframe <- rbind(results_dataframe, out)
}
saveRDS(results_dataframe, paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/RetroFun.RVS_results_seq_chr", chr, ".RDS"))

genome_results <- rbind(genome_results, results_dataframe)
}
saveRDS(genome_results, paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/RetroFun.RVS_results_seq_all_chromosomes.RDS"))
data.table::fwrite(data.table::data.table(missing_TADs), "/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/empty_TADs.txt")
