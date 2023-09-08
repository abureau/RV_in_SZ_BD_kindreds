#module load r/4.2 StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(stringr, quietly=T)
library(GenomicRanges, quietly=T)
library(RetroFunRVS, quietly=T)
options(scipen=999)

#Scripts pour RetroFun-RVS
#1ere étape: Du fait que certaines colonnes ne sont pas formattées correctement pour l'application de
#RetroFun-RVS, on fait la correspondance entre les fichiers .ped et les données familiales
#2ème étape: On préprocesse les fichiers .ped en utilisant la fonction agg.genos.by.fam de RetroFun-RVS
#3ème étape: On crée les fichiers d'annotations pour chaque TAD (à automatiser pour concilier les fichiers crées par Jasmin et les fichiers présents dans CRHs_by_TAD)
#4ème étape: Execution de la fonction RetroFun-RVS pour chaque TAD

pathimpu <- args[1]
output <- args[2]

print(pathimpu)
print(output)

setwd(pathimpu)

#JR 2023-08-07, Il est important de retirer les variants ▒ MAF==0 et ceux qui rencontre un probleme d'inconsistence (MAF==NA)
genome_results <- data.frame()
missing_TADs <- c()

df.ped.2021 = read.table("/lustre03/project/6033529/quebec_10x/data/freeze/GCbroad_seq_inbred.pre", header=FALSE, sep=" ")
colnames(df.ped.2021) = c("famid","id","fid","mid","sex","affected")
subset.fam = c("103","105","110","115","119","121","124","125","126","129","131","133","151","182","207","210","211","212","217","220","224","228","230","233","234","235","238","255")
df.ped.2021 = df.ped.2021[df.ped.2021$famid%in%subset.fam,]

#null.with.consanguinity = read.table("/lustre03/project/6033529/quebec_10x/data/CRHs_iPSC_neurons/launch_RetroFunRVS/null_var_seqped2021_with_consanguinity.txt", header=TRUE, sep="\t")
#attributes(null.with.consanguinity)$distinguishHomo = TRUE
load("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS_cryptique/expected.variance.consanguinity.cryptique.seq.RData")

for(chr in 1:22){
  print(paste0("Traitement du chr ", chr))
  
  list_ped_by_chr <- list.files(pathimpu,pattern=paste0("^","chr",chr,"_",".*txt"), full.names = TRUE)
  
  if(length(list_ped_by_chr)==0){
    print(paste0("Pas de CRHs pour dans le chromosome"), chr)
    next
  }
  
  for(file in list_ped_by_chr){
    
    
    CRH <- gsub(".*[_]([^.]+)[.].*", "\\1", file)
      
    print(paste0("Traitement du CRH:",CRH))
      
    results_dataframe <- data.frame()
    df.CRH <- read.table(file, header = FALSE)
    colnames(df.CRH) <- c("chrom", "startChrom", "endChrom", "name")
      
    datafile = paste0(pathimpu, "/seq_FINAL_chr_", chr, "_CRH_", CRH, "_with_mask")
      
    pedfile.path = file.path(paste0(datafile,".ped"))

    if(file.exists(pedfile.path)) pedfile = read.table(paste0(datafile, ".ped"))
    else next

    if(ncol(pedfile) <7){
      print("Pas de variant dans ce CRH")
      next
    }
      
    for(i in 7:ncol(pedfile)){
      pedfile[,i] <- dplyr::case_when(
        pedfile[,i] == 0 ~ 0,
        pedfile[,i] == 1 ~ 2,
        pedfile[,i] == 2 ~ 1)
    }
      
      
    pedfile.fam.infos = data.frame(pedfile[,1], pedfile[,2], pedfile[,2], pedfile[,2])
    colnames(pedfile.fam.infos) = c("famid", "id", "findex", "mindex")
    pedfile.fam.infos.full = merge(pedfile.fam.infos, df.ped.2021[,c("famid","id","sex","affected")], by=c("famid","id"), sort = FALSE)
      
    pedfile <- pedfile[pedfile$V2 %in% pedfile.fam.infos.full$id, ]
    all(pedfile$V2 == pedfile.fam.infos.full$id)
    pedfile[,1:6] = pedfile.fam.infos.full
    
     
    agg.genos.by.fam = tryCatch(agg.genos.by.fam(pedfile.path=NULL, pedfile=pedfile, correction="none"),
	error = function(e) e)
    if(inherits(agg.genos.by.fam, "error")) next
        
    mapfile = read.table(paste0(datafile,".map"), header=FALSE)
    variants = str_split_fixed(mapfile$V2, ":", 4)
      
    GRanges.variants =  GRanges(seqnames=variants[,1], ranges=IRanges(start=as.numeric(variants[,2]),end=as.numeric(variants[,2])))
    GRanges.CRH.by.chrom = GRanges(seqnames=df.CRH$chrom, ranges=IRanges(start=df.CRH$startChrom, end=df.CRH$endChrom),name=df.CRH$name)
      
    annotation.matrix = matrix(1, ncol=1, nrow=length(GRanges.variants))
      
    df.annotation = data.frame(annotation.matrix)
    colnames(df.annotation) = c("Burden_Original")
      
    results = RetroFun.RVS(expected.variance.consanguinity.cryptique.seq, agg.genos.by.fam, Z_annot = annotation.matrix, W = rep(1, nrow(annotation.matrix)), independence=FALSE)
    print(results)
    n_CRH <- 1
    scores <- unname(unlist(results[1]))
    out <- data.frame("chr" = chr, "CRH_name" = CRH, "score" = scores, "ACAT" = results[,2], "Fisher" = results[,3])
      
    results_dataframe <- rbind(results_dataframe, out)
    
    saveRDS(results_dataframe, paste0(output,"/RetroFun.RVS_results_seq_chr", chr,"CRH",CRH, ".RDS"))
      
    genome_results <- rbind(genome_results, results_dataframe)
  }
    
}

saveRDS(genome_results, paste0(output,"/RetroFun.RVS_results_CRHs_seq_all_chromosomes.RDS"))
  


