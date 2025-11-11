#---- RetroFun-RVS ----
#To run using sbatch. Just change the phenotype if needed.
#pheno=GCbr; sbatch --job-name=${pheno}_retrofun --export=pheno=${pheno} /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples_without_mask/call/envoi_RetroFun_RVS.sh

library(stringr); library(GenomicRanges); library(RetroFunRVS); library(dplyr); library(kinship2); library("RetroFunRVS");library(foreach); library(data.table)
options(scipen=999)

args <- commandArgs(TRUE)
(pheno <- as.character(args[[1]]))
(consanguinity <- as.logical(args[[2]]))
(exons <- as.logical(args[[3]]))
(remove_singletons <- as.logical(args[[4]]))

#Lancer RetroFun-RVS
#1ere ?tape: Du fait que certaines colonnes ne sont pas formatt?es correctement pour l'application de RetroFun-RVS, on fait la correspondance entre les fichiers .ped et les donn?es familiales 
#2?me ?tape: On pr?processe les fichiers .ped en utilisant la fonction agg.genos.by.fam de RetroFun-RVS
#3?me ?tape: On cr?e les fichiers d'annotations pour chaque TAD (? automatiser pour concilier les fichiers cr?es par Jasmin et les fichiers pr?sents dans CRHs_by_TAD)
#4?me ?tape: Execution de la fonction RetroFun-RVS pour chaque TAD

pathAB_ped <- "/lustre09/project/6033529/schizo/data_AB/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped/"
path_retrofun <- "/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag_without_mask/RetroFunRVS/"
#path_retrofun <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag_without_mask/RetroFunRVS"
setwd(paste0(path_retrofun, "/TADs"))
# Répertoire pour les fichiers temporaires
path_tmp <- "/scratch/bureau/data"

#CRH membership problem
membership_equi <- fread(paste0(path_retrofun, "/objets_ped/CRH_problem_membership_equivalence.txt"))

#Function to load pedigrees easily.
loadRData <- function(file_name){load(file_name); get(ls()[ls() != "file_name"])}

#Function for CHRs overlapping 1 TAD.
RetroFun.RVS_run <- function(pheno, with_exons, consanguinity, maxvar = 300, remove_singletons = FALSE){
  #pheno, string, add the phenotype name.
  #with_exons, logical, TRUE if we add the exons to the CRHs for the analysis.
  #consanguinity, logical, TRUE if we consider consanguinity loops in the analysis.
  #maxvar, integer, maximum normal of variants to test together.
  #remove_singletons, logical, TRUE if we need to remove singletons from the analysis.
  
  if(remove_singletons){
    #singletons <- fread("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag_without_mask/imputation_comb/merged_with_seq/freq/impute5_gigi2_combined_seq_fam_singleton.snplist", header = FALSE)
    singletons <- fread("/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag_without_mask/imputation_comb/merged_with_seq/freq/impute5_gigi2_combined_seq_fam_singleton.snplist", header = FALSE)
    out_remove_singletons <- "_without_singletons"
  } else {
    out_remove_singletons <- ""
  }

  #Load the phenotype files with consanguinity loops.
  #Affected must be coded by a 2.
  df.ped.2021 <- loadRData(paste0(path_retrofun, "/objets_ped/ped", pheno, "_orig.RData"))
  df.ped.2021 <- data.frame(df.ped.2021$famid, df.ped.2021$id, df.ped.2021$sex, df.ped.2021$affected+1)
  colnames(df.ped.2021) <- c("famid","id", "sex","affected")

  genome_results <- data.frame()
  missing_TADs <- c()
  for(chr in 1:22){
    print(paste0("chr ", chr))
    TADs <- fread(paste0(path_retrofun, "/TADs/TADs_list_chr", chr, ".bed"), header = FALSE)
    
    #For some TADs, no variant was found. We remove them. It must not be the case, but I keep this code to be sure.
    fileschr <- list.files(paste0(path_retrofun, "/TADs"))
    fileschrTAD <- fileschr[grep(pattern = paste0("chr", chr, "_TAD_.*\\.frq"), x = fileschr)]
    fileschrTAD <- stringr::str_replace(fileschrTAD, ".frq", "")
    TADs_to_keep <- sort(as.numeric(str_split_fixed(fileschrTAD, "_", 9)[,9]))
    missing <- setdiff(TADs$V4, TADs_to_keep)
    if(length(missing)!=0){missing_TADs <- c(missing_TADs, paste0("chr_", chr, "_TAD_", setdiff(TADs$V4, TADs_to_keep)))}
    
    results_dataframe <- data.frame()
    for(TAD in TADs_to_keep){
      print(paste0("TAD ", TAD))
      
      datafile <- paste0("impute5_gigi2_combined_seq_RV_FINAL_chr", chr, "_TAD_", TAD)
      if(with_exons){
        #If the addition of exons had an impact on the .ped, import the modified one.
        if(file.exists(paste0(path_retrofun, "/TADs/with_exons/", datafile, ".ped"))){datafile <- paste0("with_exons/", datafile)}
      }
      #MAF <- fread(paste0(path_retrofun, "/TADs/", datafile, ".frq"))
      
      #Load the .map and .ped of the TAD
      #Minor alleles must be coded as 2. PLINK codes minor to 1 and major to 2, so we must change them.
      mapfile <- read.table(paste0(path_retrofun, "/TADs/", datafile, ".map"), header=FALSE)
      pedfile <- read.table(paste0(path_retrofun, "/TADs/", datafile, ".ped"))
      #Remove singletons if needed
      if(remove_singletons){
        not_singletons <- which(!mapfile$V2 %in% singletons$V1)
        mapfile <- mapfile[not_singletons,]
        pedfile <- pedfile[,c(1:6, sort( c(5+(not_singletons*2), 6+(not_singletons*2)) ))]
      }

      for(i in 7:ncol(pedfile)){pedfile[,i] <- case_when(pedfile[,i] == 0 ~ 0, pedfile[,i] == 1 ~ 2, pedfile[,i] == 2 ~ 1)}
      
      #Remove 2620b which is a duplicate of 2620 in our data.
      pedfile <- pedfile[pedfile[,2] != "2620b",]
      
      #add loaded phenotypes
      pedfile.fam.infos <- setNames(data.frame(pedfile[,1], pedfile[,2], pedfile[,3], pedfile[,4]), c("famid", "id", "findex", "mindex"))
      pedfile.fam.infos.full <- merge(x = pedfile.fam.infos, y = df.ped.2021[,c("famid","id","sex","affected")], by = c("famid","id"), sort = FALSE, all.x = TRUE)
      
      #Associate the phenotype in the .ped and keep the affected families only.
      subset.fam <- pedfile.fam.infos.full %>% group_by(famid) %>% summarise(n_affected = sum(affected==2)) %>% filter(n_affected!=0) %>% select(famid) %>% as.vector()
      df.ped.2021 <- df.ped.2021[df.ped.2021$famid %in% subset.fam$famid,]
      pedfile.fam.infos.full <- pedfile.fam.infos.full[pedfile.fam.infos.full$famid %in% subset.fam$famid,]

      #Check
      pedfile <- pedfile[pedfile$V2 %in% pedfile.fam.infos.full$id, ]
      all(pedfile$V2 == pedfile.fam.infos.full$id)
      pedfile[,1:6] <- pedfile.fam.infos.full
      
      #Adjust the pedigree for the 3 problematic families in agg.genos.by.fam.
      if(consanguinity){
        null_name <- paste0("expected.variance.consanguinity.cryptique.",pheno,".rds")
        correction <- "none"
        out_consanguinity <- "with_consanguinity"
      } else {
        null_name <- paste0("expected.variance.",pheno,".rds")
        correction <- "replace"
        out_consanguinity <- "without_consanguinity"
      }
      if (pheno %in% c("GCbr","GCna"))
      {
      null <- readRDS(paste0(path_retrofun, "/objets_ped/", null_name))
      fam_split_119 <- readRDS(paste0(path_retrofun, "/objets_ped/fam119splitted.rds"))
      pedfile$V1[pedfile$V1 == "119" & pedfile$V2 %in% fam_split_119$id[fam_split_119$fam=="119-1"]] <- "119-1"
      pedfile$V1[pedfile$V1 == "119" & pedfile$V2 %in% fam_split_119$id[fam_split_119$fam=="119-2"]] <- "119-2"
      fam_split_131 <- readRDS(paste0(path_retrofun, "/objets_ped/fam131splitted.rds"))
      pedfile$V1[pedfile$V1 == "131" & pedfile$V2 %in% fam_split_131$id[fam_split_131$fam=="131-1"]] <- "131-1"
      pedfile$V1[pedfile$V1 == "131" & pedfile$V2 %in% fam_split_131$id[fam_split_131$fam=="131-2"]] <- "131-2"
      fam_split_255 <- readRDS(paste0(path_retrofun, "/objets_ped/fam255splitted.rds"))
      pedfile$V1[pedfile$V1 == "255" & pedfile$V2 %in% fam_split_255$id[fam_split_255$fam=="255-1"]] <- "255-1"
      pedfile$V1[pedfile$V1 == "255" & pedfile$V2 %in% fam_split_255$id[fam_split_255$fam=="255-2"]] <- "255-2"
      }
      else null <- readRDS(paste0(pathAB_ped, null_name))
      
      #Create annotation files using the .map.
      variants <- str_split_fixed(mapfile$V2, ":", 4)
      GRanges.variants <-  GRanges(seqnames=variants[,1], ranges=IRanges(start=as.numeric(variants[,2]),end=as.numeric(variants[,2])))
      
      #TAD and CRHs matching
      if(with_exons){
        CRHs.by.TAD <- read.table(paste0(path_retrofun, "/TADs/with_exons/TADs_for_CRHs_overlap_1_TAD/chr", chr, "/chr", chr, "_", TADs$V2[TAD]+1, "_", TADs$V3[TAD], ".bed"), header=FALSE)
        if(CRHs.by.TAD[1,1] == "chrom"){CRHs.by.TAD <- CRHs.by.TAD[2:nrow(CRHs.by.TAD),]} #It appears that sometimes, the header is part of the dataframe, remove the line.
        colnames(CRHs.by.TAD) <- c("chrom", "chromStart", "chromEnd", "name")
        out_exons <- "CRHs_with_exons"
      } else {
        #CRHs.by.TAD <- read.table(paste0("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/TADs_for_CRHs_overlap_1_TAD/chr", chr, "/chr", chr, "_", TADs$V2[TAD]+1, "_", TADs$V3[TAD], ".bed"), header=TRUE)
        CRHs.by.TAD <- read.table(paste0("/lustre09/project/6033529/schizo/data/WGS_bs_2022/liftover_hg38_executable/TADs_for_CRHs_overlap_1_TAD/chr", chr, "/chr", chr, "_", TADs$V2[TAD]+1, "_", TADs$V3[TAD], ".bed"), header=TRUE)
        out_exons <- "CRHs_only"
      }
      
      #Split by CRH
      GRanges.CRHs.by.TAD <- GRanges(seqnames=CRHs.by.TAD$chrom, ranges=IRanges(start=CRHs.by.TAD$chromStart, end=CRHs.by.TAD$chromEnd),name=CRHs.by.TAD$name)
      split.by.CRH <- split(GRanges.CRHs.by.TAD,GRanges.CRHs.by.TAD$name)
      
      #Annotation matrix, the first column is the Burden
      annotation.matrix <- matrix(0, ncol=length(split.by.CRH), nrow=length(GRanges.variants))
      for(i in 1:length(split.by.CRH)){
        c <- countOverlaps(GRanges.variants,split.by.CRH[[i]])
        c[c>1] <- 1
        annotation.matrix[,i] <- c 
      }
      
      annotation.matrix = cbind(1, annotation.matrix)
      agg.genos.by.fam <- agg.genos.by.fam(pedfile.path=NULL, pedfile=pedfile, Z_annot=annotation.matrix, correction = correction)
      
      nvarTAD <- length(agg.genos.by.fam$index_variants)
      cat(nvarTAD,"\n")
      nw <- (nvarTAD-1)%/%maxvar + 1
      
      # Si le nombre de variants > maxvar, on decoupe en fenetre de maxvar
      if (nvarTAD > maxvar) {
        # Nombre de fen?tres
        wmat = matrix(0,ncol = nw, nrow=length(GRanges.variants))
        for (w in 1:(nw-1)) {
          wmat[agg.genos.by.fam$index_variants[maxvar*(w-1)+1]:agg.genos.by.fam$index_variants[maxvar*w],w] = 1
        }
        wmat[agg.genos.by.fam$index_variants[maxvar*(nw-1)+1]:nrow(wmat),nw] = 1
        # Ici il faut retirer la colonne de 1 initiale
        annotation.matrix = cbind(wmat, annotation.matrix[,-1])
        # On recalcule le nombre de variants par fenetre
        agg.genos.by.fam <- agg.genos.by.fam(pedfile.path=NULL, pedfile=pedfile, Z_annot=annotation.matrix, exclude_annot = 1:nw, correction = correction)
      }
      
      df.annotation <- setNames(data.frame(annotation.matrix), c(paste0("Burden_",1:nw), paste0("CRH",names(split.by.CRH))))
      results <- RetroFun.RVS(null, agg.genos.by.fam, Z_annot = df.annotation, W = rep(1, nrow(df.annotation)), independence=FALSE)
      
      #Create the by TAD output.
      n_CRH <- length(results)-2
      names_CRH <- str_split_fixed(names(results)[1:n_CRH], "_", 2)[,2]
      scores <- unname(unlist(results[1:n_CRH]))
      out <- data.frame("chr" = rep(chr, n_CRH), "TAD" = rep(TAD, n_CRH), "TAD_name" = names_CRH, "score" = scores, "ACAT" = rep(results[[n_CRH+1]], n_CRH), "Fisher" = rep(results[[n_CRH+2]], n_CRH))
      results_dataframe <- rbind(results_dataframe, out)
    }
    genome_results <- rbind(genome_results, results_dataframe)
  }
  #saveRDS(genome_results, paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/without_mask/RetroFunRVS_results_seq_all_chromosomes_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
  #fwrite(data.table(missing_TADs), paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/without_mask/empty_TADs_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".txt"))
  saveRDS(genome_results, paste0("/lustre09/project/6033529/schizo/results_AB/RetroFunRVS/WGS_bs_2022_500samples/without_mask/RetroFunRVS_results_seq_all_chromosomes_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
  data.table::fwrite(data.table::data.table(missing_TADs), paste0("/lustre09/project/6033529/schizo/results_AB/RetroFunRVS/WGS_bs_2022_500samples/without_mask/empty_TADs_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".txt"))
}

#Function for CHRs overlapping 0 or 2 TADs.
RetroFun.RVS.overlap02_run <- function(pheno, with_exons, consanguinity, remove_singletons = FALSE){
  #pheno, string, add the phenotype name.
  #with_exons, logical, TRUE if we add the exons to the CRHs for the analysis.
  #consanguinity, logical, TRUE if we consider consanguinity loops in the analysis.
  #remove_singletons, logical, TRUE if we need to remove singletons from the analysis.
  
  if(remove_singletons){
    #singletons <- fread("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag_without_mask/imputation_comb/merged_with_seq/freq/impute5_gigi2_combined_seq_fam_singleton.snplist", header = FALSE)
    singletons <- fread("/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag_without_mask/imputation_comb/merged_with_seq/freq/impute5_gigi2_combined_seq_fam_singleton.snplist", header = FALSE)
    out_remove_singletons <- "_without_singletons"
  } else {
    out_remove_singletons <- ""
  }

  #Load the phenotype files with consanguinity loops.
  #Affected must be coded by a 2.
  df.ped.2021 <- loadRData(paste0(path_retrofun, "/objets_ped/ped", pheno, "_orig.RData"))
  df.ped.2021 <- data.frame(df.ped.2021$famid, df.ped.2021$id, df.ped.2021$sex, df.ped.2021$affected+1)
  colnames(df.ped.2021) <- c("famid","id", "sex","affected")

  #Path to the data
  for(overlap in c(0,2)){
    path_data <- paste0(path_retrofun, "/TADs/overlap_", overlap)
    if(with_exons){path_data <- paste0(path_data, "/with_exons")}
    
    genome_results <- data.frame()
    files_path_data <- list.files(path_data, pattern = ".txt$")
    for(chr in 1:22){
      print(paste0("chr ", chr, " for the case with ", overlap, " overlap"))
      list_ped_by_chr <- files_path_data[grep(paste0("chr", chr, "_"), files_path_data)]
      if(length(list_ped_by_chr)==0){print(paste0("No CRHs for chromosome ", chr));next}
      
      for(file in list_ped_by_chr){
        CRH <- gsub(".*[_]([^.]+)[.].*", "\\1", file)
        print(paste0("CRH ", CRH))
        df.CRH <- read.table(paste0(path_data, "/" , file), header = FALSE)
        if(df.CRH[1,1] == "chrom"){df.CRH <- df.CRH[2:nrow(df.CRH),]} #It appears that sometimes, the header is part of the dataframe, remove the line.
        df.CRH[,2] <- as.numeric(df.CRH[,2]); df.CRH[,3]<- as.numeric(df.CRH[,3])
        if(is.null(df.CRH$name)){df.CRH$name <- as.numeric(CRH)}
        colnames(df.CRH) <- c("chrom", "startChrom", "endChrom", "name")
        pedfile.path <- file.path(paste0(path_data, "/impute5_gigi2_combined_seq_RV_FINAL_", gsub(".txt", "", file), ".ped"))
        if(file.exists(pedfile.path)){pedfile = read.table(pedfile.path)} else {next}
        mapfile <- read.table(paste0(path_data, "/impute5_gigi2_combined_seq_RV_FINAL_", gsub(".txt", "", file), ".map"), header=FALSE)
        if(ncol(pedfile) <7){ print("0 variant in this CRH"); next}
        #Remove singletons if needed
        if(remove_singletons){
          not_singletons <- which(!mapfile$V2 %in% singletons$V1)
          mapfile <- mapfile[not_singletons,]
          pedfile <- pedfile[,c(1:6, sort( c(5+(not_singletons*2), 6+(not_singletons*2)) ))]
        }

        for(i in 7:ncol(pedfile)){pedfile[,i] <- case_when(pedfile[,i] == 0 ~ 0, pedfile[,i] == 1 ~ 2, pedfile[,i] == 2 ~ 1)}
        
        #Remove 2620b which is a duplicate of 2620 in our data.
        pedfile <- pedfile[pedfile[,2] != "2620b",]
        
        #Add loaded phenotypes
        pedfile.fam.infos <- setNames(data.frame(pedfile[,1], pedfile[,2], pedfile[,3], pedfile[,4]), c("famid", "id", "findex", "mindex"))
        pedfile.fam.infos.full <- merge(x = pedfile.fam.infos, y = df.ped.2021[,c("famid","id","sex","affected")], by = c("famid","id"), sort = FALSE, all.x = TRUE)
        
        #Associate the phenotype in the .ped and keep the affected families only.
        subset.fam <- pedfile.fam.infos.full %>% group_by(famid) %>% summarise(n_affected = sum(affected==2)) %>% filter(n_affected!=0) %>% select(famid) %>% as.vector()
        df.ped.2021 <- df.ped.2021[df.ped.2021$famid %in% subset.fam$famid,]
        pedfile.fam.infos.full <- pedfile.fam.infos.full[pedfile.fam.infos.full$famid %in% subset.fam$famid,]
        #Check
        pedfile <- pedfile[pedfile$V2 %in% pedfile.fam.infos.full$id, ]
        all(pedfile$V2 == pedfile.fam.infos.full$id)
        pedfile[,1:6] <- pedfile.fam.infos.full
        
        #Adjust the pedigree for the 3 problematic families in agg.genos.by.fam.
        if(consanguinity){
          null_name <- paste0("expected.variance.consanguinity.cryptique.",pheno,".rds")
          correction <- "none"
          out_consanguinity <- "with_consanguinity"
        } else {
          null_name <- paste0("expected.variance.",pheno,".rds")
          correction <- "replace"
          out_consanguinity <- "without_consanguinity"
        }
        if(with_exons){out_exons <- "CRHs_with_exons"} else {out_exons <- "CRHs_only"}
        
        if (pheno %in% c("GCbr","GCna"))
        {
        null <- readRDS(paste0(path_retrofun, "/objets_ped/", null_name))
        fam_split_119 <- readRDS(paste0(path_retrofun, "/objets_ped/fam119splitted.rds"))
        pedfile$V1[pedfile$V1 == "119" & pedfile$V2 %in% fam_split_119$id[fam_split_119$fam=="119-1"]] <- "119-1"
        pedfile$V1[pedfile$V1 == "119" & pedfile$V2 %in% fam_split_119$id[fam_split_119$fam=="119-2"]] <- "119-2"
        fam_split_131 <- readRDS(paste0(path_retrofun, "/objets_ped/fam131splitted.rds"))
        pedfile$V1[pedfile$V1 == "131" & pedfile$V2 %in% fam_split_131$id[fam_split_131$fam=="131-1"]] <- "131-1"
        pedfile$V1[pedfile$V1 == "131" & pedfile$V2 %in% fam_split_131$id[fam_split_131$fam=="131-2"]] <- "131-2"
        fam_split_255 <- readRDS(paste0(path_retrofun, "/objets_ped/fam255splitted.rds"))
        pedfile$V1[pedfile$V1 == "255" & pedfile$V2 %in% fam_split_255$id[fam_split_255$fam=="255-1"]] <- "255-1"
        pedfile$V1[pedfile$V1 == "255" & pedfile$V2 %in% fam_split_255$id[fam_split_255$fam=="255-2"]] <- "255-2"
        }
        else null <- readRDS(paste0(pathAB_ped, null_name))
        
        variants <- str_split_fixed(mapfile$V2, ":", 4)
        GRanges.variants <-  GRanges(seqnames=variants[,1], ranges=IRanges(start=as.numeric(variants[,2]),end=as.numeric(variants[,2])))
        GRanges.CRH.by.chrom <- GRanges(seqnames=df.CRH$chrom, ranges=IRanges(start=df.CRH$startChrom, end=df.CRH$endChrom),name=df.CRH$name)
        
        annotation.matrix <- matrix(1, ncol=1, nrow=length(GRanges.variants))
        agg.genos.by.fam <- tryCatch(agg.genos.by.fam(pedfile.path=NULL, pedfile=pedfile, Z_annot=annotation.matrix, correction=correction), error = function(e) e)
        if(inherits(agg.genos.by.fam, "error")) next
        df.annotation <- data.frame(annotation.matrix)
        colnames(df.annotation) <- c("Burden_Original")
        
        results <- RetroFun.RVS(null, agg.genos.by.fam, Z_annot = annotation.matrix, W = rep(1, nrow(annotation.matrix)), independence=FALSE)
        scores <- unname(unlist(results[1]))
        out <- data.frame("chr" = chr, "CRH_name" = CRH, "result" = scores)
        genome_results <- rbind(genome_results, out)
      }
    }
    #saveRDS(genome_results, paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/without_mask/RetroFunRVS_results_seq_all_chromosomes_overlap_", overlap, "_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
    saveRDS(genome_results, paste0("/lustre09/project/6033529/schizo/results_AB/RetroFunRVS/WGS_bs_2022_500samples/without_mask/RetroFunRVS_results_seq_all_chromosomes_overlap_", overlap, "_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
  }
}


#Function for genes from the litterature
RetroFun.RVS.genes.litt_run <- function(pheno, with_exons, strict = FALSE, consanguinity, remove_singletons = FALSE){
  #pheno, string, add the phenotype name.
  #with_exons, logical, TRUE if we consider exonic variants ONLY.
  #strict, logical, TRUE if we consider exonic variants ONLY without synonymous variants. Only TRUE if with_exons if TRUE.
  #consanguinity, logical, TRUE if we consider consanguinity loops in the analysis.
  #remove_singletons, logical, TRUE if we need to remove singletons from the analysis.
  
  if(remove_singletons){
    #singletons <- fread("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag_without_mask/imputation_comb/merged_with_seq/freq/impute5_gigi2_combined_seq_fam_singleton.snplist", header = FALSE)
    singletons <- fread("/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag_without_mask/imputation_comb/merged_with_seq/freq/impute5_gigi2_combined_seq_fam_singleton.snplist", header = FALSE)
    out_remove_singletons <- "_without_singletons"
  } else {
    out_remove_singletons <- ""
  }

  if(strict & !with_exons){stop("Strict cannot be TRUE if with_exons is FALSE")}

  #Path to the data
  if(with_exons){
    path_data <- paste0(path_retrofun, "/genes_litt/exons_only")
    out_exons <- "only_exonic_variants"
    if(strict){
      path_data <- paste0(path_retrofun, "/genes_litt/exons_only_strict")
      out_exons <- "only_exonic_variants_without_synonymous"
    }
  }else{
    path_data <- paste0(path_retrofun, "/genes_litt")
    out_exons <- "all_variants"
  }

  file <- "impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.ped"
  path_gene_info <- paste0(path_data, "/bp_sz_genes_effects_to_keep_seq_FINAL.txt")
  genes_list <- read.table(paste0(path_retrofun, "/genes_litt/2024-06-25_selected_genes_ID_BP_SZ.txt"), header = FALSE)
  genes_paper <- read.table(paste0(path_retrofun, "/genes_litt/2024-11-01_selected_genes_ID_and_paper.txt"), header = FALSE)

  mapfile <- read.table(paste0(path_data, "/" , gsub(".ped", "", file), ".map"), header = FALSE)
  colnames(mapfile) <- c("chrom", "ID", "genpos", "pos")
  pedfile <- read.table(paste0(path_data, "/" , file))
  #Remove singletons if needed
  if(remove_singletons){
    not_singletons <- which(!mapfile$ID %in% singletons$V1)
    mapfile <- mapfile[not_singletons,]
    pedfile <- pedfile[,c(1:6, sort( c(5+(not_singletons*2), 6+(not_singletons*2)) ))]
  }

  for(i in 7:ncol(pedfile)){pedfile[,i] <- case_when(pedfile[,i] == 0 ~ 0, pedfile[,i] == 1 ~ 2, pedfile[,i] == 2 ~ 1)}
  
  #Remove 2620b which is a duplicate of 2620 in our data.
  pedfile <- pedfile[pedfile[,2] != "2620b",]
  
  #Load the phenotype files with consanguinity loops.
  #Affected must be coded by a 2.
  df.ped.2021 <- loadRData(paste0(path_retrofun, "/objets_ped/ped", pheno, "_orig.RData"))
  df.ped.2021 <- data.frame(df.ped.2021$famid, df.ped.2021$id, df.ped.2021$sex, df.ped.2021$affected+1)
  colnames(df.ped.2021) <- c("famid","id","sex","affected")
  pedfile.fam.infos <- setNames(data.frame(pedfile[,1], pedfile[,2], pedfile[,3], pedfile[,4]), c("famid", "id", "findex", "mindex"))
  pedfile.fam.infos.full <- merge(x = pedfile.fam.infos, y = df.ped.2021[,c("famid","id","sex","affected")], by = c("famid","id"), sort = FALSE, all.x = TRUE)
  
  #Associate the phenotype in the .ped and keep the affected families only.
  subset.fam <- pedfile.fam.infos.full %>% group_by(famid) %>% summarise(n_affected = sum(affected==2)) %>% filter(n_affected!=0) %>% select(famid) %>% as.vector()
  df.ped.2021 <- df.ped.2021[df.ped.2021$famid %in% subset.fam$famid,]
  pedfile.fam.infos.full <- pedfile.fam.infos.full[pedfile.fam.infos.full$famid %in% subset.fam$famid,]
  #Check
  pedfile <- pedfile[pedfile$V2 %in% pedfile.fam.infos.full$id, ]
  all(pedfile$V2 == pedfile.fam.infos.full$id)
  pedfile[,1:6] <- pedfile.fam.infos.full
  
  #Adjust the pedigree for the 3 problematic families in agg.genos.by.fam.
  if(consanguinity){
    null_name <- paste0("expected.variance.consanguinity.cryptique.",pheno,".rds")
    correction <- "none"
    out_consanguinity <- "with_consanguinity"
  } else {
    null_name <- paste0("expected.variance.",pheno,".rds")
    correction <- "replace"
    out_consanguinity <- "without_consanguinity"
  }
  
  if (pheno %in% c("GCbr","GCna"))
  {
  null <- readRDS(paste0(path_retrofun, "/objets_ped/", null_name))
  fam_split_119 <- readRDS(paste0(path_retrofun, "/objets_ped/fam119splitted.rds"))
  pedfile$V1[pedfile$V1 == "119" & pedfile$V2 %in% fam_split_119$id[fam_split_119$fam=="119-1"]] <- "119-1"
  pedfile$V1[pedfile$V1 == "119" & pedfile$V2 %in% fam_split_119$id[fam_split_119$fam=="119-2"]] <- "119-2"
  fam_split_131 <- readRDS(paste0(path_retrofun, "/objets_ped/fam131splitted.rds"))
  pedfile$V1[pedfile$V1 == "131" & pedfile$V2 %in% fam_split_131$id[fam_split_131$fam=="131-1"]] <- "131-1"
  pedfile$V1[pedfile$V1 == "131" & pedfile$V2 %in% fam_split_131$id[fam_split_131$fam=="131-2"]] <- "131-2"
  fam_split_255 <- readRDS(paste0(path_retrofun, "/objets_ped/fam255splitted.rds"))
  pedfile$V1[pedfile$V1 == "255" & pedfile$V2 %in% fam_split_255$id[fam_split_255$fam=="255-1"]] <- "255-1"
  pedfile$V1[pedfile$V1 == "255" & pedfile$V2 %in% fam_split_255$id[fam_split_255$fam=="255-2"]] <- "255-2"
  }
  else null <- readRDS(paste0(pathAB_ped, null_name))
  
  #By gene annotation matrix and by paper annotation matrix
  annotation.matrix <- matrix(data=0, nrow = nrow(mapfile), ncol = nrow(genes_list))
  colnames(annotation.matrix) <- genes_list$V1
  annotation.matrix.paper <- matrix(data=0, nrow = nrow(mapfile), ncol = length(unique(genes_paper$V2)))
  colnames(annotation.matrix.paper) <- unique(genes_paper$V2)
  results_by_gene <- list()
  for(gene in genes_list$V1){
    gene_paper <- genes_paper$V2[genes_paper$V1 == gene]
    gene_var_ID <- system(paste0('grep -w ', gene, ' ', path_gene_info, '| cut -d " " -f 1 | sed "s/ID=//g" | sed "s/,//g"'), intern = TRUE)
    gene_position <- which(mapfile$ID %in% gene_var_ID)
    results_by_gene[[gene]][["ped"]] <- pedfile[,c(1:6, sort( c(5+(gene_position*2), 6+(gene_position*2)) ))]
    results_by_gene[[gene]][["agg.genos.by.fam"]] <- tryCatch(agg.genos.by.fam(pedfile.path=NULL, pedfile = results_by_gene[[gene]][["ped"]], correction=correction), error = function(e) e)
    results_by_gene[[gene]][["n"]] <- length(results_by_gene[[gene]][["agg.genos.by.fam"]]$index_variants)
    if(inherits(results_by_gene[[gene]][["agg.genos.by.fam"]], "error")){results_by_gene[[gene]][["results"]] <- NA; next}
    results_by_gene[[gene]][["results"]] <- RetroFun.RVS(null, results_by_gene[[gene]][["agg.genos.by.fam"]],
                                                         Z_annot = matrix(1, ncol = 1, nrow = length(gene_position)), W = rep(1, length(gene_position)), independence=FALSE)
    annotation.matrix[gene_position, gene] <- 1
    annotation.matrix.paper[gene_position, gene_paper] <- 1
  }
  annotation.matrix <- cbind(annotation.matrix.paper, annotation.matrix)
  
  #Add a total burden for every variants.
  annotation.matrix <- cbind(1, annotation.matrix)
  colnames(annotation.matrix)[1] <- "allGenes"
  agg.genos.by.fam <- tryCatch(agg.genos.by.fam(pedfile.path=NULL, pedfile=pedfile, Z_annot=annotation.matrix, exclude_annot = 1, correction=correction), error = function(e) e)
  if(inherits(agg.genos.by.fam, "error")) next
  n_by_annot <- colSums(annotation.matrix[agg.genos.by.fam$index_variants,])
  df.annotation <- data.frame(annotation.matrix)
  results <- RetroFun.RVS(null, agg.genos.by.fam, Z_annot = df.annotation, W = rep(1, nrow(df.annotation)), independence=FALSE)
  results_n <- data.frame("p" = unlist(results), "n" = c(n_by_annot, length(agg.genos.by.fam$index_variants), length(agg.genos.by.fam$index_variants))); rownames(results_n) <- names(results)
  results_n$p_analyse_gene_seul <- c(NA, NA, NA, NA, NA, unlist(sapply(results_by_gene, FUN = function(x){x[["results"]]})), NA, NA)
  results_n$p_analyse_gene_seul_n <- c(NA, NA, NA, NA, NA, unlist(sapply(results_by_gene, FUN = function(x){x[["n"]]})), NA, NA)
  #saveRDS(results_n, paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/without_mask/RetroFunRVS_results_seq_genes_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))    
  saveRDS(results_n, paste0("/lustre09/project/6033529/schizo/results_AB/RetroFunRVS/WGS_bs_2022_500samples/without_mask/RetroFunRVS_results_seq_genes_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))    
}

#Function for genes from SynGO pathways
RetroFun.RVS.genes.pathways_run <- function(pheno, with_exons, strict = FALSE, consanguinity, remove_singletons = FALSE){
  #pheno, string, add the phenotype name.
  #with_exons, logical, TRUE if we consider exonic variants ONLY.
  #strict, logical, TRUE if we consider exonic variants ONLY without synonymous variants. Only TRUE if with_exons if TRUE.
  #consanguinity, logical, TRUE if we consider consanguinity loops in the analysis.
  #remove_singletons, logical, TRUE if we need to remove singletons from the analysis.
  
  if(remove_singletons){
    #singletons <- fread("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag_without_mask/imputation_comb/merged_with_seq/freq/impute5_gigi2_combined_seq_fam_singleton.snplist", header = FALSE)
    singletons <- fread("/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag_without_mask/imputation_comb/merged_with_seq/freq/impute5_gigi2_combined_seq_fam_singleton.snplist", header = FALSE)
    out_remove_singletons <- "_without_singletons"
  } else {
    out_remove_singletons <- ""
  }

  if(strict & !with_exons){stop("Strict cannot be TRUE if with_exons is FALSE")}

  #Path to the data
  if(with_exons){
    path_data <- paste0(path_retrofun, "/pathways/exons_only")
    out_exons <- "only_exonic_variants"
    if(strict){
      path_data <- paste0(path_retrofun, "/pathways/exons_only_strict")
      out_exons <- "only_exonic_variants_without_synonymous"
    }
  }else{
    path_data <- paste0(path_retrofun, "/pathways")
    out_exons <- "all_variants"
  }
  file <- "impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.ped"
  path_gene_info <- paste0(path_data, "/pathways_genes_effects_to_keep_seq_FINAL.txt")
  onto <- openxlsx::read.xlsx(paste0(path_retrofun, "/pathways/syngo_ontologies.xlsx"))

  mapfile <- read.table(paste0(path_data, "/" , gsub(".ped", "", file), ".map"), header = FALSE)
  colnames(mapfile) <- c("chrom", "ID", "genpos", "pos")
  pedfile <- read.table(paste0(path_data, "/" , file))
  #Remove singletons if needed
  if(remove_singletons){
    not_singletons <- which(!mapfile$ID %in% singletons$V1)
    mapfile <- mapfile[not_singletons,]
    pedfile <- pedfile[,c(1:6, sort( c(5+(not_singletons*2), 6+(not_singletons*2)) ))]
  }

  for(i in 7:ncol(pedfile)){pedfile[,i] <- case_when(pedfile[,i] == 0 ~ 0, pedfile[,i] == 1 ~ 2, pedfile[,i] == 2 ~ 1)}
  
  #Remove 2620b which is a duplicate of 2620 in our data.
  pedfile <- pedfile[pedfile[,2] != "2620b",]
  
  #Load the phenotype files with consanguinity loops.
  #Affected must be coded by a 2.
  df.ped.2021 <- loadRData(paste0(path_retrofun, "/objets_ped/ped", pheno, "_orig.RData"))
  df.ped.2021 <- data.frame(df.ped.2021$famid, df.ped.2021$id, df.ped.2021$sex, df.ped.2021$affected+1)
  colnames(df.ped.2021) <- c("famid","id","sex","affected")
  pedfile.fam.infos <- setNames(data.frame(pedfile[,1], pedfile[,2], pedfile[,3], pedfile[,4]), c("famid", "id", "findex", "mindex"))
  pedfile.fam.infos.full <- merge(x = pedfile.fam.infos, y = df.ped.2021[,c("famid","id","sex","affected")], by = c("famid","id"), sort = FALSE, all.x = TRUE)
  
  #Associate the phenotype in the .ped and keep the affected families only.
  subset.fam <- pedfile.fam.infos.full %>% group_by(famid) %>% summarise(n_affected = sum(affected==2)) %>% filter(n_affected!=0) %>% select(famid) %>% as.vector()
  df.ped.2021 <- df.ped.2021[df.ped.2021$famid %in% subset.fam$famid,]
  pedfile.fam.infos.full <- pedfile.fam.infos.full[pedfile.fam.infos.full$famid %in% subset.fam$famid,]
  #Check
  pedfile <- pedfile[pedfile$V2 %in% pedfile.fam.infos.full$id, ]
  all(pedfile$V2 == pedfile.fam.infos.full$id)
  pedfile[,1:6] <- pedfile.fam.infos.full
  
  #Adjust the pedigree for the 3 problematic families in agg.genos.by.fam.
  if(consanguinity){
    null_name <- paste0("expected.variance.consanguinity.cryptique.",pheno,".rds")
    correction <- "none"
    out_consanguinity <- "with_consanguinity"
  } else {
    null_name <- paste0("expected.variance.",pheno,".rds")
    correction <- "replace"
    out_consanguinity <- "without_consanguinity"
  }
  
  if (pheno %in% c("GCbr","GCna"))
  {
  null <- readRDS(paste0(path_retrofun, "/objets_ped/", null_name))
  fam_split_119 <- readRDS(paste0(path_retrofun, "/objets_ped/fam119splitted.rds"))
  pedfile$V1[pedfile$V1 == "119" & pedfile$V2 %in% fam_split_119$id[fam_split_119$fam=="119-1"]] <- "119-1"
  pedfile$V1[pedfile$V1 == "119" & pedfile$V2 %in% fam_split_119$id[fam_split_119$fam=="119-2"]] <- "119-2"
  fam_split_131 <- readRDS(paste0(path_retrofun, "/objets_ped/fam131splitted.rds"))
  pedfile$V1[pedfile$V1 == "131" & pedfile$V2 %in% fam_split_131$id[fam_split_131$fam=="131-1"]] <- "131-1"
  pedfile$V1[pedfile$V1 == "131" & pedfile$V2 %in% fam_split_131$id[fam_split_131$fam=="131-2"]] <- "131-2"
  fam_split_255 <- readRDS(paste0(path_retrofun, "/objets_ped/fam255splitted.rds"))
  pedfile$V1[pedfile$V1 == "255" & pedfile$V2 %in% fam_split_255$id[fam_split_255$fam=="255-1"]] <- "255-1"
  pedfile$V1[pedfile$V1 == "255" & pedfile$V2 %in% fam_split_255$id[fam_split_255$fam=="255-2"]] <- "255-2"
  }
  else null <- readRDS(paste0(pathAB_ped, null_name))
  
  #By gene annotation matrix and by paper annotation matrix
  annotation.matrix <- matrix(data=0, nrow = nrow(mapfile), ncol = nrow(onto))
  colnames(annotation.matrix) <- onto$id
  results_by_onto <- list()
  for(id in onto$id){
    gene_id <- onto$hgnc_symbol[onto$id == id]
    write <- data.table(gsub("(.*)", "|\\1|", strsplit(gene_id, ", ")[[1]], fixed = FALSE))
#    fwrite(write, paste0(path_data, "/gene_i_", pheno, "_", out_exons, "_", out_consanguinity, ".txt"), col.names = FALSE, row.names = FALSE)
#    gene_var_ID <- system(paste0('grep -f ', path_data, '/gene_i_', pheno, '_', out_exons, '_', out_consanguinity, '.txt ', path_gene_info, '| cut -d " " -f 1 | sed "s/ID=//g" | sed "s/,//g"'), intern = TRUE)
    fwrite(write, paste0(path_tmp, "/gene_i_", pheno, "_", out_exons, "_", out_consanguinity, ".txt"), col.names = FALSE, row.names = FALSE)
    gene_var_ID <- system(paste0('grep -f ', path_tmp, '/gene_i_', pheno, '_', out_exons, '_', out_consanguinity, '.txt ', path_gene_info, '| cut -d " " -f 1 | sed "s/ID=//g" | sed "s/,//g"'), intern = TRUE)
    var_position <- which(mapfile$ID %in% gene_var_ID)
    results_by_onto[[id]][["ped"]] <- pedfile[,c(1:6, sort( c(5+(var_position*2), 6+(var_position*2)) ))]
    results_by_onto[[id]][["agg.genos.by.fam"]] <- tryCatch(agg.genos.by.fam(pedfile.path=NULL, pedfile = results_by_onto[[id]][["ped"]], correction=correction), error = function(e) e)
    results_by_onto[[id]][["n"]] <- length(results_by_onto[[id]][["agg.genos.by.fam"]]$index_variants)
    if(inherits(results_by_onto[[id]][["agg.genos.by.fam"]], "error")){results_by_onto[[id]][["results"]] <- NA; next}
    results_by_onto[[id]][["results"]] <- RetroFun.RVS(null, results_by_onto[[id]][["agg.genos.by.fam"]],
                                                         Z_annot = matrix(1, ncol = 1, nrow = length(var_position)), W = rep(1, length(var_position)), independence=FALSE)
    annotation.matrix[var_position, id] <- 1
  }
  #system(paste0('rm ', path_data, '/gene_i_', pheno, '_', out_exons, '_', out_consanguinity, '.txt'))
  system(paste0('rm ', path_tmp, '/gene_i_', pheno, '_', out_exons, '_', out_consanguinity, '.txt'))
  
  #Add a total burden for every variants.
  annotation.matrix <- cbind(1, annotation.matrix)
  colnames(annotation.matrix)[1] <- "Burden"
  agg.genos.by.fam <- tryCatch(agg.genos.by.fam(pedfile.path=NULL, pedfile=pedfile, Z_annot=annotation.matrix, exclude_annot = 1, correction=correction), error = function(e) e)
  if(inherits(agg.genos.by.fam, "error")) next
  n_by_annot <- colSums(annotation.matrix[agg.genos.by.fam$index_variants,])
  df.annotation <- data.frame(annotation.matrix)
  results <- RetroFun.RVS(null, agg.genos.by.fam, Z_annot = df.annotation, W = rep(1, nrow(df.annotation)), independence=FALSE)
  results_n <- data.frame("p" = unlist(results), "n" = c(n_by_annot, length(agg.genos.by.fam$index_variants), length(agg.genos.by.fam$index_variants))); rownames(results_n) <- names(results)
  results_n$p_analyse_onto_seul <- c(NA, unlist(sapply(results_by_onto, FUN = function(x){x[["results"]]})), NA, NA)
  results_n$p_analyse_onto_seul_n <- c(NA, unlist(sapply(results_by_onto, FUN = function(x){x[["n"]]})), NA, NA)
  #saveRDS(results_n, paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/without_mask/RetroFunRVS_results_seq_pathways_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))    
  saveRDS(results_n, paste0("/lustre09/project/6033529/schizo/results_AB/RetroFunRVS/WGS_bs_2022_500samples/without_mask/RetroFunRVS_results_seq_pathways_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))    
}

#Run
RetroFun.RVS_run(pheno = pheno, with_exons = exons, consanguinity = consanguinity, remove_singletons = remove_singletons)
RetroFun.RVS.overlap02_run(pheno = pheno, with_exons = exons, consanguinity = consanguinity, remove_singletons = remove_singletons)

#In this case, if exons is true, the analysis is made ONLY among the variants that are exonic.
#if it is false, then all variants in the genes are included.
if(exons){
  #RetroFun.RVS.genes.litt_run(pheno = pheno, with_exons = exons, consanguinity = consanguinity, remove_singletons = remove_singletons)
  #RetroFun.RVS.genes.litt_run(pheno = pheno, with_exons = exons, strict = TRUE, consanguinity = consanguinity, remove_singletons = remove_singletons)
  RetroFun.RVS.genes.pathways_run(pheno = pheno, with_exons = exons, consanguinity = consanguinity, remove_singletons = remove_singletons)
  RetroFun.RVS.genes.pathways_run(pheno = pheno, with_exons = exons, strict = TRUE, consanguinity = consanguinity, remove_singletons = remove_singletons)  
}

#MAF objects
#This code import the functions used to output the MAF(freq=TRUE) or the minor alleles(freq=FALSE) of the 10 most significants RetroFun-RVS results.
if(exons){out_exons <- "CRHs_with_exons"} else {out_exons <- "CRHs_only"}
if(consanguinity){out_consanguinity <- "with_consanguinity"} else {out_consanguinity <- "without_consanguinity"}
if(remove_singletons){out_remove_singletons <- "_without_singletons"} else {out_remove_singletons <- ""}
#source("/lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples_without_mask/call/analyse_MAF_n_table.R")
source("/lustre09/project/6033529/schizo/scripts_AB/without_mask/analyse_MAF_n_table.R")

MAFs_by_TAD <- sign_CRH_results(freq = TRUE, pheno = pheno, with_exons = exons, consanguinity = consanguinity, remove_singletons = remove_singletons)
n_by_TAD <- sign_CRH_results(freq = FALSE, pheno = pheno, with_exons = exons, consanguinity = consanguinity, remove_singletons = remove_singletons)
results_by_TAD <- purrr::map(Map(cbind, MAFs_by_TAD, n_by_TAD), function(x){x[sort(colnames(x))]})
#saveRDS(results_by_TAD, paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/without_mask/MAF_n_variants_10_CRH_sign_by_fam_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
saveRDS(results_by_TAD, paste0("/lustre09/project/6033529/schizo/results_AB/RetroFunRVS/WGS_bs_2022_500samples/without_mask/MAF_n_variants_10_CRH_sign_by_fam_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))

MAFs_by_TAD_0 <- sign_CRH_results_overlap(freq = TRUE, pheno = pheno, overlap = 0, with_exons = exons, consanguinity = consanguinity, remove_singletons = remove_singletons)
n_by_TAD_0 <- sign_CRH_results_overlap(freq = FALSE, pheno = pheno, overlap = 0, with_exons = exons, consanguinity = consanguinity, remove_singletons = remove_singletons)
results_by_TAD_0 <- purrr::map(Map(cbind, MAFs_by_TAD_0, n_by_TAD_0), function(x){x[sort(colnames(x))]})
#saveRDS(results_by_TAD_0, paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/without_mask/MAF_n_variants_10_CRH_sign_by_fam_overlap_0_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
saveRDS(results_by_TAD_0, paste0("/lustre09/project/6033529/schizo/results_AB/RetroFunRVS/WGS_bs_2022_500samples/without_mask/MAF_n_variants_10_CRH_sign_by_fam_overlap_0_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))

MAFs_by_TAD_2 <- sign_CRH_results_overlap(freq = TRUE, pheno = pheno, overlap = 2, with_exons = exons, consanguinity = consanguinity, remove_singletons = remove_singletons)
n_by_TAD_2 <- sign_CRH_results_overlap(freq = FALSE, pheno = pheno, overlap = 2, with_exons = exons, consanguinity = consanguinity, remove_singletons = remove_singletons)
results_by_TAD_2 <- purrr::map(Map(cbind, MAFs_by_TAD_2, n_by_TAD_2), function(x){x[sort(colnames(x))]})
#saveRDS(results_by_TAD_2, paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/without_mask/MAF_n_variants_10_CRH_sign_by_fam_overlap_2_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
saveRDS(results_by_TAD_2, paste0("/lustre09/project/6033529/schizo/results_AB/RetroFunRVS/WGS_bs_2022_500samples/without_mask/MAF_n_variants_10_CRH_sign_by_fam_overlap_2_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))

if(consanguinity){out_consanguinity <- "with_consanguinity"} else {out_consanguinity <- "without_consanguinity"}
if(exons){
  out_exons <- "only_exonic_variants"
  MAFs_by_genes_litt <- sign_genes_litt_results(freq = TRUE, pheno = pheno, with_exons = exons, consanguinity = consanguinity, remove_singletons = remove_singletons)
  n_by_genes_litt <- sign_genes_litt_results(freq = FALSE, pheno = pheno, with_exons = exons, consanguinity = consanguinity, remove_singletons = remove_singletons)
  results_by_genes_litt <- purrr::map(Map(cbind, MAFs_by_genes_litt, n_by_genes_litt), function(x){x[sort(colnames(x))]})
  #saveRDS(results_by_genes_litt, paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/without_mask/MAF_n_variants_10_genes_sign_by_fam_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
  saveRDS(results_by_genes_litt, paste0("/lustre09/project/6033529/schizo/results_AB/RetroFunRVS/WGS_bs_2022_500samples/without_mask/MAF_n_variants_10_genes_sign_by_fam_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
  
  MAFs_by_genes_pathways <- sign_genes_pathways_results(freq = TRUE, pheno = pheno, with_exons = exons, consanguinity = consanguinity, remove_singletons = remove_singletons)
  n_by_genes_pathways <- sign_genes_pathways_results(freq = FALSE, pheno = pheno, with_exons = exons, consanguinity = consanguinity, remove_singletons = remove_singletons)
  results_by_genes_pathways <- purrr::map(Map(cbind, MAFs_by_genes_pathways, n_by_genes_pathways), function(x){x[sort(colnames(x))]})
  #saveRDS(results_by_genes_pathways, paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/without_mask/MAF_n_variants_10_pathways_sign_by_fam_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
  saveRDS(results_by_genes_pathways, paste0("/lustre09/project/6033529/schizo/results_AB/RetroFunRVS/WGS_bs_2022_500samples/without_mask/MAF_n_variants_10_pathways_sign_by_fam_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
  
  out_exons <- "only_exonic_variants_without_synonymous"
  MAFs_by_genes_litt_strict <- sign_genes_litt_results(freq = TRUE, pheno = pheno, with_exons = exons, strict = TRUE, consanguinity = consanguinity, remove_singletons = remove_singletons)
  n_by_genes_litt_strict <- sign_genes_litt_results(freq = FALSE, pheno = pheno, with_exons = exons, strict = TRUE, consanguinity = consanguinity, remove_singletons = remove_singletons)
  results_by_genes_litt_strict <- purrr::map(Map(cbind, MAFs_by_genes_litt_strict, n_by_genes_litt_strict), function(x){x[sort(colnames(x))]})
  #saveRDS(results_by_genes_litt_strict, paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/without_mask/MAF_n_variants_10_genes_sign_by_fam_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
  saveRDS(results_by_genes_litt_strict, paste0("/lustre09/project/6033529/schizo/results_AB/RetroFunRVS/WGS_bs_2022_500samples/without_mask/MAF_n_variants_10_genes_sign_by_fam_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
  
  MAFs_by_genes_pathways_strict <- sign_genes_pathways_results(freq = TRUE, pheno = pheno, with_exons = exons, strict = TRUE, consanguinity = consanguinity, remove_singletons = remove_singletons)
  n_by_genes_pathways_strict <- sign_genes_pathways_results(freq = FALSE, pheno = pheno, with_exons = exons, strict = TRUE, consanguinity = consanguinity, remove_singletons = remove_singletons)
  results_by_genes_pathways_strict <- purrr::map(Map(cbind, MAFs_by_genes_pathways_strict, n_by_genes_pathways_strict), function(x){x[sort(colnames(x))]})
  #saveRDS(results_by_genes_pathways_strict, paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/without_mask/MAF_n_variants_10_pathways_sign_by_fam_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
  saveRDS(results_by_genes_pathways_strict, paste0("/lustre09/project/6033529/schizo/results_AB/RetroFunRVS/WGS_bs_2022_500samples/without_mask/MAF_n_variants_10_pathways_sign_by_fam_", pheno, "_", out_exons, "_", out_consanguinity, out_remove_singletons, ".RDS"))
}