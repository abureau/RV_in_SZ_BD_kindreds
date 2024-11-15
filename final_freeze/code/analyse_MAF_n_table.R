#module load StdEnv/2020 plink/1.9b_6.21-x86_64 gcc/9.3.0 vcftools/0.1.16 bcftools/1.16 r/4.2
path_retrofun <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS"
setwd(path_retrofun)
library(stringr); library(GenomicRanges); library(dplyr); library(RetroFunRVS)

#Import pedigree data with the phenotype
loadRData <- function(file_name){load(file_name); get(ls()[ls() != "file_name"])}
df.ped.2021 <- loadRData(paste0(path_retrofun, "/objets_ped/ped", pheno, "_orig.RData"))
df.ped.2021 <- data.frame(df.ped.2021$famid, df.ped.2021$id, df.ped.2021$sex, df.ped.2021$affected+1)
colnames(df.ped.2021) <- c("famid","id","sex","affected")

#This function is used to extract variants in certain position ranges.
is_within_any_range <- function(value, ranges_df) {
  any(value >= ranges_df$chromStart& value <= ranges_df$chromEnd)
}

#This function is used to compute the MAF(freq=T) or the number of minor alleles(freq=F).
MAF <- function(ped, var, freq = TRUE){
  ped_var <- ped[,c(5+(var*2), 6+(var*2))]
  tab <- sort(table(unlist(ped_var)))
  #Minor alleles are 1 and major alleles are 2. 
  #If the length of tab is 1 and the only allele is "2", then the MAF is 0;
  #If the length of tab is 1 and the only allele is "1", then the MAF is 1;
  if(length(tab)==1){
    if(names(tab) == "2"){tab <- c(0,tab)}else if(names(tab) == "1"){tab <- c(tab,0)}
  } 
  if(freq){return(tab[[1]]/sum(tab))}else{return(tab[[1]])}
}

#---- Overlap 1 ----
#This function is used to output the MAF(freq=TRUE) or to output the minor alleles(freq=FALSE) in a list of data.frames.
sign_CRH_results <- function(freq, pheno, with_exons, consanguinity){
  if(with_exons){out_exons <- "CRHs_with_exons"} else {out_exons <- "CRHs_only"}
  if(consanguinity){out_consanguinity <- "with_consanguinity"} else {out_consanguinity <- "without_consanguinity"}

  results <- readRDS(paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/RetroFunRVS_results_seq_all_chromosomes_", pheno, "_", out_exons, "_", out_consanguinity, ".RDS"))
  #Use the 10 most significant CRHs.
  test_var <- "score"
  results <- results[order(results[[test_var]]),]
  sign_CRH <- results[1:10, 1:3]

  #Initialize the results object.
  MAFs_by_TAD <- vector(mode = "list", length = nrow(sign_CRH))
  names(MAFs_by_TAD) <- sign_CRH$TAD_name
  if(freq){out_name <- "MAF"}else{out_name <- "n"}
  for(CRH in 1:nrow(sign_CRH)){
    chr <- sign_CRH$chr[CRH]
    TAD <- sign_CRH$TAD[CRH]
    #Path to the TAD data where the present CRH is found.
    datafile <- paste0("impute5_gigi2_combined_seq_RV_FINAL_chr", chr, "_TAD_", TAD)
    if(with_exons){
      #If the addition of exons had an impact on the .ped, import the modified one.
      if(file.exists(paste0(path_retrofun, "/TADs/with_exons/", datafile, ".ped"))){datafile <- paste0("with_exons/", datafile)}
    }
    ped <- read.table(paste0(path_retrofun, "/TADs/", datafile, ".ped"))
    map <- read.table(paste0(path_retrofun, "/TADs/", datafile, ".map"))
    TADs_chr <- data.table::fread(paste0(path_retrofun, "/TADs/TADs_list_chr", chr, ".bed"), header = FALSE)
    
    #Remove 2620b which is a duplicate of 2620 in our data.
    ped <- ped[ped[,2] != "2620b",]

    #Import and split the TAD data to keep the variants in this CRH only.
    if(with_exons){
      CRHs.by.TAD <- read.table(paste0(path_retrofun, "/TADs/with_exons/TADs_for_CRHs_overlap_1_TAD/chr", chr, "/chr", chr, "_", TADs_chr$V2[TAD]+1, "_", TADs_chr$V3[TAD], ".bed"), header=FALSE)
      if(CRHs.by.TAD[1,1] == "chrom"){CRHs.by.TAD <- CRHs.by.TAD[2:nrow(CRHs.by.TAD),]} #It appears that sometimes, the header is part of the dataframe, remove the line.
      colnames(CRHs.by.TAD) <- c("chrom", "chromStart", "chromEnd", "name")
    } else {
      CRHs.by.TAD <- read.table(paste0("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/TADs_for_CRHs_overlap_1_TAD/chr", chr, "/chr", chr, "_", TADs_chr$V2[TAD]+1, "_", TADs_chr$V3[TAD], ".bed"), header=TRUE)
      colnames(CRHs.by.TAD) <- c("chrom", "chromStart", "chromEnd", "name")
    }
    GRanges.CRHs.by.TAD <- GRanges(seqnames=CRHs.by.TAD$chrom, ranges=IRanges(start=CRHs.by.TAD$chromStart, end=CRHs.by.TAD$chromEnd),name=CRHs.by.TAD$name)
    #Add the phenotype to the ped file.
    ped <- merge(x = ped, y = df.ped.2021[,c("famid","id","affected")], by.x = c("V1", "V2"), by.y=c("famid","id"), sort = FALSE, all.x = TRUE)
    subset.fam <- ped %>% group_by(V1) %>% summarise(n_affected = sum(affected==2)) %>% filter(n_affected!=0) %>% select(V1) %>% as.vector()
    ped <- ped[ped$V1 %in% subset.fam$V1,]

    ped$V6 <- ped$affected; ped$affected <- NULL
    fam_split_119 <- readRDS(paste0(path_retrofun, "/objets_ped/fam119splitted.rds"))
    ped$V1[ped$V1 == "119" & ped$V2 %in% fam_split_119$id[fam_split_119$fam=="119-1"]] <- "119-1"
    ped$V1[ped$V1 == "119" & ped$V2 %in% fam_split_119$id[fam_split_119$fam=="119-2"]] <- "119-2"
    fam_split_131 <- readRDS(paste0(path_retrofun, "/objets_ped/fam131splitted.rds"))
    ped$V1[ped$V1 == "131" & ped$V2 %in% fam_split_131$id[fam_split_131$fam=="131-1"]] <- "131-1"
    ped$V1[ped$V1 == "131" & ped$V2 %in% fam_split_131$id[fam_split_131$fam=="131-2"]] <- "131-2"
    fam_split_255 <- readRDS(paste0(path_retrofun, "/objets_ped/fam255splitted.rds"))
    ped$V1[ped$V1 == "255" & ped$V2 %in% fam_split_255$id[fam_split_255$fam=="255-1"]] <- "255-1"
    ped$V1[ped$V1 == "255" & ped$V2 %in% fam_split_255$id[fam_split_255$fam=="255-2"]] <- "255-2"
    if(consanguinity){correction <- "none"}else{correction <- "replace"}
    ped_inv <- ped
    for(i in 7:ncol(ped_inv)){ped_inv[,i] <- case_when(ped_inv[,i] == 0 ~ 0, ped_inv[,i] == 1 ~ 2, ped_inv[,i] == 2 ~ 1)}
    agg.genos.by.fam <- agg.genos.by.fam(pedfile.path=NULL, pedfile=ped_inv, correction = correction)
    #extract_agg <- agg.genos.by.fam$index_variants
    #ped <- ped[,c(1:6, sort( c(5+(extract_agg*2), 6+(extract_agg*2)) ))]
    #map <- map[extract_agg, ]

    map_GRanges <- GRanges(seqnames = paste0("chr", map$V1), ranges = IRanges(start = map$V4, end = map$V4))
    overlap <- findOverlaps(map_GRanges, GRanges.CRHs.by.TAD)
    extract <- unique(queryHits(overlap))
    #Just in case it happens that no variants are found in the CRH (it shouldn't happen).
    if(sum(extract)==0){MAFs_by_TAD[[CRH]] <- "no variants are found in this CRH"}
    ped_CRH <- ped[,c(1:6, sort( c(5+(extract*2), 6+(extract*2)) ))]
    map_CRH <- map[extract, ]

    #Set an index for the phenotype status of each subject in each family.
    ped_CRH$fam_pheno <- paste0(ped_CRH$V1, "_", ped_CRH$V6)
    for(CRH_var in 1:nrow(map_CRH)){
      #Compute the MAF or the number of minor alleles by phenotype status in each family. 
      MAFs_var <- data.frame(do.call(rbind, by(ped_CRH, ped_CRH$fam_pheno, MAF, var = CRH_var, freq = freq, simplify = FALSE)))
      colnames(MAFs_var) <- paste0("MAF_", map_CRH$V2[CRH_var])
      MAFs_var$fam_pheno <- rownames(MAFs_var)
      if(CRH_var == 1){MAFs <- MAFs_var} else {MAFs <- merge(x = MAFs, MAFs_var, by = "fam_pheno")}
    }
    MAFs[, c("fam", "pheno")] <- str_split_fixed(MAFs$fam_pheno, "_", 2)
    #Keep the variants present in at least one affected subject.
    extract_var <- apply(select(MAFs[MAFs$pheno == 2,], -fam, -pheno, -fam_pheno), 2, function(x){ifelse(sum(x==0)!=sum(MAFs$pheno == 2), TRUE, FALSE)})
    MAFs <- MAFs[, c("fam_pheno", names(extract_var)[extract_var], "fam", "pheno")]
    #Remove the families that have none of the remaining variants
    extract_fam_which <- apply(select(MAFs[MAFs$pheno == 2,], -fam, -pheno, -fam_pheno), 1, function(x){ifelse(sum(x==0)!=sum(extract_var), TRUE, FALSE)})
    extract_fam <- gsub("-1|-2", "", MAFs$fam) %in% unique(gsub("-1|-2", "", MAFs$fam[MAFs$pheno == 2][extract_fam_which]))
    MAFs <- merge(x = MAFs, y = aggregate(extract_fam, list(fam = MAFs$fam), function(x){ifelse(any(unlist(x) == 1), TRUE, FALSE)}), by = "fam")
    MAFs <- MAFs[which(MAFs$x), c("fam_pheno", names(extract_var)[extract_var])]
    MAFs_by_TAD[[CRH]] <- data.frame(t(select(MAFs, -fam_pheno))) %>% setNames(paste0(MAFs$fam_pheno, "_", out_name))
  }
  return(MAFs_by_TAD)
}

#---- Overlap 0, 2 ----
#This function is used to output the MAF(freq=TRUE) or to output the minor alleles(freq=FALSE) in a list of data.frames.
sign_CRH_results_overlap <- function(freq, pheno, overlap = c(0,2), with_exons, consanguinity){
  if(with_exons){out_exons <- "CRHs_with_exons"} else {out_exons <- "CRHs_only"}
  if(consanguinity){out_consanguinity <- "with_consanguinity"} else {out_consanguinity <- "without_consanguinity"}

  results <- readRDS(paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/RetroFunRVS_results_seq_all_chromosomes_overlap_", overlap, "_", pheno, "_", out_exons, "_", out_consanguinity, ".RDS"))
  path_data <- paste0(path_retrofun, "/TADs/overlap_", overlap)
  if(with_exons){path_data <- paste0(path_data, "/with_exons")}
  #Use the 10 most significant CRHs.
  results <- results[order(results[["result"]]),]
  sign_CRH <- results[1:10, 1:3]

  #Initialize the results object.
  MAFs_by_TAD <- vector(mode = "list", length = nrow(sign_CRH))
  names(MAFs_by_TAD) <- sign_CRH$TAD_name
  if(freq){out_name <- "MAF"}else{out_name <- "n"}
  for(CRH in 1:nrow(sign_CRH)){
    chr <- sign_CRH$chr[CRH]
    CRH_value <- sign_CRH$CRH_name[CRH]
    #Path to the data where the present CRH is found.
    datafile <- paste0("/impute5_gigi2_combined_seq_RV_FINAL_chr", chr, "_CRH_", CRH_value)
    ped <- read.table(paste0(path_data, "/", datafile, ".ped"))
    map <- read.table(paste0(path_data, "/", datafile, ".map"))

    #Remove 2620b which is a duplicate of 2620 in our data.
    ped <- ped[ped[,2] != "2620b",]

    #Add the phenotype to the ped file.
    ped <- merge(x = ped, y = df.ped.2021[,c("famid","id","affected")], by.x = c("V1", "V2"), by.y=c("famid","id"), sort = FALSE, all.x = TRUE)
    subset.fam <- ped %>% group_by(V1) %>% summarise(n_affected = sum(affected==2)) %>% filter(n_affected!=0) %>% select(V1) %>% as.vector()
    ped <- ped[ped$V1 %in% subset.fam$V1,]

    ped$V6 <- ped$affected; ped$affected <- NULL
    fam_split_119 <- readRDS(paste0(path_retrofun, "/objets_ped/fam119splitted.rds"))
    ped$V1[ped$V1 == "119" & ped$V2 %in% fam_split_119$id[fam_split_119$fam=="119-1"]] <- "119-1"
    ped$V1[ped$V1 == "119" & ped$V2 %in% fam_split_119$id[fam_split_119$fam=="119-2"]] <- "119-2"
    fam_split_131 <- readRDS(paste0(path_retrofun, "/objets_ped/fam131splitted.rds"))
    ped$V1[ped$V1 == "131" & ped$V2 %in% fam_split_131$id[fam_split_131$fam=="131-1"]] <- "131-1"
    ped$V1[ped$V1 == "131" & ped$V2 %in% fam_split_131$id[fam_split_131$fam=="131-2"]] <- "131-2"
    fam_split_255 <- readRDS(paste0(path_retrofun, "/objets_ped/fam255splitted.rds"))
    ped$V1[ped$V1 == "255" & ped$V2 %in% fam_split_255$id[fam_split_255$fam=="255-1"]] <- "255-1"
    ped$V1[ped$V1 == "255" & ped$V2 %in% fam_split_255$id[fam_split_255$fam=="255-2"]] <- "255-2"
    if(consanguinity){correction <- "none"}else{correction <- "replace"}
    ped_inv <- ped
    for(i in 7:ncol(ped_inv)){ped_inv[,i] <- case_when(ped_inv[,i] == 0 ~ 0, ped_inv[,i] == 1 ~ 2, ped_inv[,i] == 2 ~ 1)}
    #agg.genos.by.fam <- agg.genos.by.fam(pedfile.path=NULL, pedfile=ped_inv, correction = correction)
    #extract_agg <- agg.genos.by.fam$index_variants
    #ped <- ped[,c(1:6, sort( c(5+(extract_agg*2), 6+(extract_agg*2)) ))]
    #map <- map[extract_agg, ]

    #Set an index for the phenotype status of each subject in each family.
    ped$fam_pheno <- paste0(ped$V1, "_", ped$V6)
    for(CRH_var in 1:nrow(map)){
      #Compute the MAF or the number of minor alleles by phenotype status in each family. 
      MAFs_var <- data.frame(do.call(rbind, by(ped, ped$fam_pheno, MAF, var = CRH_var, freq = freq, simplify = FALSE)))
      colnames(MAFs_var) <- paste0("MAF_", map$V2[CRH_var])
      MAFs_var$fam_pheno <- rownames(MAFs_var)
      if(CRH_var == 1){MAFs <- MAFs_var} else {MAFs <- merge(x = MAFs, MAFs_var, by = "fam_pheno")}
    }
    MAFs[, c("fam", "pheno")] <- str_split_fixed(MAFs$fam_pheno, "_", 2)
    #Keep the variants present in at least one affected subject.
    extract_var <- apply(select(MAFs[MAFs$pheno == 2,], -fam, -pheno, -fam_pheno), 2, function(x){ifelse(sum(x==0)!=sum(MAFs$pheno == 2), TRUE, FALSE)})
    MAFs <- MAFs[, c("fam_pheno", names(extract_var)[extract_var], "fam", "pheno")]
    #Remove the families that have none of the remaining variants
    extract_fam_which <- apply(select(MAFs[MAFs$pheno == 2,], -fam, -pheno, -fam_pheno), 1, function(x){ifelse(sum(x==0)!=sum(extract_var), TRUE, FALSE)})
    extract_fam <- gsub("-1|-2", "", MAFs$fam) %in% unique(gsub("-1|-2", "", MAFs$fam[MAFs$pheno == 2][extract_fam_which]))
    MAFs <- merge(x = MAFs, y = aggregate(extract_fam, list(fam = MAFs$fam), function(x){ifelse(any(unlist(x) == 1), TRUE, FALSE)}), by = "fam")
    MAFs <- MAFs[which(MAFs$x), c("fam_pheno", names(extract_var)[extract_var])]
    MAFs_by_TAD[[CRH]] <- data.frame(t(select(MAFs, -fam_pheno))) %>% setNames(paste0(MAFs$fam_pheno, "_", out_name))
  }
  return(MAFs_by_TAD)
}

#---- Genes ----
sign_genes_litt_results <- function(freq, pheno, with_exons, strict = FALSE, consanguinity){
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
  if(consanguinity){out_consanguinity <- "with_consanguinity"} else {out_consanguinity <- "without_consanguinity"}

  results <- readRDS(paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/RetroFunRVS_results_seq_genes_", pheno, "_", out_exons, "_", out_consanguinity, ".RDS"))
  results <- results[1:(nrow(results)-2),] #Remove ACAT and Fisher
  results <- results[6:nrow(results),] #Remove genes by paper
  #Use the 10 most significant CRHs.
  test_var <- "p_analyse_gene_seul"
  results <- results[order(results[[test_var]]),]
  sign_gene <- results[1:10, test_var, drop = FALSE]

  file <- "impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.ped"
  path_gene_info <- paste0(path_data, "/bp_sz_genes_effects_to_keep_seq_FINAL.txt")
  ped <- read.table(paste0(path_data, "/" , file))
  map <- read.table(paste0(path_data, "/" , gsub(".ped", "", file), ".map"), header = FALSE)
  #Remove 2620b which is a duplicate of 2620 in our data.
  ped <- ped[ped[,2] != "2620b",]

  ped <- merge(x = ped, y = df.ped.2021[,c("famid","id","affected")], by.x = c("V1", "V2"), by.y=c("famid","id"), sort = FALSE, all.x = TRUE)
  subset.fam <- ped %>% group_by(V1) %>% summarise(n_affected = sum(affected==2)) %>% filter(n_affected!=0) %>% select(V1) %>% as.vector()
  ped <- ped[ped$V1 %in% subset.fam$V1,]
  ped$V6 <- ped$affected; ped$affected <- NULL
  fam_split_119 <- readRDS(paste0(path_retrofun, "/objets_ped/fam119splitted.rds"))
  ped$V1[ped$V1 == "119" & ped$V2 %in% fam_split_119$id[fam_split_119$fam=="119-1"]] <- "119-1"
  ped$V1[ped$V1 == "119" & ped$V2 %in% fam_split_119$id[fam_split_119$fam=="119-2"]] <- "119-2"
  fam_split_131 <- readRDS(paste0(path_retrofun, "/objets_ped/fam131splitted.rds"))
  ped$V1[ped$V1 == "131" & ped$V2 %in% fam_split_131$id[fam_split_131$fam=="131-1"]] <- "131-1"
  ped$V1[ped$V1 == "131" & ped$V2 %in% fam_split_131$id[fam_split_131$fam=="131-2"]] <- "131-2"
  fam_split_255 <- readRDS(paste0(path_retrofun, "/objets_ped/fam255splitted.rds"))
  ped$V1[ped$V1 == "255" & ped$V2 %in% fam_split_255$id[fam_split_255$fam=="255-1"]] <- "255-1"
  ped$V1[ped$V1 == "255" & ped$V2 %in% fam_split_255$id[fam_split_255$fam=="255-2"]] <- "255-2"
  ped_inv <- ped
  for(i in 7:ncol(ped_inv)){ped_inv[,i] <- case_when(ped_inv[,i] == 0 ~ 0, ped_inv[,i] == 1 ~ 2, ped_inv[,i] == 2 ~ 1)}

  #Initialize the results object.
  MAFs_by_gene <- vector(mode = "list", length = nrow(sign_gene))
  names(MAFs_by_gene) <- rownames(sign_gene)
  if(freq){out_name <- "MAF"}else{out_name <- "n"}
  for(gene in 1:nrow(sign_gene)){
    gene_name <- gsub("Score_", "", rownames(sign_gene)[gene])
    gene_var_ID <- system(paste0('grep -w ', gene_name, ' ', path_gene_info, '| cut -d " " -f 1 | sed "s/ID=//g" | sed "s/,//g"'), intern = TRUE)
    extract <- which(map[[2]] %in% gene_var_ID)
    ped_gene <- ped[,c(1:6, sort( c(5+(extract*2), 6+(extract*2)) ))]
    ped_inv_gene <- ped_inv[,c(1:6, sort( c(5+(extract*2), 6+(extract*2)) ))]
    map_gene <- map[extract, ]

    if(consanguinity){correction <- "none"}else{correction <- "replace"}
    agg.genos.by.fam <- agg.genos.by.fam(pedfile.path=NULL, pedfile=ped_inv_gene, correction = correction)
    extract_agg <- agg.genos.by.fam$index_variants
    ped_gene <- ped_gene[,c(1:6, sort( c(5+(extract_agg*2), 6+(extract_agg*2)) ))]
    map_gene <- map_gene[extract_agg, ]

    #Set an index for the phenotype status of each subject in each family.
    ped_gene$fam_pheno <- paste0(ped_gene$V1, "_", ped_gene$V6)
    for(gene_var in 1:nrow(map_gene)){
      #Compute the MAF or the number of minor alleles by phenotype status in each family. 
      MAFs_var <- data.frame(do.call(rbind, by(ped_gene, ped_gene$fam_pheno, MAF, var = gene_var, freq = freq, simplify = FALSE)))
      colnames(MAFs_var) <- paste0("MAF_", map_gene$V2[gene_var])
      MAFs_var$fam_pheno <- rownames(MAFs_var)
      if(gene_var == 1){MAFs <- MAFs_var} else {MAFs <- merge(x = MAFs, MAFs_var, by = "fam_pheno")}
    }
    MAFs[, c("fam", "pheno")] <- str_split_fixed(MAFs$fam_pheno, "_", 2)
    #Keep the variants present in at least one affected subject.
    extract_var <- apply(select(MAFs[MAFs$pheno == 2,], -fam, -pheno, -fam_pheno), 2, function(x){ifelse(sum(x==0)!=sum(MAFs$pheno == 2), TRUE, FALSE)})
    MAFs <- MAFs[, c("fam_pheno", names(extract_var)[extract_var], "fam", "pheno")]
    #Remove the families that have none of the remaining variants
    extract_fam_which <- apply(select(MAFs[MAFs$pheno == 2,], -fam, -pheno, -fam_pheno), 1, function(x){ifelse(sum(x==0)!=sum(extract_var), TRUE, FALSE)})
    extract_fam <- gsub("-1|-2", "", MAFs$fam) %in% unique(gsub("-1|-2", "", MAFs$fam[MAFs$pheno == 2][extract_fam_which]))
    MAFs <- merge(x = MAFs, y = aggregate(extract_fam, list(fam = MAFs$fam), function(x){ifelse(any(unlist(x) == 1), TRUE, FALSE)}), by = "fam")
    MAFs <- MAFs[which(MAFs$x), c("fam_pheno", names(extract_var)[extract_var])]
    MAFs_by_gene[[gene]] <- data.frame(t(select(MAFs, -fam_pheno))) %>% setNames(paste0(MAFs$fam_pheno, "_", out_name))
  }
  return(MAFs_by_gene)
}

#---- Pathways ----
sign_genes_pathways_results <- function(freq, pheno, with_exons, strict = FALSE, consanguinity){
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
  if(consanguinity){out_consanguinity <- "with_consanguinity"} else {out_consanguinity <- "without_consanguinity"}

  results <- readRDS(paste0("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/WGS_bs_2022_500samples/RetroFunRVS_results_seq_pathways_", pheno, "_", out_exons, "_", out_consanguinity, ".RDS"))
  results <- results[1:(nrow(results)-2),] #Remove ACAT and Fisher
  #Use the 10 most significant CRHs.
  test_var <- "p_analyse_onto_seul"
  results <- results[order(results[[test_var]]),]
  sign_onto <- results[1:10, test_var, drop = FALSE]

  file <- "impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.ped"
  path_gene_info <- paste0(path_data, "/pathways_genes_effects_to_keep_seq_FINAL.txt")
  onto_file <- openxlsx::read.xlsx(paste0(path_retrofun, "/pathways/syngo_ontologies.xlsx"))
  ped <- read.table(paste0(path_data, "/" , file))
  map <- read.table(paste0(path_data, "/" , gsub(".ped", "", file), ".map"), header = FALSE)
  #Remove 2620b which is a duplicate of 2620 in our data.
  ped <- ped[ped[,2] != "2620b",]

  ped <- merge(x = ped, y = df.ped.2021[,c("famid","id","affected")], by.x = c("V1", "V2"), by.y=c("famid","id"), sort = FALSE, all.x = TRUE)
  subset.fam <- ped %>% group_by(V1) %>% summarise(n_affected = sum(affected==2)) %>% filter(n_affected!=0) %>% select(V1) %>% as.vector()
  ped <- ped[ped$V1 %in% subset.fam$V1,]
  ped$V6 <- ped$affected; ped$affected <- NULL
  fam_split_119 <- readRDS(paste0(path_retrofun, "/objets_ped/fam119splitted.rds"))
  ped$V1[ped$V1 == "119" & ped$V2 %in% fam_split_119$id[fam_split_119$fam=="119-1"]] <- "119-1"
  ped$V1[ped$V1 == "119" & ped$V2 %in% fam_split_119$id[fam_split_119$fam=="119-2"]] <- "119-2"
  fam_split_131 <- readRDS(paste0(path_retrofun, "/objets_ped/fam131splitted.rds"))
  ped$V1[ped$V1 == "131" & ped$V2 %in% fam_split_131$id[fam_split_131$fam=="131-1"]] <- "131-1"
  ped$V1[ped$V1 == "131" & ped$V2 %in% fam_split_131$id[fam_split_131$fam=="131-2"]] <- "131-2"
  fam_split_255 <- readRDS(paste0(path_retrofun, "/objets_ped/fam255splitted.rds"))
  ped$V1[ped$V1 == "255" & ped$V2 %in% fam_split_255$id[fam_split_255$fam=="255-1"]] <- "255-1"
  ped$V1[ped$V1 == "255" & ped$V2 %in% fam_split_255$id[fam_split_255$fam=="255-2"]] <- "255-2"
  ped_inv <- ped
  for(i in 7:ncol(ped_inv)){ped_inv[,i] <- case_when(ped_inv[,i] == 0 ~ 0, ped_inv[,i] == 1 ~ 2, ped_inv[,i] == 2 ~ 1)}

  #Initialize the results object.
  MAFs_by_onto <- vector(mode = "list", length = nrow(sign_onto))
  names(MAFs_by_onto) <- rownames(sign_onto)
  if(freq){out_name <- "MAF"}else{out_name <- "n"}
  for(onto in 1:nrow(sign_onto)){
    onto_name <- gsub("\\.", ":", gsub("Score_", "", rownames(sign_onto)[onto]))
    onto_genes <- onto_file$hgnc_symbol[onto_file$id == onto_name]
    write <- data.table::data.table(gsub("(.*)", "|\\1|", strsplit(onto_genes, ", ")[[1]], fixed = FALSE))
    data.table::fwrite(write, paste0(path_data, "/gene_i_", pheno, "_", out_exons, "_", out_consanguinity, ".txt"), col.names = FALSE, row.names = FALSE)
    onto_var_ID <- system(paste0('grep -f ', path_data, '/gene_i_', pheno, '_', out_exons, '_', out_consanguinity, '.txt ', path_gene_info, '| cut -d " " -f 1 | sed "s/ID=//g" | sed "s/,//g"'), intern = TRUE)
    system(paste0('rm ', path_data, '/gene_i_', pheno, '_', out_exons, '_', out_consanguinity, '.txt'))
    extract <- which(map[[2]] %in% onto_var_ID)
    ped_onto <- ped[,c(1:6, sort( c(5+(extract*2), 6+(extract*2)) ))]
    ped_inv_onto <- ped_inv[,c(1:6, sort( c(5+(extract*2), 6+(extract*2)) ))]
    map_onto <- map[extract, ]

    if(consanguinity){correction <- "none"}else{correction <- "replace"}
    agg.genos.by.fam <- agg.genos.by.fam(pedfile.path=NULL, pedfile=ped_inv_onto, correction = correction)
    extract_agg <- agg.genos.by.fam$index_variants
    ped_onto <- ped_onto[,c(1:6, sort( c(5+(extract_agg*2), 6+(extract_agg*2)) ))]
    map_onto <- map_onto[extract_agg, ]

    #Set an index for the phenotype status of each subject in each family.
    ped_onto$fam_pheno <- paste0(ped_onto$V1, "_", ped_onto$V6)
    for(onto_var in 1:nrow(map_onto)){
      #Compute the MAF or the number of minor alleles by phenotype status in each family. 
      MAFs_var <- data.frame(do.call(rbind, by(ped_onto, ped_onto$fam_pheno, MAF, var = onto_var, freq = freq, simplify = FALSE)))
      colnames(MAFs_var) <- paste0("MAF_", map_onto$V2[onto_var])
      MAFs_var$fam_pheno <- rownames(MAFs_var)
      if(onto_var == 1){MAFs <- MAFs_var} else {MAFs <- merge(x = MAFs, MAFs_var, by = "fam_pheno")}
    }
    MAFs[, c("fam", "pheno")] <- str_split_fixed(MAFs$fam_pheno, "_", 2)
    #Keep the variants present in at least one affected subject.
    extract_var <- apply(select(MAFs[MAFs$pheno == 2,], -fam, -pheno, -fam_pheno), 2, function(x){ifelse(sum(x==0)!=sum(MAFs$pheno == 2), TRUE, FALSE)})
    MAFs <- MAFs[, c("fam_pheno", names(extract_var)[extract_var], "fam", "pheno")]
    #Remove the families that have none of the remaining variants
    extract_fam_which <- apply(select(MAFs[MAFs$pheno == 2,], -fam, -pheno, -fam_pheno), 1, function(x){ifelse(sum(x==0)!=sum(extract_var), TRUE, FALSE)})
    extract_fam <- gsub("-1|-2", "", MAFs$fam) %in% unique(gsub("-1|-2", "", MAFs$fam[MAFs$pheno == 2][extract_fam_which]))
    MAFs <- merge(x = MAFs, y = aggregate(extract_fam, list(fam = MAFs$fam), function(x){ifelse(any(unlist(x) == 1), TRUE, FALSE)}), by = "fam")
    MAFs <- MAFs[which(MAFs$x), c("fam_pheno", names(extract_var)[extract_var])]
    MAFs_by_onto[[onto]] <- data.frame(t(select(MAFs, -fam_pheno))) %>% setNames(paste0(MAFs$fam_pheno, "_", out_name))
  }
  return(MAFs_by_onto)
}
