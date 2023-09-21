#module load r/4.2 StdEnv/2020
pathdata <- "/lustre03/project/6033529/quebec_10x/data/freeze/QC/mendel_corrected"
setwd(pathdata)
library(stringr); library(GenomicRanges); library(dplyr)

#Import pedigree data with the phenotype
df.ped.2021 = read.table("/lustre03/project/6033529/quebec_10x/data/freeze/GCbroad_seq_inbred.pre", header=FALSE, sep=" ")
colnames(df.ped.2021) = c("famid","id","fid","mid","sex","affected")
subset.fam = c("103","105","110","115","119","121","124","125","126","129","131","133","151","182","207","210","211","212","217","220","224","228","230","233","234","235","238","255")
df.ped.2021 = df.ped.2021[df.ped.2021$famid%in%subset.fam,]
#Remove subject 791 with a missing phenotype.
df.ped.2021 <- df.ped.2021[df.ped.2021$id != 791,]

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



#Are we using the results from the analysis with consanguinity? If yes, set `consanguinity` to TRUE, else, set it to FALSE.
consanguinity <- TRUE

if(consanguinity){
  results <- readRDS("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/RetroFun.RVS_results_seq_all_chromosomes_with_consanguinity.RDS")
  pathout <- "/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/MAF_n_variants_10_CRH_sign_by_fam_pheno_with_consanguinity.RDS"
  }else{
  results <- readRDS("/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/RetroFun.RVS_results_seq_all_chromosomes_without_consanguinity_replace.RDS")
  pathout <- "/lustre03/project/6033529/quebec_10x/results/RetroFunRVS/MAF_n_variants_10_CRH_sign_by_fam_pheno.RDS"
}

#Use the 10 most significant CRHs.
results <- results[order(results$score),]
sign_CRH <- results[1:10, 1:3]

#This function is used to output the MAF(freq=TRUE) or to output the minor alleles(freq=FALSE) in a list of data.frames.
sign_CRH_results <- function(freq){
  #Initialize the results object.
  MAFs_by_TAD <- vector(mode = "list", length = nrow(sign_CRH))
  names(MAFs_by_TAD) <- sign_CRH$TAD_name
  if(freq){out_name <- "MAF"}else{out_name <- "n"}
  for(CRH in 1:nrow(sign_CRH)){
    #Path to the TAD data where the present CRH is found.
    datafile <- paste0("seq_FINAL_chr_", sign_CRH$chr[CRH], "_TAD_", sign_CRH$TAD[CRH], "_with_mask")
    TADs_chr <- data.table::fread(paste0(pathdata, "/TADs_list_chr", sign_CRH$chr[CRH], ".bed"), header = FALSE)
    
    #Import and split the TAD data to keep the variants in this CRH only.
    CRHs.by.TAD <- read.table(paste0("/lustre03/project/6033529/quebec_10x/data/CRHs_iPSC_neurons/liftover_hg38/TADs_for_CRHs_overlap_1_TAD/chr", sign_CRH$chr[CRH], "/chr", sign_CRH$chr[CRH], "_", TADs_chr$V2[sign_CRH$TAD[CRH]]+1, "_", TADs_chr$V3[sign_CRH$TAD[CRH]], ".bed"), header=TRUE)
    map <- read.table(paste0(pathdata, "/", datafile, ".map"))
    ped <- read.table(paste0(pathdata, "/", datafile, ".ped"))
    extract <- sapply(map$V4, function(x) is_within_any_range(x, CRHs.by.TAD[paste0("CRH", CRHs.by.TAD$name) == sign_CRH$TAD_name[CRH], c("chromStart", "chromEnd")]))
    #Just in case it happens that no variants are found in the CRH (it shouldn't happen).
    if(sum(extract)==0){MAFs_by_TAD[[CRH]] <- "no variants are found in this CRH"}
    ped_CRH <- ped[,c(1:6, sort( c(5+(which(extract)*2), 6+(which(extract)*2)) ))]
    map_CRH <- map[extract, ]
    #Add the phenotype to the ped file.
    ped_CRH <- merge(x = ped_CRH, y = df.ped.2021[,c("famid","id","affected")], by.x = c("V1", "V2"), by.y=c("famid","id"), sort = FALSE)
    ped_CRH$V6 <- ped_CRH$affected; ped_CRH$affected <- NULL
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
    extract_fam <- apply(select(MAFs, -fam, -pheno, -fam_pheno), 1, function(x){ifelse(sum(x==0)!=sum(extract_var), TRUE, FALSE)})
    MAFs <- merge(x = MAFs, y = aggregate(extract_fam, list(fam = MAFs$fam), function(x){ifelse(any(unlist(x) == 1), TRUE, FALSE)}), by = "fam")
    MAFs <- MAFs[which(MAFs$x), c("fam_pheno", names(extract_var)[extract_var])]
    MAFs_by_TAD[[CRH]] <- data.frame(t(select(MAFs, -fam_pheno))) %>% setNames(paste0(MAFs$fam_pheno, "_", out_name))
  }
  return(MAFs_by_TAD)
}
#Output the MAF.
MAFs_by_TAD <- sign_CRH_results(freq = TRUE)
#Output the number of minor alleles.
n_by_TAD <- sign_CRH_results(freq = FALSE)
#Bind the results.
results_by_TAD <- Map(cbind, MAFs_by_TAD, n_by_TAD)
sort(colnames(results_by_TAD[[1]]))
results_by_TAD <- purrr::map(results_by_TAD, function(x){x[sort(colnames(x))]})
#Save the results.
saveRDS(results_by_TAD, pathout)

