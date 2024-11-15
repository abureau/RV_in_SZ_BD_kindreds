library(stringr); library(GenomicRanges); library(RetroFunRVS); library(RVS); library(dplyr); library(kinship2); library(RetroFunRVS);library(foreach); library(doParallel)
packageVersion("RVS")
options(scipen=999)

#For parallel computing
args <- commandArgs(TRUE)
(pheno <- as.character(args[[1]])) #Phenotype used
(fam_idx <- as.numeric(args[[2]])) #Family used
(consanguinity <- as.logical(args[[3]])) #Do we take consanguinity into account? TRUE or FALSE?
(cryptic <- as.logical(args[[4]])) #Do we take cryptic relations into account? TRUE or FALSE?
(dist <- as.logical(args[[5]])) #Do we want to output the probability distribution? TRUE or FALSE?
if(dist){out_dist <- ".prob.dist"}else{out_dist <- ""}
path_ped <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped"

loadRData <- function(file_name){load(file_name); get(ls()[ls() != "file_name"])}
seqped2021 <- loadRData(paste0(path_ped, "/ped", pheno, "_orig.RData"))
pedigree <- seqped2021[fam_idx]

#Kinship moyen entre fondateurs
region <- data.table::fread("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/regionsFamilles.txt", header = TRUE)
fam <- unique(pedigree$fam)
kinshipCoeff <- dplyr::case_when(fam %in% region$Famille[region[[2]] == "Saguenay"] ~ 0.010648164,
                                 fam %in% region$Famille[region[[2]] == "Beauce"] ~ 0.005516408,
                                 fam %in% region$Famille[region[[2]] == "NB_Iles"] ~ 0.00172)
print("*** valeur de kinshipCoeff pour la famille en question ***"); print(fam); print(region[[2]][region$Famille == fam]); print(kinshipCoeff)

#Split problematic families. We now have a list of pedigrees.
fam_prob <- fam_idx %in% c(10, 18, 47)
if (fam_prob){
  source("/lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples/call/RetroFunRVS_null_object_split_fam_prob.R")
  print(pedigree); print("*** Pedigree after splitting ***"); print(split_pedigree)
  pedigree <- split_pedigree
  saveRDS(pedigree, paste0(path_ped, "/fam", fam, "splitted.rds"))
}

#multicores.
ncores <- Sys.getenv("SLURM_CPUS_PER_TASK") 
registerDoParallel(cores=ncores)
print(ncores)

if(cryptic){
  if(!consanguinity) {stop("consanguinity needs to be TRUE when cryptic is TRUE")}
  expected.variance.seq.fam <- suppressMessages(compute.null.parallel(pedigree = pedigree, distinguishHomo = consanguinity, cryptic.relatedness = cryptic, kinshipCoeff = kinshipCoeff, out.prob.dist = dist))
  attributes(expected.variance.seq.fam)$distinguishHomo = T
  saveRDS(expected.variance.seq.fam, file = paste0(path_ped, paste0("/expected.variance.consanguinity.cryptique.seq.fam", fam_idx, out_dist, ".rds")))
}else if(consanguinity){
  expected.variance.seq.fam <- suppressMessages(compute.null.parallel(pedigree = pedigree, distinguishHomo = consanguinity, cryptic.relatedness = cryptic, out.prob.dist = dist))
  attributes(expected.variance.seq.fam)$distinguishHomo = T
  saveRDS(expected.variance.seq.fam, file = paste0(path_ped, paste0("/expected.variance.consanguinity.seq.fam", fam_idx, out_dist, ".rds")))
}else{
  expected.variance.seq.fam <- suppressMessages(compute.null.parallel(pedigree = pedigree, distinguishHomo = consanguinity, cryptic.relatedness = cryptic, out.prob.dist = dist))
  saveRDS(expected.variance.seq.fam, file = paste0(path_ped, paste0("/expected.variance.seq.fam", fam_idx, out_dist, ".rds")))
}