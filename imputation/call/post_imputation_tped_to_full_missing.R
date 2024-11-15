library(data.table)

#for parallel computing in sbatch
args <- commandArgs(TRUE)
chr <- as.numeric(args[1])

path_impute <- paste0("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/imputation_fam/chr", chr)
ped <- fread(paste0(path_impute, "/geno_impute_chr", chr, ".tped"))
nobs <- ncol(ped)-4
for(id in 1:(nobs/2)){
  ped_id <- ped[, (3+(id*2)):(4+(id*2))]
  to_change <- which((rowSums(ped_id == 0) %% 2) == 1)
  if(length(to_change) == 0){
    next
  } else {
    ped[to_change, (3+(id*2)):(4+(id*2))] <- 0
  }
}
fwrite(ped, paste0(path_impute, "/geno_impute_chr", chr, "_full_missing.tped"), row.names = FALSE, col.names = FALSE, sep = "\t")
