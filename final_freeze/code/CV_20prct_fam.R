library(data.table); library(dplyr)

path_rvs <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS"
max_fam <- ceiling(0.20*48) #We have 48 families.

#Parallel by chromosome
args <- commandArgs(TRUE)
chr <- as.numeric(args[1])

  to_exclude <- c()
  print(paste0("chr : ", chr))
  geno <- fread(paste0(path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.ped'))
  map <- fread(paste0(path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.map'))
  freq <- fread(paste0(path_rvs, '/impute5_gigi2_combined_seq_RV_chr', chr, '.frq'))
  CV <- which(freq$MAF>0.01)
  geno <- dplyr::select(geno, c(1:6, sort(c(5+(CV*2), 6+(CV*2)))))
  map <- map[CV]; freq <- freq[CV]
  n_var <- length(CV)
  for(var in 1:n_var){
    var_i <- setNames(dplyr::select(geno, c(1, 5+(var*2), 6+(var*2))), c("FID", "var_geno1", "var_geno2"))
    var_i <- data.table::melt(var_i, id.vars = "FID")
    n_fam_minor_present <- sum(aggregate(x = var_i$value, by = list(var_i$FID), FUN = function(x){any(unlist(x)==1)})$x)
    if(n_fam_minor_present >= max_fam){to_exclude <- c(to_exclude, map$V2[var])}
  }
  #Remove the concerned variants.
  fwrite(data.table(to_exclude), paste0(path_rvs, "/filtered_CV_to_exclude_chr_", chr, ".snplist"), row.names = FALSE, col.names = FALSE)
  system(paste0("plink --file ", path_rvs, "/impute5_gigi2_combined_seq_RV_chr", chr, " --exclude ", path_rvs, "/filtered_CV_to_exclude_chr_", chr, ".snplist --keep-allele-order --nonfounders --freq --recode --out ", path_rvs, "/impute5_gigi2_combined_seq_RV_FINAL_chr", chr))
  system(paste0("rm ", path_rvs, "/impute5_gigi2_combined_seq_RV_chr", chr, ".ped ", path_rvs, "/impute5_gigi2_combined_seq_RV_chr", chr, ".map"))
