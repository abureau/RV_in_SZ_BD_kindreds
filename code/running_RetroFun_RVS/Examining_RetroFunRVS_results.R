# Results with distinguishHomo

# Reading results files
RetroFun.res = readRDS("RetroFun.RVS_results_seq_all_chromosomes_with_consanguinity.RDS")
RetroFun0.res = readRDS("overlap_0/RetroFun.RVS_results_CRHs_seq_all_chromosomes.RDS")
RetroFun2.res = readRDS("overlap_2/RetroFun.RVS_results_CRHs_seq_all_chromosomes.RDS")

# Displaying top results
RetroFun.res[order(RetroFun.res$score)[1:6],]
RetroFun0.res[order(RetroFun0.res$score)[1:6],]
RetroFun2.res[order(RetroFun2.res$score)[1:6],]

# Results without distinguishHomo

# Reading results files
RetroFunOrig.res = readRDS("RetroFun.RVS_results_seq_all_chromosomes_without_consanguinity_replace.RDS")
RetroFunOrig0.res = readRDS("overlap_0_replace/RetroFun.RVS_results_CRHs_seq_all_chromosomes_without_consanguinity_replace.RDS")
RetroFunOrig2.res = readRDS("overlap_2_replace/RetroFun.RVS_results_CRHs_seq_all_chromosomes_without_consanguinity_replace.RDS")

# Displaying top results
RetroFunOrig.res[order(RetroFunOrig.res$score)[1:6],]
RetroFunOrig0.res[order(RetroFunOrig0.res$score)[1:6],]
RetroFunOrig2.res[order(RetroFunOrig2.res$score)[1:6],]
