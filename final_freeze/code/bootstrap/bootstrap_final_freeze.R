# Code for recomputing with the bootstrap the p-values of the best signal of the CRH analysis
library(stringr); library(GenomicRanges); library(RetroFunRVS); library(dplyr)

path_retrofun <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/RetroFunRVS"
#setwd(paste0(path_retrofun, "/TADs"))

#Function to load pedigrees easily.
loadRData <- function(file_name){load(file_name); get(ls()[ls() != "file_name"])}

# Sans relations cryptiques
TAD.chr21.1.orig = agg.genos.annotations.TAD(pheno="GCbr", with_exons=TRUE, consanguinity=FALSE, path_retrofun=path_retrofun, chr=21, TAD=1, maxvar = 300)
dim(TAD.chr21.1.orig$Z_annot)
null_name <- "seq.prob.dist.rds"
null <- readRDS(paste0(path_retrofun, "/objets_ped/", null_name))
exp_var_name = "expected.variance.seq.rds"
exp_var = readRDS(paste0(path_retrofun, "/objets_ped/", exp_var_name))

results <- RetroFun.RVS(exp_var, TAD.chr21.1.orig$agg.genos.by.fam, Z_annot = TAD.chr21.1.orig$Z_annot, W = rep(1, nrow(TAD.chr21.1.orig$Z_annot)), independence=FALSE)

set.seed(10)
boot.TAD.chr21.1.orig = bootRetroFunRVS(TAD.chr21.1.orig$agg.genos.by.fam, n.unique.config.by.fam=null$N, Z_annot = TAD.chr21.1.orig$Z_annot[,-4], prob.sharing.by.fam=null$pattern.prob,nullval=exp_var,nrep=1000)

Var.Stat = unlist(compute.Var.by.Annot(exp_var,TAD.chr21.1.orig$agg.genos.by.fam,TAD.chr21.1.orig$Z_annot,W= rep(1, nrow(TAD.chr21.1.orig$Z_annot)),independence = FALSE))

pmat = pnorm(boot.TAD.chr21.1.orig/matrix(sqrt(Var.Stat[-4]),1000,3,byrow = T),0,1,lower.tail = FALSE)

boot.TAD.chr21.1.orig.ACAT = apply(pmat,1,function(x) ACAT::ACAT(x[!is.nan(x)]))
mean(boot.TAD.chr21.1.orig.ACAT<=results$ACAT)
# [1] 0.001

system.time(bootRetroFunRVS(TAD.chr21.1.orig$agg.genos.by.fam, n.unique.config.by.fam=null$N, Z_annot = TAD.chr21.1.orig$Z_annot[,-4], prob.sharing.by.fam=null$pattern.prob,nullval=exp_var,nrep=1))
system.time(resample.genos.by.fam(agg.genos.by.fam = TAD.chr21.1.orig$agg.genos.by.fam, n.unique.config.by.fam=null$N, prob.sharing.by.fam=null$pattern.prob))

CRH998.agg = extract_CRH(TAD.chr21.1.orig,"CRH998")
CRH998.boot = bootRetroFunRVS(CRH998.agg, n.unique.config.by.fam=null$N, prob.sharing.by.fam=null$pattern.prob,nullval=exp_var,nrep=10000)

annotation.matrix = matrix(1, ncol=1, nrow=max(CRH998.agg$index_variants))
score_orig.CRH998 = compute.Burden.by.Annot(exp_var, CRH998.agg,  Z_annot = annotation.matrix, W = rep(1, nrow(annotation.matrix)))$B[1,1]

mean(CRH998.boot>=score_orig.CRH998)
# [1] 0.0023

CRH1006.agg = extract_CRH(TAD.chr21.1.orig,"CRH1006")

# Avec relations cryptiques
TAD.chr20.1.cryptic = agg.genos.annotations.TAD(pheno="GCbr", with_exons=TRUE, consanguinity=TRUE, path_retrofun=path_retrofun, chr=20, TAD=1, maxvar = 300)
exp_var_cryptic_name <- "expected.variance.consanguinity.cryptique.seq.rds"
exp_var_cryptic <- readRDS(paste0(path_retrofun, "/objets_ped/", exp_var_cryptic_name))

results <- RetroFun.RVS(exp_var_cryptic, TAD.chr20.1.cryptic$agg.genos.by.fam, Z_annot = TAD.chr20.1.cryptic$Z_annot, W = rep(1, nrow(TAD.chr20.1.cryptic$Z_annot)), independence=FALSE)
