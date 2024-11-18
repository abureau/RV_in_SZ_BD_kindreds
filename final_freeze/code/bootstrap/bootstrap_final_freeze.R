# Code for recomputing with the bootstrap the p-values of the best signal of the CRH analysis
library(stringr); library(GenomicRanges); library(RetroFunRVS); library(dplyr)

source("extrait_pathway.R")
source("extrait_CRH_TAD.R")

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

set.seed(10)
CRH998.agg = extract_CRH(TAD.chr21.1.orig,"CRH998")
CRH998.boot = bootRetroFunRVS(CRH998.agg, n.unique.config.by.fam=null$N, prob.sharing.by.fam=null$pattern.prob,nullval=exp_var,nrep=10000)

annotation.matrix = matrix(1, ncol=1, nrow=max(CRH998.agg$index_variants))
score_orig.CRH998 = compute.Burden.by.Annot(exp_var, CRH998.agg,  Z_annot = annotation.matrix, W = rep(1, nrow(annotation.matrix)))$B[1,1]

mean(CRH998.boot>=score_orig.CRH998)
# [1] 0.0023

CRH1006.agg = extract_CRH(TAD.chr21.1.orig,"CRH1006")

# Avec relations cryptiques
TAD.chr21.1.cryptic = agg.genos.annotations.TAD(pheno="GCbr", with_exons=TRUE, consanguinity=TRUE, path_retrofun=path_retrofun, chr=21, TAD=1, maxvar = 300)
null_cryptic_name <- "consanguinity.cryptique.seq.prob.dist.rds"
null.cryptic <- readRDS(paste0(path_retrofun, "/objets_ped/", null_name))
exp_var_cryptic_name <- "expected.variance.consanguinity.cryptique.seq.rds"
exp_var.cryptic <- readRDS(paste0(path_retrofun, "/objets_ped/", exp_var_cryptic_name))

results <- RetroFun.RVS(exp_var_cryptic, TAD.chr21.1.cryptic$agg.genos.by.fam, Z_annot = TAD.chr21.1.cryptic$Z_annot, W = rep(1, nrow(TAD.chr21.1.cryptic$Z_annot)), independence=FALSE)

set.seed(10)
boot.TAD.chr21.1.cryptic = bootRetroFunRVS(TAD.chr21.1.cryptic$agg.genos.by.fam, n.unique.config.by.fam=null.cryptic$N, Z_annot = TAD.chr21.1.cryptic$Z_annot[,-4], prob.sharing.by.fam=null.cryptic$pattern.prob,nullval=exp_var.cryptic,nrep=1000)
Var.Stat = unlist(compute.Var.by.Annot(exp_var.cryptic,TAD.chr21.1.cryptic$agg.genos.by.fam,TAD.chr21.1.cryptic$Z_annot,W= rep(1, nrow(TAD.chr21.1.cryptic$Z_annot)),independence = FALSE))

pmat = pnorm(boot.TAD.chr21.1.cryptic/matrix(sqrt(Var.Stat[-4]),1000,3,byrow = T),0,1,lower.tail = FALSE)

boot.TAD.chr21.1.cryptic.ACAT = apply(pmat,1,function(x) ACAT::ACAT(x[!is.nan(x)]))
mean(boot.TAD.chr21.1.cryptic.ACAT<=results$ACAT)
# [1] 0.001

set.seed(10)
CRH998.cryptic.agg = extract_CRH(TAD.chr21.1.cryptic,"CRH998")
CRH998.cryptic.boot = bootRetroFunRVS(CRH998.cryptic.agg, n.unique.config.by.fam=null.cryptic$N, prob.sharing.by.fam=null.cryptic$pattern.prob,nullval=exp_var.cryptic,nrep=10000)

annotation.matrix = matrix(1, ncol=1, nrow=max(CRH998.agg$index_variants))
score_cryptic.CRH998 = compute.Burden.by.Annot(exp_var, CRH998.cryptic.agg,  Z_annot = annotation.matrix, W = rep(1, nrow(annotation.matrix)))$B[1,1]
mean(CRH998.cryptic.boot>=score_cryptic.CRH998)
#[1] 0.0019

# Code for recomputing with the bootstrap the p-values of the best signal of the pathway analysis

# Sans relations cryptiques
GO.0098939.orig = agg.genos.pathways(pheno="GCbr", with_exons=TRUE, consanguinity=FALSE, path_retrofun=path_retrofun, id="GO:0098939")
set.seed(11)
GO.0098939.boot = bootRetroFunRVS(GO.0098939.orig, n.unique.config.by.fam=null$N, prob.sharing.by.fam=null$pattern.prob,nullval=exp_var,nrep=10000)

annotation.matrix = matrix(1, ncol=1, nrow=max(GO.0098939.orig$index_variants))
score_orig.GO.0098939 = compute.Burden.by.Annot(exp_var, GO.0098939.orig,  Z_annot = annotation.matrix, W = rep(1, nrow(annotation.matrix)))$B[1,1]

mean(GO.0098939.boot>=score_orig.GO.0098939)
# [1] 7e-04

# Avec relations cryptiques
GO.0098939.cryptic = agg.genos.pathways(pheno="GCbr", with_exons=TRUE, consanguinity=TRUE, path_retrofun=path_retrofun, id="GO:0098939")
set.seed(11)
GO.0098939.cryptic.boot = bootRetroFunRVS(GO.0098939.cryptic, n.unique.config.by.fam=null.cryptic$N, prob.sharing.by.fam=null.cryptic$pattern.prob,nullval=exp_var.cryptic,nrep=10000)

annotation.matrix = matrix(1, ncol=1, nrow=max(GO.0098939.cryptic$index_variants))
score_cryptic.GO.0098939 = compute.Burden.by.Annot(exp_var.cryptic, GO.0098939.cryptic,  Z_annot = annotation.matrix, W = rep(1, nrow(annotation.matrix)))$B[1,1]

mean(GO.0098939.cryptic.boot>=score_cryptic.GO.0098939)
#[1] 7e-04
