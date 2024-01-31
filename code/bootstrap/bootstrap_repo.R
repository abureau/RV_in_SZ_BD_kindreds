# Code for recomputing with the bootstrap the p-values of the five signals 
# with an ACAT-combined p-value inferior to the Bonferroni-corrected significance level of 0.05.
# All computations performed with distinguishHomo

library(RetroFunRVS)
load("expected.variance.consanguinity.cryptique.seq.RData")
load("SPAPseqprob_cryptique.RData")
load("seqinfo_pour_bootstrap.RData")


load("agg.genos.by.fam_CRH1280.RData")


set.seed(10)
boot.crypt.CRH1280 = bootRetroFunRVS(agg.genos.by.fam, n.unique.config.by.fam=seqinfo.distinguishomo.N.list, prob.sharing.by.fam=seqinfo.pattern.prob.distinguishHomo.unlist,nullval=expected.variance.consanguinity.cryptique.seq,nrep=10000)

annotation.matrix = matrix(1, ncol=1, nrow=max(agg.genos.by.fam$index_variants))
score_crypt.CRH1280 = compute.Burden.by.Annot(expected.variance.consanguinity.cryptique.seq, 
                                      agg.genos.by.fam,  Z_annot = annotation.matrix, W = rep(1, nrow(annotation.matrix)))$B[1,1]

mean(boot.crypt.CRH1280>=score_crypt.CRH1280)

load("agg.genos.by.fam_CRH1094.RData")

set.seed(10)
boot.crypt.CRH1094 = bootRetroFunRVS(agg.genos.by.fam, n.unique.config.by.fam=seqinfo.distinguishomo.N.list, prob.sharing.by.fam=seqinfo.pattern.prob.distinguishHomo.unlist,nullval=expected.variance.consanguinity.cryptique.seq,nrep=100000)

annotation.matrix = matrix(1, ncol=1, nrow=max(agg.genos.by.fam$index_variants))
score_crypt.CRH1094 = compute.Burden.by.Annot(expected.variance.consanguinity.cryptique.seq, 
                                             agg.genos.by.fam,  Z_annot = annotation.matrix, W = rep(1, nrow(annotation.matrix)))$B[1,1]

# Here the boostrap always returns values inferior to observed score because the level of cryptic relatedness is insufficient to explain homozygous sharing of variants
mean(boot.crypt.CRH1094>=score_crypt.CRH1094)


load("agg.genos.by.fam_CRH694.RData")

set.seed(10)
boot.crypt.CRH694 = bootRetroFunRVS(agg.genos.by.fam, n.unique.config.by.fam=seqinfo.distinguishomo.N.list, prob.sharing.by.fam=seqinfo.pattern.prob.distinguishHomo.unlist,nullval=expected.variance.consanguinity.cryptique.seq,nrep=10000)

annotation.matrix = matrix(1, ncol=1, nrow=max(agg.genos.by.fam$index_variants))
score_crypt.CRH694 = compute.Burden.by.Annot(expected.variance.consanguinity.cryptique.seq, 
                                      agg.genos.by.fam,  Z_annot = annotation.matrix, W = rep(1, nrow(annotation.matrix)))$B[1,1]

mean(boot.crypt.CRH694>=score_crypt.CRH694)


load("agg.genos.by.fam_CRH277.RData")

set.seed(10)
boot.crypt.CRH277 = bootRetroFunRVS(agg.genos.by.fam, n.unique.config.by.fam=seqinfo.distinguishomo.N.list, prob.sharing.by.fam=seqinfo.pattern.prob.distinguishHomo.unlist,nullval=expected.variance.consanguinity.cryptique.seq,nrep=10000)

annotation.matrix = matrix(1, ncol=1, nrow=max(agg.genos.by.fam$index_variants))
score_crypt.CRH277 = compute.Burden.by.Annot(expected.variance.consanguinity.cryptique.seq, 
                                      agg.genos.by.fam,  Z_annot = annotation.matrix, W = rep(1, nrow(annotation.matrix)))$B[1,1]

mean(boot.crypt.CRH277>=score_crypt.CRH277)


load("agg.genos.by.fam_CRH1438.RData")

set.seed(10)

boot.crypt.CRH1438 = bootRetroFunRVS(agg.genos.by.fam, n.unique.config.by.fam=seqinfo.distinguishomo.N.list, prob.sharing.by.fam=seqinfo.pattern.prob.distinguishHomo.unlist,nullval=expected.variance.consanguinity.cryptique.seq,nrep=100000)

annotation.matrix = matrix(1, ncol=1, nrow=max(agg.genos.by.fam$index_variants))
score_crypt.CRH1438 = compute.Burden.by.Annot(expected.variance.consanguinity.cryptique.seq, 
                                        agg.genos.by.fam,  Z_annot = annotation.matrix, W = rep(1, nrow(annotation.matrix)))$B[1,1]
mean(boot.crypt.CRH1438>=score_crypt.CRH1438)


# Gathering the results

pboot_crypt.vec = c(CRH1280=mean(boot.crypt.CRH1280>=score_crypt.CRH1280), CRH1094=mean(boot.crypt.CRH1094>=score_crypt.CRH1094),CRH694=mean(boot.crypt.CRH694>=score_crypt.CRH694),CRH277=mean(boot.crypt.CRH277>=score_crypt.CRH277), CRH1438=mean(boot.crypt.CRH1438>=score_crypt.CRH1438))
save(pboot_crypt.vec,file="pboot_crypt.RData")