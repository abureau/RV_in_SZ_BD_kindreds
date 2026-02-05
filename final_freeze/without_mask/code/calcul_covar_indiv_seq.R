# Calcul des variances et covariances par individu sur les distributions de probabilité avec consanguinité existantes
library(RetroFunRVS)
library(kinship2)
pheno="GCbr"

path_ped <- "/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped"
pathAB_ped <- "/lustre09/project/6033529/schizo/data_AB/WGS_bs_2022/500_samples_cag/RetroFunRVS/objets_ped"
#path_retrofun <- "/lustre09/project/6033529/schizo/data/WGS_bs_2022/500_samples_cag_without_mask/RetroFunRVS/"
if (pheno %in% c("GCbr","GCna"))
{
  null_name <- paste0("/expected.variance.consanguinity.cryptique.seq.rds")
  null <- readRDS(paste0(path_ped, null_name))
  ped = readRDS(paste0(pathAB_ped,"/pedigrees_51.RDS"))
}
else 
{
    null_name <- paste0("/expected.variance.consanguinity.cryptique.",pheno,".rds")
    null <- readRDS(paste0(pathAB_ped, null_name))
    loadRData <- function(file_name){load(file_name); get(ls()[ls() != "file_name"])}
    ped <- loadRData(paste0(path_ped, "/ped", pheno, "_orig.RData"))
    ped_48 <- vector("list", nfam)
    
    for (i in seq_along(fam.vec)) {
      # Extraction de  la famille i sous forme de pedigree
      ped_48[[i]] <- ped[fam.vec[i]]
    }
    ped = ped_48
}


# Boucle de lecture des familles
#fam.vec <- names(ped)
#fam.vec <- unique(ped$fam)
#nfam=length(fam.vec)
nfam=48
probs = list()

#  probs[[fam.vec[fam]]] = sapply(tmp$distributions[[fam.vec[fam]]],function(l) l$probs)
for (f in 1:nfam)
{
  if (pheno == "GCbr")
    tmp = readRDS(paste0(path_ped, paste0("/expected.variance.consanguinity.cryptique.seq.fam",f,".prob.dist.rds")))
  else tmp = readRDS(paste0(path_ped, paste0("/expected.variance.consanguinity.cryptique.",pheno,".fam",f,".prob.dist.rds")))
  for (fam in names(tmp$distributions))
    probs[[fam]] = sapply(tmp$distributions[[fam]],function(l) l$probs)
}
saveRDS(probs,file="probs.RDS")
probs = readRDS("probs.RDS")

# Note importante: les familles ne sont pas dans le même ordre dans ped et dans probs
# Valeurs par défaut de eta = 0.95, etap = 0.95
covar.list = compute.null.by.indiv(ped,distinguishHomo=T,eta=0.95,etap=0.95,probs.list=probs)
saveRDS(covar.list,file="covar.list.RDS")
null_update95 = replace.covar.by.fam(null,covar.list)

saveRDS(null_update95, file = paste0(pathAB_ped, paste0("/expected.variance.consanguinity.cryptique.eta95.",pheno,".rds")))

# Valeurs de eta = 0.9, etap = 0.9
covar9.list = compute.null.by.indiv(ped,distinguishHomo=T,eta=0.9,etap=0.9,probs.list=probs)
saveRDS(covar9.list,file="covar9.list.RDS")
null_update9 = replace.covar.by.fam(null,covar9.list)

saveRDS(null_update9, file = paste0(pathAB_ped, paste0("/expected.variance.consanguinity.cryptique.eta9.",pheno,".rds")))

