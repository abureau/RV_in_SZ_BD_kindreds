# Functions used to extract objects for specific TADs and CRHs

# Function for extracting the aggregated genotypes object and annotation matrix of a TAD encompassing entire CHRs, 
# i.e. the CRH overlap only that TAD.
agg.genos.annotations.TAD <- function(pheno, with_exons, consanguinity, path_retrofun, chr, TAD, maxvar = 300){
  #pheno, string, add the phenotype name.
  #with_exons, logical, TRUE if we add the exons to the CRHs for the analysis.
  #consanguinity, logical, TRUE if we consider consanguinity loops in the analysis.
  #maxvar, integer, maximum normal of variants to test together.
  
  TADs <- data.table::fread(paste0(path_retrofun, "/TADs/TADs_list_chr", chr, ".bed"), header = FALSE)
  datafile <- paste0("impute5_gigi2_combined_seq_RV_FINAL_chr", chr, "_TAD_", TAD)
  if(with_exons){
    #If the addition of exons had an impact on the .ped, import the modified one.
    if(file.exists(paste0(path_retrofun, "/TADs/with_exons/", datafile, ".ped"))){datafile <- paste0("with_exons/", datafile)}
  }

  #Load the .ped of the TAD
  #Minor alleles must be coded as 2. PLINK codes minor to 1 and major to 2, so we must change them.
  pedfile <- read.table(paste0(path_retrofun, "/TADs/", datafile, ".ped"))
  for(i in 7:ncol(pedfile)){pedfile[,i] <- case_when(pedfile[,i] == 0 ~ 0, pedfile[,i] == 1 ~ 2, pedfile[,i] == 2 ~ 1)}
  
  #Remove 2620b which is a duplicate of 2620 in our data.
  pedfile <- pedfile[pedfile[,2] != "2620b",]
  
  #Load the phenotype files with consanguinity loops.
  #Affected must be coded by a 2.
  df.ped.2021 <- loadRData(paste0(path_retrofun, "/objets_ped/ped", pheno, "_orig.RData"))
  df.ped.2021 <- data.frame(df.ped.2021$famid, df.ped.2021$id, df.ped.2021$findex, df.ped.2021$mindex, df.ped.2021$sex, df.ped.2021$affected+1)
  colnames(df.ped.2021) <- c("famid","id","findex","mindex","sex","affected")
  pedfile.fam.infos <- setNames(data.frame(pedfile[,1], pedfile[,2]), c("famid", "id"))
  pedfile.fam.infos.full <- merge(x = pedfile.fam.infos, y = df.ped.2021[,c("famid","id","findex","mindex","sex","affected")], by = c("famid","id"), sort = FALSE, all.x = TRUE)
  
  #Associate the phenotype in the .ped and keep the affected families only.
  subset.fam <- pedfile.fam.infos.full %>% group_by(famid) %>% summarise(n_affected = sum(affected==2)) %>% filter(n_affected!=0) %>% select(famid) %>% as.vector()
  df.ped.2021 <- df.ped.2021[df.ped.2021$famid %in% subset.fam$famid,]
  pedfile.fam.infos.full <- pedfile.fam.infos.full[pedfile.fam.infos.full$famid %in% subset.fam$famid,]
  #Check
  pedfile <- pedfile[pedfile$V2 %in% pedfile.fam.infos.full$id, ]
  all(pedfile$V2 == pedfile.fam.infos.full$id)
  pedfile[,1:6] <- pedfile.fam.infos.full
  
  #Adjust the pedigree for the 3 problematic families in agg.genos.by.fam.
  if(consanguinity){
    correction <- "none"
    out_consanguinity <- "with_consanguinity"
  } else {
    correction <- "replace"
    out_consanguinity <- "without_consanguinity"
  }
  fam_split_119 <- readRDS(paste0(path_retrofun, "/objets_ped/fam119splitted.rds"))
  pedfile$V1[pedfile$V1 == "119" & pedfile$V2 %in% fam_split_119$id[fam_split_119$fam=="119-1"]] <- "119-1"
  pedfile$V1[pedfile$V1 == "119" & pedfile$V2 %in% fam_split_119$id[fam_split_119$fam=="119-2"]] <- "119-2"
  fam_split_131 <- readRDS(paste0(path_retrofun, "/objets_ped/fam131splitted.rds"))
  pedfile$V1[pedfile$V1 == "131" & pedfile$V2 %in% fam_split_131$id[fam_split_131$fam=="131-1"]] <- "131-1"
  pedfile$V1[pedfile$V1 == "131" & pedfile$V2 %in% fam_split_131$id[fam_split_131$fam=="131-2"]] <- "131-2"
  fam_split_255 <- readRDS(paste0(path_retrofun, "/objets_ped/fam255splitted.rds"))
  pedfile$V1[pedfile$V1 == "255" & pedfile$V2 %in% fam_split_255$id[fam_split_255$fam=="255-1"]] <- "255-1"
  pedfile$V1[pedfile$V1 == "255" & pedfile$V2 %in% fam_split_255$id[fam_split_255$fam=="255-2"]] <- "255-2"
  
  #Create annotation files using the .map.
  mapfile <- read.table(paste0(path_retrofun, "/TADs/", datafile, ".map"), header=FALSE)
  variants <- str_split_fixed(mapfile$V2, ":", 4)
  GRanges.variants <-  GRanges(seqnames=variants[,1], ranges=IRanges(start=as.numeric(variants[,2]),end=as.numeric(variants[,2])))
  
  #TAD and CRHs matching
  if(with_exons){
    CRHs.by.TAD <- read.table(paste0(path_retrofun, "/TADs/with_exons/TADs_for_CRHs_overlap_1_TAD/chr", chr, "/chr", chr, "_", TADs$V2[TAD]+1, "_", TADs$V3[TAD], ".bed"), header=FALSE)
    if(CRHs.by.TAD[1,1] == "chrom"){CRHs.by.TAD <- CRHs.by.TAD[2:nrow(CRHs.by.TAD),]} #It appears that sometimes, the header is part of the dataframe, remove the line.
    colnames(CRHs.by.TAD) <- c("chrom", "chromStart", "chromEnd", "name")
    out_exons <- "CRHs_with_exons"
  } else {
    CRHs.by.TAD <- read.table(paste0("/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/liftover_hg38_executable/TADs_for_CRHs_overlap_1_TAD/chr", chr, "/chr", chr, "_", TADs$V2[TAD]+1, "_", TADs$V3[TAD], ".bed"), header=TRUE)
    out_exons <- "CRHs_only"
  }
  
  #Split by CRH
  GRanges.CRHs.by.TAD <- GRanges(seqnames=CRHs.by.TAD$chrom, ranges=IRanges(start=CRHs.by.TAD$chromStart, end=CRHs.by.TAD$chromEnd),name=CRHs.by.TAD$name)
  split.by.CRH <- split(GRanges.CRHs.by.TAD,GRanges.CRHs.by.TAD$name)
  
  #Annotation matrix, the first column is the Burden
  annotation.matrix <- matrix(0, ncol=length(split.by.CRH), nrow=length(GRanges.variants))
  for(i in 1:length(split.by.CRH)){
    c <- countOverlaps(GRanges.variants,split.by.CRH[[i]])
    c[c>1] <- 1
    annotation.matrix[,i] <- c 
  }
  annotation.matrix = cbind(1, annotation.matrix)
  agg.genos.by.fam <- agg.genos.by.fam(pedfile.path=NULL, pedfile=pedfile, Z_annot=annotation.matrix, correction = correction)
  
  nvarTAD <- length(agg.genos.by.fam$index_variants)
  cat(nvarTAD,"\n")
  nw <- (nvarTAD-1)%/%maxvar + 1
  
  # Si le nombre de variants > maxvar, on decoupe en fenetre de maxvar
  if (nvarTAD > maxvar) {
    # Nombre de fen?tres
    wmat = matrix(0,ncol = nw, nrow=length(GRanges.variants))
    for (w in 1:(nw-1)) {
      wmat[agg.genos.by.fam$index_variants[maxvar*(w-1)+1]:agg.genos.by.fam$index_variants[maxvar*w],w] = 1
    }
    wmat[agg.genos.by.fam$index_variants[maxvar*(nw-1)+1]:nrow(wmat),nw] = 1
    # Ici il faut retirer la colonne de 1 initiale
    annotation.matrix = cbind(wmat, annotation.matrix[,-1])
    # On recalcule le nombre de variants par fenetre
    agg.genos.by.fam <- agg.genos.by.fam(pedfile.path=NULL, pedfile=pedfile, Z_annot=annotation.matrix, exclude_annot = 1:nw, correction = correction)
  }
  
  df.annotation <- setNames(data.frame(annotation.matrix), c(paste0("Burden_",1:nw), paste0("CRH",names(split.by.CRH))))
  
  list(agg.genos.by.fam=agg.genos.by.fam, Z_annot = df.annotation)
}

# Function for extracting a CRH from an object returned by agg.genos.annotations.TAD
extract_CRH = function(TAD.obj,selec)
{
  iv = which(TAD.obj$Z_annot[TAD.obj$agg.genos.by.fam$index_variants,selec]==1)
  if (length(iv)>0) 
  {
    temp = TAD.obj$agg.genos.by.fam$ped_agg[,c(1,iv+1)]
    fam.count = apply(temp[,-1,drop=F],1,sum)
    temp = temp[fam.count>0,,drop=F]
    list(ped_agg=temp,index_variants=TAD.obj$agg.genos.by.fam$index_variants[iv])
  }
  else list(ped_agg=NULL,index_variants=NULL)
}