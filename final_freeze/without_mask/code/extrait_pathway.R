  ## Function for extracting the aggregated genotypes object for genes from a SynGO pathway
agg.genos.pathways <- function(pheno, with_exons, strict = FALSE, consanguinity, path_retrofun, id){
  #pheno, string, add the phenotype name.
  #with_exons, logical, TRUE if we consider exonic variants ONLY.
  #strict, logical, TRUE if we consider exonic variants ONLY without synonymous variants. Only TRUE if with_exons if TRUE.
  #consanguinity, logical, TRUE if we consider consanguinity loops in the analysis.
  #path_retrofun, string, Path to root folder for required data
  #id, string, Pathway ID
  
  if(strict & !with_exons){stop("Strict cannot be TRUE if with_exons is FALSE")}
  
  #Path to the data
  if(with_exons){
    path_data <- paste0(path_retrofun, "/pathways/exons_only")
    out_exons <- "only_exonic_variants"
    if(strict){
      path_data <- paste0(path_retrofun, "/pathways/exons_only_strict")
      out_exons <- "only_exonic_variants_without_synonymous"
    }
  }else{
    path_data <- paste0(path_retrofun, "/pathways")
    out_exons <- "all_variants"
  }
  file <- "impute5_gigi2_combined_seq_RV_FINAL_genome_all_gene.ped"
  path_gene_info <- paste0(path_data, "/pathways_genes_effects_to_keep_seq_FINAL.txt")
  onto <- openxlsx::read.xlsx(paste0(path_retrofun, "/pathways/syngo_ontologies.xlsx"))
  
  mapfile <- read.table(paste0(path_data, "/" , gsub(".ped", "", file), ".map"), header = FALSE)
  colnames(mapfile) <- c("chrom", "ID", "genpos", "pos")
  pedfile <- read.table(paste0(path_data, "/" , file))
  for(i in 7:ncol(pedfile)){pedfile[,i] <- case_when(pedfile[,i] == 0 ~ 0, pedfile[,i] == 1 ~ 2, pedfile[,i] == 2 ~ 1)}
  
  #Remove 2620b which is a duplicate of 2620 in our data.
  pedfile <- pedfile[pedfile[,2] != "2620b",]
  
  #Load the phenotype files with consanguinity loops.
  #Affected must be coded by a 2.
  df.ped.2021 <- loadRData(paste0(path_retrofun, "/objets_ped/ped", pheno, "_orig.RData"))
  df.ped.2021 <- data.frame(df.ped.2021$famid, df.ped.2021$id, df.ped.2021$sex, df.ped.2021$affected+1)
  colnames(df.ped.2021) <- c("famid","id","sex","affected")
  pedfile.fam.infos <- setNames(data.frame(pedfile[,1], pedfile[,2], pedfile[,3], pedfile[,4]), c("famid", "id", "findex", "mindex"))
  pedfile.fam.infos.full <- merge(x = pedfile.fam.infos, y = df.ped.2021[,c("famid","id","sex","affected")], by = c("famid","id"), sort = FALSE, all.x = TRUE)
  
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
  
    gene_id <- onto$hgnc_symbol[onto$id == id]
    write <- data.table::data.table(gsub("(.*)", "|\\1|", strsplit(gene_id, ", ")[[1]], fixed = FALSE))
    data.table::fwrite(write, paste0("gene_i_", pheno, "_", out_exons, "_", out_consanguinity, ".txt"), col.names = FALSE, row.names = FALSE)
    gene_var_ID <- system(paste0('grep -f gene_i_', pheno, '_', out_exons, '_', out_consanguinity, '.txt ', path_gene_info, '| cut -d " " -f 1 | sed "s/ID=//g" | sed "s/,//g"'), intern = TRUE)
    system(paste0('rm gene_i_', pheno, '_', out_exons, '_', out_consanguinity, '.txt'))
    var_position <- which(mapfile$ID %in% gene_var_ID)
    pedframe <- pedfile[,c(1:6, sort( c(5+(var_position*2), 6+(var_position*2)) ))]
    agg.genos.by.fam(pedfile.path=NULL, pedfile = pedframe, correction=correction)
}
