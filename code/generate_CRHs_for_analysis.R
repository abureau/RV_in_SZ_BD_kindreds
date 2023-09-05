#Script permettant de générer les CRHs et TADs utilisés dans l'analyse de données
#Il convient d'importer les fonctions du fichier CRH.R pour préprocesser les données de CRHs
#Les chemins d'accès et de sortie sont à changer
library(igraph)
library(GenomicRanges)
library(rtracklayer)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ensembldb)
library(EnsDb.Hsapiens.v86)

#Fichier nécessaire pour le liftover de hg19 vers hg38
chainObject = import.chain("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\hg19ToHg38.over.chain")

#Importation des fichiers de CRHs sous format igraph + Fichiers de sortie du score ABC
#CRHs_iPSC_NEU = readRDS("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\CRHs_iPSC_NEU.rds")
ABC_Pred = process.ABC("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\EnhancerPredictions_NEU_iPSC.txt")

#Liftover des Enhancers et Gènes de hg19 vers hg38
ABC_Pred.for.liftover = ABC_Pred

Enhancers.for.liftover = unique(ABC_Pred.for.liftover[,c("chr.x","start","end","name")])

GRanges.Enhancers.for.liftover = GRanges(seqnames = Enhancers.for.liftover$chr.x, ranges=IRanges(start=Enhancers.for.liftover$start, end=Enhancers.for.liftover$end),name=Enhancers.for.liftover$name)
GRanges.Enhancers.hg38 = liftOver(GRanges.Enhancers.for.liftover, chainObject)

Genes.for.liftover = na.omit(unique(ABC_Pred.for.liftover[,c("chr.x","startProm","endProm","TargetGene")]))

GRanges.Genes.for.liftover = GRanges(seqnames = Genes.for.liftover$chr.x, ranges=IRanges(start=Genes.for.liftover$start, end=Genes.for.liftover$end),name=Genes.for.liftover$TargetGene)
GRanges.Genes.hg38 = liftOver(GRanges.Genes.for.liftover, chainObject)

table(sapply(GRanges.Genes.hg38,length))
#0     1     2     3     4     5     6     9    13 
#3 13491    72    32    15     3     2     1     1

table(sapply(GRanges.Enhancers.hg38,length))
#0     1     2     3     4     5     7     9    18    21 
#1 30686    78    15     8     7     1     1     1     1
#La plupart des gènes et enhancers sont en correspondance 1-1 vers hg38

GRanges.Enhancers.hg38.clean = lapply(1:length(GRanges.Enhancers.for.liftover), function(x){
  chr = Enhancers.for.liftover[x,"chr.x"]
  GRanges.Enhancers.hg38[[x]][seqnames(GRanges.Enhancers.hg38[[x]])==chr]
})

GRanges.Genes.hg38.clean = lapply(1:length(GRanges.Genes.for.liftover), function(x){
  chr = Genes.for.liftover[x,"chr.x"]
  GRanges.Genes.hg38[[x]][seqnames(GRanges.Genes.hg38[[x]])==chr]
})

index.null.liftover.Enhancers.hg38 = which(sapply(GRanges.Enhancers.hg38.clean, length)==0)
index.null.liftover.Genes.hg38 = which(sapply(GRanges.Genes.hg38.clean, length)==0)

GRanges.Enhancers.hg38.clean = GRanges.Enhancers.hg38.clean[-index.null.liftover.Enhancers.hg38]
GRanges.Genes.hg38.clean = GRanges.Genes.hg38.clean[-index.null.liftover.Genes.hg38]
which(sapply(GRanges.Genes.hg38.clean,function(x) x$name)=="GSTT1")
#On concatène les cas ou le liftover correspond à plusieurs régions dans hg38: Enhancers
min_start_Enhancers_liftoverhg38 = sapply(GRanges.Enhancers.hg38.clean, function(x) {
  min(start(x))
})

max_end_Enhancers_liftoverhg38 = sapply(GRanges.Enhancers.hg38.clean, function(x) {
  max(end(x))
})


GRanges_Enhancers.hg38.final = GRanges(seqnames=Enhancers.for.liftover[-index.null.liftover.Enhancers.hg38,"chr.x"],
                                      ranges=IRanges(start=min_start_Enhancers_liftoverhg38,end=max_end_Enhancers_liftoverhg38), name=Enhancers.for.liftover[-index.null.liftover.Enhancers.hg38,"name"])

#On concatène les cas ou le liftover correspond à plusieurs régions dans hg38: Gènes
min_start_Genes_liftoverhg38 = sapply(GRanges.Genes.hg38.clean, function(x) {
  min(start(x))
})

max_end_Genes_liftoverhg38 = sapply(GRanges.Genes.hg38.clean, function(x) {
  max(end(x))
})

GRanges_Genes.hg38.final = GRanges(seqnames=Genes.for.liftover[-index.null.liftover.Genes.hg38,"chr.x"],
                                       ranges=IRanges(start=min_start_Genes_liftoverhg38,end=max_end_Genes_liftoverhg38), TargetGene=Genes.for.liftover[-index.null.liftover.Genes.hg38,"TargetGene"])
df_Enhancers.hg38.final = as.data.frame(GRanges_Enhancers.hg38.final)
colnames(df_Enhancers.hg38.final) = c("chr.x","start.hg38","end.hg38","widthEnh.hg38","strand","name")
df_Genes.hg38.final = as.data.frame(GRanges_Genes.hg38.final)
colnames(df_Genes.hg38.final) = c("chr.x","startProm.hg38","endProm.hg38","widthProm.hg38","strand","TargetGene")

ABC_Pred.for.liftover = merge(merge(ABC_Pred.for.liftover, df_Enhancers.hg38.final[,c("chr.x","start.hg38","end.hg38", "name")],by=c("chr.x", "name"), all.x=TRUE), df_Genes.hg38.final[,c("chr.x","startProm.hg38","endProm.hg38","TargetGene")], by=c("chr.x","TargetGene"),all.x=TRUE)
ABC_Pred.for.liftover = ABC_Pred.for.liftover[,c("chr.x","start.hg38","end.hg38", "name","startProm.hg38","endProm.hg38","TargetGene","TargetGeneTSS", "CellType","ABC.Score",     "typeOf", "chr.y")]
colnames(ABC_Pred.for.liftover)=c("chr.x","start", "end", "name", "startProm","endProm", "TargetGene","TargetGeneTSS", "CellType","ABC.Score",     "typeOf", "chr.y")

CRHs_iPSC_NEU.liftover.hg38=create.CRHs(ABC_Pred.for.liftover, "ABC")
CRHs_iPSC_NEU = create.CRHs(ABC_Pred, "ABC")

#On importe les données de TADs
TADs_DI_iPSC_NEU = import("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\TADs_DI_iPSC_NEU.bed", format = "bed")

#Liftover des TADs de hg19 vers hg38
TADs_DI_iPSC_NEU_hg38 = liftOver(TADs_DI_iPSC_NEU, chainObject)

length(which(sapply(TADs_DI_iPSC_NEU_hg38, length)==0))#2 TADs pour lesquels on a pas de correspondances sur hg38
TADs_DI_iPSC_NEU.sub = TADs_DI_iPSC_NEU[-which(sapply(TADs_DI_iPSC_NEU_hg38, length)==0)]

TADs_DI_iPSC_NEU_hg38 = TADs_DI_iPSC_NEU_hg38[which(sapply(TADs_DI_iPSC_NEU_hg38, length)>0)]
df_TADs_DI_iPSC_NEU = as.data.frame(TADs_DI_iPSC_NEU.sub)

#Ici on s'assure que les niveaux soient identiques entre les deux objets
df_TADs_DI_iPSC_NEU$seqnames = factor(df_TADs_DI_iPSC_NEU$seqnames, levels=levels(seqnames(TADs_DI_iPSC_NEU_hg38[[1]])))


#On retire certaines inconsistences: liftover pas sur le même chromosome
TADs_DI_iPSC_NEU_hg38.clean = lapply(1:length(TADs_DI_iPSC_NEU.sub), function(x){
  chr = df_TADs_DI_iPSC_NEU$seqnames[x]
  TADs_DI_iPSC_NEU_hg38[[x]][seqnames(TADs_DI_iPSC_NEU_hg38[[x]])==chr]
})

#On concatène les cas ou le liftover correspond à plusieurs régions dans hg38
min_start_TAD_liftoverhg38 = sapply(TADs_DI_iPSC_NEU_hg38.clean, function(x) {
  min(start(x))
})
index_inconsistencies.start = which(min_start_TAD_liftoverhg38==Inf|min_start_TAD_liftoverhg38==-Inf)

min_start_TAD_liftoverhg38 = min_start_TAD_liftoverhg38[-index_inconsistencies.start]
max_end_TAD_liftoverhg38 = sapply(TADs_DI_iPSC_NEU_hg38.clean, function(x) {
  max(end(x))
})
index_inconsistencies.end = which(max_end_TAD_liftoverhg38==Inf|max_end_TAD_liftoverhg38==-Inf)

max_end_TAD_liftoverhg38 = max_end_TAD_liftoverhg38[-index_inconsistencies.end]


TADs_DI_iPSC_NEU_hg38.final = GRanges(seqnames=df_TADs_DI_iPSC_NEU[-unique(c(index_inconsistencies.start,index_inconsistencies.end)),"seqnames"],
        ranges=IRanges(start=min_start_TAD_liftoverhg38,end=max_end_TAD_liftoverhg38))


decompose_CRHsTADs_DI_iPSC_NEU_iPSC = decompose(create.CRHs(ABC_Pred, "ABC"))

#On ajoute l'appartenance au CRH pour chaque élément du fichier ABC_pred
ABC_Pred_with_membership = add.membership(CRHs_iPSC_NEU,ABC_Pred)
ABC_Pred_with_membership.liftover.hg38 = add.membership(CRHs_iPSC_NEU.liftover.hg38,ABC_Pred.for.liftover)

#Le fichier est découpé sur la base du CRH
split_CRHs_NEU = split(ABC_Pred_with_membership, ABC_Pred_with_membership$membership)

split_CRHs_NEU.liftover.hg38 = split(ABC_Pred_with_membership.liftover.hg38, ABC_Pred_with_membership.liftover.hg38$membership)

#On crée un fichier de CRH correspondant aux éléments les plus éloignés du CRH
CRHs_clusters = coverage.By.CRH(CRHs_iPSC_NEU,ABC_Pred_with_membership ,method="ABC")$clusters
CRHs_clusters.liftover.hg38 = coverage.By.CRH(CRHs_iPSC_NEU.liftover.hg38,ABC_Pred_with_membership.liftover.hg38 ,method="ABC")$clusters

#Correspondance entre les CRHs et les TADs
index_overlaps_CRHs_TADs = findOverlaps(CRHs_clusters, TADs_DI_iPSC_NEu)
overlaps_CRHs_TADs = countOverlaps(CRHs_clusters, TADs_DI_iPSC_NEu)

index_overlaps_CRHs_TADs.liftover.hg38 = findOverlaps(CRHs_clusters.liftover.hg38, TADs_DI_iPSC_NEU_hg38.final)
overlaps_CRHs_TADs.liftover.hg38 = countOverlaps(CRHs_clusters.liftover.hg38, TADs_DI_iPSC_NEU_hg38.final)


#Nombre de TADs pour lesquels on a des CRHs qui ne chevauchent qu'1 seul TAD
unique(subjectHits(index_overlaps_CRHs_TADs[queryHits(index_overlaps_CRHs_TADs)%in%which(overlaps_CRHs_TADs%in%c(1))]))
#Nombre de TADs pour lesquels on a des CRHs qui chevauchent qu'2 TADs
unique(subjectHits(index_overlaps_CRHs_TADs[queryHits(index_overlaps_CRHs_TADs)%in%which(overlaps_CRHs_TADs%in%c(2))]))

CRHs_wich_overlap_0TAD = which(overlaps_CRHs_TADs==0)
CRHs_wich_overlap_1TAD = which(overlaps_CRHs_TADs==1)
CRHs_wich_overlap_2TADs = which(overlaps_CRHs_TADs==2)

queryHits(index_overlaps_CRHs_TADs)[queryHits(index_overlaps_CRHs_TADs)%in%CRHs_wich_overlap_0TAD]
table(table(queryHits(index_overlaps_CRHs_TADs)[queryHits(index_overlaps_CRHs_TADs)%in%CRHs_wich_overlap_1TAD]))
table(table(queryHits(index_overlaps_CRHs_TADs)[queryHits(index_overlaps_CRHs_TADs)%in%CRHs_wich_overlap_2TADs]))
#Les résultats sont cohérents avec ce qu'il est attendu

colnames.GE = c("chrom","chromStart", "chromEnd", "name")

list.Genes.Enhancers.by.CRH = lapply(1:length(split_CRHs_NEU), function(x){
  Promoters = na.omit(unique(split_CRHs_NEU[[x]][,c("chr.x", "startProm", "endProm", "membership")]))
  Enhancers = na.omit(unique(split_CRHs_NEU[[x]][,c("chr.x","start", "end", "membership")]))
  
  colnames(Promoters) = colnames(Enhancers) = colnames.GE
  PE = rbind(Promoters,Enhancers)
  PE
})

sum(sapply(1:1632, function(CRH) min(list.Genes.Enhancers.by.CRH[[CRH]]$chromStart) == start(CRHs_clusters[CRH]) & max(list.Genes.Enhancers.by.CRH[[CRH]]$chromEnd) == end(CRHs_clusters[CRH]))) #== 1632 le nombre de CRHs dans les données 
#Pour l'ensemble des CRHs il y a une parfaite correspondance entre les objets list.Genes.Enhancers.by.CRH et CRHs_clusters
#Les débuts et fin des CRHs coincident parfaitement 

#CRHs qui chevauchent 2 TADs
#Création des fichiers bed 
for(CRH in CRHs_wich_overlap_2TADs){
  chr = unique(list.Genes.Enhancers.by.CRH[[CRH]]$chrom)
  write.table(list.Genes.Enhancers.by.CRH[[CRH]][,1:3], paste0("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\CRHs_overlap_2_TADs\\CRH_NEU_iPSC_",CRH,"_",chr,".bed"), col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")
}

#CRHs qui ne chevauchent aucun TADs
#Création des fichiers bed
for(CRH in CRHs_wich_overlap_0TAD){
  chr = unique(list.Genes.Enhancers.by.CRH[[CRH]]$chrom)
  write.table(list.Genes.Enhancers.by.CRH[[CRH]][,1:3], paste0("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\CRHs_overlap_0_TAD\\CRH_NEU_iPSC_",CRH,"_",chr,".bed"), col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")
}

queryHits(index_overlaps_CRHs_TADs[subjectHits(index_overlaps_CRHs_TADs)==729])%in%queryHits(index_overlaps_CRHs_TADs[queryHits(index_overlaps_CRHs_TADs)%in%CRHs_wich_overlap_1TAD])
#correspondance pas unique, ce qui explique qu'il y a des CRHs qui chevauchents plusieurs TADs

length(unique(subjectHits(index_overlaps_CRHs_TADs[queryHits(index_overlaps_CRHs_TADs)%in%CRHs_wich_overlap_1TAD&subjectHits(index_overlaps_CRHs_TADs)%in%subset_index_TADs_CRHs_wich_overlap_1TAD])))
#=685 
subset_index_TADs_CRHs_wich_overlap_1TAD = unique(subjectHits(index_overlaps_CRHs_TADs[queryHits(index_overlaps_CRHs_TADs)%in%CRHs_wich_overlap_1TAD]))

#CRHs à retirer de l'analyse
unique(queryHits(index_overlaps_CRHs_TADs[subjectHits(index_overlaps_CRHs_TADs)%in%subset_index_TADs_CRHs_wich_overlap_1TAD]))[which(!unique(queryHits(index_overlaps_CRHs_TADs[subjectHits(index_overlaps_CRHs_TADs)%in%subset_index_TADs_CRHs_wich_overlap_1TAD]))%in%queryHits(index_overlaps_CRHs_TADs[queryHits(index_overlaps_CRHs_TADs)%in%CRHs_wich_overlap_1TAD&subjectHits(index_overlaps_CRHs_TADs)%in%subset_index_TADs_CRHs_wich_overlap_1TAD]))]
length(unique(queryHits(index_overlaps_CRHs_TADs[subjectHits(index_overlaps_CRHs_TADs)%in%subset_index_TADs_CRHs_wich_overlap_1TAD])))-length(unique(queryHits(index_overlaps_CRHs_TADs[subjectHits(index_overlaps_CRHs_TADs)%in%subset_index_TADs_CRHs_wich_overlap_1TAD]))[which(!unique(queryHits(index_overlaps_CRHs_TADs[subjectHits(index_overlaps_CRHs_TADs)%in%subset_index_TADs_CRHs_wich_overlap_1TAD]))%in%queryHits(index_overlaps_CRHs_TADs[queryHits(index_overlaps_CRHs_TADs)%in%CRHs_wich_overlap_1TAD&subjectHits(index_overlaps_CRHs_TADs)%in%subset_index_TADs_CRHs_wich_overlap_1TAD]))])

#TADs avec CRHs pour lesquels les CRHs ne chevauchent qu'un seul TAD
#à refaire rouler pour obtenir des CRHs qui ne chevauchent réellement qu' 1 seul TAD
for(i in unique(subjectHits(index_overlaps_CRHs_TADs[queryHits(index_overlaps_CRHs_TADs)%in%CRHs_wich_overlap_1TAD&subjectHits(index_overlaps_CRHs_TADs)%in%subset_index_TADs_CRHs_wich_overlap_1TAD]))){
  write.table(do.call("rbind",list.Genes.Enhancers.by.CRH[queryHits(index_overlaps_CRHs_TADs[subjectHits(index_overlaps_CRHs_TADs)==i])]),
              paste0("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\TADs_for_CRHs_overlap_1_TAD\\", seqnames(TADs_DI_iPSC_NEu[i]),"_",start(TADs_DI_iPSC_NEu[i]),"_", end(TADs_DI_iPSC_NEu[i]), ".bed"), col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")
  
}


for(TAD in subset_index_TADs_CRHs_wich_overlap_1TAD){
  chr = seqnames(TADs_DI_iPSC_NEu[TAD])
  rtracklayer::export(TADs_DI_iPSC_NEu[TAD], paste0("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\Subset_TADs\\",chr,"\\",paste0("TAD_",start(TADs_DI_iPSC_NEu[TAD]),"_", end(TADs_DI_iPSC_NEu[TAD]), ".bed")), format="bed")
}

#Analyse avec liftover

index_overlaps_CRHs_TADs.liftover.hg38 = findOverlaps(CRHs_clusters.liftover.hg38, TADs_DI_iPSC_NEU_hg38.final)
overlaps_CRHs_TADs.liftover.hg38 = countOverlaps(CRHs_clusters.liftover.hg38, TADs_DI_iPSC_NEU_hg38.final)


#Nombre de TADs pour lesquels on a des CRHs qui ne chevauchent qu'1 seul TAD
unique(subjectHits(index_overlaps_CRHs_TADs.liftover.hg38[queryHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%which(overlaps_CRHs_TADs.liftover.hg38%in%c(1))]))
#Nombre de TADs pour lesquels on a des CRHs qui chevauchent qu'2 TADs
unique(subjectHits(index_overlaps_CRHs_TADs.liftover.hg38[queryHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%which(overlaps_CRHs_TADs.liftover.hg38%in%c(2))]))

CRHs_wich_overlap_0TAD.liftover.hg38 = which(overlaps_CRHs_TADs.liftover.hg38==0)
CRHs_wich_overlap_1TAD.liftover.hg38 = which(overlaps_CRHs_TADs.liftover.hg38==1)
CRHs_wich_overlap_2TADs.liftover.hg38 = which(overlaps_CRHs_TADs.liftover.hg38==2)

sum(overlaps_CRHs_TADs.liftover.hg38%in%c(0,1,2))

queryHits(index_overlaps_CRHs_TADs.liftover.hg38)[queryHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%CRHs_wich_overlap_0TAD.liftover.hg38]
table(table(queryHits(index_overlaps_CRHs_TADs.liftover.hg38)[queryHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%CRHs_wich_overlap_1TAD.liftover.hg38]))
table(table(queryHits(index_overlaps_CRHs_TADs.liftover.hg38)[queryHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%CRHs_wich_overlap_2TADs.liftover.hg38]))
#Les résultats sont cohérents avec ce qu'il est attendu

colnames.GE = c("chrom","chromStart", "chromEnd", "name")

list.Genes.Enhancers.by.CRH.liftover.hg38 = lapply(1:length(split_CRHs_NEU.liftover.hg38), function(x){
  Promoters = na.omit(unique(split_CRHs_NEU.liftover.hg38[[x]][,c("chr.x", "startProm", "endProm", "membership")]))
  Enhancers = na.omit(unique(split_CRHs_NEU.liftover.hg38[[x]][,c("chr.x","start", "end", "membership")]))
  
  colnames(Promoters) = colnames(Enhancers) = colnames.GE
  PE = rbind(Promoters,Enhancers)
  PE
})

sum(sapply(1:1632, function(CRH) min(list.Genes.Enhancers.by.CRH.liftover.hg38[[CRH]]$chromStart) == start(CRHs_clusters.liftover.hg38[CRH]) & max(list.Genes.Enhancers.by.CRH.liftover.hg38[[CRH]]$chromEnd) == end(CRHs_clusters.liftover.hg38[CRH]))) #== 1632 le nombre de CRHs dans les données 
#Pour l'ensemble des CRHs il y a une parfaite correspondance entre les objets list.Genes.Enhancers.by.CRH et CRHs_clusters
#Les débuts et fin des CRHs coincident parfaitement 

#CRHs qui chevauchent 2 TADs
#Création des fichiers bed 
for(CRH in CRHs_wich_overlap_2TADs.liftover.hg38){
  chr = unique(list.Genes.Enhancers.by.CRH.liftover.hg38[[CRH]]$chrom)
  write.table(list.Genes.Enhancers.by.CRH.liftover.hg38[[CRH]][,1:3], paste0("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\liftover_hg38\\CRHs_overlap_2_TADs\\CRH_NEU_iPSC_",CRH,"_",chr,".bed"), col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")
}

#CRHs qui ne chevauchent aucun TADs
#Création des fichiers bed
for(CRH in CRHs_wich_overlap_0TAD.liftover.hg38){
  chr = unique(list.Genes.Enhancers.by.CRH.liftover.hg38[[CRH]]$chrom)
  write.table(list.Genes.Enhancers.by.CRH.liftover.hg38[[CRH]][,1:3], paste0("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\liftover_hg38\\CRHs_overlap_0_TAD\\CRH_NEU_iPSC_",CRH,"_",chr,".bed"), col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")
}

queryHits(index_overlaps_CRHs_TADs[subjectHits(index_overlaps_CRHs_TADs)==729])%in%queryHits(index_overlaps_CRHs_TADs[queryHits(index_overlaps_CRHs_TADs)%in%CRHs_wich_overlap_1TAD])
#correspondance pas unique, ce qui explique qu'il y a des CRHs qui chevauchents plusieurs TADs

subset_index_TADs_CRHs_wich_overlap_1TAD.liftover.hg38 = unique(subjectHits(index_overlaps_CRHs_TADs.liftover.hg38[queryHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%CRHs_wich_overlap_1TAD.liftover.hg38]))

length(unique(subjectHits(index_overlaps_CRHs_TADs.liftover.hg38[queryHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%CRHs_wich_overlap_1TAD.liftover.hg38&subjectHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%subset_index_TADs_CRHs_wich_overlap_1TAD.liftover.hg38])))
#=546

#CRHs à retirer de l'analyse
unique(queryHits(index_overlaps_CRHs_TADs.liftover.hg38[subjectHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%subset_index_TADs_CRHs_wich_overlap_1TAD.liftover.hg38]))[which(!unique(queryHits(index_overlaps_CRHs_TADs.liftover.hg38[subjectHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%subset_index_TADs_CRHs_wich_overlap_1TAD.liftover.hg38]))%in%queryHits(index_overlaps_CRHs_TADs.liftover.hg38[queryHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%CRHs_wich_overlap_1TAD.liftover.hg38&subjectHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%subset_index_TADs_CRHs_wich_overlap_1TAD.liftover.hg38]))]
length(unique(queryHits(index_overlaps_CRHs_TADs.liftover.hg38[subjectHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%subset_index_TADs_CRHs_wich_overlap_1TAD.liftover.hg38])))-length(unique(queryHits(index_overlaps_CRHs_TADs.liftover.hg38[subjectHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%subset_index_TADs_CRHs_wich_overlap_1TAD.liftover.hg38]))[which(!unique(queryHits(index_overlaps_CRHs_TADs.liftover.hg38[subjectHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%subset_index_TADs_CRHs_wich_overlap_1TAD.liftover.hg38]))%in%queryHits(index_overlaps_CRHs_TADs.liftover.hg38[queryHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%CRHs_wich_overlap_1TAD.liftover.hg38&subjectHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%subset_index_TADs_CRHs_wich_overlap_1TAD.liftover.hg38]))])

#TADs avec CRHs pour lesquels les CRHs ne chevauchent qu'un seul TAD
#à refaire rouler pour obtenir des CRHs qui ne chevauchent réellement qu' 1 seul TAD
for(i in unique(subjectHits(index_overlaps_CRHs_TADs.liftover.hg38[queryHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%CRHs_wich_overlap_1TAD.liftover.hg38&subjectHits(index_overlaps_CRHs_TADs.liftover.hg38)%in%subset_index_TADs_CRHs_wich_overlap_1TAD.liftover.hg38]))){
  write.table(do.call("rbind",list.Genes.Enhancers.by.CRH.liftover.hg38[queryHits(index_overlaps_CRHs_TADs.liftover.hg38[subjectHits(index_overlaps_CRHs_TADs.liftover.hg38)==i])]),
              paste0("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\liftover_hg38\\TADs_for_CRHs_overlap_1_TAD\\", seqnames(TADs_DI_iPSC_NEU_hg38.final[i]),"_",start(TADs_DI_iPSC_NEU_hg38.final[i]),"_", end(TADs_DI_iPSC_NEU_hg38.final[i]), ".bed"), col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")
  
}


for(TAD in subset_index_TADs_CRHs_wich_overlap_1TAD.liftover.hg38){
  chr = seqnames(TADs_DI_iPSC_NEU_hg38.final[TAD])
  rtracklayer::export(TADs_DI_iPSC_NEU_hg38.final[TAD], paste0("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\liftover_hg38\\Subset_TADs\\",chr,"\\",paste0("TAD_",start(TADs_DI_iPSC_NEU_hg38.final[TAD]),"_", end(TADs_DI_iPSC_NEU_hg38.final[TAD]), ".bed")), format="bed")
}





