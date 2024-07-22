load("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\grandtab.RData")


ABC_Pred = process.ABC("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\EnhancerPredictions_NEU_iPSC.txt")
TADs_DI_iPSC_NEU = rtracklayer::import("C:\\Users\\loicm\\Documents\\CRHs_SZ\\data\\TADs_DI_iPSC_NEU.bed", format = "bed")

#Liftover des Enhancers et Gènes de hg19 vers hg38
ABC_Pred.for.liftover = ABC_Pred

liftover_hg38_TADs = read.table("C:\\Users\\loicm\\Downloads\\hglft_genome_18652_b182c0.bed", header=FALSE)
colnames(liftover_hg38_TADs) = c("chrom_on_hg38","start_hg38","end_hg38","original_TAD_hg19","n_correspondance")

liftover_hg38_TADs$chrom_on_hg19 = sapply(strsplit(liftover_hg38_TADs$original_TAD_hg19, ":", fixed=TRUE), function(x) x[1])

sum(liftover_hg38_TADs$chrom_on_hg38 != liftover_hg38_TADs$chrom_on_hg19)
#Le nombre de TADs pour lesquels le liftover a donné des régions situées sur d'autres chromosomes
#Je retire donc ces régions 

liftover_hg38_TADs = liftover_hg38_TADs[liftover_hg38_TADs$chrom_on_hg38 == liftover_hg38_TADs$chrom_on_hg19,]
GRanges_TADs_liftover_hg38 = GRanges(seqnames = liftover_hg38_TADs$chrom_on_hg38 ,ranges=IRanges(start=liftover_hg38_TADs$start_hg38,end=liftover_hg38_TADs$end_hg38))

boxplot(width(TADs_DI_iPSC_NEU), width(GRanges_TADs_liftover_hg38)) #Les résultats sont consistants avant et après liftover

liftover_hg38_Enhancers = read.table("C:\\Users\\loicm\\Downloads\\hglft_genome_339cf_b1bdb0.bed", header=FALSE)
colnames(liftover_hg38_Enhancers) = c("chrom_on_hg38","start_hg38","end_hg38","original_Enhancer_hg19","n_correspondance", "strand")

liftover_hg38_Enhancers$chrom_on_hg19 = sapply(strsplit(stringr::str_extract(liftover_hg38_Enhancers$original_Enhancer_hg19, "(?=>|)(.+)(?=\\:)"), "|", fixed=TRUE), function(x) x[2])

liftover_hg38_Enhancers = liftover_hg38_Enhancers[liftover_hg38_Enhancers$chrom_on_hg19==liftover_hg38_Enhancers$chrom_on_hg38,]
GRanges_Enhancers_liftover_hg38 = GRanges(seqnames = liftover_hg38_Enhancers$chrom_on_hg38 ,ranges=IRanges(start=liftover_hg38_Enhancers$start_hg38,end=liftover_hg38_Enhancers$end_hg38), name=liftover_hg38_Enhancers$original_Enhancer_hg19)

liftover_hg38_Genes = read.table("C:\\Users\\loicm\\Downloads\\hglft_genome_f83_b21e10.bed", header=FALSE)
colnames(liftover_hg38_Genes) = c("chrom_on_hg38","start_hg38","end_hg38","original_Genes_hg19","n_correspondance", "strand")

liftover_hg38_Genes$chrom_on_hg19 = sapply(strsplit(liftover_hg38_Genes$original_Genes_hg19, "_", fixed=TRUE), function(x) x[1])
liftover_hg38_Genes$TargetGene = sapply(strsplit(liftover_hg38_Genes$original_Genes_hg19, "_", fixed=TRUE), function(x) x[2])

liftover_hg38_Genes = liftover_hg38_Genes[liftover_hg38_Genes$chrom_on_hg19==liftover_hg38_Genes$chrom_on_hg38,]
GRanges_Genes_liftover_hg38 = GRanges(seqnames = liftover_hg38_Genes$chrom_on_hg38 ,ranges=IRanges(start=liftover_hg38_Genes$start_hg38,end=liftover_hg38_Genes$end_hg38), TargetGene=liftover_hg38_Genes$TargetGene)

df_Enhancers.hg38.final = as.data.frame(GRanges_Enhancers_liftover_hg38)
colnames(df_Enhancers.hg38.final) = c("chr.x","start.hg38","end.hg38","widthEnh.hg38","strand","name")
df_Genes.hg38.final = as.data.frame(GRanges_Genes_liftover_hg38)
colnames(df_Genes.hg38.final) = c("chr.x","startProm.hg38","endProm.hg38","widthProm.hg38","strand","TargetGene")

ABC_Pred.for.liftover = merge(merge(ABC_Pred.for.liftover, df_Enhancers.hg38.final[,c("chr.x","start.hg38","end.hg38", "name")],by=c("chr.x", "name"), all.x=TRUE), df_Genes.hg38.final[,c("chr.x","startProm.hg38","endProm.hg38","TargetGene")], by=c("chr.x","TargetGene"),all.x=TRUE)
ABC_Pred.for.liftover.tmp = ABC_Pred.for.liftover
ABC_Pred.for.liftover = ABC_Pred.for.liftover[,c("chr.x","start.hg38","end.hg38", "name","startProm.hg38","endProm.hg38","TargetGene","TargetGeneTSS", "CellType","ABC.Score",     "typeOf", "chr.y")]
colnames(ABC_Pred.for.liftover)=c("chr.x","start", "end", "name", "startProm","endProm", "TargetGene","TargetGeneTSS", "CellType","ABC.Score",     "typeOf", "chr.y")

ABC_Pred.for.liftover$name = paste0(ABC_Pred.for.liftover$start,"-" ,ABC_Pred.for.liftover$end)
CRHs_iPSC_NEU.liftover.hg38=create.CRHs(ABC_Pred.for.liftover, "ABC")

ABC_Pred_with_membership.liftover.hg38 = add.membership(CRHs_iPSC_NEU.liftover.hg38,ABC_Pred.for.liftover)

write.table(ABC_Pred_with_membership.liftover.hg38, "C:\\Users\\loicm\\Documents\\CRHs_SZ\\EnhancerPredictions_NEU_iPSC_liftoverhg38.txt", col.names = TRUE, row.names = TRUE, sep="\t", quote=F)


po = coverage.By.CRH(CRHs_iPSC_NEU.liftover.hg38, ABC_Pred_with_membership.liftover.hg38, method="ABC")

unique(ABC_Pred_with_membership.liftover.hg38[ABC_Pred_with_membership.liftover.hg38$membership==221,"TargetGene"])
unique(ABC_Pred_with_membership.liftover.hg38[ABC_Pred_with_membership.liftover.hg38$membership==221,c("start", "end")])

unique(ABC_Pred_with_membership.liftover.hg38[ABC_Pred_with_membership.liftover.hg38$membership==277,"TargetGene"])
unique(ABC_Pred_with_membership.liftover.hg38[ABC_Pred_with_membership.liftover.hg38$membership==277,"name"])

unique(ABC_Pred_with_membership.liftover.hg38[ABC_Pred_with_membership.liftover.hg38$membership==33,"TargetGene"])
unique(ABC_Pred_with_membership.liftover.hg38[ABC_Pred_with_membership.liftover.hg38$membership==33,"name"])

unique(ABC_Pred_with_membership.liftover.hg38[ABC_Pred_with_membership.liftover.hg38$membership==1280,"TargetGene"])
unique(ABC_Pred_with_membership.liftover.hg38[ABC_Pred_with_membership.liftover.hg38$membership==1280,"name"])

unique(ABC_Pred_with_membership.liftover.hg38[ABC_Pred_with_membership.liftover.hg38$membership==1094,"TargetGene"])
unique(ABC_Pred_with_membership.liftover.hg38[ABC_Pred_with_membership.liftover.hg38$membership==1094,"name"])

unique(ABC_Pred_with_membership.liftover.hg38[ABC_Pred_with_membership.liftover.hg38$membership==1438,"TargetGene"])
unique(ABC_Pred_with_membership.liftover.hg38[ABC_Pred_with_membership.liftover.hg38$membership==1438,"name"])

unique(ABC_Pred_with_membership.liftover.hg38[ABC_Pred_with_membership.liftover.hg38$membership==694,"TargetGene"])
unique(ABC_Pred_with_membership.liftover.hg38[ABC_Pred_with_membership.liftover.hg38$membership==694,"name"])

#Plot CRHs with significant signals
V(CRHs_iPSC_NEU.liftover.hg38)$Elements = ifelse(names(V(CRHs_iPSC_NEU.liftover.hg38))%in% ABC_Pred_with_membership.liftover.hg38$TargetGene, "Promoters", "Distal")


#Scpefying the color depending on element type
V(CRHs_iPSC_NEU.liftover.hg38)$color = ifelse(names(V(CRHs_iPSC_NEU.liftover.hg38))%in% ABC_Pred_with_membership.liftover.hg38$TargetGene, "#00BFC4","#F8766D")
V(CRHs_iPSC_NEU.liftover.hg38)$Name = names(V(CRHs_iPSC_NEU.liftover.hg38))

decompose_CRHs_iPSC_NEU.liftover.hg38 = decompose(CRHs_iPSC_NEU.liftover.hg38)

ggnet2(decompose_CRHs_iPSC_NEU.liftover.hg38[[1438]], size=4,color="Elements", shape="Elements", palette = "Set1",label="Name",label.size = 3,legend.size = 20,
       legend.position = "none")
