### Creation du fichier pedigree
# modifs de Jasmin juin 2023
# Prépare les fichiers pour rentrer dans Morgan pour tous les chromosomes

# À exécuter à partir de ~/home/jricard/gl_auto_data


for(j in 1:22){

setwd("/home/jricard/gl_auto_data/corrected_pedigree")

dir.create(paste0("chr",j))

setwd(paste0("chr",j))

system(paste0("cp /home/jricard/gl_auto_data/glauto.par ./"))

# 3 premieres colonnes : subjectID fatherID motherID
# 4e colonne : sexe
# 5e colonne : dummy trait (values 0) used only by gl_auto

# Fichier pre obtenu à cet endroit : 
pre <- read.table("/home/jricard/gl_auto_data/GCbroad.pre")
pre[pre[,1] == "245" & pre[,3]==2162, 3]=0 ;pre[pre[,1] == "245" & pre[,4]==2163, 4]=0 
pre[pre[,3] %in% c(0, 3186, 2314, 2381), 3]=0 ; pre[pre[,4] %in% c(0, 3187, 2313, 2380), 4]=0

ped <- read.table(paste0("/home/jricard/gl_auto_data/corrected_pedigree/GSA_genotypes_clumped_mendel_chr",j,".ped"))

#Ne conserver l'info que pour nos sujets
#pre <- pre[pre$V2 %in% ped$V2,]

#Ordonner le pedigree
#pre2 <- pre[order(pre[,c(1)], -pre[,c(2)]),]
#colnames(pre2)[1:6] <- c("fam", "ind", "father", "mother", "sex", "pheno")
#pre2[,7:9] <- NULL
#pre2_genlib <- GENLIB::gen.genealogy(pre2)
#pre <- GENLIB::gen.genout(pre2_genlib, sorted = FALSE)
#pre$fam <- pre2[match(pre$ind, pre2$ind), "fam"]
#pre$pheno <- pre2[match(pre$ind, pre2$ind), "pheno"]
#pre <- pre[,c("fam", "ind", "father", "mother", "sex", "pheno")]

# combiner avec noms de famille avec nom de sujet pour assurer pas de doublons dans les noms de sujet
pre[,2]=paste(pre[,1],pre[,2],sep="_")
pere=pre[,3]; mere=pre[,4]
pre[,3]=paste(pre[,1],pere,sep="_"); pre[,4]=paste(pre[,1],mere,sep="_")
pre[pere==0,3]="0"; pre[mere==0,4]="0"

pedigree <- pre[,2:5]
#creation de la variable dummy trait
pedigree$V6 <- 0

# Il faut ajouter ceci au debut du fichier pour le pedigree

# On change seulement le pedigree size selon notre nb de sujets
cat("input pedigree size ",nrow(pedigree), 
"\ninput pedigree record names 3 integers 2\ninput pedigree record trait 1 integer 2\n*****\n",file="GCbroad_v2.ped")

write.table(pedigree,"GCbroad_v2.ped",quote=F,row.names=F,col.names=F,append=T)



### Extraction des noms des marqueurs

op <- options(scipen=999)

position <- read.table(paste0("/home/jricard/gl_auto_data/corrected_pedigree/GSA_genotypes_clumped_mendel_chr",j,".map"))
names(position) <- c("chr","rs","positionGenetique", "positionPhysique")
source("/mnt-biostats/gen_map_GRCh38/gen_pos_interpolation_grch38.R")
position$positionGenetique <- gen_pos_interpolation_grch38(infos.chr = position$chr, infos.pos = position$positionPhysique, dir = "/mnt-biostats/gen_map_GRCh38", method = "linear", ncores = 1)
which.keep=which(!duplicated(position$positionGenetique))
which.geno.keep=c(2*which.keep-1,2*which.keep)
which.geno.keep=which.geno.keep[order(which.geno.keep)]
position=position[which.keep,]

# On fixe les fréquences d'allèle à 0.5
nb_freq = rep("0.5 0.5",nrow(position))

# On cree le texte qui doit etre ecrit pour les frequences
text2=paste("set markers",1:length(nb_freq),"allele freq",nb_freq)

write.table(text2,"framework_freq.txt",quote=F,row.names=F,col.names=F)

## On ajoute 0.0001 a certain SNPs pour ne pas avoir des positions 'egales'
# position$pos_gen_fin <- position$positionGenetique

# pos.tmp=position$positionGenetique
# while(sum(duplicated(pos.tmp))>0) pos.tmp[duplicated(pos.tmp)]=pos.tmp[duplicated(pos.tmp)]+0.0001
# position$pos_gen_fin=pos.tmp

write.table(position$positionGenetique,"framework_gcbroad.map_af_geno",quote=F,row.names=F,col.names=F)



# # Pour avoir le fichier des positions pour les framwork markers (utilisable avec GIGI)

write.table(position$positionGenetique,"framework_pos.txt",quote=F,row.names=F,col.names=F)



# # Pour le fichier .MAP_AF_GENO utilise dans gl_auto, il suffit de combiner les 3 fichiers ci-haut :

# # 1. Positions
# # 2. Frequence des alleles
# # 3. Genotypes

### conversion des alleles: l'allele avec la premiere lettre dans l'ordre alphabetique obtient le 1 et l'autre allele obtient le 2 (ex : A=1 et G=2)
geno=ped[,-(1:6)]
geno=geno[,which.geno.keep]
geno=apply(geno,2,as.character)
geno[geno=="0"]="Z"
tmp=matrix(as.vector(as.matrix(geno)),ncol=ncol(geno)/2)
geno.new=matrix(as.vector(as.matrix(apply(tmp,2,function(x) as.numeric(as.factor(x))))),ncol=ncol(geno))
geno.new[geno.new==3]=0

text2=paste("set markers",1:length(nb_freq),"allele freq",nb_freq)

write.table(as.matrix(c("map markers position\n",position$positionGenetique,text2,paste0("\nset markers data ",length(position$positionGenetique)))),file="framework_gcbroad.map_af_geno",quote=F,col.names=F,row.names=F)
write.table(cbind(paste0(ped[,1], "_", ped[,2]),geno.new),file="framework_gcbroad.map_af_geno",quote=F,col.names=F,row.names=F,append=T)

}


# Voir fichiers dans ~/3dossier-calculs/projets/Imputation_familiale/MORGAN/ pour plus de details


# Le fichier de parametres pour lancer l'execution de gl_auto se trouve :
# ~/3dossier-calculs/projets/Imputation_familiale/MORGAN/repertoireSPAP_chr*/glauto_*.pa
# On obtient le fichier framework_*.pre$V1s (dans le meme dossier) qui sera utilise pour GIGI
# Le fichier de sortie se trouve dans :
# ~/3dossier-calculs/projets/Imputation_familiale/MORGAN/repertoireSPAP_chr*/gl_auto_*.out




