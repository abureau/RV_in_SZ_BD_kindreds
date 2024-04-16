#To generate Figure S16.
library(data.table); library(dplyr); library(kinship2); library(glue)

path_data <- "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/phasing/seq_FINAL_with_mask_gen_map_var_chr7"
pos <- "chr7:44131068:G:T"
fam <- "238"
pedigree_path <- "/lustre03/project/6033529/quebec_10x/data/freeze/GCbroad_seq_inbred.pre"
df.ped.2021 = read.table("/lustre03/project/6033529/quebec_10x/data/freeze/GCbroad_seq_inbred.pre", header=FALSE, sep=" ")
colnames(df.ped.2021) = c("famid","id","fid","mid","sex","affected")
subset.fam = c("103","105","110","115","119","121","124","125","126","129","131","133","151","182","207","210","211","212","217","220","224","228","230","233","234","235","238","255")
path_out <- "/home/jasmric"
width <- 1400
height <-  900

#Create the data
system(glue("grep {fam} {path_data}.fam | cut -d ' ' -f-2 > {path_out}/fam.keep"))
system(glue("echo {pos} > {path_out}/var.extract"))
system(glue("plink --file {path_data} --extract {path_out}/var.extract --keep {path_out}/fam.keep --recode --out  {path_out}/data"))

#Read the data and remove a possible "chr" header on the chromosome number.
sample_fam <- fread(glue("{path_out}/data.ped"), header = FALSE)

#Prepare whole pedigree for the plot.
pedigree <- df.ped.2021[df.ped.2021$famid %in% fam,]
sample_fam_n <- nrow(sample_fam); pedigree_n <- nrow(pedigree)
sample_fam <- merge(x = sample_fam, y = pedigree, by.x = "V2", by.y = "id", all.y = TRUE, sort = FALSE)
str_haplo <- apply(sample_fam[!is.na(sample_fam$V7) & !is.na(sample_fam$V8), c("V7", "V8")], 1, function(x){paste0(x[1], "|", x[2])})
sample_fam <- sample_fam[, c("famid", "V2", "fid", "mid", "sex", "affected")] %>% setNames(., c("famid", "id", "dadid", "momid", "sex", "affected"))
sample_fam$affected <- sample_fam$affected-1; sample_fam$affected[sample_fam$affected==-1] <- NA

#Plot.
ped <- pedigree(id=sample_fam$id, dadid=sample_fam$dadid, momid=sample_fam$momid, sex=sample_fam$sex, famid=sample_fam$famid, affected = sample_fam$affected)
str_haplo <- c(str_haplo, rep("", pedigree_n-sample_fam_n))
str_haplo <- paste(sample_fam$id, str_haplo, sep="\n")

png(glue("{path_out}/seq_var_{pos}_fam_{fam}.png"), width=width, height=height)
plot(ped[1], cex = 1, id = str_haplo,
     symbolsize = 1, branch = 1, packed = TRUE)
title(main = paste0("Family ", fam), sub = pos, cex.main = 1.5)
dev.off()
