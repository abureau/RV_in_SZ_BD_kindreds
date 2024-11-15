library(argparse); library(data.table); options(encoding = "latin1")
parser <- ArgumentParser(description = "Modify a .fam file using a reference pedifree")
parser$add_argument("-ref", type = "character", help = "reference pedigree")
parser$add_argument("-fam", type = "character", help = ".fam file to modify")

#To test
#args <- list("fam" = "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/QC/seq_multiallelic.fam",
#	     "ref" = "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/manip/seq_pedigree.txt")

args <- parser$parse_args()
fam <- setNames(fread(args$fam, header = FALSE, colClasses = c(rep("character", 6)) ), paste0("fam", c("V1", "V2", "V3", "V4", "V5", "V6")))
ref <- setNames(fread(args$ref, header = FALSE, colClasses = c(rep("character", 6)) ), paste0("ref", c("V1", "V2", "V3", "V4", "V5", "V6")))
match <- merge(x = fam, y = ref, by.y = "refV2", by.x = "famV2", all.x = TRUE, sort = FALSE)
match <- match[,c("famV1", "famV2", "refV3", "refV4", "refV5", "refV6")]
match[is.na(match)] <- 0
fwrite(fam, paste0(args$fam, "~"), col.names = FALSE, row.names = FALSE, sep = "\t")
fwrite(match, args$fam, col.names = FALSE, row.names = FALSE, sep = "\t")