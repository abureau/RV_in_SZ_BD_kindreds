library(argparse); library(data.table); library(dplyr); options(encoding = "latin1")
#To test
#args <- list("")

parser <- ArgumentParser(description = "Transpose a dosage file column by column")
parser$add_argument("-dosage", type = "character", help = "Path to the dosages file to tranpose")
parser$add_argument("-out", type = "character", help = "Output file name")
args <- parser$parse_args()

dosage <- fread(args$dosage, header = FALSE)
n_subject <- nrow(dosage)
n_var <- ncol(dosage)

for(var in 1:n_var){
  dosage_i <- dosage[[var]]
  cat(paste(dosage_i, collapse = " "), file = args$out, append = TRUE, sep = "\n")
}