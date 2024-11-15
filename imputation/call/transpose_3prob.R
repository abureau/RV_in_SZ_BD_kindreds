library(argparse); library(data.table); library(dplyr); options(encoding = "latin1")
#To test
#args <- list("")

parser <- ArgumentParser(description = "Transpose a probabilities file column by column")
parser$add_argument("-prob", type = "character", help = "Path to the probabilities file to tranpose")
parser$add_argument("-prob_type", type = "character", help = "Transposition type based on the format of the prob file. 
                                                              var_subject = Rows are variants and columns are subject, we wish to swithc that. or 
                                                              subject_var = Rows are subjects and columns are variants, we wish to swithc that")
parser$add_argument("-out", type = "character", help = "Output file name")
args <- parser$parse_args()

if(args$prob_type == "var_subject"){
  prob <- fread(args$prob, header = FALSE)
  n_subject <- ncol(prob)/3
  n_var <- nrow(prob)

  for(subject in 1:n_subject){
    i <- c(subject*3-2, subject*3-1, subject*3)
    prob_i <- select(prob, all_of(i))
    out_i_space <- as.vector(t(prob_i))
    cat(paste(out_i_space, collapse = " "), file = args$out, append = TRUE, sep = "\n")
  }
}

if(args$prob_type == "subject_var"){
  prob <- fread(args$prob, header = FALSE)
  n_subject <- nrow(prob)
  n_var <- ncol(prob)/3
  out_matrix <- matrix(nrow = n_var, ncol = n_subject*3)
  for(var in 1:n_var){
    i <- c(var*3-2, var*3-1, var*3)
    prob_i <- select(prob, all_of(i))
    out_matrix[var, ] <- as.vector(t(prob_i))
  }
  fwrite(data.table(out_matrix), file = args$out, sep = " ", col.names = FALSE, row.names = FALSE)
}