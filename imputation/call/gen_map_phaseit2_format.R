library(argparse);library(data.table); library(dplyr); options(encoding = "latin1")
parser <- ArgumentParser(description = "Format the genetic map to be used in ShapeIt2 using the Hapmap reference")
parser$add_argument("-ref_dir", type = "character", default = "/lustre03/project/6033529/IBD_epilepsy/data/genetic_map_grch38", help = "Directory where the Hapmap PLINK files are found [default \"%(default)s\"]")
parser$add_argument("-out", type = "character", help = "Directory where the formatted genetic map will be saved")

#To test
#args <- list("out" = "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/phasing/gen_map",
#	       "ref_dir" = "/lustre03/project/6033529/IBD_epilepsy/data/genetic_map_grch38")

args <- parser$parse_args()
for(chr in 1:22){
    basename <- paste0("plink.chr", chr, ".GRCh38.map")
    mapfile <- file.path(args$ref_dir, basename)
    map_chr <- fread(mapfile) %>%
        mutate(V2=1000000*(diff(c(0,V3)) / diff(c(0,V4)))) %>%
        select(V4, V2, V3) %>%
        setNames(c("pposition", "rrate", "gposition"))
    #Header (pposition rrate gposition); space delimited
    #pposition=The physical position (bp) [integer]
    #rrate=The recombination rate (cM/Mb=1000000*cM/b) [float]
    #gposition=The genetic position (cM) [float]
    mapfile_out <- file.path(args$out, basename)
    fwrite(map_chr, mapfile_out, col.names = TRUE, row.names = FALSE, sep = "\t")
}
