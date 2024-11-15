library(argparse); options(encoding = "latin1")
parser <- ArgumentParser(description = "Interpolate genetic position using a GRCh38 reference Hapmap reference")
parser$add_argument("-bim", type = "character", help = ".bim file to generate genetic positions in")
parser$add_argument("-ref_dir", type = "character", default = "/lustre03/project/6033529/IBD_epilepsy/data/genetic_map_grch38", help = "Directory where the Hapmap PLINK files are found [default \"%(default)s\"]")
parser$add_argument("-cores", type = "integer", default = 1, help = "Number of cores to be uses for computation [default \"%(default)s\"]")
parser$add_argument("-method", type = "character", default = "linear", help = "Method to be used for interpolation (spline or linear) [default \"%(default)s\"]")

#To test
#args <- list("bim" = "/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/QC/seq_FINAL_with_mask.bim",
#	       "ref_dir" = "/lustre03/project/6033529/IBD_epilepsy/data/genetic_map_grch38", "cores" = 1)

library(bigsnpr); library(data.table)

args <- parser$parse_args()
if(!args$method %in% c("spline", "linear")){error("Specified method is not implemented")}

gen_pos_interpolation_grch38 <- function (infos.chr, infos.pos, dir, method, ncores = 1)
{
    #Inspired by the bigsnpr package by Florien Prive

    snp_split(infos.chr, function(ind.chr, pos, dir, method) {
        chr <- attr(ind.chr, "chr")
        basename <- paste0("plink.chr", chr, ".GRCh38.map")
        mapfile <- file.path(dir, basename)
        map.chr <- bigreadr::fread2(mapfile, showProgress = FALSE, nThread = 1)
        pos.chr <- pos[ind.chr]
        ind <- match(pos.chr, map.chr$V4) #matching the physical pos found in the reference
        new_pos <- map.chr$V3[ind] #genetic pos is extracted for physical pos that are found in the reference
        indNA <- which(is.na(ind)) #which variants need to be interpolated?
        if (length(indNA) > 0) {
            if(method == "spline"){
                new_pos[indNA] <- suppressWarnings(stats::spline(map.chr$V4,
                    map.chr$V3, xout = pos.chr[indNA], method = "hyman")$y)
            } else {
                new_pos[indNA] <- suppressWarnings(stats::approx(map.chr$V4,
                    map.chr$V3, xout = pos.chr[indNA], yleft = 0, yright = max(map.chr$V3))$y)
            }
        }
        new_pos
    }, combine = "c", pos = infos.pos, dir = dir, method = method, ncores = ncores)
}

bim <- fread(args$bim, header = FALSE)
fwrite(bim, paste0(args$bim, "~"), col.names = FALSE, row.names = FALSE, sep = "\t")
gen_pos <- gen_pos_interpolation_grch38(infos.chr = bim$V1, infos.pos = bim$V4, dir = args$ref_dir, method = args$method, ncores = args$cores)
bim$V3 <- gen_pos
fwrite(bim, args$bim, col.names = FALSE, row.names = FALSE, sep = "\t")