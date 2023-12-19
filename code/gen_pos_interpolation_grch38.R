library(bigsnpr)
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
