#This code is used to split the 10th, 18th and 47th families when running the code to get the NULL objects for RetroFun-RVS

if(fam_idx == 10){
  new_founders <- list(c(2395,2562), c(2560,2561)); last_founders <- c(2646,  2647)
}else if(fam_idx == 18){
  new_founders <- list(c(2753, 2754), c(2795, 2794, 2889, 2890)); last_founders <- c(3066,  3067)
}else if(fam_idx == 47){
  new_founders <- list(c(843, 844, 845, 851), c(841, 842, 829, 827)); last_founders <- c(839, 840)
}

to_split <- data.frame(pedigree)
split_ped <- function(pedigree, parents){
  kids <- 0; family <- c(parents)
  while(length(kids)!=0){
    kids <- pedigree$id[pedigree$dadid %in% parents | pedigree$momid %in% parents]
    family <- unique(c(family, kids))
    parents <- kids 
  }
  out <- unique(unlist(pedigree[pedigree$id %in% family, c("id", "dadid", "momid"),]))
  return(out)
}

sep1 <- split_ped(pedigree = to_split, parents = new_founders[[1]])
if(length(new_founders[[1]])<=2){
  sep1 <- sep1[!sep1 %in% last_founders]
  sep1 <- to_split[to_split$id %in% sep1,]
  sep1$dadid[sep1$dadid %in% last_founders] <- 0; sep1$momid[sep1$momid %in% last_founders] <- 0
} else {
  sep1 <- to_split[to_split$id %in% sep1,]
}
#sep1 <- pedigree(id=sep1$id, dadid=sep1$dadid, momid=sep1$momid, sex=sep1$sex, affected=sep1$affected)

sep2 <- split_ped(pedigree = to_split, parents = new_founders[[2]])
if(length(new_founders[[2]])<=2){
  sep2 <- sep2[!sep2 %in% last_founders]
  sep2 <- to_split[to_split$id %in% sep2,]
  sep2$dadid[sep2$dadid %in% last_founders] <- 0; sep2$momid[sep2$momid %in% last_founders] <- 0
} else {
  sep2 <- to_split[to_split$id %in% sep2,]
}
#sep2 <- pedigree(id=sep2$id, dadid=sep2$dadid, momid=sep2$momid, sex=sep2$sex, affected=sep2$affected)
#split_pedigree <- list(sep1, sep2)

split_pedigree <- pedigree(id = c(sep1$id, sep2$id),
         dadid = c(sep1$dadid, sep2$dadid),
         momid = c(sep1$momid, sep2$momid), 
         sex = c(sep1$sex, sep2$sex), 
         affected = c(sep1$affected, sep2$affected),
         famid = c(rep(paste0(fam, "-1"), nrow(sep1)), rep(paste0(fam, "-2"), nrow(sep2))))