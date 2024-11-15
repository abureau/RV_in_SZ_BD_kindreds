#Phasing using ShapeIt5 2024-04-04

#The softwares are here
shapeit5_common=/lustre03/project/6033529/SOFT/shapeit5/phase_common/bin/phase_common
duohmm=/lustre03/project/6033529/SOFT/duohmm_v0.1.7/duohmm

#The hg38 genetic map references from the HapMap are downloaded here:
gen_map=/lustre03/project/6033529/IBD_epilepsy/data/genetic_map_grch38

#GCbroad phenotype file, will be used to complete the pedigree.
pheno=/lustre03/project/6033529/quebec_10x/data/phenotypes/GCbroad.pre
cut -f-6 ${pheno} > ${path_data}/seq_pedigree.fam

#Genotyping data
path_data=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/genotyping
mkdir ${path_data}/gen_map

#To compare with sequencing data
path_phasing=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/phasing
mask=with_mask
mkdir ${path_data}/rm_flip

#Verify that the pedigree is well formated. We now need the full pedigree: photo2021. Set mendelien errors to missing.
seq_ped=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples/Photo2021.pre
cut -f-6 ${seq_ped} > ${path_data}/Photo2021.pre

for puce in omni gsa 
do
  Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples/call/update_pheno.R -fam ${path_data}/${puce}_grch38.fam -ref ${path_data}/Photo2021.pre
  #To correct mendel error, use plink1, might need to play with the loaded modules...
  plink --bfile ${path_data}/${puce}_grch38 --allow-extra-chr --set-me-missing --mendel-duos --mendel-multigen --make-bed --out ${path_data}/${puce}_grch38
  mv ${path_data}/${puce}_grch38.log ${path_data}/${puce}_grch38_mendel.log

  for chr in {22..1}
  do
    #Break the dataset by chromosome.
    plink --bfile ${path_data}/${puce}_grch38 --chr ${chr} --allow-extra-chr --make-bed --out ${path_data}/${puce}_grch38_chr${chr}

    #Add genetic positions. Be sure that every markers are in the genetic map, if not, trimming is needed to avoid errors.
    Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/gen_pos_interpolation_grch38.R -bim ${path_data}/${puce}_grch38_chr${chr}.bim -ref_dir ${gen_map} -method "linear"

    #Rearrange genetic maps for ShapeIt5
    cat <(echo -e "pos\tchr\tcM") <(awk 'BEGIN{OFS="\t"} {print $4, $1, $3}' ${gen_map}/plink.chr${chr}.GRCh38.map) > ${path_data}/gen_map/plink.chr${chr}.GRCh38.map

    #We wish to phase on a dataset where the variants that are not included in the genetic map are removed.
    cat ${gen_map}/plink.chr${chr}.GRCh38.map | cut -d " " -f 1,4 | sed -n '1p;$p' > ${path_data}/gen_map/gen_pos_ranges_chr${chr}_tmp.txt
    awk '{a[$1] = a[$1] " " $2} END {for (i in a) print i a[i]}' ${path_data}/gen_map/gen_pos_ranges_chr${chr}_tmp.txt | awk 'BEGIN {OFS="\t"} {print $0, $1}' > ${path_data}/gen_map/gen_pos_ranges_chr${chr}.txt; rm ${path_data}/gen_map/gen_pos_ranges_chr${chr}_tmp.txt
    plink --bfile ${path_data}/${puce}_grch38_chr${chr} --extract range ${path_data}/gen_map/gen_pos_ranges_chr${chr}.txt --make-bed -out ${path_data}/${puce}_grch38_in_genmap_chr${chr}
    rm ${path_data}/${puce}_grch38_chr${chr}.*

    #We need to keep variants that are present in the sequences too. Also, we need to flip stand and switch allele order if needed.
    awk 'BEGIN {FS=OFS="\t"} {print $1, $4, $5, $6, $2}' ${path_data}/${puce}_grch38_in_genmap_chr${chr}.bim > ${path_data}/rm_flip/${puce}_grch38_in_genmap_chr${chr}_comp.bim
    bcftools query -f '%CHROM\t%POS\t%ALT\t%REF\n'  ${path_phasing}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf > ${path_data}/rm_flip/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}_comp.bim
    Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples/call/2023-12-05_bim_matching.R -bim1 ${path_data}/rm_flip/${puce}_grch38_in_genmap_chr${chr}_comp -bim2 ${path_data}/rm_flip/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}_comp

    #Generate the .vcf versions without the variants not in the sequences
    plink --bfile ${path_data}/${puce}_grch38_in_genmap_chr${chr} --extract ${path_data}/rm_flip/${puce}_grch38_in_genmap_chr${chr}_comp_SNP_to_keep.txt --reference-allele ${path_data}/rm_flip/${puce}_grch38_in_genmap_chr${chr}_comp_SNP_to_switch.txt --flip ${path_data}/rm_flip/${puce}_grch38_in_genmap_chr${chr}_comp_SNP_to_flip.txt --recode vcf -out ${path_data}/${puce}_grch38_in_genmap_chr${chr}

    #Create fam and map files in the format of shapeit5. In the .fam file, the fam needs to be appended to the IDs
    awk '{print $1"_"$2, $1"_"$3, $1"_"$4}' ${path_data}/${puce}_grch38_in_genmap_chr${chr}.fam | awk -v OFS='\t'  -F' ' '{for(i=1;i<=NF;i++) if($i~/_0$/) $i=0}1' > ${path_data}/${puce}_grch38_in_genmap_chr${chr}_shapeit5.fam

    #Phase common variants
    bcftools +fill-tags ${path_data}/${puce}_grch38_in_genmap_chr${chr}.vcf | bgzip -c > ${path_data}/${puce}_grch38_in_genmap_chr${chr}.vcf.gz
    bcftools index -t -f ${path_data}/${puce}_grch38_in_genmap_chr${chr}.vcf.gz
    ${shapeit5_common} --input ${path_data}/${puce}_grch38_in_genmap_chr${chr}.vcf.gz \
                       --pedigree ${path_data}/${puce}_grch38_in_genmap_chr${chr}_shapeit5.fam \
                       --region ${chr} \
                       --map ${path_data}/gen_map/plink.chr${chr}.GRCh38.map \
                       --output ${path_data}/${puce}_PHASED_grch38_in_genmap_chr${chr}.bcf \
                       --thread 8

    #Convert to haps format
    plink2 --bcf ${path_data}/${puce}_PHASED_grch38_in_genmap_chr${chr}.bcf --export haps --out ${path_data}/${puce}_PHASED_grch38_in_genmap_chr${chr}

    #the .sample file needs to contain the parents info
    mv ${path_data}/${puce}_PHASED_grch38_in_genmap_chr${chr}.sample ${path_data}/${puce}_PHASED_grch38_in_genmap_chr${chr}.sample~
    cat <(printf "ID_1 ID_2 missing father mother sex\n0 0 0 D D D\n") <(awk -F'\t' '{print $0 "\t0"}' ${path_data}/${puce}_grch38_in_genmap_chr${chr}.fam | awk -F'\t' '{print $1 "\t" $2 "\t" $7 "\t" $3 "\t" $4 "\t" $5}' | tr "\t" " ") > ${path_data}/${puce}_PHASED_grch38_in_genmap_chr${chr}.sample

    #Infering variants with likely genotyping errors based on the duoHMM algorith
    ${duohmm} -H ${path_data}/${puce}_PHASED_grch38_in_genmap_chr${chr} \
    	       -M ${path_data}/gen_map/plink.chr${chr}.GRCh38.map \
      	       -G ${path_data}/${puce}_PHASED_grch38_in_genmap_chr${chr}.generr

    #Remove files what wont be used 
    rm ${path_data}/${puce}_grch38_in_genmap_chr${chr}.fam ${path_data}/${puce}_grch38_in_genmap_chr${chr}.bed ${path_data}/${puce}_grch38_in_genmap_chr${chr}.bim ${path_data}/${puce}_grch38_in_genmap_chr${chr}.vcf* ${path_data}/${puce}_grch38_in_genmap_chr${chr}.log
    rm ${path_data}/${puce}_grch38_in_genmap_chr${chr}_shapeit5.fam
    rm ${path_data}/${puce}_PHASED_grch38_in_genmap_chr${chr}.haps ${path_data}/${puce}_PHASED_grch38_in_genmap_chr${chr}.sample ${path_data}/${puce}_PHASED_grch38_in_genmap_chr${chr}.sample~ ${path_data}/${puce}_PHASED_grch38_in_genmap_chr${chr}.log
  done
done