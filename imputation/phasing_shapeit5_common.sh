#Phasing using ShapeIt5 2023-13-22

#The software is here
shapeit5=/lustre03/project/6033529/SOFT/shapeit5/phase_common/bin/phase_common
duohmm=/lustre03/project/6033529/SOFT/duohmm_v0.1.7/duohmm

#We will use the QCed data with mask.
mask=with_mask

#Results will be exported here:
path_data=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/phasing

#The hg38 genetic map references from the HapMap are downloaded here:
gen_map=/lustre03/project/6033529/IBD_epilepsy/data/genetic_map_grch38
path_gen_map=${path_data}/gen_map

#Compute the recombination rate and generate the genetic map in the correct format
Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/gen_map_phaseit2_format.R -ref_dir ${gen_map} -out ${path_gen_map}

#GCbroad phenotype file, will be used to complete the pedigree.
pheno=/lustre03/project/6033529/quebec_10x/data/phenotypes/GCbroad.pre
cut -f-6 ${pheno} > ${path_data}/seq_pedigree.fam

for chr in {1.22}
do
  #Our QCed sequencing data are found here (we used the data with mask):
  seq_data=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/QC/seq_FINAL_chr${chr}_with_mask.vcf.gz

  #Generate bfiles.
  plink2 --vcf ${seq_data} --make-bed --out ${path_data}/seq_FINAL_chr${chr}_with_mask

  #Add genetic positions. Be sure that every markers are in the genetic map, if not, trimming is needed to avoid errors.
  Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/gen_pos_interpolation_grch38.R -bim ${path_data}/seq_FINAL_chr${chr}_with_mask.bim -ref_dir ${gen_map} -method "linear"

  #We wish to phase on a dataset where the variants that are not included in the genetic map are removed.
  cat ${gen_map}/plink.chr${chr}.GRCh38.map | cut -d " " -f 1,4 | sed -n '1p;$p' > ${path_data}/gen_pos_ranges_chr${chr}_tmp.txt
  awk '{a[$1] = a[$1] " " $2} END {for (i in a) print i a[i]}' ${path_data}/gen_pos_ranges_chr${chr}_tmp.txt | awk 'BEGIN {OFS="\t"} {print $0, $1}' > ${path_data}/gen_pos_ranges_chr${chr}.txt; rm ${path_data}/gen_pos_ranges_chr${chr}_tmp.txt
  plink2 --bfile ${path_data}/seq_FINAL_chr${chr}_with_mask --extract range ${path_data}/gen_pos_ranges_chr${chr}.txt --make-bed -out ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}
  rm ${path_data}/seq_FINAL_chr${chr}_with_mask.*

  #Generate the .vcf versions without the variants
  plink2 --bfile ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr} --recode vcf -out ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}

  #Create fam and map files in the format of shapeit5. In the .fam file, the fam needs to be appended to the IDs
  cut -f2- ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam > ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam~
  paste <(cut -d "_" -f1 ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam~) <(cut -d "_" -f2- ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam~) > ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam
  Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples/call/update_pheno.R -fam ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam -ref ${path_data}/seq_pedigree.fam
  awk -v OFS='\t' '{print $1"_"$2, $1"_"$3, $1"_"$4}' ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam | awk -v OFS='\t'  -F' ' '{for(i=1;i<=NF;i++) if($i~/_0$/) $i=0}1' > ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}_shapeit5.fam
  awk -v OFS='\t' '{print $4, $1, $3}' ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.bim > ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}_shapeit5.bim 

  #Phase common variants
  bcftools +fill-tags ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.vcf | bgzip -c > ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.vcf.gz
  bcftools index -t -f ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.vcf.gz
  ${shapeit5} --input ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.vcf.gz \
            --pedigree ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}_shapeit5.fam \
            --region ${chr} \
            --map ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}_shapeit5.bim \
            --output ${path_data}/seq_FINAL_PHASED_${mask}_chr${chr}.bcf \
            --thread 8

  #Convert to haps format
  plink2 --bcf ${path_data}/seq_FINAL_PHASED_${mask}_chr${chr}.bcf --export haps --out ${path_data}/seq_FINAL_PHASED_${mask}_chr${chr}

  #the .sample file needs to contain the parents info
  mv ${path_data}/seq_FINAL_PHASED_${mask}_chr${chr}.sample ${path_data}/seq_FINAL_PHASED_${mask}_chr${chr}.sample~
  cat <(printf "ID_1 ID_2 missing father mother sex\n0 0 0 D D D\n") <(awk -F'\t' '{print $0 "\t0"}' ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam | awk -F'\t' '{print $1 "\t" $2 "\t" $7 "\t" $3 "\t" $4 "\t" $5}' | tr "\t" " ") > ${path_data}/seq_FINAL_PHASED_${mask}_chr${chr}.sample

  #Infering variants with likely genotyping errors based on the duoHMM algorith
  ${duohmm} -H ${path_data}/seq_FINAL_PHASED_${mask}_chr${chr} \
  	    -M ${path_gen_map}/plink.chr${chr}.GRCh38.map \
    	    -G ${path_data}/seq_FINAL_PHASED_${mask}_chr${chr}.generr

  #Remove files what wont be used
  rm ${path_data}/seq_FINAL_with_mask_in_genmap_chr${chr}.fam ${path_data}/seq_FINAL_with_mask_in_genmap_chr${chr}.fam~ ${path_data}/seq_FINAL_with_mask_in_genmap_chr${chr}.bed ${path_data}/seq_FINAL_with_mask_in_genmap_chr${chr}.bim ${path_data}/seq_FINAL_with_mask_in_genmap_chr${chr}.vcf
  rm ${path_data}/seq_FINAL_with_mask_in_genmap_chr${chr}_shapeit5.fam ${path_data}/seq_FINAL_with_mask_in_genmap_chr${chr}_shapeit5.bim
  rm ${path_data}/seq_FINAL_PHASED_with_mask_chr${chr}.haps ${path_data}/seq_FINAL_PHASED_with_mask_chr${chr}.sample ${path_data}/seq_FINAL_PHASED_with_mask_chr${chr}.sample~ ${path_data}/seq_FINAL_PHASED_with_mask_chr${chr}.log
done