#Phasing using ShapeIt5 2023-13-22

#The softwares are here
shapeit5_common=/lustre03/project/6033529/SOFT/shapeit5/phase_common/bin/phase_common
shapeit5_rare=/lustre03/project/6033529/SOFT/shapeit5/phase_rare/bin/phase_rare
duohmm=/lustre03/project/6033529/SOFT/duohmm_v0.1.7/duohmm

#The chunks that we use are provided by the authors of ShapeIt5. It is based on the UKBB data
#These chunks only help to parallelize the phasing of rare variants...
chunks=/lustre03/project/6033529/SOFT/shapeit5/resources/chunks/b38/UKB_WGS_200k/

#We will use the QCed data with mask.
mask=with_mask

#Results will be exported here:
path_data=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/phasing

#The hg38 genetic map references from the HapMap are downloaded here:
gen_map=/lustre03/project/6033529/IBD_epilepsy/data/genetic_map_grch38

#GCbroad phenotype file, will be used to complete the pedigree.
pheno=/lustre03/project/6033529/quebec_10x/data/phenotypes/GCbroad.pre
cut -f-6 ${pheno} > ${path_data}/seq_pedigree.fam

for chr in {1..22}
do
  #Our QCed sequencing data are found here (we used the data with mask):
  seq_data=/lustre03/project/6033529/quebec_10x/data/WGS_bs_2022/500_samples_cag/QC/seq_FINAL_chr${chr}_with_mask.vcf.gz

  #Generate bfiles.
  plink2 --vcf ${seq_data} --make-bed --out ${path_data}/seq_FINAL_chr${chr}_with_mask

  #Add genetic positions. Be sure that every markers are in the genetic map, if not, trimming is needed to avoid errors.
  Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022/call/gen_pos_interpolation_grch38.R -bim ${path_data}/seq_FINAL_chr${chr}_with_mask.bim -ref_dir ${gen_map} -method "linear"

  #Rearrange genetic maps for ShapeIt5
  cat <(echo -e "pos\tchr\tcM") <(awk 'BEGIN{OFS="\t"} {print $4, $1, $3}' ${gen_map}/plink.chr${chr}.GRCh38.map) > ${path_data}/gen_map/plink.chr${chr}.GRCh38.map

  #We wish to phase on a dataset where the variants that are not included in the genetic map are removed.
  cat ${gen_map}/plink.chr${chr}.GRCh38.map | cut -d " " -f 1,4 | sed -n '1p;$p' > ${path_data}/gen_map/gen_pos_ranges_chr${chr}_tmp.txt
  awk '{a[$1] = a[$1] " " $2} END {for (i in a) print i a[i]}' ${path_data}/gen_map/gen_pos_ranges_chr${chr}_tmp.txt | awk 'BEGIN {OFS="\t"} {print $0, $1}' > ${path_data}/gen_map/gen_pos_ranges_chr${chr}.txt; rm ${path_data}/gen_map/gen_pos_ranges_chr${chr}_tmp.txt
  plink2 --bfile ${path_data}/seq_FINAL_chr${chr}_with_mask --extract range ${path_data}/gen_map/gen_pos_ranges_chr${chr}.txt --make-bed -out ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}
  rm ${path_data}/seq_FINAL_chr${chr}_with_mask.*

  #Generate the .vcf versions without the variants
  plink2 --bfile ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr} --recode vcf -out ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}

  #Create fam and map files in the format of shapeit5. In the .fam file, the fam needs to be appended to the IDs
  cut -f2- ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam > ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam~
  paste <(cut -d "_" -f1 ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam~) <(cut -d "_" -f2- ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam~) > ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam
  Rscript /lustre03/project/6033529/quebec_10x/scripts/WGS_bs_2022_500samples/call/update_pheno.R -fam ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam -ref ${path_data}/seq_pedigree.fam
  awk -v OFS='\t' '{print $1"_"$2, $1"_"$3, $1"_"$4}' ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam | awk -v OFS='\t'  -F' ' '{for(i=1;i<=NF;i++) if($i~/_0$/) $i=0}1' > ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}_shapeit5.fam

  #Phase common variants
  bcftools +fill-tags ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.vcf | bgzip -c > ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.vcf.gz
  bcftools index -t -f ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.vcf.gz
  ${shapeit5_common} --input ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.vcf.gz \
                     --filter-maf 0.001 \
                     --pedigree ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}_shapeit5.fam \
                     --region ${chr} \
                     --map ${path_data}/gen_map/plink.chr${chr}.GRCh38.map \
                     --output ${path_data}/seq_FINAL_PHASED_${mask}_chr${chr}_scaffold.bcf \
                     --thread 8

  #Run the phasing of rare variants by chunks, as proposed by the authors.
  sed "s/chr//g" ${chunks}/large_chunks_25cM/chunks_chr${chr}.txt > ${path_data}/chunks/chunks_chr${chr}.txt
  while read LINE; do
	  CHK=$(echo $LINE | awk '{ print $1; }')
	  SRG=$(echo $LINE | awk '{ print $3; }')
	  IRG=$(echo $LINE | awk '{ print $4; }')
	  ${shapeit5_rare} --input ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.vcf.gz \
                     --scaffold ${path_data}/seq_FINAL_PHASED_${mask}_chr${chr}_scaffold.bcf \
                     --pedigree ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}_shapeit5.fam \
                     --map ${path_data}/gen_map/plink.chr${chr}.GRCh38.map \
                     --input-region ${IRG} \
                     --scaffold-region ${SRG} \
                     --output ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}_chunk${CHK}.bcf \
                     --thread 8
  done < ${path_data}/chunks/chunks_chr${chr}.txt

  #Merge the phased files.
  chunks_file=$(cut -f 1 ${chunks}/large_chunks_25cM/chunks_chr${chr}.txt)
  for CHK in ${chunks_file}; do echo "${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}_chunk${CHK}.bcf" >> ${path_data}/chunks.files_chr${chr}.txt; done
  bcftools concat -n -Ob -o ${path_data}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf -f ${path_data}/chunks.files_chr${chr}.txt
  bcftools index ${path_data}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf

  #Convert to haps format
  plink2 --bcf ${path_data}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.bcf --export haps --out ${path_data}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}

  #the .sample file needs to contain the parents info
  mv ${path_data}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.sample ${path_data}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.sample~
  cat <(printf "ID_1 ID_2 missing father mother sex\n0 0 0 D D D\n") <(awk -F'\t' '{print $0 "\t0"}' ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam | awk -F'\t' '{print $1 "\t" $2 "\t" $7 "\t" $3 "\t" $4 "\t" $5}' | tr "\t" " ") > ${path_data}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.sample

  #Infering variants with likely genotyping errors based on the duoHMM algorith
  ${duohmm} -H ${path_data}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr} \
  	    -M ${path_data}/gen_map/plink.chr${chr}.GRCh38.map \
    	    -G ${path_data}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.generr

  #Remove files what wont be used 
  rm ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.fam~ ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.bed ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.bim ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.vcf* ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}.log
  rm ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}_shapeit5.fam
  rm ${path_data}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.haps ${path_data}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.sample ${path_data}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.sample~ ${path_data}/seq_FINAL_PHASED_${mask}_in_genmap_chr${chr}.log
  rm ${path_data}/chunks.files_chr${chr}.txt ${path_data}/seq_FINAL_PHASED_${mask}_chr${chr}_scaffold.* ${path_data}/seq_FINAL_${mask}_in_genmap_chr${chr}_chunk*
done