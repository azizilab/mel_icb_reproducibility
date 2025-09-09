#!/bin/bash

### Numbat preprocessing
### Author: Jana Biermann, PhD

declare -a samples
samples=(
F01_pre
F01_on
F02_pre
F02_on
F03_post1_pre2
F03_post1_on2
F04_pre
F05_pre
F06_post1_pre2
F07_post1_pre2
F08_post
F09_post
F10_post
F12_pre
F12_post
F15_pre
F16_pre
F16_post1_pre2
F17_post
F18_post
F20_post1_pre2
F22_post
F23_post
F25_post
F26_post
F27_post1_pre2
F28_post1_pre2
F29_post1_pre2
F30_post
F31_post
R310_pre
R310_on1
R204_pre
R294_on
R308_pre
R310_on2
R319_pre
R319_on
R328_on
R329_on
R334_pre
R354_pre)

for pat in ${samples[@]}; do
  aws s3 cp s3://melanoma-ribas/cellranger/v7.0.0/${pat}/possorted_genome_bam.bam data/${pat}/ --quiet
  aws s3 cp s3://melanoma-ribas/cellranger/v7.0.0/${pat}/possorted_genome_bam.bam.bai data/${pat}/ --quiet
  
  if [ -f "data/${pat}/barcodes_${pat}_filtered.csv" ]; then
    echo "BC file for ${pat} exists."
  else 
    Rscript brain_mets/melanoma/numbat/numbat_1_get_filtered_bc.R
  fi

  brain_mets/melanoma/numbat/pileup_and_phase.R --label ${pat} \
  --samples ${pat} \
  --bams data/${pat}/possorted_genome_bam.bam \
  --barcodes data/${pat}/barcodes_${pat}_filtered.csv \
  --gmap Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
  --snpvcf numbat/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf \
  --paneldir numbat/1000G_hg38 \
  --outdir data/${pat} \
  --ncores 32

  rm data/${pat}/possorted_genome_bam.*

  aws s3 sync data/${pat}/ s3://melanoma-ribas/numbat/${pat}/ --quiet
done
