#!/bin/bash
# path to the genome reference:
# gzip -d hg38_genome.gz
# awk '$3 == "CDS"' "/home/alu/netlandes/MHCpan/hg38_genome" > hg38_genome_only_cds.gtf
# awk -F'\t' -v OFS="\t" '{print $1,$4,$5,$9,$7,$8}' "/home/alu/netlandes/MHCpan/hg38_genome_only_cds" > bed6_of_genom_38.bed
genome="/home/alu/netlandes/MHCpan/bed6_of_ensemble.bed"
chr="$1"
mut_position="$2"
k="$3"
start_seq=$((mut_position-15))
end_seq=$((mut_position+15))
# path to create temp file of BED format for query
out_temp="/home/alu/netlandes/MHCpan/temp_folder/TEMP${k}.bed"
# create the temp file with the BED entry  $genome
printf "%s\t%s\t%s\n" "$chr" "$start_seq" "$end_seq" > "$out_temp"
bedtools intersect -wa -a "$genome" -b "$out_temp" 
