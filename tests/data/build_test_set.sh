#!/bin/bash
NGS_base=/home/sqreb/nfs/rna4_data/cotton/new_genome/gosSpp2_FL_NGS_rep2_star
three_GS=/home/sqreb/nfs/rna4_data/cotton/new_genome/gosSpp2_iso_star/all_fixed/all_fixed_star.sort.bam
TSS=/home/sqreb/nfs/rna4_data/cotton/new_genome/identify_element/code/TSS.csv
TES=/home/sqreb/nfs/rna4_data/cotton/new_genome/identify_element/code/TES.csv
loc=test_loc.csv
python3 build_test_set.py --bam ${NGS_base}/*.all.sort.bam ${three_GS} --TXS ${TSS} ${TES} -o . --loc ${loc}
