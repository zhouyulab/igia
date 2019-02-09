#!/bin/bash

# to speedup test, you could use smaller number of --ngs data.
igia  -g data/genome.fa --size data/chrom.size \
    --tgs data/all_fixed_star.sort.bam \
    --tss data/TSS.csv \
    --tes data/TES.csv \
    --ngs data/Anther.sort.bam data/Bract.sort.bam data/Cotyledon.sort.bam \
    -r "1+-,1-+,2++,2--" \
    -v \
    -o igia_demo

