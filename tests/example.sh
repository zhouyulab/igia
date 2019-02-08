#!/bin/bash

igia  -g data/genome.fa --size data/chrom.size \
    --tgs data/all_fixed_star.sort.bam \
    --tss data/TSS.csv \
    --tes data/TES.csv \
    --ngs data/Anther.sort.bam data/Bract.sort.bam data/Cotyledon.sort.bam \
    data/Flower-0DPA.sort.bam data/Ovule-0DPA.sort.bam \
    data/Ovule-10DPA.sort.bam data/Ovule-20DPA.sort.bam \
    data/Ovule-5DPA.sort.bam data/Petal.sort.bam data/Petiole.sort.bam \
    data/Phloem.sort.bam data/Seedling_root.sort.bam \
    data/Seedling_stem.sort.bam data/Seed.sort.bam \
    data/Sepal.sort.bam data/Stamen.sort.bam \
    -r "1+-,1-+,2++,2--" \
    -v \
    -o igia_res
