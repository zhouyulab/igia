import pysam
import argparse
import os

def pick_bam(input, output, loc_list):
    bam_file = pysam.AlignmentFile(input, "rb")
    with pysam.AlignmentFile(output, "wb", template=bam_file) as f:
        for chrom, start, end in loc_list:
            reads = list(bam_file.fetch(chrom, start, end))
            reads = filter(lambda x: not x.reference_end is None, reads)
            reads = filter(lambda x: not x.reference_start is None, reads)
            reads = list(filter(lambda x: (x.reference_start >= start) and (x.reference_end <= end), reads))
            for read in reads:
                f.write(read)
    bam_file.close()

def pick_TXS(input, output, loc_list):
    with open(input, "r") as f:
        TXS_lines = f.readlines()
    TXS_dict = dict()
    for line in TXS_lines:
        chrom, loc, strand = line.rstrip("\n").split("\t")
        loc = int(loc)
        if chrom in TXS_dict.keys():
            TXS_dict[chrom][loc] = line
        else:
            TXS_dict[chrom] = {loc: line}

    with open(output, "w") as f:
        for chrom, start, end in loc_list:
            if chrom in TXS_dict.keys():
                locs = TXS_dict[chrom].keys()
                locs_list = list(filter(lambda x: (x >= start) and (x <= end), locs))
                for loc in locs_list:
                    f.write(TXS_dict[chrom][loc])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", nargs="*", type=str, dest="bam", metavar="file.bam", required=True)
    parser.add_argument("-o", type=str, dest="out_folder", metavar="out_folder", required=True)
    parser.add_argument("--TXS", nargs="*", type=str, dest="TXS", metavar="TXS.csv", required=True)
    parser.add_argument("--loc", type=str, dest="loc", metavar="loc.csv", required=True)
    args = parser.parse_args()
    with open(args.loc, "r") as f:
        loc_list = list(map(lambda x: x.rstrip("\n").split("\t"), f.readlines()))
        loc_list = list(map(lambda x: [x[0], int(x[1]), int(x[2])], loc_list))
    for bam_input in args.bam:
        file_name = os.path.split(bam_input)[1]
        out_file_name = file_name.split(".")[0] + ".bam"
        bam_output = os.path.join(args.out_folder, out_file_name)
        pick_bam(bam_input, bam_output, loc_list)
        sorted_file_name = file_name.split(".")[0] + ".sort.bam"
        pysam.sort("-o",sorted_file_name,out_file_name)
        pysam.index(sorted_file_name)
    for TXS_input in args.TXS:
        file_name = os.path.split(TXS_input)[1]
        out_file_name = file_name.split(".")[0] + ".csv"
        TXS_output = os.path.join(args.out_folder, out_file_name)
        pick_TXS(TXS_input, TXS_output, loc_list)
