from igia.utils import iterline, Bed12
import sys
import argparse
import numpy as np
from bx.intervals.cluster import ClusterTree
from bx.intervals.intersection import IntervalTree


class Cluster(object):
    """ Gene cluster """
    def __init__(self, cluster_name):
        self.name = cluster_name
        self.gene_list = list()
        self.exon_blocks = list()
        self.chrom = None
        self.start = None
        self.end = None
        self.strand = None

    def add_gene(self, gene):
        """
        Add gene to this cluster

        Args:
            gene ():

        Returns:

        """
        if self.strand is None:
            self.chrom = gene.chrom
            self.start = gene.chromStart
            self.end = gene.chromEnd
            self.strand = gene.strand
        else:
            self.start = min(self.start, gene.chromStart)
            self.end = max(self.end, gene.chromEnd)
            assert self.chrom == gene.chrom
            assert self.strand == gene.strand
        self.gene_list.append(gene)

    def build_exon_block(self):
        """Merge exons in all genes"""
        block_tree = ClusterTree(0, 0)
        for gene in self.gene_list:
            for start, end in gene.iterblock():
                block_tree.insert(start, end, 0)

        self.exon_blocks = [(x[0], x[1]) for x in block_tree.getregions()]

    @classmethod
    def block2position(cls, blocks):
        position_list = list()
        for start, end in blocks:
            position_list += list(range(start, end))
        return position_list

    def compute_exon_overlap_len(self, gene):
        if self.chrom != gene.chrom:
            return 0
        if self.strand != gene.strand:
            return 0
        cluster_position = set(self.block2position(self.exon_blocks))
        gene_position = set(self.block2position(gene.iterblock()))
        overlap_len = len(cluster_position.intersection(gene_position))
        return overlap_len

    def write_mapping(self, f):
        for gene in self.gene_list:
            f.write("{0}\t{1}\n".format(self.name, gene.name))

    def write2bed6(self, f):
        text = "{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n".format(
            chrom=self.chrom, start=self.start, end=self.end,
            name=self.name, score=len(self.gene_list), strand=self.strand
        )
        f.write(text)


def load_gene(f_gene):
    """
    Load gene information from Bed12 format file

    Args:
        f_gene (str): Gene annotation

    Returns:
        dict: Gene annotation
    """
    gene_dict = dict()
    for line in iterline(f_gene):
        gene = Bed12(line)
        gene_dict[gene.name] = gene
    return gene_dict


def load_gl(f_gl):
    """Generate cluster:gene from resulted file of clusterGenes"""
    for line in iterline(f_gl):
        data = line.split("\t")
        cluster = data[0]
        gene = data[2]
        yield (cluster, gene)


def assign_gene_to_cluster(gene_list, cluster_tree, f_out_no_overlap):
    """Assign gene to most overlaped cluster

    :param gene_list: gene list
    :param cluster_tree:
    :param f_out_no_overlap: Output file to write
    """
    for gene in gene_list:
        key = (gene.chrom, gene.strand)
        overlaped_clusters = cluster_tree[key].find(gene.chromStart, gene.chromEnd)
        if not overlaped_clusters:
            #sys.stderr.write("{0} not overlapped with any cluster\n".format(gene.name))
            gene.write(f_out_no_overlap)
            continue

        exon_overlap_len = [
            cluster.compute_exon_overlap_len(gene) 
            for cluster in overlaped_clusters]

        if max(exon_overlap_len) == 0:
            #sys.stderr.write("{0} not overlapped with any cluster on exon\n".format(gene.name))
            gene.write(f_out_no_overlap)
            continue

        # assign gene to the most overlapped cluster
        overlaped_clusters[np.argmax(exon_overlap_len)].add_gene(gene)


def merge_gene_into_cluster(args):
    """
    Merge external genes into clusters of genes from clusterGenes
    """
    args = parse_args(args)
    f_gl = args.f_gl
    f_gl_gene = args.f_gl_gene
    f_ext_gene = args.f_ext_gene
    f_out = args.f_out
    f_out_no_overlap = args.f_out_no_overlap

    print("Loading gl gene ...")
    gl_gene_dict = load_gene(f_gl_gene)

    print("Loading gl ...")
    cluster_dict = dict()
    for cluster, gene in load_gl(f_gl):
        if cluster not in cluster_dict.keys():
            new_cluster = Cluster(cluster)
            cluster_dict[cluster] = new_cluster
        assert gene in gl_gene_dict.keys(), "Cannot find {0} in {1}".format(gene, f_gl_gene)
        cluster_dict[cluster].add_gene(gl_gene_dict[gene])

    cluster_list = list(cluster_dict.values())
    # Build Chrom:Strand IntervalTree    
    ctree = dict()
    for cluster in cluster_list:
        cluster.build_exon_block()
        key = (cluster.chrom, cluster.strand)
        if key not in ctree:
            ctree[key] = IntervalTree()
        ctree[key].insert(cluster.start, cluster.end, cluster)

    print("Loading external gene ...")
    ext_gene = list(load_gene(f_ext_gene).values())

    print("Assigning gene into clusters ...")
    with open(f_out_no_overlap, "w") as f:
        assign_gene_to_cluster(ext_gene, ctree, f)

    with open(f_out, "w") as f:
        for cluster in cluster_list:
            cluster.write_mapping(f)


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--gl", type=str, dest="f_gl", metavar="cluster.gl", required=True)
    parser.add_argument("--gl-gene", type=str, dest="f_gl_gene", metavar="gl_gene.bed12", required=True)
    parser.add_argument("--ins-gene", type=str, dest="f_ext_gene", metavar="ext_gene.bed12", required=True)
    parser.add_argument("-o", type=str, dest="f_out", metavar="output.csv", required=True)
    parser.add_argument("--no-overlap", type=str, dest="f_out_no_overlap", metavar="no_overlap.bed12", required=True)
    return parser.parse_args(args)


def run():
    merge_gene_into_cluster(sys.argv[1:])


if __name__ == "__main__":
    run()
