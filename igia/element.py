from igia.utils import Interval, GenomeFile, AlignReadMethod
import igia.GVAR as GVAR
import networkx
import numpy as np
import sys, os
import time
import copy


class GeneLinkage(object):
    """ Linkage of read cluster. """
    def __init__(self, chrom, start, end, ngs_bam_list, tgs_bam_list, ann=None):
        self.ival = Interval(chrom, start, end, ".")
        self.ngs_bam_list = ngs_bam_list
        self.tgs_bam_list = tgs_bam_list
        self.ann = ann
        self.chrsize = self.ngs_bam_list[0].chrom2size()[chrom]
        self.tgs_read_list = list()
        self.ann_read_list = list()

    def __str__(self):
        return "GeneLinkage: {0}:{1}-{2}".format(
            self.ival.chrom, self.ival.start, self.ival.end)

    __repr__ = __str__

    def filter_txs(self, txs_list):
        """
        Select txs in the linkage
        Args:
            txs_list (list): TXS list to filter

        Returns:
            list: List of TXS
        """
        txs_list = list(filter(
            lambda x: (self.ival.start <= x[1] <= self.ival.end) and (x[0] == self.ival.chrom), txs_list))
        return txs_list


    def fetch_reads(self):
        """
        Fetch reads in this linkage

        Returns:
            list: List of reads
        """
        # Fetch annotation reads
        if self.ann is not None:
            self.ann_read_list = self.ann.fetch_reads_in_ival(self.ival)

        # Fetch TGS reads
        for tgs_bam_obj in self.tgs_bam_list:
            self.tgs_read_list += tgs_bam_obj.fetch_reads_in_ival(self.ival)

    def find_clusters(self):
        """ Compute gene cluster """
        clusters = GeneLinkageFinder.build_tgs_cluster(
            self.tgs_read_list + self.ann_read_list, self.ngs_bam_list, self.ival.chrom,
            GVAR.TGS_OVERLAP_RATIO, GVAR.TXS_DIFF,
            GVAR.TGS_STRAND_FRACTION, GVAR.SHARE_INTRON_CUTOFF, GVAR.NGS_LIMIT_READ_COVER_PER_BASE)
        return clusters

    def cluster2gene(self, cluster, f_genome=None):
        """
        Convert Cluster object to Gene object

        Args:
            cluster (list): Reads cluster
            f_genome (str): Genome file

        Returns:
            Gene: Gene object
        """
        if cluster[0].is_reverse:
            cluster_strand = "-"
        else:
            cluster_strand = "+"

        cluster_start = max(min([x.reference_start for x in cluster]) - 2 * GVAR.TXS_DIFF, 1)
        cluster_end = min(max([x.reference_end for x in cluster]) + 2 * GVAR.TXS_DIFF, self.chrsize)

        ann_reads = list(filter(lambda x: x.get_tag('TP') == "ANN", cluster))
        tgs_read_cluster = list(filter(lambda x: x.get_tag('TP') == "TGS", cluster))
        return Gene(self.ival.chrom, cluster_start, cluster_end, cluster_strand,
                    tgs_read_cluster, self.ngs_bam_list, ann_reads, f_genome)

    def split2gene(self, f_genome=None):
        """
        Split GeneLinkage into Gene cluster.

        Args:
            f_genome (str): Genome file

        Returns:
            list: List of gene
        """
        # Fetch TGS and ANN reads
        self.fetch_reads()

        # Clustering reads
        clusters = self.find_clusters()

        # Cluster2gene
        gene_list = list()
        for cluster in clusters:
            gene = self.cluster2gene(cluster, f_genome)
            gene_list.append(gene)
        gene_list.sort(key=lambda x: x.ival.start)
        return gene_list


class GeneLinkageFinder(object):
    """ Abstract class for gene linkage """
    @classmethod
    def compute_read_overlap_ratio(cls, read1, read2, strand=True):
        """
        Compute overlap ratio between two reads.
        Args:
            read1 (AlignedSegment): Read1
            read2 (AlignedSegment): Read2
            strand (bool): strand specific

        Returns:
            tuple: (overlap_len/read1, overlap_len/read2)
        """
        if strand:
            if read1.is_reverse != read2.is_reverse:
                return (0, 0)
        if read1.reference_start > read2.reference_end:
            return (0, 0)
        if read1.reference_end < read2.reference_start:
            return (0, 0)
        overlap_len = len(set(read1.positions) & set(read2.positions))
        return (overlap_len / len(read1.positions), overlap_len / len(read2.positions))

    @classmethod
    def build_exon_overlap_cluster(cls, reads, overlap_ratio=0.5):
        """
        Build exon overlap cluster
        Args:
            reads (list): List of reads
            overlap_ratio (float): Overlap ratio to label two read in one cluster

        Returns:
            list: List of reads cluster
        """
        # Build exon overlap based graph.
        exon_ovarlap_ratio_graph = networkx.Graph()
        for indx1 in range(len(reads)):
            exon_ovarlap_ratio_graph.add_edge(indx1, indx1)
            for indx2 in range(indx1):
                if max(cls.compute_read_overlap_ratio(
                        reads[indx1], reads[indx2])) > overlap_ratio:
                    exon_ovarlap_ratio_graph.add_edge(indx1, indx2)
        clusters = list()
        for connected_component_indxs in networkx.connected_components(exon_ovarlap_ratio_graph):
            clusters.append([reads[x] for x in connected_component_indxs])
        return clusters

    @classmethod
    def split_cluster_by_overlap(cls, cluster):
        """ Split cluster by if overlap statement """
        read_len = len(cluster)
        ovarlap_graph = networkx.Graph()

        for indx1 in range(read_len):
            ovarlap_graph.add_edge(indx1, indx1)
            for indx2 in range(indx1):
                read1 = cluster[indx1]
                read2 = cluster[indx2]
                if (read1.reference_start >= read2.reference_end) or (read2.reference_start >= read1.reference_end):
                    continue
                ovarlap_graph.add_edge(indx1, indx2)

        # The relationship_tuple represents the overlap relationship of one isoform to others.
        # Using the isoforms which have same relationship_tuple to form a cluster
        cliques = list(networkx.find_cliques(ovarlap_graph))
        cliques_dict = dict()
        for indx in range(read_len):
            relationship_tuple = tuple([indx in x for x in cliques])
            if relationship_tuple in cliques_dict.keys():
                cliques_dict[relationship_tuple].append(indx)
            else:
                cliques_dict[relationship_tuple] = [indx]

        new_clusters = list()
        for indxs in cliques_dict.values():
            sub_cluster = [cluster[indx] for indx in indxs]
            assert all([(read.is_reverse == sub_cluster[0].is_reverse) for read in sub_cluster])
            new_clusters.append(sub_cluster)
        return new_clusters

    @classmethod
    def merge_clusters(cls, clusters):
        """Merge clusters into reads"""
        reads = list()
        for cluster in clusters:
            reads += cluster
        return reads

    @classmethod
    def find_intron(cls, read):
        """
        Find intron in a read
        Args:
            read (AlignedSegment): Read

        Returns:
            list: List of intron
        """
        blocks = read.blocks
        if len(blocks) < 2:
            return list()
        block_start = [x[0] for x in blocks]
        block_end = [x[1] for x in blocks]
        return list(zip(block_end[:-1], block_start[1:]))

    @classmethod
    def find_intron_in_cluster(cls, cluster):
        """
        Find all intron in a cluster

        Args:
            cluster (list): Reads cluster

        Returns:
            set: Set of intron
        """
        introns = set()
        for read in cluster:
            introns |= set(cls.find_intron(read))
        return introns

    @classmethod
    def filter_cluster_by_strand(
            cls, cluster_list, strand_fraction=0.5, intron_cutoff=0.1, exon_cutoff=0.5):
        """Clusters sharing same introns are considered as one gene"""
        cluster_num = len(cluster_list)
        intron_list = list()
        read_num_list = list()
        strand_list = list()
        for cluster in cluster_list:
            intron_list.append(cls.find_intron_in_cluster(cluster))
            read_num_list.append(sum([read.get_tag('FL') for read in cluster]))
            strand_list.append(cluster[0].is_reverse)

        loc_list = [
            (min([read.reference_start for read in cluster]),
             max([read.reference_end for read in cluster]))
            for cluster in cluster_list]

        # Build sharing intron radio matrix
        sharing_intron_matrix = np.zeros([cluster_num, cluster_num]) + np.diag(np.ones(cluster_num))
        for indx1 in range(cluster_num):
            for indx2 in range(indx1):
                first_loc = loc_list[indx1]
                second_loc = loc_list[indx2]
                # If two clusters are not overlapped, continue as independent
                if (first_loc[0] >= second_loc[1]) or (second_loc[0] >= first_loc[1]):
                    continue

                fitst_cluster_introns = intron_list[indx1]
                second_cluster_introns = intron_list[indx2]
                same_introns = fitst_cluster_introns.intersection(second_cluster_introns)
                first_intron_num = len(fitst_cluster_introns)
                second_intron_num = len(second_cluster_introns)
                same_intron_num = len(same_introns)
                # At least one cluster is non-spliced
                if (first_intron_num == 0) or (second_intron_num == 0):
                    exon_overlap_ratio = cls.compute_read_overlap_ratio(
                        cluster_list[indx1][0], cluster_list[indx2][0], False)
                    if max(exon_overlap_ratio) > exon_cutoff:
                        sharing_intron_matrix[indx1, indx2] = 1
                        sharing_intron_matrix[indx2, indx1] = 1
                    else:
                        sharing_intron_matrix[indx1, indx2] = 0
                        sharing_intron_matrix[indx2, indx1] = 0
                # Two cluster are both spliced
                else:
                    sharing_intron_matrix[indx1, indx2] = same_intron_num / first_intron_num
                    sharing_intron_matrix[indx2, indx1] = same_intron_num / second_intron_num

        # Filter cluster by sharing intron clusters' strand
        clean_cluster_list = list()
        for indx in range(cluster_num):
            overlapped_cluster_indxs = np.where(sharing_intron_matrix[indx, :] > intron_cutoff)[0]
            same_strand_num = 0
            diff_strand_num = 0
            for overlap_indx in overlapped_cluster_indxs:
                if strand_list[indx] == strand_list[overlap_indx]:
                    same_strand_num += read_num_list[overlap_indx]
                else:
                    diff_strand_num += read_num_list[overlap_indx]
            if same_strand_num / (same_strand_num + diff_strand_num) >= strand_fraction:
                clean_cluster_list.append(cluster_list[indx])
        return clean_cluster_list

    @classmethod
    def filter_nonspliced_cluster_by_cov(cls, cluster, bam_list, chrom, cov_cutoff, min_tgs_num=2):
        """ Non-spliced cluster should have enough FPKM """
        if cls.find_intron_in_cluster(cluster) and len(cluster) >= min_tgs_num:
            return True
        positions = set()
        for read in cluster:
            positions |= set(read.positions)
        boundary_start = min([x.reference_start for x in cluster])
        boundary_end = max([x.reference_end for x in cluster])
        positions = np.array(sorted(positions)) - boundary_start
        if cluster[0].is_reverse:
            strand = "-"
        else:
            strand = "+"
        boundary_ival = Interval(chrom, boundary_start, boundary_end, strand)
        boundary_ival.build_cov(bam_list)
        if np.any(boundary_ival.cov.sig[:, positions].mean(1).reshape([len(bam_list), 1]) > cov_cutoff):
            return True
        else:
            return False

    @classmethod
    def build_tgs_cluster(cls, reads, bam_list, chrom, overlap_ratio=0.5,
                          txs_diff=400, strand_fraction=0.5, intron_cutoff=0.1, cov_cutoff=None):
        """Use strand filtered tgs reads to build exonic overlapped gene cluster."""
        if cov_cutoff is None:
            cov_cutoff = np.array([1] * len(bam_list)).reshape([len(bam_list), 1])

        start_read_num = len(reads)

        # Build clusters by exon overlap
        exon_overlap_clusters = cls.build_exon_overlap_cluster(reads, overlap_ratio)

        # Filter cluster by strand
        clean_clusters = cls.filter_cluster_by_strand(
            exon_overlap_clusters, strand_fraction, intron_cutoff, overlap_ratio)

        # Filter cluster by cov
        clean_clusters = list(
            filter(
                lambda x: cls.filter_nonspliced_cluster_by_cov(
                    x, bam_list, chrom, cov_cutoff), clean_clusters))

        # If some reads are removed, build cluster again
        passed_reads = cls.merge_clusters(clean_clusters)
        passed_read_num = len(passed_reads)
        if passed_read_num < start_read_num:
            clean_clusters = cls.build_tgs_cluster(
                passed_reads, bam_list, chrom, overlap_ratio,
                txs_diff, strand_fraction, intron_cutoff)
        return clean_clusters


class Gene(object):
    """ Linkage of gene cluster. """
    def __init__(self, chrom, start, end, strand, tgs_read_list, ngs_bam_list, ann_read_list=None, f_genome=None):
        if ann_read_list is None:
            ann_read_list = list()

        self.ival = Interval(chrom, start, end, strand)
        self.tgs_read_list = tgs_read_list
        self.ngs_bam_list = ngs_bam_list
        self.ann_read_list = ann_read_list
        self.intron_list = None
        self.internal_exon_list = None
        self.tss_site_list = None
        self.tes_site_list = None
        self.tss_exon_list = None
        self.tes_exon_list = None
        self.build_cov()
        if f_genome is not None:
            self.genome = GenomeFile(f_genome)
        else:
            self.genome = None

    def __str__(self):
        return "Gene: {0}:{1}-{2} {3}".format(
            self.ival.chrom, self.ival.start, self.ival.end, self.ival.strand)

    __repr__ = __str__

    def build_cov(self):
        """ Create wiggle like signal coverage from bam files """
        self.ival.build_cov(self.ngs_bam_list)

    def identify_intron(self):
        """ Identify intron glue code. """
        # Due to NGS intron obj have spliced info, use NGS intron first.
        assert self.ival.cov is not None
        intron_set_ngs = set(NgsElementDiscover.identify_intron(self.ival, self.ngs_bam_list))
        intron_set_tgs = set(TgsElementDiscover.identify_intron(self.ival, self.tgs_read_list))
        self.intron_list = ElementDiscover.fix_mapping_error(
            intron_set_ngs, intron_set_tgs, self.genome, self.tgs_read_list, self.ngs_bam_list)

    def identify_tss_site(self, ext_tss_site_list):
        """ Identify tss site glue code. """
        assert self.intron_list is not None

        tss_site_list = copy.deepcopy(ext_tss_site_list)

        tgs_upsmost_site = -1
        if self.tgs_read_list:
            if self.ival.strand == "+":
                tgs_upsmost_site = np.min([x.reference_start for x in self.tgs_read_list])
            else:
                tgs_upsmost_site = np.max([x.reference_end for x in self.tgs_read_list])

        # If TGS end is farther than known TSS site, add it.
        if tgs_upsmost_site > 0:
            if tss_site_list:
                if self.ival.strand == "+":
                    tgs_ext_dist = min(tss_site_list) - tgs_upsmost_site
                else:
                    tgs_ext_dist = tgs_upsmost_site - max(tss_site_list)

                if tgs_ext_dist >= GVAR.TXS_DIFF:
                    tss_site_list.append(tgs_upsmost_site)
            else:
                tss_site_list.append(tgs_upsmost_site)

        if not tss_site_list:
            if self.ival.strand == "+":
                tss_site_list = list({x.reference_start for x in self.ann_read_list})
            else:
                tss_site_list = list({x.reference_end for x in self.ann_read_list})
        self.tss_site_list = tss_site_list

    def identify_tes_site(self, ext_tes_site_list):
        """ Identify tes site glue code. """
        assert self.intron_list is not None
        tes_site_list = copy.deepcopy(ext_tes_site_list)

        tgs_downsmost_site = -1
        if self.ival.strand == "+":
            if self.tgs_read_list:
                tgs_downsmost_site = np.max(
                    [x.reference_end for x in self.tgs_read_list])
        else:
            if self.tgs_read_list:
                tgs_downsmost_site = np.min(
                    [x.reference_start for x in self.tgs_read_list])

        # If TGS end is farther than known TES site, add it.
        if tgs_downsmost_site > 0:
            if tes_site_list:
                if self.ival.strand == "+":
                    tgs_ext_dist = tgs_downsmost_site - max(tes_site_list)
                else:
                    tgs_ext_dist = min(tes_site_list) - tgs_downsmost_site
                if tgs_ext_dist >= GVAR.TXS_DIFF:
                    tes_site_list.append(tgs_downsmost_site)
            else:
                tes_site_list.append(tgs_downsmost_site)

        # If not TES site, try to use annotation tes.
        if not tes_site_list:
            if self.ival.strand == "+":
                tes_site_list = list(
                    {x.reference_end for x in self.ann_read_list})
            else:
                tes_site_list = list(
                    {x.reference_start for x in self.ann_read_list})
        self.tes_site_list = tes_site_list

    def identify_tss_exon(self):
        """ Identify tss exon glue code. """
        assert self.tss_site_list is not None
        assert self.tes_site_list is not None
        assert self.intron_list is not None
        self.tss_exon_list = ElementDiscover.identify_tss_exon(
            self.tss_site_list, self.tes_site_list, self.intron_list, self.ival)

    def identify_tes_exon(self):
        """ Identify tes exon glue code. """
        assert self.tss_site_list is not None
        assert self.tes_site_list is not None
        assert self.intron_list is not None
        self.tes_exon_list = ElementDiscover.identify_tes_exon(
            self.tss_site_list, self.tes_site_list, self.intron_list, self.ival)

    def identify_internal_exon(self):
        """ Identify internal exon glue code. """
        assert self.intron_list is not None
        self.internal_exon_list, _, _ = ElementDiscover.identify_internal_exon(
            self.intron_list, self.ival)

    def identify_element(self, ext_tss_site_list, ext_tes_site_list, paraclu_path=None):
        """ Identify transcript elements by TGS and NGS data. """
        if (not ext_tss_site_list) and (paraclu_path is not None):
            ext_tss_site_list = ElementDiscover.find_txs_by_tgs(self.tgs_read_list, "tss", paraclu_path)
        if (not ext_tes_site_list) and (paraclu_path is not None):
            ext_tes_site_list = ElementDiscover.find_txs_by_tgs(self.tgs_read_list, "tes", paraclu_path)
        self.identify_intron()
        self.identify_internal_exon()
        self.identify_tss_site(ext_tss_site_list)
        self.identify_tes_site(ext_tes_site_list)
        self.identify_tss_exon()
        self.identify_tes_exon()

    def write_element2bed6(self, f_intron, f_internal_exon, f_tss_exon, f_tes_exon, gene_name):
        """ Write elements into file with Bed6 format """
        assert all(map(lambda x: x is not None,
                       [self.intron_list, self.tss_exon_list,
                        self.tes_exon_list, self.internal_exon_list]))

        for indx, intron in enumerate(self.intron_list):
            intron.write2bed6("{0}_intron_{1}".format(gene_name, indx+1), f_intron)

        for indx, internal_exon in enumerate(self.internal_exon_list):
            internal_exon.write2bed6("{0}_internal_exon_{1}".format(gene_name, indx + 1), f_internal_exon)

        for indx, tss_exon in enumerate(self.tss_exon_list):
            tss_exon.write2bed6("{0}_tss_exon_{1}".format(gene_name, indx + 1), f_tss_exon)

        for indx, tes_exon in enumerate(self.tes_exon_list):
            tes_exon.write2bed6("{0}_tes_exon_{1}".format(gene_name, indx + 1), f_tes_exon)

    def ready4assembly(self):
        return all(map(
            lambda x: x is not None, [
            self.ival.cov.sig, self.intron_list, self.tss_exon_list, self.tes_exon_list, self.internal_exon_list]))

    def has_element(self):
        return len(self.tes_exon_list) > 0 and len(self.tss_exon_list) > 0


class Exon(Interval):
    """ A Structure to store exon data. """
    def __init__(self, chrom, start, end, strand):
        super().__init__(chrom, start, end, strand)
        self.tss_exon = False
        self.tes_exon = False

    def __str__(self):
        return "Exon: {0}:{1}-{2} {3}".format(self.chrom, self.start, self.end, self.strand)

    __repr__ = __str__

    def set_tss_exon(self, is_tss_exon):
        """ Label this exon as TSS exon or not """
        self.tss_exon = is_tss_exon

    def set_tes_exon(self, is_tes_exon):
        """ Label this exon as TES exon or not """
        self.tes_exon = is_tes_exon


class ElementDiscover(object):
    """ General method in identify transcript elements. """

    @classmethod
    def detect_txs_by_tgs(cls, txs_list, paraclu_path, maxlen=20, minDensityRise=4, minReads=5):
        """Detect txs by paraclu"""
        tmp_floder = "/tmp/igia/{0}".format(hash(time.time()))
        if not os.path.exists(tmp_floder):
            os.makedirs(tmp_floder)
        f_paracludata = os.path.join(tmp_floder, "paracludata.csv")

        with open(f_paracludata, "w") as foh:
            for i in sorted(set(txs_list)):
                foh.write("\t".join(map(str, ["chr1", '+', i, txs_list.count(i)])) + "\n")

        f_paraclu = os.path.join(tmp_floder, "paraclu.csv")
        cmdline = "%s %f %s > %s" % (
            os.path.join(paraclu_path, "paraclu"),
            minReads,
            f_paracludata,
            f_paraclu)
        os.system(cmdline)

        f_paraclucut = os.path.join(tmp_floder, "paraclu-cut.csv")
        cmdline = "%s -l %d -d %d %s > %s" % (
            os.path.join(paraclu_path, "paraclu-cut.sh"),
            maxlen,
            minDensityRise,
            f_paraclu,
            f_paraclucut)
        os.system(cmdline)

        paraclu_res = list()
        for line in open(f_paraclucut):
            if line.startswith("chr"):
                f = line.rstrip().split("\t")
                paraclu_res.append((int((int(f[2])+int(f[3]))/2), int(f[4])))
        os.remove(os.path.join(tmp_floder, "paracludata.csv"))
        os.remove(os.path.join(tmp_floder, "paraclu.csv"))
        os.remove(os.path.join(tmp_floder, "paraclu-cut.csv"))
        os.rmdir(tmp_floder)
        res = list()
        for tmp_txs, cnt in sorted(paraclu_res, key=lambda x: x[0]):
            overlap_txs_cnts = [x[1] for x in paraclu_res if abs(x[0] - tmp_txs) <= GVAR.TXS_DIFF]
            if max(overlap_txs_cnts) > cnt:
                continue
            same_cnt_txs = [x[0] for x in paraclu_res if (abs(x[0] - tmp_txs) <= GVAR.TXS_DIFF) and (x[1] == cnt)]
            if (len(same_cnt_txs) > 1) and (sorted(same_cnt_txs).index(tmp_txs) != 0):
                continue
            res.append(tmp_txs)
        return res

    @classmethod
    def find_txs_by_tgs(cls, tgs_list, txs_type, paraclu_path):
        """
        Detect TXS by TGS data

        Args:
            tgs_list (list): List of TGS reads
            txs_type (str): TXS type in [tss, tes]
            paraclu_path (str): Path to paraclu

        Returns:
            list: List of TXS
        """
        assert txs_type in ["tss", "tes"]
        if not tgs_list:
            return list()
        txs_list = list()
        for tgs in tgs_list:
            fl_num = tgs.get_tag('FL')
            if (txs_type == "tss") ^ tgs.is_reverse:
                txs_list += [tgs.reference_start] * fl_num
            else:
                txs_list += [tgs.reference_end] * fl_num
        txs_res = cls.detect_txs_by_tgs(txs_list, paraclu_path)
        return txs_res

    @classmethod
    def identify_tss_exon(cls, tss_site_list, tes_site_list, intron_list, gene_ival):
        """ Predict TSS exon. """
        tss_site_list.sort()
        tss_exon_list = list()
        if gene_ival.strand == "+":
            # [TSS, TES] or [TSS, intron_boundary]
            intron_junction = [x.start for x in intron_list]
            boundary_list = sorted(tes_site_list + intron_junction)
            for tss_site in tss_site_list:
                tss_exon_boundary_list = list(filter(
                    lambda x: (GVAR.MINEXON_LEN <= x - tss_site <= GVAR.MAXEXON_LEN), boundary_list))
                if tss_exon_boundary_list:
                    tss_exon_list.append(Exon(
                        gene_ival.chrom, tss_site, tss_exon_boundary_list[0], gene_ival.strand))
        else:
            # [TES, TSS] or [intron_boundary, TSS]
            intron_junction = [x.end for x in intron_list]
            boundary_list = sorted(tes_site_list + intron_junction)
            for tss_site in tss_site_list:
                tss_exon_boundary_list = list(filter(
                    lambda x: (GVAR.MINEXON_LEN <= tss_site - x <= GVAR.MAXEXON_LEN), boundary_list))
                if tss_exon_boundary_list:
                    tss_exon_list.append(Exon(
                        gene_ival.chrom, tss_exon_boundary_list[-1], tss_site, gene_ival.strand))

        for tss_exon in tss_exon_list:
            tss_exon.inherit_cov_from(gene_ival)
            tss_exon.set_tss_exon(True)
        return tss_exon_list

    @classmethod
    def identify_tes_exon(cls, tss_site_list, tes_site_list, intron_list, gene_ival):
        """ Predict TES exon. """
        tes_site_list.sort()
        tes_exon_list = list()
        if gene_ival.strand == "+":
            # [TSS, TES] or [intron_boundary, TES]
            intron_junction = [x.end for x in intron_list]
            boundary_list = sorted(tss_site_list + intron_junction)
            for tes_site in tes_site_list:
                tes_exon_boundary_list = list(filter(
                    lambda x: (GVAR.MINEXON_LEN <= tes_site - x <= GVAR.MAXEXON_LEN), boundary_list))
                if tes_exon_boundary_list:
                    tes_exon_list.append(Exon(
                        gene_ival.chrom, tes_exon_boundary_list[-1], tes_site, gene_ival.strand))
        else:
            # [TES, TSS] or [TES, intron_boundary]
            intron_junction = [x.start for x in intron_list]
            boundary_list = sorted(tss_site_list + intron_junction)
            for tes_site in tes_site_list:
                tes_exon_boundary_list = list(filter(
                    lambda x: (GVAR.MINEXON_LEN <= x - tes_site <= GVAR.MAXEXON_LEN), boundary_list))
                if tes_exon_boundary_list:
                    tes_exon_list.append(Exon(
                        gene_ival.chrom, tes_site, tes_exon_boundary_list[0], gene_ival.strand))
        for tes_exon in tes_exon_list:
            tes_exon.inherit_cov_from(gene_ival)
            tes_exon.set_tes_exon(True)
        return tes_exon_list

    @classmethod
    def enumerate_exon(cls, intron_list):
        """ Enumerate possible internal exon by intron. """
        starts = {x.start for x in intron_list}
        ends = {x.end for x in intron_list}
        
        ss = set() # splice site
        for spliced_intron in filter(lambda x: x.spliced, intron_list):
            ss.add(spliced_intron.start)
            ss.add(spliced_intron.end)
        ss = [0] + sorted(ss) + [sys.maxsize] 

        exon_set = set()
        for i in range(len(ss)-1):
            start = ss[i]
            end = ss[i+1]
            tmp_starts = [p for p in starts if start <= p <= end]
            tmp_ends = [p for p in ends if start <= p <= end]
            exon_set |= {
                (x, y) for x in tmp_ends for y in tmp_starts
                if (GVAR.MINEXON_LEN <= (y - x) <= GVAR.MAXEXON_LEN)}
   
        unused_EI_junction = starts.difference({x[1] for x in exon_set})
        unused_IE_junction = ends.difference({x[0] for x in exon_set})
        return sorted(exon_set), unused_EI_junction, unused_IE_junction

    @classmethod
    def has_gap(cls, exon):
        """ Judge if a region have too many gap to by a exon. """
        gaplen = np.sum(exon.cov.sig < GVAR.NGS_LIMIT_READ_COVER_PER_BASE, 1).min()
        return gaplen <= GVAR.MAX_GAP_INEXON

    @classmethod
    def identify_internal_exon(cls, intron_list, gene_ival):
        """ Identify internal exon by intron. """
        exon_loc, unused_EI_junction, unused_IE_junction = cls.enumerate_exon(intron_list)
        exon_list = [Exon(gene_ival.chrom, x[0], x[1], gene_ival.strand) for x in exon_loc]
        
        for exon in exon_list:
            exon.inherit_cov_from(gene_ival)
        exons = list(filter(cls.has_gap, exon_list))
        
        return exons, unused_EI_junction, unused_IE_junction

    @classmethod
    def fetch_neighbor_seq(cls, start, end, ival, genome):
        """Extraction sequence near a interval"""
        assert start < ival.start < ival.end < end
        return genome.find_sequence(ival.chrom, start, ival.start) + genome.find_sequence(ival.chrom, ival.end, end)

    rev_code = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    @classmethod
    def fetch_splice_site(cls, intron, genome):
        """Extraction splice site"""
        first_ss = genome.find_sequence(intron.chrom, intron.start, intron.start + 2).upper()
        second_ss = genome.find_sequence(intron.chrom, intron.end - 2, intron.end).upper()
        if intron.strand == "+":
            return (first_ss, second_ss)
        else:
            first_ss = "".join(map(lambda x: cls.rev_code[x], first_ss))[::-1]
            second_ss = "".join(map(lambda x: cls.rev_code[x], second_ss))[::-1]
            return (second_ss, first_ss)

    SS = {("GT", "AG"): 3, ("GC", "AG"): 2, ("AT", "AC"): 1}
    @classmethod
    def ss2pri(cls, ss):
        """Compute splice site priority"""
        if ss in cls.SS.keys():
            return cls.SS[ss]
        else:
            return 0

    @classmethod
    def blocks2cigar(cls, blocks):
        """
        Produce cigar tuple by exon blocks

        Args:
            blocks (list): list of exon block

        Returns:
            tuple: cigar tuple

        """
        cigar = list()
        for indx, (block_start, block_end) in enumerate(blocks):
            cigar.append((0, block_end - block_start))
            if indx < (len(blocks) - 1):
                cigar.append((3, blocks[indx+1][0] - block_end))
        cigar = tuple(cigar)
        return cigar

    @classmethod
    def adjust_intron_position(cls, tgs_read, intron):
        """Adjust intron position if possible"""
        if tgs_read.reference_name != intron.chrom:
            return None
        if tgs_read.is_reverse ^ (intron.strand == "-"):
            return None
        if (tgs_read.reference_start > intron.end) or (tgs_read.reference_end < intron.start):
            return None
        new_blocks = list()
        blocks = tgs_read.blocks
        for indx, block in enumerate(blocks):
            this_block = list(block)
            if indx == (len(blocks) - 1):
                if this_block[0] < intron.start < this_block[1]:
                    this_block = [this_block[0], intron.start]
            else:
                if this_block[0] < intron.start < blocks[indx + 1][0]:
                    this_block = [this_block[0], intron.start]

            if indx == 0:
                if this_block[0] < intron.end < this_block[1]:
                    this_block = [intron.end, this_block[1]]
            else:
                if blocks[indx - 1][1] < intron.end < this_block[1]:
                    this_block = [intron.end, this_block[1]]
            new_blocks.append(tuple(this_block))
        tgs_read.cigar = cls.blocks2cigar(new_blocks)
        return None

    @classmethod
    def fix_mapping_error(cls, ngs_intron_set, tgs_intron_set, genome, tgs_reads, ngs_bam_list):
        """Unified NGS and TGS intron and adjust tgs reads' intron position"""
        if genome is None:
            return list(ngs_intron_set) + list(tgs_intron_set.difference(ngs_intron_set))
        ngs_intron_list = list(ngs_intron_set)
        tgs_intron_list = list(tgs_intron_set)
        rm_ngs_indx_set = set()
        rm_tgs_indx_set = set()

        for ngs_indx, ngs_intron in enumerate(ngs_intron_list):
            matched_tgs_indx_list = list()
            same_intron = False
            ngs_neibour_start = ngs_intron.start - GVAR.MAPPING_SLIDING_ERROR
            ngs_neibour_end = ngs_intron.end + GVAR.MAPPING_SLIDING_ERROR
            ngs_neibour_seq = cls.fetch_neighbor_seq(ngs_neibour_start, ngs_neibour_end, ngs_intron, genome)
            overlap_same_len_tgs_indx = list()
            # Only select the TGS intron which have adjacent position and same local mRNA sequence with NGS intron
            for tgs_indx, tgs_intron in enumerate(tgs_intron_list):
                if tgs_intron.chrom != ngs_intron.chrom:
                    continue
                if tgs_intron.strand != ngs_intron.strand:
                    continue
                if max(tgs_intron.start, ngs_intron.start) >= min(tgs_intron.end, ngs_intron.end):
                    continue
                if (tgs_intron.end - tgs_intron.start) != (ngs_intron.end - ngs_intron.start):
                    continue
                if (tgs_intron.start <= ngs_neibour_start) or (tgs_intron.end >= ngs_neibour_end):
                    continue
                if tgs_intron.start == ngs_intron.start:
                    same_intron = True
                    matched_tgs_indx_list = list()
                    break
                tgs_neibour_seq = cls.fetch_neighbor_seq(ngs_neibour_start, ngs_neibour_end, tgs_intron, genome)
                if ngs_neibour_seq != tgs_neibour_seq:
                    overlap_same_len_tgs_indx.append(tgs_indx)
                    continue
                matched_tgs_indx_list.append(tgs_indx)

            if same_intron:
                continue
            ngs_reads = list()
            for bam in ngs_bam_list:
                ngs_reads += list(bam.smart_fetch(ngs_intron.chrom, ngs_intron.start, ngs_intron.end))

            fix = False
            tmp_ngs_reads = list(filter(
                lambda read: AlignReadMethod.has_intron(read, ngs_intron.start, ngs_intron.end), ngs_reads))

            for tgs_indx in overlap_same_len_tgs_indx:
                tgs_intron = tgs_intron_list[tgs_indx]
                tmp_tgs_reads = list(filter(
                lambda read: AlignReadMethod.has_intron(read, tgs_intron.start, tgs_intron.end), tgs_reads))
                left_genome_seq = genome.find_sequence(
                    ngs_intron.chrom,
                    min(ngs_intron.start, tgs_intron.start),
                    max(ngs_intron.start, tgs_intron.start)).upper()
                right_genome_seq = genome.find_sequence(
                    ngs_intron.chrom,
                    min(ngs_intron.end, tgs_intron.end),
                    max(ngs_intron.end, tgs_intron.end)).upper()
                ngs_cnt = 0
                tgs_cnt = 0
                if ngs_intron.start < tgs_intron.start:
                    for ngs_read in tmp_ngs_reads:
                        ngs_seq = AlignReadMethod.fetch_seq_by_ref_loc(ngs_read, ngs_intron.end, tgs_intron.end)
                        if (ngs_seq == left_genome_seq) and (ngs_seq != right_genome_seq):
                            ngs_cnt += 1
                    for tgs_read in tmp_tgs_reads:
                        tgs_seq = AlignReadMethod.fetch_seq_by_ref_loc(tgs_read, ngs_intron.start, tgs_intron.start)
                        if (tgs_seq != left_genome_seq) and (tgs_seq == right_genome_seq):
                            tgs_cnt += 1
                else:
                    for ngs_read in tmp_ngs_reads:
                        ngs_seq = AlignReadMethod.fetch_seq_by_ref_loc(ngs_read, tgs_intron.start, ngs_intron.start)
                        if (ngs_seq != left_genome_seq) and (ngs_seq == right_genome_seq):
                            ngs_cnt += 1
                    for tgs_read in tmp_tgs_reads:
                        tgs_seq = AlignReadMethod.fetch_seq_by_ref_loc(tgs_read, tgs_intron.end, ngs_intron.end)
                        if (tgs_seq == left_genome_seq) and (tgs_seq != right_genome_seq):
                            tgs_cnt += 1
                if (ngs_cnt > 0 and ngs_cnt/len(tmp_ngs_reads)>0.5) and (tgs_cnt == 0 or tgs_cnt/len(tmp_tgs_reads)<0.5):
                    rm_ngs_indx_set.add(ngs_indx)
                    fix = True
                if (tgs_cnt > 0 and tgs_cnt/len(tmp_ngs_reads)>0.5) and (ngs_cnt == 0 or ngs_cnt/len(tmp_tgs_reads)<0.5):
                    rm_tgs_indx_set.add(tgs_indx)
                    for read in tmp_tgs_reads:
                        cls.adjust_intron_position(read, ngs_intron)
                    fix = True

            if fix:
                continue

            # Using the highest priority of the intron
            if matched_tgs_indx_list:
                ngs_ss = cls.fetch_splice_site(ngs_intron, genome)
                tgs_ss_list = list(map(
                    lambda x: cls.fetch_splice_site(tgs_intron_list[x], genome), matched_tgs_indx_list))

                ngs_pri = cls.ss2pri(ngs_ss)
                tgs_pri_list = list(map(cls.ss2pri, tgs_ss_list))

                if max(tgs_pri_list) >= ngs_pri:
                    # Use TGS intron
                    rm_ngs_indx_set.add(ngs_indx)
                    highest_pri = max(tgs_pri_list)
                    for indx, pri in enumerate(tgs_pri_list):
                        if pri < highest_pri:
                            rm_tgs_indx_set.add(matched_tgs_indx_list[indx])
                    highest_pri_introns = [tgs_intron_list[matched_tgs_indx_list[indx]]
                                           for indx in range(len(tgs_pri_list)) if tgs_pri_list[indx] == highest_pri]
                    if len(highest_pri_introns) == 1:
                        highest_pri_introns[0].set_spliced(ngs_intron.spliced)
                        highest_pri_introns[0].set_spliced_readnum(ngs_intron.spliced_readnum)

                else:
                    # Use NGS intron
                    rm_tgs_indx_set |= set(matched_tgs_indx_list)
                    for tgs_read in tgs_reads:
                        cls.adjust_intron_position(tgs_read, ngs_intron)

        clean_ngs_intron_set = set()
        for indx, ngs_intron in enumerate(ngs_intron_list):
            if indx not in rm_ngs_indx_set:
                clean_ngs_intron_set.add(ngs_intron)
        clean_tgs_intron_set = set()
        for indx, tgs_intron in enumerate(tgs_intron_list):
            if indx not in rm_tgs_indx_set:
                clean_tgs_intron_set.add(tgs_intron)
        return list(clean_ngs_intron_set) + list(clean_tgs_intron_set.difference(clean_ngs_intron_set))


class NgsElementDiscover(ElementDiscover):
    """ Method class to identify transcript elements by ngs data. """
    @classmethod
    def identify_intron(cls, gene_ival, ngs_bams):
        """ Predict intron by NGS data. """
        ngs_junction_graph = JunctionGraphNgs(gene_ival)
        intron_list = ngs_junction_graph.identify_intron(ngs_bams)
        return intron_list


class TgsElementDiscover(ElementDiscover):
    """ Method class to identify transcript elements by tgs data. """
    @classmethod
    def identify_intron(cls, gene_ival, tgs_read_list):
        """ Predict intron by TGS data. """
        tgs_junction_graph = JunctionGraphTgs(gene_ival)
        intron_list = tgs_junction_graph.identify_intron(tgs_read_list)
        return intron_list


class Intron(Interval):
    """ A Structure to store intron data """
    def __init__(self, chrom, start, end, strand, spliced_readnum=0):
        super().__init__(chrom, start, end, strand)
        self.psi = None
        self.spliced = False
        self.spliced_readnum = spliced_readnum

    def __str__(self):
        return "Intron: {0}:{1}-{2} {3}".format(self.chrom, self.start, self.end, self.strand)

    __repr__ = __str__

    def is_spliced(self, gene_ival, pircutoff=GVAR.SPLICED_INTRON_PIR_CUTOFF):
        """
        Compute if this intron is almost spliced

        Args:
            gene_ival (Interval): Gene interval
            pircutoff: PIR cutoff

        Returns:
            bool

        """
        self.inherit_cov_from(gene_ival)
        ss1_unspliced_num = self.cov[self.start].sum()
        ss2_unspliced_num = self.cov[self.end - 1].sum()
        psi = self.spliced_readnum / (
            self.spliced_readnum + np.mean([ss1_unspliced_num, ss2_unspliced_num]))
        self.psi = psi
        if 1 - psi < pircutoff:
            self.spliced = True

    def set_spliced(self, spliced):
        """ Set splice information """
        self.spliced = spliced

    def set_spliced_readnum(self, spliced_readnum):
        """ Set spliced read number """
        self.spliced_readnum = spliced_readnum


class JunctionGraph(object):
    """ Use junction reads to build junction graph to identify intron. """
    edge_degree_cutoff = 0
    node_degree_cutoff = 0
    edge_degree_ratio = GVAR.JUNCTION_BOUNDARY_RATIO

    def __init__(self, ival):
        self.ival = ival
        self.junction_list = list()
        self.junction_graph = networkx.MultiGraph()

    @classmethod
    def _find_junction(cls, read, check_len=GVAR.JUNCTION_READ_ENDEXON_LEN):
        """ Find junction in read. """
        junction_list = list()
        chrom = read.reference_name
        first_type = "EI"
        second_type = "IE"
        if read.is_reverse:
            strand = "-"
        else:
            strand = "+"

        start = read.reference_start
        block_start = start
        cigar = read.cigartuples
        cigar_num = len(cigar)

        #  Fetch junction and check if junction two ends are valid
        for indx, (block_type, block_len) in enumerate(cigar):
            if (block_type == 3) and (indx not in (0, cigar_num - 1)):
                # Check left side
                len_checked = 0
                is_junction = True
                for check_indx in range(indx-1, -1, -1):
                    if len_checked >= check_len:
                        break
                    check_block_type, check_block_len = cigar[check_indx]
                    if check_block_type not in (0, 3):
                        is_junction = False
                        break
                    if check_block_type == 0:
                        len_checked += check_block_len
                if not is_junction:
                    continue

                # Check right side
                len_checked = 0
                for check_indx in range(indx + 1, cigar_num):
                    if len_checked >= check_len:
                        break
                    check_block_type, check_block_len = cigar[check_indx]
                    if check_block_type not in (0, 3):
                        is_junction = False
                        break
                    if check_block_type == 0:
                        len_checked += check_block_len

                if is_junction:
                    junction_list.append((
                        (chrom, block_start, strand, first_type),
                        (chrom, block_start + block_len, strand, second_type)))

            if block_type in (0, 2, 3):  # Match, Deletion, Gap
                block_start += block_len

        return junction_list

    def _build(self, bam_list):
        raise NotImplementedError()

    def _filter(self):
        """ Filter nodes and edges in junction graph. """
        # Remove junction with low expression level.
        edge_list = self.junction_graph.edges()
        edges = set(edge_list)
        for edge in edges:
            if self.junction_graph.number_of_edges(*edge) < self.edge_degree_cutoff:
                self.junction_graph.remove_edges_from([edge] * self.junction_graph.number_of_edges(*edge))

        # Remove junction boundary with low expression level.
        remove_node_list = list()
        for node in self.junction_graph.nodes():
            if self.junction_graph.degree(node) < self.node_degree_cutoff:
                remove_node_list.append(node)

        while remove_node_list:
            self.junction_graph.remove_nodes_from(remove_node_list)
            remove_node_list = list()
            for node in self.junction_graph.nodes():
                if self.junction_graph.degree(node) < self.node_degree_cutoff:
                    remove_node_list.append(node)

        # Remove junction with low expression fraction
        tv = False
        for node in self.junction_graph.nodes():
            neighbors_edge_list = self.junction_graph.edges(node)
            all_neighbors_edge_count = len(neighbors_edge_list)
            neighbors_edge_set = set(neighbors_edge_list)
            for neighbors_edge in neighbors_edge_set:
                neighbors_edge_count = self.junction_graph.number_of_edges(*neighbors_edge)
                if (neighbors_edge_count / all_neighbors_edge_count) < self.edge_degree_ratio:
                    self.junction_graph.remove_edges_from([neighbors_edge] * neighbors_edge_count)
                    tv = True

        # When junction changed, we need to compute junction boundary expression level again.
        if tv:
            self._filter()

    def _graph2intron(self):
        """ Make introns by using junction graph information. """
        junction_list = list(self.junction_graph.edges())
        intron_list = [
            Intron(self.ival.chrom, min(x[0][1], x[1][1]), max(x[0][1], x[1][1]), 
                self.ival.strand, junction_list.count(x))
            for x in set(junction_list)]
        return intron_list


class JunctionGraphNgs(JunctionGraph):
    """ Use ngs junction reads to build junction graph to identify intron. """
    edge_degree_cutoff = GVAR.MIN_JUNCTION_READ_NGS
    node_degree_cutoff = GVAR.MIN_JUNCTION_BOUNDARY_READ_NGS

    def _is_junction_read(self, read):
        """
        Judge if a ngs read is junction read.

        Args:
            read (AlignedSegment): Read

        Returns:
            bool: if this read is junction read
        """
        try:
            if len(read.cigartuples) < 3: # has >=1 intron
                return False
            if read.is_secondary:
                return False
            return True
        except:
            return False

    def _build(self, bam_list):
        """ Build junction graph. """
        junction_read_list = list()
        for bam in bam_list:
            junction_read_list += list(filter(
                lambda x: self._is_junction_read(x), bam.fetch_reads_in_ival(self.ival)))
        if GVAR.RULE == "single_end":
            for read in junction_read_list:
                self.junction_graph.add_edges_from(self._find_junction(read))
        else:  # For pair end read, same junction only count one time.
            read_dict = dict()
            for read in junction_read_list:
                if read.query_name in read_dict.keys():
                    read_dict[read.query_name].append(read)
                else:
                    read_dict[read.query_name] = [read]
            for pair_list in read_dict.values():
                junction_set = set()
                for read in pair_list:
                    junction_set |= set(self._find_junction(read))

                self.junction_graph.add_edges_from(junction_set)
        #print(sorted(set(self.junction_graph.edges())))

    def identify_intron(self, bam_list):
        """ Reads to intron. """
        self._build(bam_list)
        self._filter()
        intron_list = self._graph2intron()
        for intron in intron_list:
            intron.is_spliced(self.ival)
        return intron_list


class JunctionGraphTgs(JunctionGraph):
    """ Use tgs reads to build junction graph to identify intron. """
    edge_degree_cutoff = GVAR.MIN_JUNCTION_READ_TGS
    node_degree_cutoff = GVAR.MIN_JUNCTION_BOUNDARY_READ_TGS

    def _build(self, tgs_read_list):
        """ Find junction in read. """
        for read in tgs_read_list:
            self.junction_graph.add_edges_from(self._find_junction(read))
        
    def _filter(self):
        """ Filter nodes and edges in junction graph. """
        if self.ival.cov is not None and self.ival.cov.sig.sum(0).max() > GVAR.MIN_NGS_COV_FOR_ALL_JUNCTION_BY_TGS:
            super()._filter()

    def identify_intron(self, tgs_read_list):
        """ Reads to intron. """
        self._build(tgs_read_list)
        self._filter()
        return self._graph2intron()


def identify_element(
        chrom, start, end, ngs_bam_list, tgs_bam_list, ext_tss_site_list, ext_tes_site_list,
        ann=None, f_genome=None, paraclu_path=None):
    """
    Identify transcript elements by TGS and NGS data.

    Args:
        chrom (str): Chrom
        start (int): Start position
        end (int): End position
        ngs_bam_list (list): List of NGS bam files
        tgs_bam_list (list): List of TGS bam files
        ext_tss_site_list (list): List of annotated TSS site (chrom, loc, strand)
        ext_tes_site_list (list): List of annotated TES site (chrom, loc, strand)
        ann (SeqFile): Annotation object
        f_genome (str): Genome file
        paraclu_path (str): Path of paraclu

    Returns:
        list: List of gene object
    """
    gene_linkage = GeneLinkage(chrom, start, end, ngs_bam_list, tgs_bam_list, ann)
    gene_list = gene_linkage.split2gene(f_genome)
    for gene in gene_list:
        tmp_ext_tss_site_list = [
            x[1] for x in ext_tss_site_list
            if (x[0] == gene.ival.chrom) and (x[2] == gene.ival.strand) and (gene.ival.start <= x[1] <= gene.ival.end)]
        tmp_ext_tes_site_list = [
            x[1] for x in ext_tes_site_list
            if (x[0] == gene.ival.chrom) and (x[2] == gene.ival.strand) and (gene.ival.start <= x[1] <= gene.ival.end)]
        gene.identify_element(tmp_ext_tss_site_list, tmp_ext_tes_site_list, paraclu_path)
    return gene_list
