import sys
import numpy as np
import pybedtools
import pysam
import scipy.stats
import os
import re
import time
import copy
from igia import GVAR


class SeqFile(object):
    '''
    Bam format file
    '''
    def __init__(self, filename, file_type):
        """
        Create a bam object

        Args:
            filename (str): Bam format file, str
            file_type (str): File type in [NGS, TGS, ANN]
        """
        assert file_type in ('NGS', 'TGS', "ANN")
        assert os.path.isfile(filename), filename
        self.bamfile = filename
        self.bam = pysam.AlignmentFile(self.bamfile, "rb")
        self.type = file_type

    def close(self):
        try:
            self.bam.close()
        except:
            pass

    def __copy__(self):
        newone = SeqFile(self.bamfile, self.type)
        return newone

    def __deepcopy__(self, memodict={}):
        return self.__copy__()

    def copy(self):
        return self.__copy__()

    def __str__(self):
        return "seq_file object\nfile dir:{filename}\nfile type:{file_type}".format(
            filename=self.bamfile, file_type=self.type)

    __repr__ = __str__

    def chromsize(self):
        """
        Compute chrom size.

        Returns:
            list: (chrom, size)
        """
        return list(zip(self.bam.references, self.bam.lengths))

    def chrom2size(self):
        """
        Compute chrom size.

        Returns:
            dict: {chrom: size}
        """
        return dict(zip(self.bam.references, self.bam.lengths))

    def genomesize(self):
        """ Compute genome size. """
        return sum(self.bam.lengths)

    def mapped_number(self):
        """ Compute read number. """
        return self.bam.mapped

    def readlen(self, top=100):
        """ Compute read length. """
        read_len_list = list()
        reads = self.bam.fetch()
        for i in range(top):
            try:
                read = next(reads)
                read_len_list.append(read.query_length)
            except:
                break
        assert len(read_len_list) > 0
        return scipy.stats.mode(read_len_list).mode[0]

    def poisbg(self, binsize, alpha=0.01):
        """
        Compute NGS signal background by poisson distribution.

        Args:
            binsize (int): Window size
            alpha (float): Confidence coefficient

        Returns:
            float: Background cutoff
        """

        cutoff = 0
        if self.type == 'NGS':
            cutoff = poiscut(self.mapped_number(), self.genomesize(), binsize, alpha)
        return cutoff

    def filter_clean_reads(self, iter):
        """ Filter read by rules. """
        iter = filter(lambda x: not x.is_secondary, iter)
        iter = filter(lambda x: x.reference_end is not None, iter)
        iter = filter(lambda x: x.reference_start is not None, iter)
        iter = filter(lambda x: (x.reference_end - x.reference_start - x.query_length) < GVAR.MAX_INTRON_LEN, iter)
        if (GVAR.RULE != "single_end") and (self.type == "NGS"):
            iter = filter(lambda x: x.is_reverse ^ x.mate_is_reverse, iter)
        return iter

    def count(self, chrom, start, end):  
        """ Compute read number in ival. """
        num = self.bam.count(chrom, start, end)
        return num

    def smart_fetch(self, chrom, start, end):  
        """ Fetch reads in ival. """
        reads_num = self.count(chrom, start, end)
        try:
            read_ram_size = sys.getsizeof(next(self.bam.fetch(chrom, start, end)))
        except:
            # Have no read
            read_ram_size = 0
        if reads_num * read_ram_size < GVAR.MAX_RAM_PER_ival:
            reads = list(map(self.pretreat, self.bam.fetch(chrom, start, end)))
            # reads = filter(lambda x: x.reference_start >= start and x.reference_end <= end, reads)
            return reads
        else:
            reads = (self.pretreat(x) for x in self.bam.fetch(chrom, start, end))
            # reads = filter(lambda x: x.reference_start >= start and x.reference_end <= end, read_generator)
            return reads

    def fetch_reads_in_ival(self, ival, skip_boundary_span=True): 
        """ Fetch reads in interval. """
        chrom = ival.chrom
        start = ival.start
        end = ival.end
        strand = ival.strand
        ival_read_iter = self.smart_fetch(chrom, start, end)
        if skip_boundary_span:
            ival_read_iter = filter(
                lambda x: x.reference_start >= start and x.reference_end <= end, ival_read_iter)

        clean_ival_read_iter = self.filter_clean_reads(ival_read_iter)
        if strand == "+":
            clean_ival_read_iter = filter(lambda x: not x.is_reverse, clean_ival_read_iter)
        elif strand == "-":
            clean_ival_read_iter = filter(lambda x: x.is_reverse, clean_ival_read_iter)

        return list(clean_ival_read_iter)

    def pretreat(self, read):
        """
        Read pretreatment

        Args:
            read (AlignedSegment): read to pretreatment

        Returns:
            AlignedSegment: read
        """
        if self.type == "NGS":
            if GVAR.RULE == "1++,1--,2+-,2-+":
                if read.is_read2:
                    read.is_reverse = not read.is_reverse
                    read.mate_is_reverse = not read.mate_is_reverse
            elif GVAR.RULE == "1+-,1-+,2++,2--":
                if read.is_read1:
                    read.is_reverse = not read.is_reverse
                    read.mate_is_reverse = not read.mate_is_reverse

        elif self.type in ("TGS", "ANN"):
            try:
                FL_num = int(re.findall(re.compile("[\+\w]+/f(\d+)p\d+/\d+"), read.query_name)[0])
                read.setTag("FL", FL_num)
            except:
                read.setTag("FL", 1)

        read.setTag("TP", self.type)
        return read


class GenomeFile(object):
    '''
    FASTA format file
    '''
    def __init__(self, file_dir):
        """
        Create a genome object

        Args:
            file_dir (str): Fasta format genome file
        """
        self.file_dir = file_dir
        self.file = pysam.FastaFile(self.file_dir)

    def find_sequence(self, chrom, start, end):
        """
        Fetch sequence from genome

        Args:
            chrom (str): Chromosome to fetch
            start (str): Start position
            end (str): End position

        Returns:
            str: Sequence
        """
        return self.file.fetch(chrom, start, end)


class Coverage(object):
    """ Wig like coverage format. """
    def __init__(self, ival):
        self.ival = ival
        self.sig = np.array([], dtype=np.int)

    def __getitem__(self, pos):
        assert (pos >= self.ival.start) and (pos <= self.ival.end)
        return self.sig[:, pos - self.ival.start]

    def __len__(self):
        """left and right are both inclusive"""
        return self.ival.end - self.ival.start + 1

    def build(self, ngs_list):
        """ Make coverage by NGS data. """
        sig = np.zeros((len(ngs_list), self.ival.end - self.ival.start + 1), np.int)
        for indx, ngs in enumerate(ngs_list):
            for read in ngs.fetch_reads_in_ival(self.ival, skip_boundary_span=False):
                pos = np.array([x for x in read.positions if self.ival.start <= x <= self.ival.end], np.int)
                idx = pos - self.ival.start
                sig[indx, idx] += 1
        self.sig = sig

    def slice(self, ival):
        """ Fetch subregion coverage. """
        assert ival.chrom == self.ival.chrom, (str(self.ival), str(ival))
        assert ival.start >= self.ival.start, (str(self.ival), str(ival))
        assert ival.end <= self.ival.end, (str(self.ival), str(ival))
        assert ival.strand == self.ival.strand, (str(self.ival), str(ival))
        new_cov = Coverage(ival)
        new_cov.sig = self.sig[:, (ival.start - self.ival.start):(ival.end - self.ival.start + 1)]
        return new_cov


class Interval(object):
    """ Interval information and operation. """
    def __init__(self, chrom, start, end, strand):
        assert strand in "+-."
        assert start < end
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.fpkm = None
        self.cov = None

    def __str__(self):
        return "Interval: {0}:{1}-{2} {3}".format(self.chrom, self.start, self.end, self.strand)

    __repr__ = __str__

    def __hash__(self):
        return hash((self.chrom, self.start, self.end, self.strand))

    def __eq__(self, other):
        return (self.chrom, self.start, self.end, self.strand) == (other.chrom, other.start, other.end, other.strand)

    def __len__(self):
        return self.end - self.start

    def build_cov(self, ngs_list):
        """ Make coverage by NGS data. """
        self.cov = Coverage(self)
        self.cov.build(ngs_list)

    def slice_cov(self, ival):
        """ Fetch subregion coverage. """
        assert self.cov is not None
        return self.cov.slice(ival)

    def inherit_cov_from(self, ival):
        """ Slice coverage signal from another interval """
        self.cov = ival.slice_cov(self)

    @classmethod
    def cov2fpkm(cls, cnt, ivlen, readnum, alpha):
        """Formular from GRIT"""
        return 1e9 / ivlen * scipy.stats.beta.ppf([alpha, 0.5, 1-alpha], cnt, readnum-cnt)

    def compute_fpkm(self, readlen_list, readnum_list, alpha=GVAR.ALPHA):
        """Compute FPKM from coverage"""
        assert self.cov is not None
        assert len(readlen_list) == len(readnum_list) == self.cov.sig.shape[0]
        fpkm_list = list()
        for indx in range(self.cov.sig.shape[0]):
            cnt = self.cov.sig[indx, :].sum() / readlen_list[indx] + 1e-6
            fpkm = self.cov2fpkm(cnt, len(self), readnum_list[indx], alpha)
            fpkm_list.append(list(fpkm))
        fpkm = np.array(fpkm_list)
        self.fpkm = fpkm

    def write2bed6(self, name, f):
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
            self.chrom, self.start, self.end, name, 0, self.strand))
        f.flush()


def load_seqinfo(ngs_bams):
    '''
    NGS_FILE_INFO_list: (file_name, read_num, read_len)
    '''
    file_names = [x.bamfile for x in ngs_bams]
    read_nums = [x.mapped_number() for x in ngs_bams]
    read_len = [x.readlen() for x in ngs_bams]
    GVAR.NGS_FILE_INFO_list = list(zip(file_names, read_nums, read_len))
    chrom_size_list = [x.chromsize() for x in ngs_bams]
    GVAR.CHROM_SIZE = chrom_size_list[0]
    assert chrom_size_list.count(GVAR.CHROM_SIZE) == len(ngs_bams)
    GVAR.GENOME_SIZE = sum([chrom_size[1] for chrom_size in GVAR.CHROM_SIZE])
    GVAR.NGS_LIMIT_READ_COVER_PER_BASE = poiscut(
        np.array(read_nums), GVAR.GENOME_SIZE, 2 * np.array(read_len), GVAR.ALPHA).reshape([len(ngs_bams), 1])


def poiscut(totread, genomesize, binsize, alpha=0.01):
    """Compute the cutoff of count for given binsize with totread number"""
    return scipy.stats.poisson.ppf(1-alpha, totread/(genomesize/binsize))


def load_txs(file):
    """
    Load tsv format TXS annotation
    Args:
        file (str): TXS file

    Returns:
        list: TXS list
    """
    if not file:
        return list()
    with open(file, "r") as f:
        txs_text = f.readlines()
        txs_data = list(map(lambda x: x.rstrip("\n").split("\t"), txs_text))
        txs_list = list(map(lambda x: (x[0], int(x[1]), x[2]), txs_data))
    return txs_list


def make_read(ref_id, name, start, cigar):
    """ Build sam format read by position """
    header = pysam.AlignmentHeader()
    header_dict = { 'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1e6, 'SN': 'chr1'}, {'LN': 1e6, 'SN': 'chr2'}] }
    header = header.from_dict(header_dict)
    tgs_read = pysam.AlignedSegment(header)
    tgs_read.query_name = name
    tgs_read.reference_id = ref_id
    tgs_read.reference_start = start
    tgs_read.cigar = cigar
    tgs_read.setTag('FL', 1)
    return tgs_read


def bed2bam(f_bed, f_chromsize, outdir):
    """
    Transfer Bed format file to Bam format file

    Args:
        f_bed (str): Bed format file
        f_chromsize (str): Chrom size file
        outdir (str): Outout folder

    Returns:
        str: Bam format file
    """
    ann = pybedtools.BedTool(f_bed)
    ann_bam = ann.to_bam(g=f_chromsize, bed12=True)
    name = os.path.basename(f_bed)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    f_bam = os.path.join(outdir, "{0}.bam".format(name))
    ann_bam.saveas(f_bam)
    sorted_bam = os.path.join(outdir, "{0}.sorted.bam".format(name))
    pysam.sort("-o", sorted_bam, f_bam)
    pysam.index(sorted_bam)
    return sorted_bam


def iterline(file):
    """Read file with cache"""
    with open(file, "r") as f:
        lines = f.readlines(1000000)
        while lines:
            for line in lines:
                if line == "\n":
                    continue
                elif line[0] == "#":
                    continue
                else:
                    yield line
            lines = f.readlines(1000000)


class Bed12(object):
    """ Bed12 format """
    def __init__(self, line):
        line = line.replace("\n", "")
        line_data = line.split("\t")
        self.chrom, chromStart, chromEnd, self.name, self.score, self.strand, thickStart, thickEnd, self.itemRgb, blockCount, blockSizes, blockStarts = line_data
        self.chromStart = int(chromStart)
        self.chromEnd = int(chromEnd)
        self.transcript_len = self.chromEnd - self.chromStart
        self.thickStart = int(thickStart)
        self.thickEnd = int(thickEnd)
        self.blockCount = int(blockCount)
        blockSizes = blockSizes.rstrip(",").split(",")
        self.blockSizes = list(map(int, blockSizes))
        blockStarts = blockStarts.rstrip(",").split(",")
        self.blockStarts = list(map(int, blockStarts))
        self.blockEnd = [self.blockSizes[i] + self.blockStarts[i] for i in range(self.blockCount)]
        self.cum_block_size = np.cumsum(self.blockSizes)

    def __str__(self):
        blockSizes = ",".join(list(map(str, self.blockSizes)))
        blockStarts = ",".join(list(map(str, self.blockStarts)))
        text = "\t".join(list(map(str, [
            self.chrom, self.chromStart, self.chromEnd, self.name, self.score, self.strand,
            self.thickStart, self.thickEnd, self.itemRgb, self.blockCount, blockSizes, blockStarts
        ]))) + "\n"
        return text

    __repr__ = __str__

    def rel2abs(self, rel_pos):
        """Transform relative mRNA position to absolute genomic position"""
        block_indx = np.where(self.cum_block_size >= rel_pos)[0].min()
        if block_indx == 0:
            cum_pos = 0
        else:
            cum_pos = self.cum_block_size[block_indx - 1]
        return self.chromStart + self.blockStarts[block_indx] + rel_pos - cum_pos

    def iterblock(self):
        """ Iterate exon """
        for indx in range(self.blockCount):
            block_start = self.chromStart + self.blockStarts[indx]
            block_end = block_start + self.blockSizes[indx]
            yield (block_start, block_end)

    def write(self, f):
        """ Write record to file """
        f.write(self.__str__())

    def find_intron(self):
        """ Find intron from the record """
        if self.blockCount == 1:
            return list()
        blocks = list(self.iterblock())
        introns = list(zip([x[1] for x in blocks[:-1]], [x[0] for x in blocks[1:]]))
        return introns

def load_ann(f_ann, size, out_dir, type):
    """
    Load transcript annotation

    Args:
        f_ann (str): Bed12 format annotation file
        size (str): Chrom size file
        out_dir (str): Output folder to storage bam format file
        type (str): Annotation file type

    Returns:
        SeqFile: Annotation sequence object
    """
    ann = None
    if f_ann:
        if not size:
            raise ValueError("Missing chromsize file")

        ann_bam = bed2bam(f_ann, size, out_dir)
        ann = SeqFile(ann_bam, type)
    return ann

class AlignReadMethod(object):
    """ Abstract class for a read """
    @classmethod
    def ref_loc2query_loc(cls, read, start, end):
        """
        Transfer the reference position to query position

        Args:
            read (AlignedSegment): Query read
            start (int): Reference start site
            end (int): Reference end site

        Returns
            tuple: query_start, query_end
        """
        tmp_ref_start = copy.deepcopy(read.reference_start)
        tmp_query_start = 0
        query_start = None
        query_end = None

        for cigar, length in read.cigar:
            ref_block_end = tmp_ref_start + length
            if cigar == 0:
                if tmp_ref_start <= start <= ref_block_end:
                    query_start = tmp_query_start + start - tmp_ref_start
                    query_end = query_start
            if (start is None) and (tmp_ref_start > start):
                return None

            if query_start is not None:
                if cigar == 0:
                    if ref_block_end <= end:
                        if tmp_ref_start < start < ref_block_end:
                            query_end += ref_block_end - start
                        else:
                            query_end += length
                    else:
                        query_end += max(end - max(start, tmp_ref_start), 0)
                        return query_start, query_end
                elif cigar in [1, 4, 5]:
                    query_end += length

            if cigar in [0, 1, 4, 5]:
                tmp_query_start += length
            if cigar in [0, 2, 3, 6]:
                tmp_ref_start += length

    @classmethod
    def has_intron(cls, read, start, end):
        """
        Determine if there is an intron inside a region

        Args:
            read (AlignedSegment): Read to scanning
            start (int): Start position
            end (int): End position

        Returns:
            bool: has intron
        """
        blocks = read.blocks
        if len(blocks) == 1:
            return False
        for indx in range(len(blocks) - 1):
            if (blocks[indx][1] == start) and (blocks[indx+1][0] == end):
                return True
        return False

    @classmethod
    def fetch_seq_by_ref_loc(cls, read, start, end):
        """
        Fetch sequence for a read

        Args:
            read (AlignedSegment): Read
            start (int): Start position
            end (int): End position

        Returns:
            str: RNA-seq sequence
        """
        qloc = cls.ref_loc2query_loc(read, start, end)
        if qloc is not None:
            query_start, query_end = cls.ref_loc2query_loc(read, start, end)
            return read.query_sequence[query_start: query_end]
