import numpy as np
from bx.intervals.cluster import ClusterTree
from igia.coverage import CountReadsPerBinWithIntron
from igia import GVAR
import logging
_logger = logging.getLogger(__name__)


class Linkage(object):
    """List of blocks by chrom"""
    def __init__(self):
        self.region = dict()

    def __str__(self):
        info = ["{chrom}:{length}".format(
            chrom=key, length=len(self.region[key])) for key in sorted(self.region)]
        return "\n".join(info)

    __repr__ = __str__

    def add_linkage(self, linkage):
        """
        Merge the information from another linkage
        Args:
            linkage (Linkage): Another linkage

        """
        for chrom, linkage_list in linkage.region.items():
            if chrom in self.region.keys():
                self.add_chr_linkage(chrom, linkage_list)
            else:
                self.region[chrom] = linkage_list

    def add_chr_linkage(self, chrom, region_list):
        """
        Merge the information from another linkage for certain chrom
        Args:
            chrom (str): Chrom to merge
            region_list (list): List of block
        """
        if chrom not in self.region.keys():
            self.region[chrom] = list()
        linkage_tree = ClusterTree(0, 0)
        for (start, end) in region_list + self.region[chrom]:
            linkage_tree.insert(start, end, 0)
        self.region[chrom] = self.__get_linkage_list(linkage_tree)

    def __get_linkage_list(self, tree):
        return [(x[0], x[1]) for x in tree.getregions()]

    def getregions(self, chrom):
        """
        Find linkage for a certain chrom.

        Args:
            chrom (str): Chrom

        Returns:
            list: List of block region
        """
        if chrom in self.region.keys():
            return self.region[chrom]
        else:
            return list()

    def iterlinkage(self):
        """ Iter linkage in whole genome """
        for chrom, linkage_list in self.region.items():
            for start, end in linkage_list:
                yield (chrom, start, end)


def find_linkage_worker(chrom_size, seq_obj):
    '''Find RNA linkage from NGS and TGS data for one BAM file'''
    linkage = Linkage()
    chrom, size = chrom_size
    bincutoff = seq_obj.poisbg(GVAR.LINKAGE_SEARCH_WINDOW_SIZE, GVAR.ALPHA)
    cr = CountReadsPerBinWithIntron(
        [seq_obj.bamfile], binLength=GVAR.LINKAGE_SEARCH_WINDOW_SIZE,
        stepSize=GVAR.LINKAGE_SEARCH_WINDOW_SIZE, samFlag_exclude=3852,
        ignoreDuplicates=True, minMappingQuality=20)
    bin_count = cr.count_reads_in_region_with_intron(chrom, 1, size)[0]
    bin_count = bin_count.reshape(bin_count.size)
    rich_indx = np.where(bin_count > bincutoff)[0]
    start_loc = 1 + rich_indx * GVAR.LINKAGE_SEARCH_WINDOW_SIZE
    end_loc = start_loc + GVAR.LINKAGE_SEARCH_WINDOW_SIZE
    rich_bins = list(zip(start_loc, end_loc))
    linkage.add_chr_linkage(chrom, rich_bins)

    return linkage


def find_linkage(bamlist):
    '''Find RNA linkage from NGS and TGS data'''
    assert len(bamlist) > 0
    linkage = Linkage()
    chrom_size_list = bamlist[0].chromsize()
    for chrom_size in chrom_size_list:
        for seq_obj in bamlist:
            _logger.debug("Start scan linkage in {0}:{1}".format(chrom_size[0], seq_obj.bamfile))
            sub_linkage = find_linkage_worker(chrom_size, seq_obj)
            _logger.debug("Finish scan linkage in {0}:{1}".format(chrom_size[0], seq_obj.bamfile))
            linkage.add_linkage(sub_linkage)

    return linkage



