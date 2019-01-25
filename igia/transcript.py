from igia.utils import Interval
from igia import GVAR
from igia.element import GeneLinkageFinder
import numpy as np
import networkx
import copy
from bx.intervals.cluster import ClusterTree
from collections import defaultdict


class Segment(Interval):
    """ Segment in identify transcript. """

    def __init__(self, chrom, start, end, strand):
        super().__init__(chrom, start, end, strand)
        self.tss_seg = False
        self.tes_seg = False
        self.spliced_seg = False
        self._fpkm = None

    def __str__(self):
        return "Segment: {0}:{1}-{2} {3}; tss={4}; tes={5}; spliced={6}".format(
            self.chrom, self.start, self.end, self.strand, self.tss_seg,
            self.tes_seg, self.spliced_seg)

    __repr__ = __str__

    def set_tss_seg(self, is_tss=False):
        """ Label the segment as TSS segment or not """
        self.tss_seg = is_tss

    def set_tes_seg(self, is_tes=False):
        """ Label the segment as TES segment or not """
        self.tes_seg = is_tes

    def set_spliced_seg(self, is_spliced=False):
        """ Label the segment as spliced segment or not """
        self.spliced_seg = is_spliced


class Isoform(object):
    """ Isoform record """

    def __init__(self, segment_array, trans_ival):
        self.ival = trans_ival
        self.segary = np.array(segment_array)
        self.seglen = None
        self.sidx = None  # index of the first included segment for this isoform
        self.eidx = None  # index of the last included segment
        self.starts_with_tss = None
        self.ends_with_tes = None
        self.invalid = None
        self.label = ""
        self.flnum = 0
        self.nflnum = 0

    def __str__(self):
        return "Isoform: segment array: {0}; tss={1}; tes={2}; invalid={3}".format(
            list(self.segary), self.starts_with_tss, self.ends_with_tes, self.invalid)

    __repr__ = __str__

    def __hash__(self):
        return hash(tuple(self.segary))

    def __eq__(self, other):
        return isinstance(other, Isoform) and (hash(self) == hash(other))

    def set_tag(self, tss_indxs, tes_indxs, label):
        """
        Set tag of this isoform if this isoform is full-length isoform

        Args:
            tss_indxs (list): List of TSS boundary index
            tes_indxs (list): List of TES boundary index
            label (str): tag to set
        """
        incl_indx = np.where(self.segary == 1)[0]
        self.sidx = incl_indx.min()
        self.eidx = incl_indx.max()
        self.seglen = self.eidx - self.sidx + 1

        if self.ival.strand == "+":
            self.starts_with_tss = self.sidx in tss_indxs
            self.ends_with_tes = self.eidx in tes_indxs
        else:
            self.ends_with_tes = self.sidx in tes_indxs
            self.starts_with_tss = self.eidx in tss_indxs

        self.invalid = np.where(self.segary < 0)[0].size > 0

        # Begin with TSS site, End with TES site, and No conflict with known elements
        if self.starts_with_tss and self.ends_with_tes and (not self.invalid):
            self.set_label(label)

    def set_label(self, label):
        """ Set label """
        self.label = label

    def iso2bed12(self, seglist, name):
        """
        Transfer this isoform to Bed12 format string
        Args:
            seglist (list): List of segment
            name (str): Isoform name

        Returns:
            str: Bed12 format string
        """
        assert len(seglist) == self.segary.size
        incl_indx = np.where(self.segary == 1)[0]
        assert incl_indx.size
        exon_tree = ClusterTree(0, 0)
        for indx in incl_indx:
            incl_seg_start = seglist[indx].start
            incl_seg_end = seglist[indx].end
            exon_tree.insert(incl_seg_start, incl_seg_end, 0)
        exon_list = sorted([(x[0], x[1]) for x in exon_tree.getregions()])
        chromStart = exon_list[0][0]
        chromEnd = exon_list[-1][1]
        blockStarts = [exon[0] - chromStart for exon in exon_list]
        blockSizes = [exon[1] - exon[0] for exon in exon_list]
        blockCount = len(exon_list)
        return "\t".join(map(str, [
            self.ival.chrom, chromStart, chromEnd, name, 0, self.ival.strand,
            chromStart, chromStart, "0,0,0", blockCount,
            ",".join(list(map(str, blockSizes))),
            ",".join(list(map(str, blockStarts)))
        ])) + "\n"

    def write2bed12(self, seglist, name, f):
        """ Write this isoform to a file with Bed12 format """
        f.write(self.iso2bed12(seglist, name))
        f.flush()


class TransAssembler(object):
    """ Identify transcript. """

    def __init__(self, gene, ann=None):
        assert gene.ready4assembly()
        self.ival = gene.ival
        self.tgs_read_list = gene.tgs_read_list
        self.ngs_bam_list = gene.ngs_bam_list
        self.tissue_num = len(self.ngs_bam_list)
        self.exon_list = gene.tss_exon_list + gene.tes_exon_list + gene.internal_exon_list
        # Filter introns outside gene composed of exons
        exon_start = min([exon.start for exon in self.exon_list])
        exon_end = max([exon.end for exon in self.exon_list])
        self.intron_list = [intron for intron in gene.intron_list
                            if intron.start > exon_start and intron.end < exon_end]
        self.ann = ann
        self._intron_path = dict()
        self._seglink_cnt = None
        self.segment_list = None
        self.tss_indx = None
        self.tes_indx = None
        self.intron_indx = None
        self.exon_indx = None
        self.iso_list = None
        self.compatible_matrix = None
        self.overlap_matrix = None
        self.subpath_matrix = None
        self.first_order_path = None
        self.isoF = None
        self.isoM = None
        self.isoC = None
        self.isoR = None  # Rescued isoform
        self.isoP = None
        self.isoA = list()  # Annotation isoform
        self.iso_cluster = None

    def __str__(self):
        return "TransAssembler: isoF={0} isoA={1} isoM={2} isoC={3} isoP={4}".format(
            *self.get_isonum())

    __repr__ = __str__

    def get_isonum(self):
        """ Compute full-length isoform number in this assembler """
        def numiso(li):
            return len(li) if li is not None else 0

        return list(map(numiso, [self.isoF, self.isoA, self.isoM, self.isoC, self.isoP]))

    def make_segment(self):
        """ Build segment list. """
        self.segment_list, spliced_segment_pair = TransDiscover.element2segment(
            self.intron_list, self.exon_list, self.ival)
        spliced_segment_pair.sort(key=lambda x: (x[0], -x[1]))
        atom_spliced_segment_pair = list()
        for s_indx, e_indx in spliced_segment_pair:
            if s_indx == e_indx:
                atom_spliced_segment_pair.append((s_indx, e_indx))
        self.spliced_segment_pair = atom_spliced_segment_pair

    def make_seg_indx(self):
        """ Compute the boundary of functional elements """
        assert self.segment_list is not None
        self.intron_indx = TransDiscover.build_element_seg_indx(self.intron_list, self.segment_list)
        self.exon_indx = TransDiscover.build_element_seg_indx(self.exon_list, self.segment_list)
        self.tss_indx = TransDiscover.build_tss_seg_indx(self.segment_list)
        self.tes_indx = TransDiscover.build_tes_seg_indx(self.segment_list)

    def rescude_isoforms(self, invalid_iso_list):
        """ Rescude the isoforms with errors """
        assert self.isoF is not None
        isoR = list()
        if not self.isoF:
            self.isoR = list()
            return None

        for iso in invalid_iso_list:
            rescude_iso_list = TransDiscover.rescue_isoform(iso, self.isoF, self.intron_indx, self.exon_indx,
                                                            self.segment_list)
            isoR += rescude_iso_list

        for iso in isoR:
            iso.set_tag(self.tss_indx, self.tes_indx, "R")

        isoR = set(filter(lambda x: (x.label == "R") and TgsFilterRule.filter_iso(
            x, self.segment_list, self.intron_indx, self.exon_indx), isoR))
        self.isoR = sorted(isoR.difference(set(self.isoF)), key=lambda x: x.sidx)

    def init_segiso(self):
        """ Init isoform  """
        assert all(map(lambda x: x is not None, [
            self.tgs_read_list, self.segment_list, self.intron_indx, self.tss_indx, self.tes_indx]))
        self.iso_list = TgsTransDiscover.refine(
            self.tgs_read_list, self.segment_list, self.intron_indx,
            self.exon_indx, self.spliced_segment_pair, GVAR.SEGMENT_TGS_OVERLAP_FRACTION, self.ival)

        for iso in self.iso_list:
            iso.set_tag(self.tss_indx, self.tes_indx, "F")

        self.isoF = sorted(set(filter(lambda x: x.label == "F", self.iso_list)), key=lambda x: x.sidx)

        invalid_iso_list = list(filter(lambda x: x.invalid, self.iso_list))
        self.rescude_isoforms(invalid_iso_list)

        # If have annotation data, add full length annotation isoform into iso list (help to merge nfl isoform).
        if self.isoF:
            self.isoA = list()
        else:
            if self.ann is not None:
                ann_read_list = self.ann.fetch_reads_in_ival(self.ival)
                ann_iso_list = TgsTransDiscover.refine(
                    ann_read_list, self.segment_list, self.intron_indx,
                    self.exon_indx, self.spliced_segment_pair, GVAR.SEGMENT_TGS_OVERLAP_FRACTION, self.ival)
                for iso in ann_iso_list:
                    iso.set_tag(self.tss_indx, self.tes_indx, "A")
                self.isoA = sorted(set(filter(lambda x: x.label == "A", ann_iso_list)), key=lambda x: x.sidx)
                self.iso_list += self.isoA

    def build_compatible_matrix(self):
        """Build the matrix that if two isoform are compatible"""
        assert self.iso_list is not None
        self.compatible_matrix = TransDiscover.build_compatible_matrix(self.iso_list)

    def build_overlap_matrix(self):
        """Build the matrix that if two isoforms are overlapped on exons"""
        assert self.iso_list is not None
        self.overlap_matrix = TransDiscover.build_overlap_matrix(self.iso_list)

    def build_subpath_matrix(self):
        """Build the matrix [i, j] that if isoform i contains isoform j"""
        assert self.compatible_matrix is not None
        self.subpath_matrix = TransDiscover.build_subpath_matrix(self.compatible_matrix, self.iso_list)

    def build_first_order_path(self):
        assert self.iso_list is not None
        pass
        # self.first_order_path = TgsTransDiscover.build_first_order_path(self.iso_list)

    def init_assembly(self):
        """Prepare isoform information for assebmly"""
        self.build_compatible_matrix()
        self.build_overlap_matrix()
        self.build_subpath_matrix()
        self.build_first_order_path()

    def cluster_iso(self):
        """ Compute the isoform cluster """
        assert all(map(lambda x: x is not None, [
            self.isoF, self.isoR, self.isoA, self.isoM, self.isoC, self.isoP]))
        all_isos = self.isoF + self.isoR + self.isoA + self.isoM + self.isoC + self.isoP
        iso_clusters = TransDiscover.cluster_iso(all_isos)
        self.iso_cluster = iso_clusters

    def intron_path(self, start_seg_indx, end_seg_indx):
        """ Build the intron path """
        key = (start_seg_indx, end_seg_indx)
        if key not in self._intron_path.keys():
            self._intron_path[key] = TransDiscover.enum_intron_path(self.intron_indx, start_seg_indx, end_seg_indx,
                                                                    self.segment_list, self.seglink_cnt)
        return self._intron_path[key]

    @property
    def seglink_cnt(self):
        """Count reads number supported the given segment pair"""
        if self._seglink_cnt is None:
            self._seglink_cnt = [TransDiscover.compute_seg_link_cnt(f, self.segment_list) for f in self.ngs_bam_list]
        return self._seglink_cnt

    def identify_isoform(self):
        """Engine for assembly"""
        self.make_segment()
        self.make_seg_indx()
        self.init_segiso()
        self.init_assembly()

        # Collect non-full and valid isoforms
        nonF_iso_indxs = np.where([(x.label not in ("F", "A")) and (not x.invalid) for x in self.iso_list])[0]
        # Drop isoforms which are contained in other isoforms 
        nonF_nonsub_iso_indxs = list(filter(
            lambda x: not TransDiscover.is_subpath(x, self.subpath_matrix), nonF_iso_indxs))

        # Merge non-full-length isoform into cluster
        merged_iso_list = TransDiscover.merge_nfl_isoforms(
            nonF_nonsub_iso_indxs, self.compatible_matrix, self.overlap_matrix,
            self.iso_list, self.tss_indx, self.tes_indx)
        isoM = set(filter(
            lambda x: (x.label == "M") and TgsFilterRule.filter_iso(
                x, self.segment_list, self.intron_indx, self.exon_indx), merged_iso_list))
        self.isoM = sorted(isoM.difference(set(self.isoF + self.isoA + self.isoR)), key=lambda x: x.sidx)

        # Complete non-full-length merged isoforms with information from full-length isoforms
        nfl_isos = set(merged_iso_list).difference(isoM)
        fl_isos = set(self.isoF + self.isoA + self.isoR + self.isoM)
        nfl_isos = list(filter(lambda x: not TransDiscover.invalid_iso_is_subpath(x, fl_isos), nfl_isos))
        isoC = set()
        for fl in [self.isoF, self.isoA, self.isoR, self.isoM]:
            fl_segmat = np.array([x.segary for x in fl])
            completed_iso_list = list()
            for merged_nfl_iso in nfl_isos:
                completed_iso_list += TransDiscover.complete_iso_by_fl_iso(
                    merged_nfl_iso, fl_segmat, self.tss_indx, self.tes_indx)
            tmp_isoC = set(filter(lambda x: (x.label == "C") and TgsFilterRule.filter_iso(
                x, self.segment_list, self.intron_indx, self.exon_indx), completed_iso_list))
            isoC |= tmp_isoC
            nfl_isos = set(completed_iso_list).difference(tmp_isoC)

        self.isoC = sorted(isoC.difference(set(self.isoF + self.isoA + self.isoR + self.isoM)), key=lambda x: x.sidx)

        # Complete partial isoforms by enumerate feasible intron path
        partial_isos = set(nfl_isos).difference(isoC)
        fl_isos = set(self.isoF + self.isoA + self.isoR + self.isoM + self.isoC)
        partial_isos = filter(lambda x: not TransDiscover.invalid_iso_is_subpath(x, fl_isos), partial_isos)
        complete_partial_isos = list()
        for iso in partial_isos:
            complete_partial_isos += TransDiscover.complete_partial_isoform(iso, self)
        isoP = set(filter(lambda x: (x.label == "P") and TgsFilterRule.filter_iso(
            x, self.segment_list, self.intron_indx, self.exon_indx), complete_partial_isos))
        self.isoP = sorted(isoP.difference(set(self.isoF + self.isoA + self.isoR + self.isoM + self.isoC)),
                           key=lambda x: x.sidx)

        # Cluster isoforms
        self.cluster_iso()

    def write_iso2bed12(self, iso_type, gene_name, isolist, f):
        """Write isoforms in same type to a file with Bed12 format"""
        if isolist is not None:
            isoindx = 0
            for iso in isolist:
                isoindx += 1
                iso_name = "{0}_{1}_{2}".format(gene_name, iso_type, isoindx)
                iso.write2bed12(self.segment_list, iso_name, f)

    def write2bed12(self, cluster_name, f_isoF, f_isoA, f_isoR, f_isoM, f_isoC, f_isoP):
        """Write all isoforms in this gene cluster to a file with Bed12 format"""
        assert self.iso_cluster is not None
        for gid, cluster in enumerate(self.iso_cluster):
            iso_types = ["F", "A", "R", "M", "C", "P"]
            iso_indxs = dict(zip(iso_types, [0] * len(iso_types)))
            f_iso = dict(zip(iso_types, [f_isoF, f_isoA, f_isoR, f_isoM, f_isoC, f_isoP]))
            for iso in cluster:
                iso_type = iso.label
                iso_indxs[iso_type] += 1
                iso_name = "{cid}_g_{gid}_iso{iso_type}_{iso_indx}".format(
                    cid=cluster_name, gid=gid + 1, iso_type=iso_type, iso_indx=iso_indxs[iso_type])
                iso.write2bed12(self.segment_list, iso_name, f_iso[iso_type])


class TgsFilterRule(object):
    """ TGS read filter rules. """

    @classmethod
    def indx2ival(cls, indx_array):
        """Convert index array to index range"""
        ival_list = list()
        for i, indx in enumerate(indx_array):
            if i == 0:
                start_indx = indx
            if i == (indx_array.size - 1):
                end_indx = indx
                ival_list.append((start_indx, end_indx))
            elif indx_array[i + 1] - indx > 1:
                end_indx = indx
                ival_list.append((start_indx, end_indx))
                start_indx = indx_array[i + 1]

        return ival_list

    @classmethod
    def is_intron(cls, intron_indx, intron_seg_indx_list):
        """
        Judge if TGS intron is in valid intron list

        Args:
            intron_indx (tuple): (start_indx, end_indx)
            intron_seg_indx_list (list): Segment intron index for all valid introns

        Returns:
            bool: Is intron or not

        """
        if intron_indx not in intron_seg_indx_list:
            return False
        return True

    @classmethod
    def is_internal_exon(cls, exon_indx, exon_seg_indx_list):
        """
        Judge if TGS internal exon is valid

        Args:
            exon_indx (tuple): (start_indx, end_indx)
            exon_seg_indx_list (list): Segment internal exon index for all valid internal exon

        Returns:
            bool: Is internal exon or not

        """
        if exon_indx not in exon_seg_indx_list:
            return False
        return True

    @classmethod
    def is_spliced_seg_in_internal_exon(cls, exon_indx, spliced_segment_indx):
        """Judge if an internal exon contains spliced segment"""
        start, end = exon_indx
        if list(filter(lambda x: start <= x <= end, spliced_segment_indx)):
            return True
        return False

    @classmethod
    def is_skip_segment(cls, read, segment):
        """
        Check if skip-segment boundaries are not included
        """
        if (segment.start in read.positions) or ((segment.end - 1) in read.positions):
            return False
        return True

    @classmethod
    def is_incl_segment(cls, read, segment, boundary='both'):
        """
        Check if incl-segment boundaries are included
        """
        assert boundary in ('both', 'left', 'right')
        if boundary == 'both':
            return (segment.start in read.positions) and ((segment.end - 1) in read.positions)
        elif boundary == 'left':
            return (segment.start in read.positions)
        else:
            return ((segment.end - 1) in read.positions)

    @classmethod
    def filter(cls, read, isoseg, segment_list, intron_seg_indx_list, exon_seg_indx_list, spliced_segment_pair):
        """
        Label the error region for this segment array

        Args:
            read (AlignedSegment): TGS read
            isoseg (narray or list): Projected segment for this TGS read
            segment_list (list): Valid segment list
            intron_seg_indx_list (list): List of intron segment index
            exon_seg_indx_list (list): List of exon segment index
            spliced_segment_pair (list):  List of spliced segment pair

        Returns:
            narray: Isoform segment array
        """
        # Fetch segment information
        spliced_segment_indx = np.where([segment.spliced_seg for segment in segment_list])[0]

        # Compute exon (start, end) in isoform record
        isoseg = np.array(isoseg)
        iso_exon_seg_indx = np.where(isoseg == 1)[0]
        iso_exon_indx_list = cls.indx2ival(iso_exon_seg_indx)

        # Compute intron (start, end) in isoform record
        start_indx = min(iso_exon_seg_indx)
        end_indx = max(iso_exon_seg_indx)
        iso_intron_seg_indx = np.where(isoseg == 0)[0]
        iso_intron_seg_indx = np.array(list(filter(
            lambda x: (x >= start_indx) and (x <= end_indx), iso_intron_seg_indx)))
        iso_intron_indx_list = cls.indx2ival(iso_intron_seg_indx)

        # Filter
        fail_indx_set = set()
        # Check if isoform intron is in valid intron list
        for intron_indx in iso_intron_indx_list:
            if not cls.is_intron(intron_indx, intron_seg_indx_list):
                fail_indx_set |= set(range(intron_indx[0], intron_indx[1] + 1))
                # print(read.qname, intron_indx, "no intron")

        # Check if internal exon contain spliced segment.
        # If contain spliced segment, forced to interrupt IR.
        IR_file_indx_set = set()
        if len(iso_exon_indx_list) > 2:
            for s_indx, e_indx in spliced_segment_pair:
                isoseg[s_indx: e_indx + 1] = 0
            # for exon_indx in iso_exon_indx_list[1:-1]:
            #     if cls.is_spliced_seg_in_internal_exon(exon_indx, spliced_segment_indx):
            #         fail_indx_set |= set(range(exon_indx[0], exon_indx[1] + 1))
            # print(read.qname, exon_indx, "spliced exon")
            # if not cls.is_internal_exon(exon_indx, exon_seg_indx_list):
            #     fail_indx_set |= set(range(exon_indx[0], exon_indx[1] + 1))
        elif len(iso_exon_indx_list) == 1:
            exon_indx = iso_exon_indx_list[0]
            if cls.is_spliced_seg_in_internal_exon(exon_indx, spliced_segment_indx):
                fail_indx_set |= set(range(exon_indx[0], exon_indx[1] + 1))
            # print(read.qname, exon_indx, "spliced exon")

        if len(iso_exon_indx_list) >= 2:
            # First exon
            first_exon = iso_exon_indx_list[0]
            if first_exon[1] - first_exon[0] > 1:
                check_exon = (first_exon[0] + 2, first_exon[1])
                if cls.is_spliced_seg_in_internal_exon(check_exon, spliced_segment_indx):
                    fail_indx_set |= set(range(check_exon[0], check_exon[1] + 1))
            # Last exon
            last_exon = iso_exon_indx_list[-1]
            if last_exon[1] - last_exon[0] > 1:
                check_exon = (last_exon[0], last_exon[1] - 2)
                if cls.is_spliced_seg_in_internal_exon(check_exon, spliced_segment_indx):
                    fail_indx_set |= set(range(check_exon[0], check_exon[1] + 1))

            if not cls.is_incl_segment(read, segment_list[start_indx], 'right'):
                fail_indx_set.add(start_indx)
                # print(read.qname, start_indx, "start exon right")
            if not cls.is_incl_segment(read, segment_list[end_indx], 'left'):
                fail_indx_set.add(end_indx)
                # print(read.qname, end_indx, "end exon left")

        # Intron segment boundary should be skipped
        for indx in iso_intron_seg_indx:
            if not cls.is_skip_segment(read, segment_list[indx]):
                fail_indx_set.add(indx)
                # print(read.qname, indx, "intron seg")

        # Exon segment boundary should be included
        if len(iso_exon_seg_indx) > 2:
            for indx in iso_exon_seg_indx[1:-1]:
                if not cls.is_incl_segment(read, segment_list[indx]):
                    fail_indx_set.add(indx)
                    # print(read.qname, indx, "exon seg")

        isoseg[sorted(fail_indx_set)] -= 2
        return isoseg

    @classmethod
    def filter_iso(cls, iso, segment_list, intron_seg_indx_list, exon_seg_indx_list):
        """
        Find if the isoform is corrected

        Args:
            iso (Isoform): Isoform to task
            segment_list (list): Valid segment list
            intron_seg_indx_list (list): List of intron segment index
            exon_seg_indx_list (list): List of exon segment index

        Returns:
            bool: Is corrected or not
        """
        '''
        Filter TGS isoform by rules
        :param isoseg: Projected segment for this TGS read
        :param segment_list: Valid segment list
        :return:
        '''
        # Fetch segment information
        spliced_segment_indx = np.where([segment.spliced_seg for segment in segment_list])[0]
        segary = iso.segary

        # Compute exon (start, end) in isoform record
        seg = np.array(segary)
        iso_exon_seg_indx = np.where(seg == 1)[0]
        iso_exon_indx_list = cls.indx2ival(iso_exon_seg_indx)

        # Compute intron (start, end) in isoform record
        start_indx = min(iso_exon_seg_indx)
        end_indx = max(iso_exon_seg_indx)
        iso_intron_seg_indx = np.where(seg == 0)[0]
        iso_intron_seg_indx = np.array(list(filter(
            lambda x: (x >= start_indx) and (x <= end_indx), iso_intron_seg_indx)))
        iso_intron_indx_list = cls.indx2ival(iso_intron_seg_indx)

        # Filter
        # Check if isoform intron is in valid intron list
        if not all([cls.is_intron(
                intron_indx, intron_seg_indx_list) for intron_indx in iso_intron_indx_list]):
            return False

        # Check if internal exon is valid
        if len(iso_exon_indx_list) > 2:
            if any([cls.is_spliced_seg_in_internal_exon(
                    exon_indx, spliced_segment_indx) for exon_indx in iso_exon_indx_list[1:-1]]):
                return False
            # if not all([cls.is_internal_exon(
            #         exon_indx, exon_seg_indx_list) for exon_indx in iso_exon_indx_list[1:-1]]):
            #     return False
        elif len(iso_exon_indx_list) == 1:
            if cls.is_spliced_seg_in_internal_exon(iso_exon_indx_list[0], spliced_segment_indx):
                return False

        if len(iso_exon_indx_list) >= 2:
            # First exon
            first_exon = iso_exon_indx_list[0]
            if first_exon[1] - first_exon[0] > 1:
                check_exon = (first_exon[0] + 2, first_exon[1])
                if cls.is_spliced_seg_in_internal_exon(check_exon, spliced_segment_indx):
                    return False

            # Last exon
            last_exon = iso_exon_indx_list[-1]
            if last_exon[1] - last_exon[0] > 1:
                check_exon = (last_exon[0], last_exon[1] - 2)
                if cls.is_spliced_seg_in_internal_exon(check_exon, spliced_segment_indx):
                    return False
        return True


class TransDiscover(object):
    """ General method in identify transcript. """

    @classmethod
    def element2segment(cls, intron_list, exon_list, trans_ival):
        """ Build transcript segments by transcript elements information. """
        # Init transcript segments
        all_element = intron_list + exon_list
        all_element_boundary = sorted(
            {element.start for element in all_element} | {element.end for element in all_element})
        strand = trans_ival.strand
        segment_list = list()
        for indx in range(len(all_element_boundary) - 1):
            segment_list.append(Segment(
                trans_ival.chrom, all_element_boundary[indx], all_element_boundary[indx + 1], strand))

        # Set TSS tag
        tss_exon_list = list(filter(lambda x: x.tss_exon, exon_list))
        for exon in tss_exon_list:
            if strand == "+":
                segment_indx = all_element_boundary.index(exon.start)
            else:
                segment_indx = all_element_boundary.index(exon.end) - 1
            segment_list[segment_indx].set_tss_seg(True)

        # Set TES tag
        tes_exon_list = list(filter(lambda x: x.tes_exon, exon_list))
        for exon in tes_exon_list:
            if strand == "+":
                segment_indx = all_element_boundary.index(exon.end) - 1
            else:
                segment_indx = all_element_boundary.index(exon.start)
            segment_list[segment_indx].set_tes_seg(True)

        # Set spliced tag
        spliced_intron_list = list(filter(lambda x: x.spliced, intron_list))
        spliced_segment_paired = list()
        for intron in spliced_intron_list:
            spliced_indx1 = all_element_boundary.index(intron.start)
            spliced_indx2 = all_element_boundary.index(intron.end) - 1
            segment_list[spliced_indx1].set_spliced_seg(True)
            segment_list[spliced_indx2].set_spliced_seg(True)
            spliced_segment_paired.append((spliced_indx1, spliced_indx2))

        # Compute FPKM
        ngs_read_num = [x[1] for x in GVAR.NGS_FILE_INFO_list]
        ngs_read_len = [x[2] for x in GVAR.NGS_FILE_INFO_list]
        for segment in segment_list:
            segment.inherit_cov_from(trans_ival)
            segment.compute_fpkm(ngs_read_len, ngs_read_num, GVAR.ALPHA)

        return segment_list, spliced_segment_paired

    @classmethod
    def build_element_seg_indx(cls, element_list, segment_list):
        """ Find intron index in segment list. """
        all_sig_boundary = sorted(
            {seg.start for seg in segment_list} | {seg.end for seg in segment_list})
        element_seg_indx_list = list()
        for element in element_list:
            start_seg_indx = all_sig_boundary.index(element.start)
            end_seg_indx = all_sig_boundary.index(element.end) - 1
            element_seg_indx_list.append((start_seg_indx, end_seg_indx))
        return element_seg_indx_list

    @classmethod
    def build_tss_seg_indx(cls, segment_list):
        """ Find TSS segment index. """
        return np.where([segment.tss_seg for segment in segment_list])[0]

    @classmethod
    def build_tes_seg_indx(cls, segment_list):
        """ Find TES segment index. """
        return np.where([segment.tes_seg for segment in segment_list])[0]

    @classmethod
    def is_compatible(cls, iso1, iso2):
        """Judge if two isoform are compatibility"""
        start_indx = max(iso1.sidx, iso2.sidx)
        end_indx = min(iso1.eidx, iso2.eidx)
        if np.where(iso1.segary == -1)[0].any() or np.where(iso2.segary == -1)[0].any():
            return False
        if start_indx > end_indx:
            return True
        for indx in range(start_indx, end_indx + 1):
            if iso1.segary[indx] != iso2.segary[indx]:
                return False
        return True

    @classmethod
    def build_compatible_matrix(cls, isolist):
        """Build isoform compatibility matrix"""
        iso_num = len(isolist)
        compatible_matrix = np.zeros([iso_num, iso_num], dtype=np.bool) + np.diag(np.ones(iso_num, dtype=np.bool))
        for indx1 in range(iso_num):
            for indx2 in range(indx1):
                tv = cls.is_compatible(isolist[indx1], isolist[indx2])
                compatible_matrix[indx1, indx2] = tv
                compatible_matrix[indx2, indx1] = tv
        return compatible_matrix

    @classmethod
    def is_overlap(cls, iso1, iso2):
        """Judge if two isoform are overlapped"""
        start_indx = max(iso1.sidx, iso2.sidx)
        end_indx = min(iso1.eidx, iso2.eidx)
        return start_indx <= end_indx

    @classmethod
    def build_overlap_matrix(cls, isolist):
        """Build isoform overlap matrix"""
        iso_num = len(isolist)
        overlap_matrix = np.zeros([iso_num, iso_num], dtype=np.bool) + np.diag(np.ones(iso_num, dtype=np.bool))
        for indx1 in range(iso_num):
            for indx2 in range(indx1):
                tv = cls.is_overlap(isolist[indx1], isolist[indx2])
                overlap_matrix[indx1, indx2] = tv
                overlap_matrix[indx2, indx1] = tv
        return overlap_matrix

    @classmethod
    def build_subpath_matrix(cls, compatible_matrix, isolist):
        """ Build isoform subpath matrix, [i, j] = T => i contains j. """
        iso_num = len(isolist)
        subpath_matrix = np.zeros_like(compatible_matrix, dtype=np.bool)
        for indx1 in range(iso_num):
            for indx2 in range(iso_num):
                if not compatible_matrix[indx1, indx2]:
                    continue
                if (isolist[indx1].sidx <= isolist[indx2].sidx) and (isolist[indx1].eidx >= isolist[indx2].eidx):
                    subpath_matrix[indx1, indx2] = True
        return subpath_matrix

    @classmethod
    def is_subpath(cls, iso_indx, subpath_matrix):
        """ Judge whether a isoform is a subpath of another isoform or not. """
        # A isoform is always a subpath of itself.
        return subpath_matrix[:, iso_indx].sum() > 1

    @classmethod
    def search_nfl_cluster(cls, iso_indxs, compatible_matrix):
        """
        The extremely complete subgraph of compatibility matrix can be seem as a cluster of non-full-length isoforms

        Args:
            iso_indxs (narray): Non-full length isoform index array
            compatible_matrix (narray): Isoform compatibility matrix

        Returns:
            set: Set of isoform index in one cluster

        """
        clusters = set()
        G = networkx.Graph()
        for indx1 in iso_indxs:
            for indx2 in np.where(compatible_matrix[indx1, :])[0]:
                G.add_edge(indx1, indx2)

        for indx in iso_indxs:
            clique_list = networkx.cliques_containing_node(G, indx)
            if clique_list:
                for clique in clique_list:
                    clusters.add(tuple(sorted(clique)))
            else:
                clusters.add((indx,))

        return clusters

    @classmethod
    def filter_compatible_with_overlap(cls, indxs, clusters, overlap_matrix):
        """Only overlapped and compatible isoforms can be merged"""
        filter_clusters = set()
        for indx in indxs:
            indx_clusters = list(filter(lambda x: indx in x, clusters))
            for cluster in indx_clusters:
                G = networkx.Graph()
                overlap_pair = [(indx1, indx2) for indx1 in cluster for indx2 in cluster
                                if overlap_matrix[indx1, indx2]]
                G.add_edges_from(overlap_pair)
                component_with_indx = list(filter(
                    lambda x: indx in x, networkx.connected_components(G)))[0]
                connected_indxs = sorted(filter(lambda x: x in component_with_indx, cluster))
                if connected_indxs:
                    filter_clusters.add(tuple(connected_indxs))
        return filter_clusters

    @classmethod
    def merge_cluster_isoforms(cls, indxs, isolist, tss_indxs, tes_indxs):
        """
        Merge multiple isoforms from same cluster

        Args:
            indxs (list): Isoform indx to be merged
            isolist (list): Isoform list
            tss_indxs (list): List of TSS segment index
            tes_indxs (list): List of TES segment index

        Returns:
            Isoform: IsoM
        """

        indxs = np.array(indxs)
        segary = np.zeros_like(isolist[0].segary)
        ival = isolist[0].ival
        for indx in indxs:
            start_indx = isolist[indx].sidx
            end_indx = isolist[indx].eidx
            segary[start_indx: (end_indx + 1)] = isolist[indx].segary[start_indx: (end_indx + 1)]
        merged_iso = Isoform(segary, ival)
        merged_iso.set_tag(tss_indxs, tes_indxs, "M")
        return merged_iso

    @classmethod
    def split_iso_by_subpath(cls, test_iso_list, parent_iso_list):
        """
        Split isoform list into subpath and not subpath.

        Args:
            test_iso_list (list): Isoform list to split
            parent_iso_list (list): Isoform list to reference

        Returns:
            tuple: (subpath, not_subpath)

        """
        all_iso_list = test_iso_list + parent_iso_list
        compatible_matrix = cls.build_compatible_matrix(all_iso_list)
        subpath_matrix = cls.build_subpath_matrix(compatible_matrix, all_iso_list)
        subpath = list()
        not_subpath = list()
        for indx, iso in enumerate(test_iso_list):
            if cls.is_subpath(indx, subpath_matrix):
                subpath.append(iso)
            else:
                not_subpath.append(iso)
        return subpath, not_subpath

    @classmethod
    def merge_nfl_isoforms(cls, iso_indxs, compatible_matrix, overlap_matrix, isolist, tss_indxs, tes_indxs):
        """Clustering non-full-length isoforms into clusters and then merge isoforms in clusters"""
        assert all(map(lambda x: x is not None, [
            iso_indxs, compatible_matrix, isolist, tss_indxs, tes_indxs]))
        # Clustering isoform by compatibility
        clusters = cls.search_nfl_cluster(iso_indxs, compatible_matrix)

        # Drop isoforms that are not overlapped
        clusters_filter_by_overlap = cls.filter_compatible_with_overlap(
            iso_indxs, clusters, overlap_matrix)

        # Merge isoforms in a cluster
        merged_iso_list = list()
        for cluster in clusters_filter_by_overlap:
            merged_iso_list.append(
                cls.merge_cluster_isoforms(cluster, isolist, tss_indxs, tes_indxs))

        # Drop merged isoforms that are subpath of other merged isoforms
        fl_merged_iso = list(filter(lambda x: x.label == "M", merged_iso_list))
        nfl_merged_iso = list(filter(lambda x: x.label != "M", merged_iso_list))
        _, rest_nfl_merged_iso = cls.split_iso_by_subpath(nfl_merged_iso, fl_merged_iso)
        return fl_merged_iso + rest_nfl_merged_iso

    @classmethod
    def complete_iso_by_fl_iso(cls, isoform, fl_segmat, tss_indxs, tes_indxs):
        """Complete isoform with information in full isoform"""
        iso_list = [isoform]
        if fl_segmat.size == 0:
            return iso_list

        strand = isoform.ival.strand
        # Complete left side
        completed_left_isos = list()
        for iso in iso_list:
            if strand == "+":
                if iso.starts_with_tss:
                    completed_left_isos.append(iso)
                    continue
                txs_indx_li = list(filter(lambda x: x <= iso.sidx, tss_indxs))
            else:
                if iso.ends_with_tes:
                    completed_left_isos.append(iso)
                    continue
                txs_indx_li = list(filter(lambda x: x <= iso.sidx, tes_indxs))

            if len(txs_indx_li) == 0:
                return []
            txs_indx = max(txs_indx_li)

            fl_paths_indxs = np.where(
                (fl_segmat[:, iso.sidx] == 1) & (fl_segmat[:, txs_indx] == 1))[0]
            fl_paths = set(map(tuple, fl_segmat[fl_paths_indxs, txs_indx: iso.sidx + 1]))
            if fl_paths:
                for path in fl_paths:
                    tmp_seg = copy.deepcopy(iso.segary)
                    tmp_seg[txs_indx: iso.sidx + 1] = path
                    new_iso = Isoform(tmp_seg, isoform.ival)
                    new_iso.set_tag(tss_indxs, tes_indxs, "C")
                    completed_left_isos.append(new_iso)
            else:
                completed_left_isos.append(iso)

        # Complete right side
        completed_isos = list()
        for iso in completed_left_isos:
            if strand == "+":
                if iso.ends_with_tes:
                    completed_isos.append(iso)
                    continue
                txs_indx_li = list(filter(lambda x: x >= iso.eidx, tes_indxs))
            else:
                if iso.starts_with_tss:
                    completed_isos.append(iso)
                    continue
                txs_indx_li = list(filter(lambda x: x >= iso.eidx, tss_indxs))

            if len(txs_indx_li) == 0:
                return []
            txs_indx = min(txs_indx_li)

            fl_paths_indxs = np.where(
                (fl_segmat[:, iso.eidx] == 1) & (fl_segmat[:, txs_indx] == 1))[0]
            fl_paths = set(map(tuple, fl_segmat[fl_paths_indxs, iso.eidx: txs_indx + 1]))
            if fl_paths:
                for path in fl_paths:
                    tmp_seg = copy.deepcopy(iso.segary)
                    tmp_seg[iso.eidx: txs_indx + 1] = path
                    new_iso = Isoform(tmp_seg, isoform.ival)
                    new_iso.set_tag(tss_indxs, tes_indxs, "C")
                    completed_isos.append(new_iso)
            else:
                completed_isos.append(iso)
        return completed_isos

    @classmethod
    def cluster_iso(cls, iso_list):
        """ Clustering isoforms by overlap """
        iso_num = len(iso_list)

        # Build overlap graph
        overlap_graph = networkx.Graph()
        for indx1 in range(iso_num):
            for indx2 in range(indx1):
                if TransDiscover.is_overlap(iso_list[indx1], iso_list[indx2]):
                    overlap_graph.add_edge(indx1, indx2)

        # The relationship_tuple represents the overlap relationship of one isoform to others.
        cliques = list(networkx.find_cliques(overlap_graph))
        cliques_dict = dict()
        for indx in range(iso_num):
            relationship_tuple = tuple([indx in x for x in cliques])
            if relationship_tuple in cliques_dict.keys():
                cliques_dict[relationship_tuple].append(indx)
            else:
                cliques_dict[relationship_tuple] = [indx]

        clusters = list()
        for indxs in cliques_dict.values():
            sub_cluster = [iso_list[indx] for indx in indxs]
            clusters.append(sub_cluster)

        return clusters

    @classmethod
    def similar_score(cls, segary1, segary2):
        """ Compute similar score between two segment array """
        return (np.array(segary1) == np.array(segary2)).sum()

    @classmethod
    def rescue_junction(cls, invalid_iso, intron_seg_indx_list, exon_seg_indx_list, segment_list):
        """ Rescue junction by validated information """
        segary = list()
        for x in invalid_iso.segary:
            if x < 0:
                segary.append(x + 2)
            else:
                segary.append(x)
        segary = np.array(segary)
        fail_indx_set = set()

        intron_seg_indx = np.where(segary == 0)[0]
        for intron_indx in TgsFilterRule.indx2ival(intron_seg_indx):
            if not TgsFilterRule.is_intron(intron_indx, intron_seg_indx_list):
                fail_indx_set |= set(range(intron_indx[0], intron_indx[1] + 1))

        splices_seg_indx = np.where([x.spliced_seg for x in segment_list])[0]
        exon_seg_indx = np.where(segary == 1)[0]
        for s, e in TgsFilterRule.indx2ival(exon_seg_indx):
            if any([(s <= x <= e) for x in splices_seg_indx]):
                fail_indx_set |= set(range(s, e + 1))
        segary[sorted(fail_indx_set)] -= 2

        fail_indx_set = set()
        exon_seg_indx = np.where(segary == 1)[0]
        for exon_indx in TgsFilterRule.indx2ival(exon_seg_indx):
            if not TgsFilterRule.is_internal_exon(exon_indx, exon_seg_indx_list):
                fail_indx_set |= set(range(exon_indx[0], exon_indx[1] + 1))
        segary[sorted(fail_indx_set)] -= 2
        invalid_iso.segary = segary

    @classmethod
    def rescue_isoform(cls, invalid_iso, isoF, intron_seg_indx_list=None, exon_seg_indx_list=None, segment_list=None):
        """Rescue invalid full-length TGS isoform by fixing errors from most similar isoform"""
        if (intron_seg_indx_list is not None) and (segment_list is not None) and (exon_seg_indx_list is not None):
            cls.rescue_junction(invalid_iso, intron_seg_indx_list, exon_seg_indx_list, segment_list)
        invalid_indx = np.where(invalid_iso.segary < 0)[0]
        if (not invalid_indx.size) and invalid_iso.starts_with_tss and invalid_iso.ends_with_tes:
            return [invalid_iso]

        # Just rescued isoforms with TSS and TES
        if not (invalid_iso.starts_with_tss and invalid_iso.ends_with_tes):
            return list()
        if np.where((invalid_indx <= invalid_iso.sidx) | (invalid_indx >= invalid_iso.eidx))[0].size:
            return list()

        overlaped_isoF = list(filter(lambda x: (x.sidx <= invalid_iso.sidx) and (x.eidx >= invalid_iso.eidx), isoF))
        if not overlaped_isoF:
            return list()

        segary_set = set()
        new_segary_set = set()
        segary = tuple(invalid_iso.segary)
        segary_set.add(segary)
        invalid_regions = TgsFilterRule.indx2ival(invalid_indx)

        for start_indx, end_indx in invalid_regions:
            matched_isoF = list(filter(
                lambda x: (x.segary[start_indx - 1] == segary[start_indx - 1]) and
                          (x.segary[end_indx + 1] == segary[end_indx + 1]) and
                          np.all(x.segary[start_indx: (end_indx + 1)] == (
                                  2 + np.array(segary[start_indx: (end_indx + 1)]))),  # 2, number_of_status
                overlaped_isoF))
            if not matched_isoF:
                segary_set = set()
                break
            else:
                for invalid_segary in segary_set:
                    similar_score = [cls.similar_score(invalid_segary, x) for x in matched_isoF]
                    best_iso = matched_isoF[np.argmax(similar_score)]
                    tmp_invalid_segary = list(invalid_segary)
                    tmp_invalid_segary[start_indx:end_indx + 1] = best_iso.segary[start_indx:end_indx + 1]
                    new_segary_set.add(tuple(tmp_invalid_segary))
            segary_set = new_segary_set
            new_segary_set = set()
        rescued_iso_list = [Isoform(x, invalid_iso.ival) for x in segary_set]
        return rescued_iso_list

    @classmethod
    def create_intron_path(cls, intron_nodes, start_seg_indx, end_seg_indx):
        """ Create intron path by intron index information """
        intron_path = [1] * (end_seg_indx - start_seg_indx + 1)
        for t, s, e in intron_nodes:
            if t != "intron":
                continue
            for i in range(s, e + 1):
                assert i - start_seg_indx >= 0
                intron_path[i - start_seg_indx] = 0
        return tuple(intron_path)

    @classmethod
    def compute_seg_link_cnt(cls, ngs_file, segli):
        """Count reads number supported the given segment pair"""
        boundarys = [x.start for x in segli] + [segli[-1].end]
        chrom = segli[0].chrom
        start = min(boundarys)
        end = max(boundarys)
        strand = segli[0].strand
        ival = Interval(chrom, start, end, strand)
        reads = ngs_file.fetch_reads_in_ival(ival)
        intron_link_cnt = defaultdict(int)
        exon_link_cnt = defaultdict(int)
        for read in reads:
            for s, e in GeneLinkageFinder.find_intron(read):
                if e < start:
                    continue
                if s >= end:
                    break
                if (s not in boundarys) or (e not in boundarys):
                    continue
                sidx = boundarys.index(s)
                eidx = boundarys.index(e) - 1
                key = (sidx, eidx)
                intron_link_cnt[key] += 1

            for s, e in read.blocks:
                for b in filter(lambda x: s < x < e, boundarys):
                    bindx = boundarys.index(b)
                    key = (bindx - 1, bindx)
                    exon_link_cnt[key] += 1
        seglink_cnt = {"intron": intron_link_cnt, "exon": exon_link_cnt}
        return seglink_cnt

    @classmethod
    def invalid_iso_is_subpath(cls, iso, fl_isos):
        """If a invalid isoform is a subpath for a full length isoform"""
        iso_seg = iso.segary
        exon_indx = np.where(iso_seg == 1)[0]
        if not exon_indx.size:
            return True
        sidx = min(exon_indx)
        eidx = max(exon_indx)
        for fl_iso in fl_isos:
            missmatch = False
            for i in range(sidx, eidx + 1):
                if iso_seg[i] < 0:
                    continue
                if iso_seg[i] != fl_iso.segary[i]:
                    missmatch = True
                    break
            if not missmatch:
                return True
        return False

    @classmethod
    def determin_intron_type(cls, intron, seglink_exon_cnt, seglink_intron_cnt, segment_list):
        """ Find if a intron could be an IR """
        s, e = intron
        tissue_num = len(seglink_exon_cnt)
        if tissue_num:
            exon_cnts = [(exon_cnt_dict[(s - 1, s)] + exon_cnt_dict[(e, e + 1)]) / 2 for exon_cnt_dict in
                         seglink_exon_cnt]
            intron_cnts = [intron_cnt_dict[(s, e)] for intron_cnt_dict in seglink_intron_cnt]
            phi = np.array([(i_cnt + 0.5) / (i_cnt + e_cnt + 1) for (i_cnt, e_cnt) in zip(intron_cnts, exon_cnts)])
            intron_tis_num = (phi > GVAR.MAJOR_AS_PHT_CUTOFF).sum()
            ir_tis_num = (phi < GVAR.MINOR_AS_PHI_CUTOFF).sum()
            both_tis_num = len(exon_cnts) - intron_tis_num - ir_tis_num
            selected_indx = np.argmax([intron_tis_num, ir_tis_num, both_tis_num])
        else:
            selected_indx = 2

        if selected_indx == 0:
            return ["intron"]
        elif selected_indx == 1:
            return ["IR"]
        else:
            if any([segment_list[x].spliced_seg for x in range(s, e + 1)]):
                return ["intron"]
            else:
                return ["intron", "IR"]

    @classmethod
    def compute_as_cnt(cls, intron, seglink_exon_cnt, seglink_intron_cnt):
        """ Compute read count for an intron """
        t, s, e = intron
        if t == "intron":
            return np.array([intron_cnt_dict[(s, e)] for intron_cnt_dict in seglink_intron_cnt])
        else:
            return np.array([(exon_cnt_dict[(s - 1, s)] + exon_cnt_dict[(e, e + 1)]) / 2 for exon_cnt_dict in
                             seglink_exon_cnt])

    @classmethod
    def enum_intron_path(cls, intron_seg_indx_list, start_seg_indx, end_seg_indx, segment_list, seglink_cnt):
        """Enumerate intron graph in certain region."""

        seglink_exon_cnt = [x["exon"] for x in seglink_cnt]
        seglink_intron_cnt = [x["intron"] for x in seglink_cnt]
        tissue_num = len(seglink_exon_cnt)

        # Filter intron in certain region
        targer_intron_indx_list = [x for x in intron_seg_indx_list if x[0] >= start_seg_indx and x[1] <= end_seg_indx]
        targer_intron_indx_list = sorted(targer_intron_indx_list, key=lambda x: (x[0], x[1]))
        if not targer_intron_indx_list:
            return set()

        # Build intron graph by position (link the nearest intron pairs)
        intron_graph = networkx.DiGraph()
        intron_edges = list()
        intron_nodes = set()

        for intron_indx1 in range(0, len(targer_intron_indx_list)):
            intron1 = targer_intron_indx_list[intron_indx1]
            intron1_li = cls.determin_intron_type(intron1, seglink_exon_cnt, seglink_intron_cnt, segment_list)
            for intron1_tv in intron1_li:
                intron_nodes.add((intron1_tv, *intron1))
            for intron_indx2 in range(intron_indx1 - 1, -1, -1):
                intron2 = targer_intron_indx_list[intron_indx2]
                intron2_li = cls.determin_intron_type(intron2, seglink_exon_cnt, seglink_intron_cnt, segment_list)
                for intron2_tv in intron2_li:
                    intron_nodes.add((intron2_tv, *intron2))
                if intron2[1] >= intron1[0]:
                    continue
                tv_break = False
                for tmp_indx in range(intron_indx2 + 1, intron_indx1):
                    tmp_intron = targer_intron_indx_list[tmp_indx]
                    if (intron2[1] < tmp_intron[0]) and (intron1[0] > tmp_intron[1]):
                        tv_break = True
                        break
                if not tv_break:
                    intron_edges += [((intron2_tv, *intron2), (intron1_tv, *intron1)) for intron2_tv in intron2_li for
                                     intron1_tv in intron1_li]

        if len(intron_nodes) > GVAR.MAX_AS_NUM_FOR_ISOP:
            #  Too complex to compute, need large memory
            return None

        intron_graph.add_nodes_from(list(intron_nodes))
        intron_graph.add_edges_from(intron_edges)
        root_nodes = {x[0] for x in dict(intron_graph.in_degree()).items() if x[1] == 0}
        leaf_nodes = {x[0] for x in dict(intron_graph.out_degree()).items() if x[1] == 0}

        # Find small probability as event
        minor_dict = dict()
        for n in intron_graph.nodes:
            t, s, e = n
            if tissue_num:
                target_as_cnt = cls.compute_as_cnt(n, seglink_exon_cnt, seglink_intron_cnt)
                overlap_nodes = list(filter(lambda n: (n[1] <= e) and (n[2] >= s), intron_graph.nodes))
                overlap_cnts = np.array(
                    [cls.compute_as_cnt(x, seglink_exon_cnt, seglink_intron_cnt) for x in overlap_nodes])
                minor_dict[n] = np.array((target_as_cnt / overlap_cnts.max(0)) <= GVAR.MINOR_AS_CNT_RATIO, dtype=np.int)
            else:
                minor_dict[n] = np.array([], dtype=np.int)
        intron_path_set = set()
        
        # Single intron node
        for n in root_nodes.intersection(leaf_nodes):
            p = cls.create_intron_path([n], start_seg_indx, end_seg_indx)
            minor = minor_dict[n]
            intron_path_set.add((p, tuple(minor)))

        # Enumerate intron path
        intron_path_num = 0
        for r in root_nodes:
            for l in leaf_nodes:
                for p in networkx.all_simple_paths(intron_graph, r, l):
                    minor = np.zeros(tissue_num)
                    for n in p:
                        minor += minor_dict[n]
                    if minor.size and (minor > 1).all():
                        continue
                    intron_path_set.add((cls.create_intron_path(p, start_seg_indx, end_seg_indx), tuple(minor)))
                    intron_path_num += 1
                    if intron_path_num > GVAR.MAX_ISO_NUM_FOR_ISOP:
                        return None
        return intron_path_set

    @classmethod
    def complete_partial_isoform_left(cls, segli, minor_list, ta):
        """ Complete the left site missing region of a partial isoform """
        if not segli:
            return list(), list()
        seg_end_indx = min(np.where(np.array(segli[0]) != 0)[0]) - 1
        if ta.ival.strand == "+":
            target_start_indx = [s for s in ta.tss_indx if
                                 (s <= seg_end_indx) and (not any([s < e < seg_end_indx for e in ta.tes_indx]))]
        else:
            target_start_indx = [s for s in ta.tes_indx if
                                 (s <= seg_end_indx) and (not any([s < e < seg_end_indx for e in ta.tss_indx]))]
        left_complete_segli = list()
        left_complete_minor_list = list()
        for s_indx in target_start_indx:
            intron_paths = ta.intron_path(s_indx, seg_end_indx)
            if intron_paths is None:
                return list(), list()
            for seg, minor in zip(segli, minor_list):
                tmp_seg_li = [list(seg) for _ in range(max(1, len(intron_paths)))]
                tmp_minor_mat = np.array([minor for _ in range(max(1, len(intron_paths)))])
                if intron_paths:
                    for indx, (path, path_minor) in enumerate(intron_paths):
                        tmp_seg_li[indx][s_indx: (seg_end_indx + 1)] = path
                        tmp_minor_mat[indx, :] += np.array(path_minor)
                else:
                    tmp_seg_li[0][s_indx: (seg_end_indx + 1)] = [1] * (seg_end_indx + 1 - s_indx)
                for indx, tmp_seg in enumerate(tmp_seg_li):
                    tmp_minor = tmp_minor_mat[indx, :]
                    if tmp_minor.size and (tmp_minor > 1).all():
                        continue
                    left_complete_segli.append(tmp_seg)
                    left_complete_minor_list.append(tmp_minor)
        return left_complete_segli, left_complete_minor_list

    @classmethod
    def complete_partial_isoform_right(cls, segli, minor_list, ta):
        """ Complete the missing region of a partial isoform """
        if not segli:
            return list(), list()
        seg_start_indx = max(np.where(np.array(segli[0]) != 0)[0]) + 1
        if ta.ival.strand == "+":
            target_end_indx = [e for e in ta.tes_indx if
                               (e >= seg_start_indx) and (not any([e > s > seg_start_indx for s in ta.tss_indx]))]
        else:
            target_end_indx = [e for e in ta.tss_indx if
                               (e >= seg_start_indx) and (not any([e > s > seg_start_indx for s in ta.tes_indx]))]
        right_complete_segli = list()
        right_complete_minor_list = list()
        for e_indx in target_end_indx:
            intron_paths = ta.intron_path(seg_start_indx, e_indx)
            if intron_paths is None:
                return list(), list()
            for seg, minor in zip(segli, minor_list):
                tmp_seg_li = [list(seg) for _ in range(max(1, len(intron_paths)))]
                tmp_minor_mat = np.array([minor for _ in range(max(1, len(intron_paths)))])
                if intron_paths:
                    for indx, (path, path_minor) in enumerate(intron_paths):
                        tmp_seg_li[indx][seg_start_indx: (e_indx + 1)] = path
                        tmp_minor_mat[indx, :] += np.array(path_minor)
                else:
                    tmp_seg_li[0][seg_start_indx: (e_indx + 1)] = [1] * (e_indx + 1 - seg_start_indx)
                for indx, tmp_seg in enumerate(tmp_seg_li):
                    tmp_minor = tmp_minor_mat[indx, :]
                    if tmp_minor.size and (tmp_minor > 1).all():
                        continue
                    right_complete_segli.append(tmp_seg)
                    right_complete_minor_list.append(tmp_minor)
        return right_complete_segli, right_complete_minor_list

    @classmethod
    def complete_partial_isoform_error(cls, segli, minor_list, ta, s_indx, e_indx):
        """ Try to corrected the errors in a isoform """
        if not segli:
            return list(), list()
        error_complete_segli = list()
        error_complete_minor_list = list()
        intron_paths = ta.intron_path(s_indx, e_indx)
        if intron_paths is None:
            return list(), list()
        for seg, minor in zip(segli, minor_list):
            tmp_seg_li = [seg] * len(intron_paths)
            tmp_minor_mat = np.array([minor for _ in range(max(1, len(intron_paths)))])
            for indx, (path, path_minor) in enumerate(intron_paths):
                tmp_seg_li[indx][s_indx: (e_indx + 1)] = path
                tmp_minor_mat[indx, :] += np.array(path_minor)
            for indx, tmp_seg in enumerate(tmp_seg_li):
                tmp_minor = tmp_minor_mat[indx, :]
                if tmp_minor.size and (tmp_minor > 1).all():
                    continue
                error_complete_segli.append(tmp_seg)
                error_complete_minor_list.append(tmp_minor)
        return error_complete_segli, error_complete_minor_list

    @classmethod
    def complete_partial_isoform(cls, isoform, ta):
        """ Complete a partial isoform """
        segli = [list(isoform.segary)]
        minor_list = [np.zeros(ta.tissue_num)]
        if ta.ival.strand == "+":
            if not isoform.starts_with_tss:
                segli, minor_list = cls.complete_partial_isoform_left(segli, minor_list, ta)
            if len(segli) > GVAR.MAX_ISO_NUM_FOR_ISOP:
                return list()
            if not isoform.ends_with_tes:
                segli, minor_list = cls.complete_partial_isoform_right(segli, minor_list, ta)
            if len(segli) > GVAR.MAX_ISO_NUM_FOR_ISOP:
                return list()
        else:
            if not isoform.ends_with_tes:
                segli, minor_list = cls.complete_partial_isoform_left(segli, minor_list, ta)
            if len(segli) > GVAR.MAX_ISO_NUM_FOR_ISOP:
                return list()
            if not isoform.starts_with_tss:
                segli, minor_list = cls.complete_partial_isoform_right(segli, minor_list, ta)
            if len(segli) > GVAR.MAX_ISO_NUM_FOR_ISOP:
                return list()
        err_seg_indx = np.where(isoform.segary < 0)[0]
        err_seg_ival = TgsFilterRule.indx2ival(err_seg_indx)
        for s, e in err_seg_ival:
            segli, minor_list = cls.complete_partial_isoform_error(segli, minor_list, ta, s, e)
        complete_iso_li = list()
        for seg in segli:
            new_iso = Isoform(seg, isoform.ival)
            new_iso.set_tag(ta.tss_indx, ta.tes_indx, "P")
            complete_iso_li.append(new_iso)
        return complete_iso_li


class TgsTransDiscover(TransDiscover):
    """Identify transcript with TGS data"""

    @classmethod
    def compute_iso_overlap_fraction(cls, read, segment):
        """ Compute the overlap fraction between a read and a segment """
        pos_array = np.array(read.positions)
        overlap_len = ((pos_array >= segment.start) & (pos_array < segment.end)).sum()
        return overlap_len / len(segment)

    @classmethod
    def refine(
            cls, tgs_read_list, segment_list, intron_seg_indx_list, exon_seg_indx_list, spliced_segment_pair,
            tgs_overlap_fraction_threshold, trans_ival):
        """Project TGS read to segment list, refine and create Isoform"""
        iso2num = {}
        for iso_read in tgs_read_list:
            isoseg = np.array([
                (cls.compute_iso_overlap_fraction(iso_read, segment) >= tgs_overlap_fraction_threshold)
                for segment in segment_list], dtype=np.int)
            if not (isoseg == 1).any():
                continue
            isoseg = TgsFilterRule.filter(
                iso_read, isoseg, segment_list, intron_seg_indx_list, exon_seg_indx_list, spliced_segment_pair)

            if any(isoseg == 1):
                segkey = tuple(isoseg)
                if segkey not in iso2num:
                    iso2num[segkey] = 0
                iso2num[segkey] += iso_read.get_tag('FL')

        iso_list = []
        for isoseg in sorted(iso2num):
            iso = Isoform(isoseg, trans_ival)
            iso.flnum = iso2num[isoseg]
            iso_list.append(iso)

        return iso_list


def identify_transcript(gene, ann=None):
    """ Method to identity transcript """
    trans = TransAssembler(gene, ann)
    trans.identify_isoform()
    return trans
