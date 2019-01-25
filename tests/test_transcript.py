import unittest
import os
from collections import defaultdict
from igia.utils import SeqFile, make_read
from igia.element import identify_element
from igia.transcript import *


class TA(object):
    def __init__(self, ival, segment_list, intron_indx, exon_indx, tss_indx, tes_indx):
        self.ival = ival
        self.segment_list = segment_list
        self.intron_indx = intron_indx
        self.exon_indx = exon_indx
        self.tss_indx = tss_indx
        self.tes_indx = tes_indx
        self._intron_path = dict()
        self.tissue_num = 0
        self.seglink_cnt = list()

    def intron_path(self, start_seg_indx, end_seg_indx):
        key = (start_seg_indx, end_seg_indx)
        if key not in self._intron_path.keys():
            self._intron_path[key] = TransDiscover.enum_intron_path(self.intron_indx, start_seg_indx, end_seg_indx,
                                                                    self.segment_list, self.seglink_cnt)
        return self._intron_path[key]


class TestElement2Segment(unittest.TestCase):
    def setUp(self):
        names = ('Bract', 'Cotyledon')
        ngs = [os.path.join('tests/data', "%s.sort.bam" % name) for name in names]
        tgs = [os.path.join('tests/data', "all_fixed_star.sort.bam")]
        ngs = [SeqFile(fname, 'NGS') for fname in ngs]
        tgs = [SeqFile(fname, 'TGS') for fname in tgs]

        # data for test element2segment
        tss_site_list = [("Chr01", x, "-") for x in [95907, 96077]]
        tes_site_list = [("Chr01", x, "-") for x in [93802]]
        self.rev_gene = identify_element("Chr01", 90000, 98000, ngs, tgs, tss_site_list, tes_site_list)[0]

        tss_site_list = [("Chr01", x, "+") for x in [485839]]
        tes_site_list = [("Chr01", x, "+") for x in [491380]]
        self.fwd_gene = identify_element("Chr01", 485500, 493000, ngs, tgs, tss_site_list, tes_site_list)[0]

    def test_element2segment_rev(self):
        exon_list = self.rev_gene.internal_exon_list + self.rev_gene.tss_exon_list + self.rev_gene.tes_exon_list
        seglist, _ = TransDiscover.element2segment(self.rev_gene.intron_list, exon_list, self.rev_gene.ival)
        self.assertEqual(str(seglist),
                         "[Segment: Chr01:93802-94148 -; tss=False; tes=True; spliced=False, " +
                         "Segment: Chr01:94148-94263 -; tss=False; tes=False; spliced=True, " +
                         "Segment: Chr01:94263-94877 -; tss=False; tes=False; spliced=False, " +
                         "Segment: Chr01:94877-94959 -; tss=False; tes=False; spliced=True, " +
                         "Segment: Chr01:94959-95353 -; tss=False; tes=False; spliced=False, " +
                         "Segment: Chr01:95353-95459 -; tss=False; tes=False; spliced=True, " +
                         "Segment: Chr01:95459-95527 -; tss=False; tes=False; spliced=False, " +
                         "Segment: Chr01:95527-95801 -; tss=False; tes=False; spliced=True, " +
                         "Segment: Chr01:95801-95874 -; tss=False; tes=False; spliced=False, " +
                         "Segment: Chr01:95874-95907 -; tss=True; tes=False; spliced=False, " +
                         "Segment: Chr01:95907-96076 -; tss=False; tes=False; spliced=False, " +
                         "Segment: Chr01:96076-96077 -; tss=True; tes=False; spliced=False]")

    def test_element2segment_fwd(self):
        exon_list = self.fwd_gene.internal_exon_list + self.fwd_gene.tss_exon_list + self.fwd_gene.tes_exon_list
        seglist, _ = TransDiscover.element2segment(self.fwd_gene.intron_list, exon_list, self.fwd_gene.ival)
        self.assertEqual(str(seglist),
                         "[Segment: Chr01:485839-486149 +; tss=True; tes=False; spliced=False, " +
                         "Segment: Chr01:486149-488629 +; tss=False; tes=False; spliced=True, " +
                         "Segment: Chr01:488629-488687 +; tss=False; tes=False; spliced=False, " +
                         "Segment: Chr01:488687-488793 +; tss=False; tes=False; spliced=True, " +
                         "Segment: Chr01:488793-488874 +; tss=False; tes=False; spliced=False, " +
                         "Segment: Chr01:488874-488989 +; tss=False; tes=False; spliced=True, " +
                         "Segment: Chr01:488989-489025 +; tss=False; tes=False; spliced=False, " +
                         "Segment: Chr01:489025-489151 +; tss=False; tes=False; spliced=True, " +
                         "Segment: Chr01:489151-489388 +; tss=False; tes=False; spliced=False, " +
                         "Segment: Chr01:489388-489475 +; tss=False; tes=False; spliced=True, " +
                         "Segment: Chr01:489475-489617 +; tss=False; tes=False; spliced=False, " +
                         "Segment: Chr01:489617-490544 +; tss=False; tes=False; spliced=True, " +
                         "Segment: Chr01:490544-491380 +; tss=False; tes=True; spliced=False, " +
                         "Segment: Chr01:491380-492288 +; tss=False; tes=True; spliced=False]")


class TestTransDiscover(unittest.TestCase):
    def setUp(self):
        names = ('Bract', 'Cotyledon')
        ngs = [os.path.join('tests/data', "%s.sort.bam" % name) for name in names]
        tgs = [os.path.join('tests/data', "all_fixed_star.sort.bam")]
        ngs = [SeqFile(fname, 'NGS') for fname in ngs]
        tgs = [SeqFile(fname, 'TGS') for fname in tgs]

        # data for test element2segment
        tss_site_list = [("Chr01", x, "-") for x in [95907, 96077]]
        tes_site_list = [("Chr01", x, "-") for x in [93802]]
        self.rev_gene = identify_element("Chr01", 90000, 98000, ngs, tgs, tss_site_list, tes_site_list)[0]

        tss_site_list = [("Chr01", x, "+") for x in [485839]]
        tes_site_list = [("Chr01", x, "+") for x in [491380]]
        self.fwd_gene = identify_element("Chr01", 485500, 493000, ngs, tgs, tss_site_list, tes_site_list)[0]

        # data for test build segment index
        exon_list = self.fwd_gene.internal_exon_list + self.fwd_gene.tss_exon_list + self.fwd_gene.tes_exon_list
        self.fwd_intron_list = self.fwd_gene.intron_list
        self.fwd_seglist, _ = TransDiscover.element2segment(self.fwd_intron_list, exon_list, self.fwd_gene.ival)

        exon_list = self.rev_gene.internal_exon_list + self.rev_gene.tss_exon_list + self.rev_gene.tes_exon_list
        self.rev_intron_list = self.rev_gene.intron_list
        self.rev_seglist, _ = TransDiscover.element2segment(self.rev_intron_list, exon_list, self.rev_gene.ival)

        # data for test isoform compatibility
        ival = Interval("test", 100, 200, "+")
        tssseg_indx = np.array([0, 1])
        tesseg_indx = np.array([9])

        """
        0123456789
        -=-=-=-==-
        -=-=---==-
        -=-=-x-==-
        -=-=-=----
        -------==-
        ---=-=-==-
        """
        iso_obj1 = Isoform(np.array([0, 1, 0, 1, 0, 1, 0, 1, 1, 0]), ival)
        iso_obj2 = Isoform(np.array([0, 1, 0, 1, 0, 0, 0, 1, 1, 0]), ival)
        iso_obj3 = Isoform(np.array([0, 1, 0, 1, 0, -1, 0, 1, 1, 0]), ival)
        iso_obj4 = Isoform(np.array([0, 1, 0, 1, 0, 1, 0, 0, 0, 0]), ival)
        iso_obj5 = Isoform(np.array([0, 0, 0, 0, 0, 0, 0, 1, 1, 0]), ival)
        iso_obj6 = Isoform(np.array([0, 0, 0, 1, 0, 1, 0, 1, 1, 0]), ival)
        self.iso_list = [iso_obj1, iso_obj2, iso_obj3, iso_obj4, iso_obj5, iso_obj6]
        for iso in self.iso_list:
            iso.set_tag(tssseg_indx, tesseg_indx, "test")

    def test_build_element_seg_indx_rev(self):
        intron_seg_indx = TransDiscover.build_element_seg_indx(self.rev_intron_list, self.rev_seglist)
        self.assertListEqual(sorted(intron_seg_indx), [(1, 1), (3, 3), (5, 5), (7, 7), (7, 10), (9, 10)])

    def test_build_element_seg_indx_fwd(self):
        intron_seg_indx = TransDiscover.build_element_seg_indx(self.fwd_intron_list, self.fwd_seglist)
        self.assertListEqual(sorted(intron_seg_indx), [(1, 1), (3, 3), (5, 5), (7, 7), (9, 9), (11, 11)])

    def test_build_tss_seg_indx_rev(self):
        tss_seg_indx = TransDiscover.build_tss_seg_indx(self.rev_seglist)
        np.testing.assert_array_equal(tss_seg_indx, np.array([9, 11], dtype=np.int))

    def test_build_tss_seg_indx_fwd(self):
        tss_seg_indx = TransDiscover.build_tss_seg_indx(self.fwd_seglist)
        np.testing.assert_array_equal(tss_seg_indx, np.array([0], dtype=np.int))

    def test_build_tes_seg_indx_rev(self):
        tes_seg_indx = TransDiscover.build_tes_seg_indx(self.rev_seglist)
        np.testing.assert_array_equal(tes_seg_indx, np.array([0], dtype=np.int))

    def test_build_tes_seg_indx_fwd(self):
        tes_seg_indx = TransDiscover.build_tes_seg_indx(self.fwd_seglist)
        np.testing.assert_array_equal(tes_seg_indx, np.array([12, 13], dtype=np.int))

    def test_is_compatible(self):
        # compatibible with itself
        self.assertTrue(TransDiscover.is_compatible(self.iso_list[0], self.iso_list[0]))
        # conflict
        self.assertFalse(TransDiscover.is_compatible(self.iso_list[0], self.iso_list[1]))
        # invalid and thus not compatible with others
        self.assertFalse(TransDiscover.is_compatible(self.iso_list[0], self.iso_list[2]))
        # non-overlap and thus compatible with each other
        self.assertTrue(TransDiscover.is_compatible(self.iso_list[3], self.iso_list[4]))
        # overlap partial reads are compatible
        self.assertTrue(TransDiscover.is_compatible(self.iso_list[3], self.iso_list[5]))

    def test_build_compatible_matrix(self):
        np.testing.assert_array_equal(
            TransDiscover.build_compatible_matrix(self.iso_list),
            np.array([[True, False, False, True, True, True],
                      [False, True, False, False, True, False],
                      [False, False, True, False, False, False],
                      [True, False, False, True, True, True],
                      [True, True, False, True, True, True],
                      [True, False, False, True, True, True]])
        )

    def test_is_overlap(self):
        self.assertTrue(TransDiscover.is_overlap(self.iso_list[0], self.iso_list[1]))
        self.assertFalse(TransDiscover.is_overlap(self.iso_list[3], self.iso_list[4]))

    def test_build_overlap_matrix(self):
        np.testing.assert_array_equal(
            TransDiscover.build_overlap_matrix(self.iso_list),
            np.array([[True, True, True, True, True, True],
                      [True, True, True, True, True, True],
                      [True, True, True, True, True, True],
                      [True, True, True, True, False, True],
                      [True, True, True, False, True, True],
                      [True, True, True, True, True, True]])
        )

    def test_is_subpath(self):
        compatible_matrix = TransDiscover.build_compatible_matrix(self.iso_list)
        subpath_matrix = TransDiscover.build_subpath_matrix(compatible_matrix, self.iso_list)
        self.assertFalse(TransDiscover.is_subpath(0, subpath_matrix))
        self.assertFalse(TransDiscover.is_subpath(1, subpath_matrix))
        self.assertFalse(TransDiscover.is_subpath(2, subpath_matrix))
        self.assertTrue(TransDiscover.is_subpath(3, subpath_matrix))
        self.assertTrue(TransDiscover.is_subpath(4, subpath_matrix))
        self.assertTrue(TransDiscover.is_subpath(5, subpath_matrix))

    def test_build_subpath_matrix(self):
        compatible_matrix = TransDiscover.build_compatible_matrix(self.iso_list)
        subpath_matrix = TransDiscover.build_subpath_matrix(compatible_matrix, self.iso_list)
        np.testing.assert_array_equal(
            subpath_matrix,
            np.array([[True, False, False, True, True, True],
                      [False, True, False, False, True, False],
                      [False, False, True, False, False, False],
                      [False, False, False, True, False, False],
                      [False, False, False, False, True, False],
                      [False, False, False, False, True, True]])
        )

    def test_search_nfl_cluster(self):
        compatible_matrix = np.array(
            [[True, True, True, False, False, True, False],
             [True, True, True, False, False, False, False],
             [True, True, True, True, False, False, False],
             [False, False, True, True, True, False, False],
             [False, False, False, True, True, True, False],
             [True, False, False, False, True, True, False],
             [False, False, False, False, False, False, True]])
        cluster = sorted(TransDiscover.search_nfl_cluster(list(range(7)), compatible_matrix))
        self.assertListEqual(cluster, [(0, 1, 2), (0, 5), (2, 3), (3, 4), (4, 5), (6,)])

    def test_filter_compatible_with_overlap(self):
        clusters = [(0, 1, 2), (3, 4, 5), (5,)]
        overlap_matrix = np.array(
            [[True, True, True, False, False, False],
             [True, True, True, False, False, False],
             [True, True, True, False, False, False],
             [False, False, False, True, True, False],
             [False, False, False, True, True, False],
             [False, False, False, False, False, True]])
        self.assertListEqual(sorted(TransDiscover.filter_compatible_with_overlap(
            list(range(6)), clusters, overlap_matrix
        )), [(0, 1, 2), (3, 4), (5,)])

    def test_merge_cluster_isoforms(self):
        merged_iso = TransDiscover.merge_cluster_isoforms(
            np.array([3, 5]), self.iso_list, np.array([0, 1]), np.array([9]))
        np.testing.assert_array_equal(merged_iso.segary, np.array([0, 1, 0, 1, 0, 1, 0, 1, 1, 0]))

    def test_merge_nfl_isoforms(self):
        ival = Interval("test", 100, 200, "+")
        tssseg_indx = np.array([0, 1])
        tesseg_indx = np.array([9])
        iso_obj1 = Isoform(np.array([0, 0, 0, 0, 0, 1, 0, 1, 1, 0]), ival)
        iso_obj2 = Isoform(np.array([0, 1, 0, 1, 0, 1, 0, 0, 0, 0]), ival)
        iso_obj3 = Isoform(np.array([0, 1, 1, 1, 0, 1, 1, 0, 0, 0]), ival)
        iso_obj4 = Isoform(np.array([1, 1, 0, 1, 0, 0, 0, 0, 0, 0]), ival)
        isolist = [iso_obj1, iso_obj2, iso_obj3, iso_obj4]
        for iso in isolist:
            iso.set_tag(tssseg_indx, tesseg_indx, "test")

        compatible_matrix = np.array(
            [[True, True, False, True],
             [True, True, False, True],
             [False, False, True, False],
             [True, True, False, True]]
        )
        overlap_matrix = np.array(
            [[True, True, True, False],
             [True, True, True, True],
             [True, True, True, True],
             [False, True, True, True]]
        )
        merged_iso_list = TransDiscover.merge_nfl_isoforms(
            np.array([0]), compatible_matrix, overlap_matrix, isolist, tssseg_indx, tesseg_indx)
        self.assertEqual(len(merged_iso_list), 1)
        np.testing.assert_array_equal(merged_iso_list[0].segary,
                                      np.array([0, 1, 0, 1, 0, 1, 0, 1, 1, 0]))

    def test_complete_iso_by_fl_iso(self):
        fl_segmat = np.array(
            [[True, True, False, True, False, True, True, False, True, False, True, True],
             [False, False, False, True, False, True, True, False, True, False, True, True],
             [True, False, True, True, False, True, True, True, False, False, False, True]])
        ival = Interval("test", 100, 200, "+")
        tssseg_indx = np.array([0, 3])
        tesseg_indx = np.array([11])
        # not have tss, not have tes
        iso_obj1 = Isoform(np.array([0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0]), ival)
        # have tss, not have tes
        iso_obj2 = Isoform(np.array([0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0]), ival)
        # have tes, not have tss
        iso_obj3 = Isoform(np.array([0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1]), ival)
        iso_obj4 = Isoform(np.array([0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1]), ival)
        # left side can not complete
        iso_obj5 = Isoform(np.array([0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0]), ival)
        # right side can not complete
        iso_obj6 = Isoform(np.array([0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0]), ival)
        # both side can not complete
        iso_obj7 = Isoform(np.array([0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0]), ival)
        # multiple complete result
        iso_obj8 = Isoform(np.array([0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0]), ival)
        isolist = [iso_obj1, iso_obj2, iso_obj3, iso_obj4, iso_obj5, iso_obj6, iso_obj7, iso_obj8]
        for iso in isolist:
            iso.set_tag(tssseg_indx, tesseg_indx, "test")

        # iso_obj1
        c_list = sorted(TransDiscover.complete_iso_by_fl_iso(
            isolist[0], fl_segmat, tssseg_indx, tesseg_indx), key=lambda x: list(x.segary))
        self.assertEqual(len(c_list), 1)
        self.assertEqual(str(c_list),
                         "[Isoform: segment array: [0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1]; tss=True; tes=True; invalid=False]")

        # iso_obj2
        c_list = sorted(TransDiscover.complete_iso_by_fl_iso(
            isolist[1], fl_segmat, tssseg_indx, tesseg_indx), key=lambda x: list(x.segary))
        self.assertEqual(len(c_list), 1)
        self.assertEqual(str(c_list),
                         "[Isoform: segment array: [0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1]; tss=True; tes=True; invalid=False]")

        # iso_obj3
        c_list = sorted(TransDiscover.complete_iso_by_fl_iso(
            isolist[2], fl_segmat, tssseg_indx, tesseg_indx), key=lambda x: list(x.segary))
        self.assertEqual(len(c_list), 1)
        self.assertEqual(str(c_list),
                         "[Isoform: segment array: [1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1]; tss=True; tes=True; invalid=False]")

        # iso_obj4
        c_list = sorted(TransDiscover.complete_iso_by_fl_iso(
            isolist[3], fl_segmat, tssseg_indx, tesseg_indx), key=lambda x: list(x.segary))
        self.assertEqual(len(c_list), 1)
        self.assertEqual(str(c_list),
                         "[Isoform: segment array: [0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1]; tss=True; tes=True; invalid=False]")

        # iso_obj5
        c_list = sorted(TransDiscover.complete_iso_by_fl_iso(
            isolist[4], fl_segmat, tssseg_indx, tesseg_indx), key=lambda x: list(x.segary))
        self.assertEqual(len(c_list), 1)
        self.assertEqual(str(c_list),
                         "[Isoform: segment array: [0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1]; tss=False; tes=True; invalid=False]")

        # iso_obj6
        c_list = sorted(TransDiscover.complete_iso_by_fl_iso(
            isolist[5], fl_segmat, tssseg_indx, tesseg_indx), key=lambda x: list(x.segary))
        self.assertEqual(len(c_list), 1)
        self.assertEqual(str(c_list),
                         "[Isoform: segment array: [0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0]; tss=True; tes=False; invalid=False]")

        # iso_obj7
        c_list = sorted(TransDiscover.complete_iso_by_fl_iso(
            isolist[6], fl_segmat, tssseg_indx, tesseg_indx), key=lambda x: list(x.segary))
        self.assertEqual(len(c_list), 1)
        self.assertEqual(str(c_list),
                         "[Isoform: segment array: [0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0]; tss=False; tes=False; invalid=False]")

        # iso_obj8
        c_list = sorted(TransDiscover.complete_iso_by_fl_iso(
            isolist[7], fl_segmat, tssseg_indx, tesseg_indx), key=lambda x: list(x.segary))
        self.assertEqual(len(c_list), 2)
        self.assertEqual(str(c_list),
                         "[Isoform: segment array: [0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1]; tss=True; tes=True; invalid=False, " +
                         "Isoform: segment array: [0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1]; tss=True; tes=True; invalid=False]")

    def test_cluster_iso(self):
        fl_segmat = np.array(
            [[True, True, False, True, False, True, True, False, True, False, True, True],
             [False, False, False, True, False, True, True, False, True, False, True, True],
             [True, False, True, True, False, True, True, True, False, False, False, True]])
        ival = Interval("test", 100, 200, "+")
        tssseg_indx = np.array([0, 3])
        tesseg_indx = np.array([11])

        iso_obj1 = Isoform(np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]), ival)
        iso_obj2 = Isoform(np.array([0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]), ival)
        iso_obj3 = Isoform(np.array([1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]), ival)
        iso_obj4 = Isoform(np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]), ival)
        iso_obj5 = Isoform(np.array([0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]), ival)
        iso_obj6 = Isoform(np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]), ival)

        isolist = [iso_obj1, iso_obj2, iso_obj3, iso_obj4, iso_obj5, iso_obj6]
        for iso in isolist:
            iso.set_tag(tssseg_indx, tesseg_indx, "test")

        # Test1
        # ++++++++++
        # ++++
        #       ++++
        test_list1 = [iso_obj1, iso_obj3, iso_obj5]
        clusters = TransDiscover.cluster_iso(test_list1)
        self.assertEqual(len(clusters), 3)

        # Test2
        # ++++++++++
        # ++++
        #  +++
        test_list2 = [iso_obj1, iso_obj3, iso_obj4]
        clusters = TransDiscover.cluster_iso(test_list2)
        self.assertEqual(len(clusters), 1)

        # Test3
        # ++++++++++
        #       ++++
        #       +++
        test_list3 = [iso_obj1, iso_obj5, iso_obj6]
        clusters = TransDiscover.cluster_iso(test_list3)
        self.assertEqual(len(clusters), 1)

        # Test4
        # ++++++++++
        #  ++++++++
        # ++++
        #  +++
        #       ++++
        #       +++
        clusters = TransDiscover.cluster_iso(isolist)
        self.assertEqual(len(clusters), 3)
        for cluster in clusters:
            self.assertEqual(len(cluster), 2)

    def test_rescue_junction(self):
        # -: intron, #: exon, X: exon error, x: intron error
        # intron :   --  --
        # iso1   : ##--##--## => ##--##--##
        # iso2   : ##xx##--## => ##--##--##
        # iso3   : ##--XX--## => ##--##--##
        # iso4   : ##xXXX--## => ##xXXX--##
        # iso5   : ##xxxX--## => ##xxxX--##
        ival = Interval("test", 100, 200, "+")
        iso1 = Isoform(np.array([1, 1, 0, 0, 1, 1, 0, 0, 1, 1]), ival)
        iso2 = Isoform(np.array([1, 1, -2, -2, 1, 1, 0, 0, 1, 1]), ival)
        iso3 = Isoform(np.array([1, 1, 0, 0, -1, -1, 0, 0, 1, 1]), ival)
        iso4 = Isoform(np.array([1, 1, -2, -1, -1, -1, 0, 0, 1, 1]), ival)
        iso5 = Isoform(np.array([1, 1, -2, -2, -2, -1, 0, 0, 1, 1]), ival)
        intron_seg = [(2, 3), (6, 7)]
        exon_seg = [(0, 1), (4, 5), (8, 9)]
        segli = [Segment("test", 100 + 10 * x, 200 + 10 * x + 10, "+") for x in range(10)]
        for i in [2, 3, 6, 7]:
            segli[i].set_spliced_seg(True)

        TransDiscover.rescue_junction(iso1, intron_seg, exon_seg, segli)
        TransDiscover.rescue_junction(iso2, intron_seg, exon_seg, segli)
        TransDiscover.rescue_junction(iso3, intron_seg, exon_seg, segli)
        TransDiscover.rescue_junction(iso4, intron_seg, exon_seg, segli)
        TransDiscover.rescue_junction(iso5, intron_seg, exon_seg, segli)

        self.assertListEqual(list(iso1.segary), [1, 1, 0, 0, 1, 1, 0, 0, 1, 1])
        self.assertListEqual(list(iso2.segary), [1, 1, 0, 0, 1, 1, 0, 0, 1, 1])
        self.assertListEqual(list(iso3.segary), [1, 1, 0, 0, 1, 1, 0, 0, 1, 1])
        self.assertListEqual(list(iso4.segary), [1, 1, -2, -1, -1, -1, 0, 0, 1, 1])
        self.assertListEqual(list(iso5.segary), [1, 1, -2, -2, -2, -1, 0, 0, 1, 1])

    def test_rescue_isoform(self):
        # -: intron, #: exon, X: exon error, x: intron error
        # isoF1: --##-#--#-#-#-##
        # isoF2: --##-#-##-#-#-##
        # TSS  : ^ ^
        # TES  :                ^
        # inv1 : --####XX#-#-#-##  One mismatch
        # inv2 : --####XX#-#X#-##  Two mismatch
        # inv3 : -#####XX#-#-#-##  No TSS
        # inv4 : --####XX#-#X#---  No TES
        # inv5 : --##-#--#-X##-##  Conflict with isoF
        # inv6 : -X####XX#-#-#-##  First index is invalid
        # inv7 : --####Xx#-#X#-##  Only isoF2 can match
        # inv8 : --####xX#-#X#-##  No isoF can match

        ival = Interval("test", 100, 200, "+")
        tssseg_indx = np.array([0, 2])
        tesseg_indx = np.array([15])
        isoF1 = Isoform(np.array([0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1]), ival)
        isoF2 = Isoform(np.array([0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1]), ival)
        inv1 = Isoform(np.array([0, 0, 1, 1, 1, 1, -2, -2, 1, 0, 1, 0, 1, 0, 1, 1]), ival)
        inv2 = Isoform(np.array([0, 0, 1, 1, 1, 1, -2, -2, 1, 0, 1, -2, 1, 0, 1, 1]), ival)
        inv3 = Isoform(np.array([0, 1, 1, 1, 1, 1, -2, -2, 1, 0, 1, 0, 1, 0, 1, 1]), ival)
        inv4 = Isoform(np.array([0, 0, 1, 1, 1, 1, -2, -2, 1, 0, 1, -2, 1, 0, 0, 0]), ival)
        inv5 = Isoform(np.array([0, 0, 1, 1, 0, 1, 0, 0, 1, 0, -2, 1, 1, 0, 1, 1]), ival)
        inv6 = Isoform(np.array([0, -2, 1, 1, 1, 1, -2, -2, 1, 0, 1, 0, 1, 0, 1, 1]), ival)
        inv7 = Isoform(np.array([0, 0, 1, 1, 1, 1, -2, -1, 1, 0, 1, 0, 1, 0, 1, 1]), ival)
        inv8 = Isoform(np.array([0, 0, 1, 1, 1, 1, -1, -2, 1, 0, 1, 0, 1, 0, 1, 1]), ival)
        isolist = [isoF1, isoF2, inv1, inv2, inv3, inv4, inv5, inv6, inv7]
        for iso in isolist:
            iso.set_tag(tssseg_indx, tesseg_indx, "test")

        # No isoF
        isoR_test1 = TransDiscover.rescue_isoform(inv1, list())
        self.assertEqual(len(isoR_test1), 0)

        # With one matched isoF
        isoR_test2 = TransDiscover.rescue_isoform(inv1, [isoF1])
        self.assertEqual(len(isoR_test2), 1)
        np.testing.assert_array_equal(
            isoR_test2[0].segary,
            np.array([0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1]))

        # With multi isoF
        isoR_test3 = TransDiscover.rescue_isoform(inv1, [isoF1, isoF2])
        self.assertEqual(len(isoR_test3), 1)
        np.testing.assert_array_equal(
            isoR_test3[0].segary,
            np.array([0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1]))

        # With multi invalid
        isoR_test4 = TransDiscover.rescue_isoform(inv2, [isoF1, isoF2])
        self.assertEqual(len(isoR_test4), 1)
        np.testing.assert_array_equal(
            isoR_test4[0].segary,
            np.array([0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1]))

        # Without TSS
        isoR_test5 = TransDiscover.rescue_isoform(inv3, [isoF1, isoF2])
        self.assertEqual(len(isoR_test5), 0)

        # Without TES
        isoR_test6 = TransDiscover.rescue_isoform(inv4, [isoF1, isoF2])
        self.assertEqual(len(isoR_test6), 0)

        # Conflict with isoF
        isoR_test7 = TransDiscover.rescue_isoform(inv5, [isoF1, isoF2])
        self.assertEqual(len(isoR_test7), 0)

        # First index is invalid
        isoR_test8 = TransDiscover.rescue_isoform(inv6, [isoF1, isoF2])
        self.assertEqual(len(isoR_test8), 0)

        # Only isoF2 can match
        isoR_test9 = TransDiscover.rescue_isoform(inv7, [isoF1, isoF2])
        self.assertEqual(len(isoR_test9), 1)
        np.testing.assert_array_equal(
            isoR_test9[0].segary,
            np.array([0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1]))

        # No isoF can match
        isoR_test10 = TransDiscover.rescue_isoform(inv8, [isoF1, isoF2])
        self.assertEqual(len(isoR_test10), 0)

    def test_invalid_iso_is_subpath(self):
        ival = Interval("test", 100, 200, "+")
        fl_iso = Isoform(np.array([0, 1, 0, 1, 0, 1, 0, 1, 1, 0]), ival)
        subpath_iso = Isoform(np.array([0, 1, 0, -1, -2, -1, 0, 1, 1, 0]), ival)
        self.assertTrue(TransDiscover.invalid_iso_is_subpath(subpath_iso, [fl_iso]))
        not_subpath_iso = Isoform(np.array([0, 1, 1, -1, -2, -1, 0, 1, 1, 0]), ival)
        self.assertFalse(TransDiscover.invalid_iso_is_subpath(not_subpath_iso, [fl_iso]))

    def test_create_intron_path(self):
        introns = [("intron", 3, 5), ("IR", 7, 7), ("intron", 9, 10)]
        self.assertEqual(TransDiscover.create_intron_path(introns, 2, 11), (1, 0, 0, 0, 1, 1, 1, 0, 0, 1))
        introns = [("intron", 3, 5), ("intron", 7, 7), ("intron", 9, 10)]
        self.assertEqual(TransDiscover.create_intron_path(introns, 2, 11), (1, 0, 0, 0, 1, 0, 1, 0, 0, 1))

    def test_enum_intron_path(self):
        introns1 = [(2, 4), (7, 7), (9, 10)]
        introns2 = [(2, 4), (2, 5), (7, 7), (9, 10), (9, 11)]
        segli1 = [Segment("test", 100 + 10 * x, 200 + 10 * x + 10, "+") for x in range(15)]
        for i in [3, 7, 9]:
            segli1[i].set_spliced_seg(True)
        segli2 = [Segment("test", 100 + 10 * x, 200 + 10 * x + 10, "+") for x in range(15)]
        for i in [3, 9]:
            segli2[i].set_spliced_seg(True)
        self.assertListEqual(sorted([x[0] for x in TransDiscover.enum_intron_path(introns1, 1, 12, segli1, [])]),
                             [(1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1)])
        self.assertListEqual(sorted([x[0] for x in TransDiscover.enum_intron_path(introns2, 1, 12, segli1, [])]),
                             sorted([(1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1),
                                     (1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1),
                                     (1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1),
                                     (1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1)
                                     ]))
        self.assertListEqual(sorted([x[0] for x in TransDiscover.enum_intron_path(introns1, 1, 12, segli2, [])]),
                             sorted([(1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1),
                                     (1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1)]))
        self.assertListEqual(sorted([x[0] for x in TransDiscover.enum_intron_path([(3, 3)], 2, 4, segli2, [])]),
                             sorted([(1, 0, 1)]))
        
        exon_dict = defaultdict(int)
        intron_dict = defaultdict(int)
        intron_dict.update({(7, 7): 100, (2, 4): 100, (9, 10): 100})
        self.assertListEqual(sorted([x[0] for x in TransDiscover.enum_intron_path(
            introns1, 1, 12,
            [Segment("test", 100 + 10 * x, 200 + 10 * x + 10, "+") for x in range(15)],
            [{"exon": exon_dict, "intron": intron_dict}])]),
            sorted([(1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1)]))
        exon_dict.update({(6, 7): 100, (7, 8): 100, (1, 2): 100, (4, 5): 100, (8, 9): 100, (9, 10): 100})
        intron_dict.update({(6, 7): 1000, (2, 3): 1000, (9, 9): 1000})
        self.assertListEqual(sorted([x[0] for x in TransDiscover.enum_intron_path(
            [(2, 4), (7, 7), (9, 10), (6, 7), (2, 3), (9, 9)], 1, 12,
            [Segment("test", 100 + 10 * x, 200 + 10 * x + 10, "+") for x in range(15)],
            [{"exon": defaultdict(int), "intron": intron_dict}])]),
                             sorted([(1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1), 
                                     (1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1), 
                                     (1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1), 
                                     (1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1) 
                                     ]))

    def test_complete_partial_isoform_left(self):
        # ^ tss, $ tes, # exon, - intron
        # index  : 012345678901234
        # TSS    : ^     ^
        # TES    :             $ $
        # Exon   : ##   ## #  ##
        #        :       #     #
        #        :            ####
        #        :             ###
        # Intron :   ---  - --
        #        :   ----   ---
        # index  : 012345678901234
        # Seg1   :        -#
        # Seg2   :      ##-#--###

        ival = Interval("test", 100, 1000, "+")
        segli = [Segment("test", 100 + 10 * x, 200 + 10 * x + 10, "+") for x in range(15)]
        for i in [3, 7, 9]:
            segli[i].set_spliced_seg(True)
        introns = [(2, 4), (2, 5), (7, 7), (9, 10), (9, 11)]
        exons = [(0, 1), (5, 6), (6, 6), (8, 8), (11, 12), (12, 12), (11, 14), (12, 14)]
        tss = [0, 6]
        tes = [12, 14]
        ta = TA(ival, segli, introns, exons, tss, tes)
        segli1 = [[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]]
        segli2 = [[0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0]]
        res1 = sorted(TransDiscover.complete_partial_isoform_left(segli1, [[]], ta)[0])
        self.assertListEqual(res1, sorted([
            [1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0]
        ]))
        res2 = sorted(TransDiscover.complete_partial_isoform_left(segli2, [[]], ta)[0])
        self.assertListEqual(res2, sorted([
            [1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0]
        ]))

        # index  : 012345678901234
        # TSS    : ^     ^
        # TES    :    $
        # Exon   : ##   ## #  ##
        #        :       #     #
        #        :            ####
        #        :             ###
        # Intron :   ---  - --
        #        :   ----   ---
        # index  : 012345678901234
        # Seg1   :        -#
        ta = TA(ival, segli, introns, exons, tss, [3])
        res3 = sorted(TransDiscover.complete_partial_isoform_left(segli1, [[]], ta)[0])
        self.assertListEqual(res3, sorted([
            [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0]
        ]))

    def test_complete_partial_isoform_right(self):
        # ^ tss, $ tes, # exon, - intron
        # index  : 012345678901234
        # TSS    : ^     ^
        # TES    :             $ $
        # Exon   : ##   ## #  ##
        #        :       #     #
        #        :            ####
        #        :             ###
        # Intron :   ---  - --
        #        :   ----   ---
        # index  : 012345678901234
        # Seg1   :        -#
        # Seg2   :      ##-#--###

        ival = Interval("test", 100, 1000, "+")
        segli = [Segment("test", 100 + 10 * x, 200 + 10 * x + 10, "+") for x in range(15)]
        for i in [3, 7, 9]:
            segli[i].set_spliced_seg(True)
        introns = [(2, 4), (2, 5), (7, 7), (9, 10), (9, 11)]
        exons = [(0, 1), (5, 6), (6, 6), (8, 8), (11, 12), (12, 12), (11, 14), (12, 14)]
        tss = [0, 6]
        tes = [12, 14]
        ta = TA(ival, segli, introns, exons, tss, tes)
        segli1 = [[0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]]
        segli2 = [[0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0]]
        res1 = sorted(TransDiscover.complete_partial_isoform_right(segli1, [[]], ta)[0])
        self.assertListEqual(res1, sorted([
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1]
        ]))
        res2 = sorted(TransDiscover.complete_partial_isoform_right(segli2, [[]], ta)[0])
        self.assertListEqual(res2, sorted([
            [0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1]
        ]))

        # index  : 012345678901234
        # TSS    :              ^
        # TES    :             $ $
        # Exon   : ##   ## #  ##
        #        :       #     #
        #        :            ####
        #        :             ###
        # Intron :   ---  - --
        #        :   ----   ---
        # index  : 012345678901234
        # Seg1   :        -#
        # Seg2   :      ##-#--###
        ta = TA(ival, segli, introns, exons, [13], tes)
        res3 = sorted(TransDiscover.complete_partial_isoform_right(segli1, [[]], ta)[0])
        self.assertListEqual(res3, sorted([
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0]
        ]))

    def test_complete_partial_isoform_error(self):
        # ^ tss, $ tes, # exon, - intron
        # index  : 012345678901234
        # Exon   : ##   ## #  ##
        #        :       #     #
        #        :            ####
        #        :             ###
        # Intron :   ---  - --
        #        :   ----   ---
        # index  : 012345678901234
        # Seg1   :      XXXX

        ival = Interval("test", 100, 1000, "+")
        segli = [Segment("test", 100 + 10 * x, 200 + 10 * x + 10, "+") for x in range(15)]
        for i in [3, 7, 9]:
            segli[i].set_spliced_seg(True)
        introns = [(2, 4), (2, 5), (7, 7), (9, 10), (9, 11)]
        exons = [(0, 1), (5, 6), (6, 6), (8, 8), (11, 12), (12, 12), (11, 14), (12, 14)]
        tss = [0, 6]
        tes = [12, 14]
        ta = TA(ival, segli, introns, exons, tss, tes)
        segli1 = [[0, 0, 0, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0]]
        res = sorted(TransDiscover.complete_partial_isoform_error(segli1, [[]], ta, 5, 8)[0])
        self.assertListEqual(res, sorted([
            [0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0]
        ]))

    def test_complete_partial_isoform(self):
        # ^ tss, $ tes, # exon, - intron
        # index  : 012345678901234
        # TSS    : ^     ^
        # TES    :             $ $
        # Exon   : ##   ## #  ##
        #        :       #     #
        #        :            ####
        #        :             ###
        # Intron :   ---  - --
        #        :   ----   ---
        # index  : 012345678901234
        # Seg1   :      XXx#

        ival = Interval("test", 100, 1000, "+")
        segli = [Segment("test", 100 + 10 * x, 200 + 10 * x + 10, "+") for x in range(15)]
        for i in [3, 7, 9]:
            segli[i].set_spliced_seg(True)
        introns = [(2, 4), (2, 5), (7, 7), (9, 10), (9, 11)]
        exons = [(0, 1), (5, 6), (6, 6), (8, 8), (11, 12), (12, 12), (11, 14), (12, 14)]
        tss = [0, 6]
        tes = [12, 14]
        ta = TA(ival, segli, introns, exons, tss, tes)
        seg = np.array([0, 0, 0, 0, 0, -1, -1, -2, 1, 0, 0, 0, 0, 0, 0])
        iso = Isoform(seg, ival)
        iso.set_tag(ta.tss_indx, ta.tes_indx, "C")
        res = TransDiscover.complete_partial_isoform(iso, ta)
        res_seg = sorted([list(x.segary) for x in res])
        self.assertListEqual(res_seg, sorted([
            [1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0],
            [1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1],
            [1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0],
            [1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1]
        ]))


class TestTgsFilterRule(unittest.TestCase):
    def setUp(self):
        self.tgs_read = make_read(0, "test", 100, ((0, 20), (3, 20), (0, 20), (3, 20), (0, 20)))

    def test_indx2ival(self):
        self.assertEqual(TgsFilterRule.indx2ival(np.array([1, 2, 3, 5, 6, 8])),
                         [(1, 3), (5, 6), (8, 8)])

    def test_is_intron(self):
        self.assertTrue(TgsFilterRule.is_intron((3, 5), [(3, 5), (6, 8)]))
        self.assertFalse(TgsFilterRule.is_intron((3, 5), [(6, 8)]))

    def test_is_internal_exon(self):
        self.assertTrue(TgsFilterRule.is_internal_exon((3, 5), [(3, 5), (6, 8)]))
        self.assertFalse(TgsFilterRule.is_internal_exon((3, 5), [(6, 8)]))

    def test_is_spliced_seg_in_internal_exon(self):
        self.assertTrue(TgsFilterRule.is_spliced_seg_in_internal_exon((3, 5), [1, 3, 6]))
        self.assertTrue(TgsFilterRule.is_spliced_seg_in_internal_exon((3, 5), [1, 4, 6]))
        self.assertTrue(TgsFilterRule.is_spliced_seg_in_internal_exon((3, 5), [1, 5, 6]))
        self.assertFalse(TgsFilterRule.is_spliced_seg_in_internal_exon((3, 5), [1, 6]))

    def test_is_skip_segment(self):
        seg = Segment("test", 120, 140, "+")
        self.assertTrue(TgsFilterRule.is_skip_segment(self.tgs_read, seg))
        seg = Segment("test", 119, 140, "+")
        self.assertFalse(TgsFilterRule.is_skip_segment(self.tgs_read, seg))
        seg = Segment("test", 120, 139, "+")
        self.assertTrue(TgsFilterRule.is_skip_segment(self.tgs_read, seg))
        seg = Segment("test", 120, 141, "+")
        self.assertFalse(TgsFilterRule.is_skip_segment(self.tgs_read, seg))

    def test_is_incl_segment(self):
        seg = Segment("test", 140, 160, "+")
        self.assertTrue(TgsFilterRule.is_incl_segment(self.tgs_read, seg))
        seg = Segment("test", 139, 160, "+")
        self.assertFalse(TgsFilterRule.is_incl_segment(self.tgs_read, seg))
        seg = Segment("test", 140, 161, "+")
        self.assertFalse(TgsFilterRule.is_incl_segment(self.tgs_read, seg))
        seg = Segment("test", 139, 160, "+")
        self.assertFalse(TgsFilterRule.is_incl_segment(self.tgs_read, seg, 'left'))
        seg = Segment("test", 140, 161, "+")
        self.assertTrue(TgsFilterRule.is_incl_segment(self.tgs_read, seg, 'left'))
        seg = Segment("test", 140, 161, "+")
        self.assertFalse(TgsFilterRule.is_incl_segment(self.tgs_read, seg, 'right'))
        seg = Segment("test", 139, 160, "+")
        self.assertTrue(TgsFilterRule.is_incl_segment(self.tgs_read, seg, 'right'))

    def test_filter_iso(self):
        seglist = [Segment("test", x, x + 20, "+") for x in range(100, 240, 20)]
        iso1 = Isoform([True, True, True, False, True, False, True], Interval("test", 100, 260, "."))
        iso2 = Isoform([True, True, True, True, True, False, True], Interval("test", 100, 260, "."))
        for indx, seg in enumerate(seglist):
            if indx % 2 == 1:
                seg.set_spliced_seg(True)
        self.assertTrue(TgsFilterRule.filter_iso(iso1, seglist, [(1, 1), (3, 3), (5, 5)], []))
        self.assertFalse(TgsFilterRule.filter_iso(iso2, seglist, [(1, 1), (3, 3), (5, 5)], []))


class TestTgsProjection2segment(unittest.TestCase):
    """ Test projection and filter """

    def setUp(self):
        self.tgs_read = make_read(0, "test", 100, ((0, 20), (3, 20), (0, 20), (3, 20), (0, 20)))
        self.seglist = [Segment("test", x[0], x[1], "+") for x in [
            (80, 100), (100, 120), (120, 140), (140, 160), (160, 180), (180, 200), (200, 220)]]
        self.intron_indx_list = [(2, 2), (4, 4)]
        self.exon_indx_list = [(1, 1), (3, 3), (5, 5)]

    def test_valid_read(self):
        iso_list = TgsTransDiscover.refine(
            [self.tgs_read], self.seglist, self.intron_indx_list, self.exon_indx_list, list(), 0.5, None)
        np.testing.assert_array_equal(iso_list[0].segary, np.array([0, 1, 0, 1, 0, 1, 0]))

    def test_gap_in_exon(self):
        # 4bp gap in the second exon
        tgs_read = make_read(0, "test", 100, ((0, 20), (3, 20), (0, 8), (3, 4), (0, 8), (3, 20), (0, 20)))
        iso_list = TgsTransDiscover.refine(
            [tgs_read], self.seglist, self.intron_indx_list, self.exon_indx_list, list(), 0.5, None)
        np.testing.assert_array_equal(iso_list[0].segary, np.array([0, 1, 0, 1, 0, 1, 0]))

    def test_micro_insertion_in_intron(self):
        # 4 bp insert in the second intron
        tgs_read = make_read(0, "test", 100, ((0, 20), (3, 8), (0, 4), (3, 8), (0, 20), (3, 20), (0, 20)))
        iso_list = TgsTransDiscover.refine(
            [tgs_read], self.seglist, self.intron_indx_list, self.exon_indx_list, list(), 0.5, None)
        np.testing.assert_array_equal(iso_list[0].segary, np.array([0, 1, 0, 1, 0, 1, 0]))
        tgs_read = make_read(0, "test", 100, ((0, 20), (3, 1), (0, 18), (3, 1), (0, 20), (3, 20), (0, 20)))
        iso_list = TgsTransDiscover.refine(
            [tgs_read], self.seglist, self.intron_indx_list, self.exon_indx_list, list(), 0.5, None)
        np.testing.assert_array_equal(iso_list[0].segary, np.array([0, 1, -1, 1, 0, 1, 0]))

    def test_no_valid_intron(self):
        intron_indx_list = [(4, 4)]
        iso_list = TgsTransDiscover.refine(
            [self.tgs_read], self.seglist, intron_indx_list, self.exon_indx_list, list(), 0.5, None)
        np.testing.assert_array_equal(iso_list[0].segary, np.array([0, 1, -2, 1, 0, 1, 0]))

    def test_exon_boundary(self):
        # left boundary of the second exon is not included
        tgs_read = make_read(0, "test", 100, ((0, 20), (3, 22), (0, 18), (3, 20), (0, 20)))
        iso_list = TgsTransDiscover.refine(
            [tgs_read], self.seglist, self.intron_indx_list, self.exon_indx_list, list(), 0.5, None)
        np.testing.assert_array_equal(iso_list[0].segary, np.array([0, 1, 0, -1, 0, 1, 0]))

        tgs_read = make_read(0, "test", 100, ((0, 18), (3, 22), (0, 20), (3, 20), (0, 20)))
        iso_list = TgsTransDiscover.refine(
            [tgs_read], self.seglist, self.intron_indx_list, self.exon_indx_list, list(), 0.5, None)
        np.testing.assert_array_equal(iso_list[0].segary, np.array([0, -1, 0, 1, 0, 1, 0]))

        tgs_read = make_read(0, "test", 100, ((0, 20), (3, 20), (0, 20), (3, 22), (0, 18)))
        iso_list = TgsTransDiscover.refine(
            [tgs_read], self.seglist, self.intron_indx_list, self.exon_indx_list, list(), 0.5, None)
        np.testing.assert_array_equal(iso_list[0].segary, np.array([0, 1, 0, 1, 0, -1, 0]))

    def test_intron_boundary(self):
        # left boundary of the second intron is not skipped
        tgs_read = make_read(0, "test", 100, ((0, 20), (3, 20), (0, 22), (3, 18), (0, 20)))
        iso_list = TgsTransDiscover.refine(
            [tgs_read], self.seglist, self.intron_indx_list, self.exon_indx_list, list(), 0.5, None)
        np.testing.assert_array_equal(iso_list[0].segary, np.array([0, 1, 0, 1, -2, 1, 0]))

        tgs_read = make_read(0, "test", 100, ((0, 20), (3, 18), (0, 22), (3, 20), (0, 20)))
        iso_list = TgsTransDiscover.refine(
            [tgs_read], self.seglist, self.intron_indx_list, self.exon_indx_list, list(), 0.5, None)
        np.testing.assert_array_equal(iso_list[0].segary, np.array([0, 1, -2, 1, 0, 1, 0]))


if __name__ == '__main__':
    unittest.main()
