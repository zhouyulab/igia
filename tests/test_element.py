#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest
import os
from igia.utils import SeqFile, load_seqinfo, make_read, GenomeFile
from igia.element import *
import pysam


class TestGeneLinkage(unittest.TestCase):
    def setUp(self):
        names = ('Bract', 'Cotyledon')
        ngs = [os.path.join('tests/data', "%s.sort.bam" % name) for name in names]
        tgs = [os.path.join('tests/data', "all_fixed_star.sort.bam")]
        self.ngs = [SeqFile(fname, 'NGS') for fname in ngs]
        self.tgs = [SeqFile(fname, 'TGS') for fname in tgs]
        self.glk = GeneLinkage("Chr01", 601500, 605000, self.ngs, self.tgs)
        
    def test_split2gene(self):
        genes = self.glk.split2gene()
        self.assertEqual(len(genes), 1)
        self.assertEqual(str(genes[0]), "Gene: Chr01:600840-605627 -")

    def test_split2gene_2(self):
        glk = GeneLinkage("Chr01", 850000, 870000, self.ngs, self.tgs)
        genes = glk.split2gene()
        self.assertEqual(len(genes), 2)
        self.assertEqual(str(genes[0]), "Gene: Chr01:852695-858604 +")
        self.assertEqual(str(genes[1]), "Gene: Chr01:861999-868481 +")


class TestGeneLinkageFinder(unittest.TestCase):
    def test_build_exon_overlap_cluster_with_same_cluster(self):
        # read1: ++++....++++
        # read2: ++++....+
        # read3: ++++....++++....++++
        read1 = make_read(0, "read1", 100, [(0, 100), (3, 100), (0, 100)])
        read2 = make_read(0, "read2", 100, [(0, 100), (3, 100), (0, 10)])
        read3 = make_read(0, "read3", 100, [(0, 100), (3, 100), (0, 100), (3, 100), (0, 100)])
        read_list = [read1, read2, read3]
        clusterlist = GeneLinkageFinder.build_exon_overlap_cluster(read_list)
        self.assertEqual(len(clusterlist), 1)

    def test_build_exon_overlap_cluster_with_no_overlap_iso(self):
        # read1: ++++....++++
        # read2:                     ++++....+
        # read3: ++++....++++....++++
        read1 = make_read(0, "read1", 100, [(0, 100), (3, 100), (0, 100)])
        read2 = make_read(0, "read2", 600, [(0, 100), (3, 100), (0, 10)])
        read3 = make_read(0, "read3", 100, [(0, 100), (3, 100), (0, 100), (3, 100), (0, 100)])
        read_list = [read1, read2, read3]
        clusterlist = GeneLinkageFinder.build_exon_overlap_cluster(read_list)
        self.assertEqual(len(clusterlist), 2)

    def test_build_exon_overlap_cluster_with_lnc_in_intron(self):
        # read1: ++++....................++++
        # read2:         ++++....++++
        read1 = make_read(0, "read1", 100, [(0, 100), (3, 500), (0, 100)])
        read2 = make_read(0, "read2", 300, [(0, 100), (3, 100), (0, 100)])
        read_list = [read1, read2]
        clusterlist = GeneLinkageFinder.build_exon_overlap_cluster(read_list)
        self.assertEqual(len(clusterlist), 2)

    def test_split_cluster_by_overlap_with_same_cluster(self):
        # read1: ++++....++++
        # read2: ++++....+
        # read3: ++++....++++....++++
        read1 = make_read(0, "read1", 100, [(0, 100), (3, 100), (0, 100)])
        read2 = make_read(0, "read2", 100, [(0, 100), (3, 100), (0, 10)])
        read3 = make_read(0, "read3", 100, [(0, 100), (3, 100), (0, 100), (3, 100), (0, 100)])
        cluster = [read1, read2, read3]
        clusterlist = GeneLinkageFinder.split_cluster_by_overlap(cluster)
        self.assertEqual(len(clusterlist), 1)

    def test_split_cluster_by_overlap_with_polycistron(self):
        # read1: ++++....++++
        # read2:                 ++++....+
        # read3: ++++....++++....++++....+
        read1 = make_read(0, "read1", 100, [(0, 100), (3, 100), (0, 100)])
        read2 = make_read(0, "read2", 500, [(0, 100), (3, 100), (0, 10)])
        read3 = make_read(0, "read3", 100, [(0, 100), (3, 100), (0, 100), (3, 100), (0, 100), (3, 100), (0, 10)])
        cluster = [read1, read2, read3]
        clusterlist = GeneLinkageFinder.split_cluster_by_overlap(cluster)
        self.assertEqual(len(clusterlist), 3)

    def test_build_exon_overlap_cluster_with_diff_strand(self):
        # read1: ++++....++++
        # read2: ----....-
        # read3: ++++....++++....++++
        read1 = make_read(0, "read1", 100, [(0, 100), (3, 100), (0, 100)])
        read2 = make_read(0, "read2", 100, [(0, 100), (3, 100), (0, 10)])
        read3 = make_read(0, "read3", 100, [(0, 100), (3, 100), (0, 100), (3, 100), (0, 100)])
        read2.is_reverse = True
        read_list = [read1, read2, read3]
        clusterlist = GeneLinkageFinder.build_exon_overlap_cluster(read_list)
        self.assertEqual(len(clusterlist), 2)

    def test_filter_cluster_by_strand_with_same_cov(self):
        # read1: ++++....++++
        # read2: ----....-
        read1 = make_read(0, "read1", 100, [(0, 100), (3, 100), (0, 100)])
        read2 = make_read(0, "read2", 100, [(0, 100), (3, 100), (0, 10)])
        read2.is_reverse = True
        read1.setTag("FL", 100)
        read2.setTag("FL", 100)
        clusterlist = [[read1], [read2]]
        res = GeneLinkageFinder.filter_cluster_by_strand(clusterlist)
        self.assertEqual(len(res), 2)

    def test_filter_cluster_by_strand_almost_pos(self):
        # read1: ++++....++++
        # read2: ----....-
        read1 = make_read(0, "read1", 100, [(0, 100), (3, 100), (0, 100)])
        read2 = make_read(0, "read2", 100, [(0, 100), (3, 100), (0, 10)])
        read2.is_reverse = True
        read1.setTag("FL", 1000)
        read2.setTag("FL", 100)
        clusterlist = [[read1], [read2]]
        res = GeneLinkageFinder.filter_cluster_by_strand(clusterlist)
        self.assertEqual(len(res), 1)
        self.assertFalse(res[0][0].is_reverse)

    def test_filter_cluster_by_strand_almost_neg(self):
        # read1: ++++....++++
        # read2: ----....-
        read1 = make_read(0, "read1", 100, [(0, 100), (3, 100), (0, 100)])
        read2 = make_read(0, "read2", 100, [(0, 100), (3, 100), (0, 10)])
        read2.is_reverse = True
        read1.setTag("FL", 100)
        read2.setTag("FL", 1000)
        clusterlist = [[read1], [read2]]
        res = GeneLinkageFinder.filter_cluster_by_strand(clusterlist)
        self.assertEqual(len(res), 1)
        self.assertTrue(res[0][0].is_reverse)

    def test_filter_cluster_by_strand_one_without_intron(self):
        # read1: ++++....++++
        # read2: ----
        read1 = make_read(0, "read1", 100, [(0, 100), (3, 100), (0, 100)])
        read2 = make_read(0, "read2", 100, [(0, 100)])
        read2.is_reverse = True
        read1.setTag("FL", 1000)
        read2.setTag("FL", 100)
        clusterlist = [[read1], [read2]]
        res = GeneLinkageFinder.filter_cluster_by_strand(clusterlist)
        self.assertEqual(len(res), 1)
        self.assertFalse(res[0][0].is_reverse)

    def test_filter_cluster_by_strand_one_without_intron_2(self):
        # read1: ++++....++++
        # read2:      --
        read1 = make_read(0, "read1", 100, [(0, 100), (3, 100), (0, 100)])
        read2 = make_read(0, "read2", 225, [(0, 50)])
        read2.is_reverse = True
        read1.setTag("FL", 1000)
        read2.setTag("FL", 100)
        clusterlist = [[read1], [read2]]
        res = GeneLinkageFinder.filter_cluster_by_strand(clusterlist)
        self.assertEqual(len(res), 2)

    def test_filter_cluster_by_strand_both_without_intron_1(self):
        # read1: ++++
        # read2: ----
        read1 = make_read(0, "read1", 100, [(0, 100)])
        read2 = make_read(0, "read2", 100, [(0, 100)])
        read2.is_reverse = True
        read1.setTag("FL", 1000)
        read2.setTag("FL", 100)
        clusterlist = [[read1], [read2]]
        res = GeneLinkageFinder.filter_cluster_by_strand(clusterlist)
        self.assertEqual(len(res), 1)
        self.assertFalse(res[0][0].is_reverse)

    def test_filter_cluster_by_strand_both_without_intron_2(self):
        # read1: ++++
        # read2:     ----
        read1 = make_read(0, "read1", 100, [(0, 100)])
        read2 = make_read(0, "read2", 200, [(0, 100)])
        read2.is_reverse = True
        read1.setTag("FL", 1000)
        read2.setTag("FL", 100)
        clusterlist = [[read1], [read2]]
        res = GeneLinkageFinder.filter_cluster_by_strand(clusterlist)
        self.assertEqual(len(res), 2)

    def test_filter_cluster_by_strand_with_diff_intron(self):
        # read1: ++++....++++
        # read2:   ----....-
        read1 = make_read(0, "read1", 100, [(0, 100), (3, 100), (0, 100)])
        read2 = make_read(0, "read2", 150, [(0, 100), (3, 100), (0, 10)])
        read2.is_reverse = True
        read1.setTag("FL", 100)
        read2.setTag("FL", 1000)
        clusterlist = [[read1], [read2]]
        res = GeneLinkageFinder.filter_cluster_by_strand(clusterlist)
        self.assertEqual(len(res), 2)


class TestJunctionGraph(unittest.TestCase):
    def setUp(self):
        names = ('Bract', 'Cotyledon')
        ngs = [os.path.join('tests/data', "%s.sort.bam" % name) for name in names]
        tgs = [os.path.join('tests/data', "all_fixed_star.sort.bam")]
        self.ngs = [SeqFile(fname, 'NGS') for fname in ngs]
        self.tgs = [SeqFile(fname, 'TGS') for fname in tgs]
        self.ival_fwd = Interval("Chr01", 485500, 493000, "+")
        self.ival_fwd.build_cov(self.ngs)
        self.ival_rev = Interval("Chr01", 601500, 605000, "-")
        self.ival_rev.build_cov(self.ngs)
        
    def test_identify_intron_NGS_rev(self):
        jgn = JunctionGraphNgs(self.ival_rev)
        introns = jgn.identify_intron(self.ngs)
        self.assertEqual(len(introns), 1)
        self.assertEqual(str(introns[0]), "Intron: Chr01:603134-603953 -")

    def test_identify_intron_TGS_rev(self):
        jgn = JunctionGraphTgs(self.ival_rev)
        tgs_read_list = self.tgs[0].fetch_reads_in_ival(self.ival_rev)
        introns = sorted(jgn.identify_intron(tgs_read_list), key=lambda x: (x.start, x.end))
        self.assertEqual(len(introns), 1)
        text_introns = "[Intron: Chr01:603134-603953 -]"
        self.assertEqual(str(introns), text_introns)

    def test_identify_intron_NGS_fwd(self):
        jgn = JunctionGraphNgs(self.ival_fwd)
        introns = sorted(jgn.identify_intron(self.ngs), key=lambda x: (x.start, x.end))
        self.assertEqual(len(introns), 6)
        text_introns = "[Intron: Chr01:486149-488629 +, Intron: Chr01:488687-488793 +, " + \
                       "Intron: Chr01:488874-488989 +, Intron: Chr01:489025-489151 +, " + \
                       "Intron: Chr01:489388-489475 +, Intron: Chr01:489617-490544 +]"
        self.assertEqual(str(introns), text_introns)

    def test_identify_intron_TGS_fwd(self):
        jgn = JunctionGraphTgs(self.ival_fwd)
        tgs_read_list = self.tgs[0].fetch_reads_in_ival(self.ival_fwd)
        introns = sorted(jgn.identify_intron(tgs_read_list), key=lambda x: (x.start, x.end))
        self.assertEqual(len(introns), 6)
        text_introns = "[Intron: Chr01:486149-488629 +, Intron: Chr01:488687-488793 +, " + \
                       "Intron: Chr01:488874-488989 +, Intron: Chr01:489025-489151 +, " + \
                       "Intron: Chr01:489388-489475 +, Intron: Chr01:489617-490544 +]"
        self.assertEqual(str(introns), text_introns)

    def test_find_junction(self):
        read = make_read(0, "test", 100,
                         [(0, 20), (3, 20), (0, 20), (1, 1), (3, 20),
                          (0, 20), (2, 10), (0, 5), (3, 20), (0, 20)])
        junctions = JunctionGraph._find_junction(read, 10)
        self.assertEqual(len(junctions), 1)
        self.assertEqual(junctions[0], (("chr1", 120, '+', 'EI'), ("chr1", 140, '+', 'IE')))


class TestIntron(unittest.TestCase):
    def setUp(self):
        names = ('Bract', 'Cotyledon')
        ngs = [os.path.join('tests/data', "%s.sort.bam" % name) for name in names]
        self.ngs = [SeqFile(fname, 'NGS') for fname in ngs]
        self.ival_rev = Interval("Chr01", 93000, 96500, "-")
        self.ival_rev.build_cov(self.ngs)
        self.ival_fwd = Interval("Chr01", 485500, 493000, "+")
        self.ival_fwd.build_cov(self.ngs)

    def test_identify_intron_NGS_rev(self):
        jgn = JunctionGraphNgs(self.ival_rev)
        introns = jgn.identify_intron(self.ngs)
        introns = sorted(introns, key=lambda x: (x.start, x.end))
        self.assertTrue(introns[3].spliced)
        self.assertFalse(introns[5].spliced)

    def test_identify_intron_NGS_fwd(self):
        jgn = JunctionGraphNgs(self.ival_fwd)
        introns = jgn.identify_intron(self.ngs)
        introns = sorted(introns, key=lambda x: (x.start, x.end))
        self.assertTrue(introns[3].spliced)
        self.assertTrue(introns[1].spliced)


class TestElementDiscover(unittest.TestCase):
    def setUp(self):
        names = ('Bract', 'Cotyledon')
        ngs = [os.path.join('tests/data', "%s.sort.bam" % name) for name in names]
        self.ngs = [SeqFile(fname, 'NGS') for fname in ngs]
        load_seqinfo(self.ngs)
        self.tgs = [SeqFile('tests/data/all_fixed_star.sort.bam', 'TGS')]
        self.ival_rev = Interval("Chr01", 93000, 96500, "-")
        self.ival_rev.build_cov(self.ngs)
        jgn_rev = JunctionGraphNgs(self.ival_rev)
        self.rev_intron = jgn_rev.identify_intron(self.ngs)
        self.ival_fwd = Interval("Chr01", 485500, 493000, "+")
        self.ival_fwd.build_cov(self.ngs)
        jgn_fwd = JunctionGraphNgs(self.ival_fwd)
        self.fwd_intron = jgn_fwd.identify_intron(self.ngs)
        self.genome = GenomeFile(os.path.join("tests", "data", "genome.fa"))

    def test_enumerate_exon(self):
        it1 = Intron('chr1', 100, 200, '+')
        it2 = Intron('chr1', 300, 400, '+')
        it3 = Intron('chr1', 500, 600, '+')
        it4 = Intron('chr1', 350, 400, '+')
        self.assertEqual(str(ElementDiscover.enumerate_exon([it1, it2])[0]), '[(200, 300)]')
        self.assertEqual(str(ElementDiscover.enumerate_exon([it1, it2, it3])[0]),
            '[(200, 300), (200, 500), (400, 500)]')
        self.assertEqual(str(ElementDiscover.enumerate_exon([it1, it2, it3, it4])[0]),
            '[(200, 300), (200, 350), (200, 500), (400, 500)]')
        it4.set_spliced(True)
        self.assertEqual(str(ElementDiscover.enumerate_exon([it1, it2, it3, it4])[0]),
            '[(200, 300), (200, 350), (400, 500)]')

        it11 = Intron('chr1', 50, 150, '+')
        it21 = Intron('chr1', 180, 190, '+')
        for it in [it1, it2, it11, it21]:
            it.set_spliced(True)
        self.assertEqual(str(ElementDiscover.enumerate_exon([it1, it2, it11, it21])[0]),
            '[(150, 180), (200, 300)]')

        # No internal exon
        self.assertEqual(str(ElementDiscover.enumerate_exon([it1])[0]), '[]')
        self.assertEqual(str(ElementDiscover.enumerate_exon([])[0]), '[]')

    def test_identify_internal_exon_rev(self):
        exons, _, _ = ElementDiscover.identify_internal_exon(self.rev_intron, self.ival_rev)
        self.assertEqual(str(exons),
            '[Exon: Chr01:94263-94877 -, Exon: Chr01:94959-95353 -, Exon: Chr01:95459-95527 -, Exon: Chr01:95801-95874 -]')

    def test_identify_internal_exon_fwd(self):
        exons, _, _ = ElementDiscover.identify_internal_exon(self.fwd_intron, self.ival_fwd)
        self.assertEqual(str(exons),
            "[Exon: Chr01:488629-488687 +, Exon: Chr01:488793-488874 +, Exon: Chr01:488989-489025 +, " +
            "Exon: Chr01:489151-489388 +, Exon: Chr01:489475-489617 +]")

    def test_identify_tss_exon_rev(self):
        tss_site_list = [95907, 96077, 96200]
        tes_site_list = [93500, 93802]
        tss_exon = ElementDiscover.identify_tss_exon(tss_site_list, tes_site_list, self.rev_intron, self.ival_rev)
        self.assertEqual(str(tss_exon),
            "[Exon: Chr01:95801-95907 -, Exon: Chr01:95801-96077 -, Exon: Chr01:96076-96200 -]")

    def test_identify_tss_exon_fwd(self):
        tss_site_list = [485839]
        tes_site_list = [491380]
        tss_exon = ElementDiscover.identify_tss_exon(tss_site_list, tes_site_list, self.fwd_intron, self.ival_fwd)
        self.assertEqual(str(tss_exon), "[Exon: Chr01:485839-486149 +]")

    def test_identify_tes_exon_rev(self):
        tss_site_list = [95907, 96077, 96200]
        tes_site_list = [93500, 93802]
        tes_exon = ElementDiscover.identify_tes_exon(tss_site_list, tes_site_list, self.rev_intron, self.ival_rev)
        self.assertEqual(str(tes_exon),
            "[Exon: Chr01:93500-94148 -, Exon: Chr01:93802-94148 -]")

    def test_identify_tes_exon_fwd(self):
        tss_site_list = [485839]
        tes_site_list = [491380]
        tes_exon = ElementDiscover.identify_tes_exon(tss_site_list, tes_site_list, self.fwd_intron, self.ival_fwd)
        self.assertEqual(str(tes_exon), "[Exon: Chr01:490544-491380 +]")

    def test_identify_tss_tes_exon(self):
        tss_site_list = [485839]
        tes_site_list = [491380]
        tes_exon = ElementDiscover.identify_tes_exon(tss_site_list, tes_site_list, list(), self.ival_fwd)
        self.assertEqual(str(tes_exon), "[Exon: Chr01:485839-491380 +]")

    def test_fetch_neighbor_seq(self):
        pos_ival = Interval("Chr01", 100, 200, "+")
        neighbor_seq = ElementDiscover.fetch_neighbor_seq(90, 210, pos_ival, self.genome)
        self.assertEqual(neighbor_seq, "GTATGTGCAATTGTGTTATG")

        neg_ival = Interval("Chr01", 100, 200, "-")
        neighbor_seq = ElementDiscover.fetch_neighbor_seq(90, 210, neg_ival, self.genome)
        self.assertEqual(neighbor_seq, "GTATGTGCAATTGTGTTATG")

    def test_fetch_splice_site(self):
        pos_intron = Intron("Chr01", 100, 150, "+")
        ss = ElementDiscover.fetch_splice_site(pos_intron, self.genome)
        self.assertEqual(ss, ("GC", "TC"))

        neg_intron = Intron("Chr01", 100, 150, "-")
        ss = ElementDiscover.fetch_splice_site(neg_intron, self.genome)
        self.assertEqual(ss, ("GA", "GC"))

    def test_ss2pri(self):
        self.assertEqual(ElementDiscover.ss2pri(("GT", "AG")), 3)
        self.assertEqual(ElementDiscover.ss2pri(("GC", "AG")), 2)
        self.assertEqual(ElementDiscover.ss2pri(("AT", "AC")), 1)
        self.assertEqual(ElementDiscover.ss2pri(("CT", "AG")), 0)

    def test_blocks2cigar(self):
        blocks = ((0, 20), (3, 80), (0, 75), (3, 42), (0, 30))
        read = make_read(0, "test", 100, blocks)
        self.assertEqual(ElementDiscover.blocks2cigar(read.blocks), blocks)

    def test_adjust_intron_position(self):
        read = make_read(0, "test", 100, [(0, 100), (3, 100), (0, 100)])
        intron = Intron(read.reference_name, 190, 310, "+")
        ElementDiscover.adjust_intron_position(read, intron)
        self.assertListEqual(read.blocks, [(100, 190), (310, 400)])

        read = make_read(0, "test", 100, [(0, 100), (3, 100), (0, 100)])
        intron = Intron(read.reference_name, 210, 290, "+")
        ElementDiscover.adjust_intron_position(read, intron)
        self.assertListEqual(read.blocks, [(100, 210), (290, 400)])

        # Difference chrom
        read = make_read(0, "test", 100, [(0, 100), (3, 100), (0, 100)])
        intron = Intron("DiffChrom", 190, 310, "+")
        ElementDiscover.adjust_intron_position(read, intron)
        self.assertListEqual(read.blocks, [(100, 200), (300, 400)])

        # Difference strand
        read = make_read(0, "test", 100, [(0, 100), (3, 100), (0, 100)])
        intron = Intron("DiffChrom", 190, 310, "-")
        ElementDiscover.adjust_intron_position(read, intron)
        self.assertListEqual(read.blocks, [(100, 200), (300, 400)])

    def test_fix_mapping_error(self):
        # pos:
        #  AG(GT...-- AG)
        # (AG --...GT)AG
        pos_intron1 = Intron("Chr01", 4*50+30, 13*50+31, "+")
        ss1 = ElementDiscover.fetch_splice_site(pos_intron1, self.genome)
        self.assertEqual(ss1, ("GT", "AG"))
        pos_intron2 = Intron("Chr01", 4*50+28, 13*50+29, "+")
        ss2 = ElementDiscover.fetch_splice_site(pos_intron2, self.genome)
        self.assertEqual(ss2, ("AG", "GT"))
        intron_list = ElementDiscover.fix_mapping_error([pos_intron1], [pos_intron2], self.genome, [], [])
        self.assertEqual(len(intron_list), 1)
        self.assertEqual(str(intron_list[0]), "Intron: Chr01:230-681 +")
        intron_list = ElementDiscover.fix_mapping_error([pos_intron2], [pos_intron1], self.genome, [], [])
        self.assertEqual(len(intron_list), 1)
        self.assertEqual(str(intron_list[0]), "Intron: Chr01:230-681 +")

        # neg:
        #  AC(CT...-- AC)
        # (AC --...CT)AC
        neg_intron1 = Intron("Chr01", 3 * 50 + 20, 8 * 50 + 11, "-")
        ss1 = ElementDiscover.fetch_splice_site(neg_intron1, self.genome)
        self.assertEqual(ss1, ("GT", "AG"))
        neg_intron2 = Intron("Chr01", 3 * 50 + 18, 8 * 50 + 9, "-")
        ss2 = ElementDiscover.fetch_splice_site(neg_intron2, self.genome)
        self.assertEqual(ss2, ("AG", "GT"))
        intron_list = ElementDiscover.fix_mapping_error([neg_intron1], [neg_intron2], self.genome, [], [])
        self.assertEqual(len(intron_list), 1)
        self.assertEqual(str(intron_list[0]), "Intron: Chr01:170-411 -")
        intron_list = ElementDiscover.fix_mapping_error([neg_intron2], [neg_intron1], self.genome, [], [])
        self.assertEqual(len(intron_list), 1)
        self.assertEqual(str(intron_list[0]), "Intron: Chr01:170-411 -")

class TestIdentifyElement(unittest.TestCase):
    def setUp(self):
        names = ('Bract', 'Cotyledon')
        ngs = [os.path.join('tests/data', "%s.sort.bam" % name) for name in names]
        tgs = [os.path.join('tests/data', "all_fixed_star.sort.bam")]
        self.ngs = [SeqFile(fname, 'NGS') for fname in ngs]
        self.tgs = [SeqFile(fname, 'TGS') for fname in tgs]

    def test_rev(self):
        tss_site_list = [("Chr01", x, "-") for x in [95907, 96077]]
        tes_site_list = [("Chr01", x, "-") for x in [93802]]
        gene_list = identify_element("Chr01", 90000, 98000, self.ngs, self.tgs, tss_site_list, tes_site_list)
        self.assertEqual(len(gene_list), 1)
        
        gene = gene_list[0]
        intron_list = sorted(gene.intron_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(intron_list),
                         "[Intron: Chr01:94148-94263 -, Intron: Chr01:94877-94959 -, Intron: Chr01:95353-95459 -, " +
                         "Intron: Chr01:95527-95801 -, Intron: Chr01:95527-96076 -, Intron: Chr01:95874-96076 -]")
        internal_exon_list = sorted(gene.internal_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(internal_exon_list),
                         "[Exon: Chr01:94263-94877 -, Exon: Chr01:94959-95353 -, " +
                         "Exon: Chr01:95459-95527 -, Exon: Chr01:95801-95874 -]")
        tss_exon_list = sorted(gene.tss_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(tss_exon_list),
                         "[Exon: Chr01:95801-95907 -, Exon: Chr01:95801-96077 -]")
        tes_exon_list = sorted(gene.tes_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(tes_exon_list), "[Exon: Chr01:93802-94148 -]")

    def test_fwd(self):
        tss_site_list = [("Chr01", x, "+") for x in [485839]]
        tes_site_list = [("Chr01", x, "+") for x in [491380]]
        gene_list = identify_element("Chr01", 485500, 493000, self.ngs, self.tgs, tss_site_list, tes_site_list)
        self.assertEqual(len(gene_list), 1)

        gene = gene_list[0]
        intron_list = sorted(gene.intron_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(intron_list),
                         "[Intron: Chr01:486149-488629 +, Intron: Chr01:488687-488793 +, " +
                         "Intron: Chr01:488874-488989 +, Intron: Chr01:489025-489151 +, " +
                         "Intron: Chr01:489388-489475 +, Intron: Chr01:489617-490544 +]")
        internal_exon_list = sorted(gene.internal_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(internal_exon_list),
                         "[Exon: Chr01:488629-488687 +, Exon: Chr01:488793-488874 +, " +
                         "Exon: Chr01:488989-489025 +, Exon: Chr01:489151-489388 +, " +
                         "Exon: Chr01:489475-489617 +]")
        tss_exon_list = sorted(gene.tss_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(tss_exon_list), "[Exon: Chr01:485839-486149 +]")
        tes_exon_list = sorted(gene.tes_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(tes_exon_list),
                         "[Exon: Chr01:490544-491380 +, Exon: Chr01:490544-492288 +]")


if __name__ == '__main__':
    unittest.main()
