#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest
from igia.linkage import *
from igia.utils import SeqFile


class TestLinkage(unittest.TestCase):
    def test_add_chr_linkage(self):
        linkage_obj = Linkage()
        linkage_obj.add_chr_linkage("test_overlap", [(100, 200), (150, 250)])
        self.assertListEqual(linkage_obj.getregions("test_overlap"), [(100, 250)])
        linkage_obj.add_chr_linkage("test_gap", [(100, 200), (300, 400)])
        self.assertListEqual(linkage_obj.getregions("test_gap"), [(100, 200), (300, 400)])
        linkage_obj.add_chr_linkage("test_edge_case1", [(100, 200), (200, 300)])
        self.assertListEqual(linkage_obj.getregions("test_edge_case1"), [(100, 300)])
        linkage_obj.add_chr_linkage("test_edge_case2", [(100, 200), (201, 300)])
        self.assertListEqual(linkage_obj.getregions("test_edge_case2"), [(100, 200), (201, 300)])
        linkage_obj.add_chr_linkage("test_continuous_adding", [(100, 200), (300, 400)])
        self.assertListEqual(linkage_obj.getregions("test_continuous_adding"), [(100, 200), (300, 400)])
        linkage_obj.add_chr_linkage("test_continuous_adding", [(150, 250), (200, 300)])
        self.assertListEqual(linkage_obj.getregions("test_continuous_adding"), [(100, 400)])

    def test_add_linkage(self):
        linkage_obj1 = Linkage()
        linkage_obj1.add_chr_linkage("test_merge_linkage", [(100, 200), (150, 250)])
        linkage_obj2 = Linkage()
        linkage_obj2.add_chr_linkage("test_merge_linkage", [(200, 300), (400, 500)])
        linkage_obj1.add_linkage(linkage_obj2)
        self.assertListEqual(linkage_obj1.getregions("test_merge_linkage"), [(100, 300), (400, 500)])


class TestFindLinkage(unittest.TestCase):
    def setUp(self):
        self.NGS_obj = SeqFile("tests/data/part_NGS.bam", "NGS")
        self.TGS_obj = SeqFile("tests/data/part_TGS.bam", "TGS")
        self.seq_list = [self.NGS_obj, self.TGS_obj]

    def test_find_linkage_NGS(self):
        lkg = find_linkage([self.NGS_obj])
        assert lkg.getregions('Chr01') == [(3901, 6101)]
        
    def test_find_linkage_TGS(self):
        lkg = find_linkage([self.TGS_obj])
        assert lkg.getregions('Chr01') == [(3501, 6201)]
        
    def test_find_linkage_merge(self):
        lkg = find_linkage(self.seq_list)
        assert lkg.getregions('Chr01') == [(3501, 6201)]


if __name__ == '__main__':
    unittest.main()
