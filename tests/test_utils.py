#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import unittest
from igia.utils import *
import numpy as np

class TestSeqFile(unittest.TestCase):
    def setUp(self):
        self.NGS_obj = SeqFile("tests/data/part_NGS.bam", "NGS")
        self.TGS_obj = SeqFile("tests/data/part_TGS.bam", "TGS")
        self.TGS2 = SeqFile("tests/data/all_fixed_star.sort.bam", "TGS")
        self.seq_list = [self.NGS_obj, self.TGS_obj]
        self.ival = Interval("Chr01", 601500, 605000, '.')

    def test_chromsize(self):
        assert self.NGS_obj.chromsize() == [('Chr01', 7000), ('Chr02', 8000)]
        assert self.TGS_obj.chromsize() == [('Chr01', 7000), ('Chr02', 8000)]

    def test_genomesize(self):
        assert self.NGS_obj.genomesize() == 15000
        assert self.TGS_obj.genomesize() == 15000

    def test_mapped_number(self):
        assert self.NGS_obj.mapped_number() == 982
        assert self.TGS_obj.mapped_number() == 74

    def test_readlen(self):
        assert self.NGS_obj.readlen() == 150

    def test_poisbg(self):
        assert self.NGS_obj.poisbg(100) == 13
        assert self.TGS_obj.poisbg(100) == 0

    def test_count(self):
        assert self.TGS2.count("Chr01", 601500, 605000) == 18

    def test_smart_fetch(self):
        assert len(list(self.TGS2.smart_fetch("Chr01", 601500, 605000))) == 18

    def test_fetch_reads_in_ival(self):
        assert len(self.TGS2.fetch_reads_in_ival(self.ival)) == 18


class TestCoverage(unittest.TestCase):
    def setUp(self):
        names = ('Bract', 'Cotyledon')
        ngs = [os.path.join('tests/data', "%s.sort.bam" % name) for name in names]
        self.ngs = [SeqFile(fname, 'NGS') for fname in ngs]
        self.cov_fwd = Coverage(Interval("Chr01", 489388-5, 489388+5, "+"))
        self.cov_rev = Coverage(Interval("Chr01", 603953-5, 603953+5, "-"))

    def test_build_fwd(self):
        self.cov_fwd.build(self.ngs[:1])
        np.testing.assert_array_equal(self.cov_fwd.sig,
                                      np.array([16, 16, 19, 19, 19, 3, 3, 3, 3, 5, 5]).reshape([1, 11]))

    def test_build_rev(self):
        self.cov_rev.build(self.ngs[:1])
        np.testing.assert_array_equal(self.cov_rev.sig,
                                      np.array([1, 1, 1, 4, 8, 47, 47, 47, 46, 47, 47]).reshape([1, 11]))

    def test_build_multi_ngs(self):
        self.cov_fwd.build(self.ngs)
        np.testing.assert_array_equal(self.cov_fwd.sig,
                                      np.array([[16, 16, 19, 19, 19, 3, 3, 3, 3, 5, 5],
                                                [19, 19, 19, 19, 19, 2, 2, 2, 2, 2, 2]]))

    def test_slice(self):
        self.cov_fwd.build(self.ngs)
        newcov = self.cov_fwd.slice(Interval("Chr01", 489388-3, 489388+3, "+"))
        np.testing.assert_array_equal(newcov.sig,
                                      np.array([[19, 19, 19, 3, 3, 3, 3],
                                                [19, 19, 19, 2, 2, 2, 2]]))

    def test_getitem(self):
        self.cov_fwd.build(self.ngs)
        np.testing.assert_array_equal(self.cov_fwd[489388], np.array([3, 2]))

    def test_len(self):
        self.assertEqual(len(self.cov_fwd), 11)


class TestInterval(unittest.TestCase):
    def setUp(self):
        names = ('Bract', 'Cotyledon')
        ngs = [os.path.join('tests/data', "%s.sort.bam" % name) for name in names]
        self.ngs = [SeqFile(fname, 'NGS') for fname in ngs]
        self.iv_fwd = Interval("Chr01", 489388-5, 489388+5, "+")
        self.iv_rev = Interval("Chr01", 603953-5, 603953+5, "-")


    def test_build_fwd(self):
        self.iv_fwd.build_cov(self.ngs[:1])
        np.testing.assert_array_equal(self.iv_fwd.cov.sig,
                                      np.array([16, 16, 19, 19, 19, 3, 3, 3, 3, 5, 5]).reshape([1, 11]))
        
    def test_build_rev(self):
        self.iv_rev.build_cov(self.ngs[:1])
        np.testing.assert_array_equal(self.iv_rev.cov.sig,
                                      np.array([1, 1, 1, 4, 8, 47, 47, 47, 46, 47, 47]).reshape([1, 11]))

    def test_compute_fpkm(self):
        self.iv_fwd.build_cov(self.ngs)
        self.iv_fwd.compute_fpkm([150, 150], [1e8, 1e8])
        np.testing.assert_almost_equal(self.iv_fwd.fpkm,
                                       np.array([[0.0017651, 0.4447810, 3.9789872],
                                                 [0.0013804, 0.4198411, 3.9104092]]))

    def test_len(self):
        self.assertEqual(len(self.iv_fwd), 10)


if __name__ == '__main__':
    unittest.main()
