from igia.utils import Bed12
from igia.cluster import Cluster
import unittest


class TestBed12(unittest.TestCase):
    def test_iterblock(self):
        line = "test	32699	40408	test1	0	+	32699	40408	0,0,0	3	431,225,470	0,1509,7239\n"
        bed = Bed12(line)
        blocks = list(bed.iterblock())
        self.assertListEqual(blocks,
                             [(32699 + 0, 32699 + 0 + 431),
                              (32699 + 1509, 32699 + 1509 + 225),
                              (32699 + 7239, 32699 + 7239 + 470)])

class TestCluster(unittest.TestCase):
    def setUp(self):
        # test1: ++..++..++
        # test2: ++......++
        # test3:   ++..++
        # test4: --......--

        lines = ["test	0	1000	test1	0	+	0	0	0,0,0	3	200,200,200	0,400,800\n",
                 "test	0	1000	test2	0	+	0	0	0,0,0	2	200,200	0,800\n",
                 "test	0	1000	test3	0	+	0	0	0,0,0	2	200,200	200,600\n",
                 "test	0	1000	test4	0	-	0	0	0,0,0	2	200,200	0,800\n"]
        self.gene_list = [Bed12(line) for line in lines]

    def test_build_exon_block(self):
        cluster = Cluster("test")
        cluster.add_gene(self.gene_list[1])
        cluster.build_exon_block()
        self.assertEqual(cluster.exon_blocks, [(0, 200), (800, 1000)])
        cluster.add_gene(self.gene_list[0])
        cluster.build_exon_block()
        self.assertEqual(cluster.exon_blocks, [(0, 200), (400, 600), (800, 1000)])

    def test_compute_exon_overlap_len(self):
        cluster = Cluster("test")
        cluster.add_gene(self.gene_list[0])
        cluster.build_exon_block()
        self.assertEqual(cluster.compute_exon_overlap_len(self.gene_list[0]), 600)
        self.assertEqual(cluster.compute_exon_overlap_len(self.gene_list[1]), 400)
        self.assertEqual(cluster.compute_exon_overlap_len(self.gene_list[2]), 0)
        self.assertEqual(cluster.compute_exon_overlap_len(self.gene_list[3]), 0)
