from igia.utils import bed2bam, SeqFile, load_seqinfo
from igia.element import identify_element
from igia.transcript import identify_transcript
import os
import unittest


class TestAnnotation(unittest.TestCase):
    def setUp(self):
        self.ann_dir = "tests/data/ann.bed12"
        self.size = "tests/data/chrom.sizes"
        names = ('Bract', 'Cotyledon')
        ngs = [os.path.join('tests/data', "%s.sort.bam" % name) for name in names]
        self.ngs = [SeqFile(fname, 'NGS') for fname in ngs]
        self.ann = SeqFile(bed2bam(self.ann_dir, self.size, "/tmp/igia"), 'ANN')
        load_seqinfo(self.ngs)

    def test_load_ann(self):
        ann = SeqFile(bed2bam(self.ann_dir, self.size, "/tmp/igia"), 'ANN')
        ann_read_list = list(ann.bam.fetch("Chr01"))
        self.assertEqual(len(ann_read_list), 3)
        self.assertListEqual(
            sorted([(read.reference_name, read.reference_start, read.reference_end) for read in ann_read_list]),
            [('Chr01', 798113, 800959), ('Chr01', 798113, 801059), ('Chr01', 4649477, 4652385)])

    def test_identify_element_singal_ann_no_txs(self):
        gene_list = identify_element("Chr01", 4648477, 4653385, self.ngs, list(), list(), list(), self.ann)
        self.assertEqual(len(gene_list), 1)

        gene = gene_list[0]
        self.assertEqual(str(gene), "Gene: Chr01:4648477-4653385 -")
        intron_list = sorted(gene.intron_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(intron_list),
                         "[Intron: Chr01:4650078-4650291 -, Intron: Chr01:4650579-4650664 -, " +
                         "Intron: Chr01:4650748-4651282 -, Intron: Chr01:4651337-4651429 -, " +
                         "Intron: Chr01:4651463-4651566 -, Intron: Chr01:4651660-4651841 -]")
        internal_exon_list = sorted(gene.internal_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(internal_exon_list),
                         "[Exon: Chr01:4650291-4650579 -, Exon: Chr01:4650664-4650748 -, " +
                         "Exon: Chr01:4651282-4651337 -, Exon: Chr01:4651429-4651463 -, " +
                         "Exon: Chr01:4651566-4651660 -]")
        tss_exon_list = sorted(gene.tss_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(tss_exon_list), "[Exon: Chr01:4651841-4652385 -]")
        tes_exon_list = sorted(gene.tes_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(tes_exon_list), "[Exon: Chr01:4649477-4650078 -]")

    def test_identify_element_singal_ann_with_tss(self):
        gene_list = identify_element(
            "Chr01", 4648477, 4653385, self.ngs, list(), [("Chr01", 4653000, "-")], list(), self.ann)
        self.assertEqual(len(gene_list), 1)

        gene = gene_list[0]
        self.assertEqual(str(gene), "Gene: Chr01:4648477-4653385 -")
        intron_list = sorted(gene.intron_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(intron_list),
                         "[Intron: Chr01:4650078-4650291 -, Intron: Chr01:4650579-4650664 -, " +
                         "Intron: Chr01:4650748-4651282 -, Intron: Chr01:4651337-4651429 -, " +
                         "Intron: Chr01:4651463-4651566 -, Intron: Chr01:4651660-4651841 -]")
        internal_exon_list = sorted(gene.internal_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(internal_exon_list),
                         "[Exon: Chr01:4650291-4650579 -, Exon: Chr01:4650664-4650748 -, " +
                         "Exon: Chr01:4651282-4651337 -, Exon: Chr01:4651429-4651463 -, " +
                         "Exon: Chr01:4651566-4651660 -]")
        tss_exon_list = sorted(gene.tss_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(tss_exon_list), "[Exon: Chr01:4651841-4653000 -]")
        tes_exon_list = sorted(gene.tes_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(tes_exon_list), "[Exon: Chr01:4649477-4650078 -]")

    def test_identify_element_singal_ann_with_tes(self):
        gene_list = identify_element(
            "Chr01", 4648477, 4653385, self.ngs, list(), list(), [("Chr01", 4649000, "-")], self.ann)
        self.assertEqual(len(gene_list), 1)

        gene = gene_list[0]
        self.assertEqual(str(gene), "Gene: Chr01:4648477-4653385 -")
        intron_list = sorted(gene.intron_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(intron_list),
                         "[Intron: Chr01:4650078-4650291 -, Intron: Chr01:4650579-4650664 -, " +
                         "Intron: Chr01:4650748-4651282 -, Intron: Chr01:4651337-4651429 -, " +
                         "Intron: Chr01:4651463-4651566 -, Intron: Chr01:4651660-4651841 -]")
        internal_exon_list = sorted(gene.internal_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(internal_exon_list),
                         "[Exon: Chr01:4650291-4650579 -, Exon: Chr01:4650664-4650748 -, " +
                         "Exon: Chr01:4651282-4651337 -, Exon: Chr01:4651429-4651463 -, " +
                         "Exon: Chr01:4651566-4651660 -]")
        tss_exon_list = sorted(gene.tss_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(tss_exon_list), "[Exon: Chr01:4651841-4652385 -]")
        tes_exon_list = sorted(gene.tes_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(tes_exon_list), "[Exon: Chr01:4649000-4650078 -]")

    def test_identify_element_singal_ann_with_txs(self):
        gene_list = identify_element(
            "Chr01", 4648477, 4653385, self.ngs, list(), [("Chr01", 4653000, "-")], [("Chr01", 4649000, "-")], self.ann)
        self.assertEqual(len(gene_list), 1)

        gene = gene_list[0]
        self.assertEqual(str(gene), "Gene: Chr01:4648477-4653385 -")
        intron_list = sorted(gene.intron_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(intron_list),
                         "[Intron: Chr01:4650078-4650291 -, Intron: Chr01:4650579-4650664 -, " +
                         "Intron: Chr01:4650748-4651282 -, Intron: Chr01:4651337-4651429 -, " +
                         "Intron: Chr01:4651463-4651566 -, Intron: Chr01:4651660-4651841 -]")
        internal_exon_list = sorted(gene.internal_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(internal_exon_list),
                         "[Exon: Chr01:4650291-4650579 -, Exon: Chr01:4650664-4650748 -, " +
                         "Exon: Chr01:4651282-4651337 -, Exon: Chr01:4651429-4651463 -, " +
                         "Exon: Chr01:4651566-4651660 -]")
        tss_exon_list = sorted(gene.tss_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(tss_exon_list), "[Exon: Chr01:4651841-4653000 -]")
        tes_exon_list = sorted(gene.tes_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(tes_exon_list), "[Exon: Chr01:4649000-4650078 -]")

    def test_identify_element_mutil_anns(self):
        gene_list = identify_element("Chr01", 795113, 808059, self.ngs, list(), list(), list(), self.ann)
        self.assertEqual(len(gene_list), 1)

        gene = gene_list[0]
        self.assertEqual(str(gene), "Gene: Chr01:797113-802059 -")
        intron_list = sorted(gene.intron_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(intron_list), "[Intron: Chr01:799412-799981 -, Intron: Chr01:800328-800876 -]")
        internal_exon_list = sorted(gene.internal_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(internal_exon_list), "[Exon: Chr01:799981-800328 -]")
        tss_exon_list = sorted(gene.tss_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(tss_exon_list), "[Exon: Chr01:800876-800959 -, Exon: Chr01:800876-801059 -]")
        tes_exon_list = sorted(gene.tes_exon_list, key=lambda x: (x.start, x.end))
        self.assertEqual(str(tes_exon_list), "[Exon: Chr01:798113-799412 -]")

    def test_identify_transcript_singal(self):
        gene_list = identify_element("Chr01", 4648477, 4653385, self.ngs, list(), list(), list(), self.ann)
        gene = gene_list[0]
        trans = identify_transcript(gene, self.ann)
        self.assertEqual(str(trans.isoA),
            "[Isoform: segment array: [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]; tss=True; tes=True; invalid=False]")

    def test_identify_transcript_multi(self):
        gene_list = identify_element("Chr01", 795113, 808059, self.ngs, list(), list(), list(), self.ann)
        gene = gene_list[0]
        trans = identify_transcript(gene, self.ann)
        self.assertEqual(str(trans.isoA),
            "[Isoform: segment array: [1, 0, 1, 0, 1, 0]; tss=True; tes=True; invalid=False, " +
            "Isoform: segment array: [1, 0, 1, 0, 1, 1]; tss=True; tes=True; invalid=False]")