from igia import __version__
from igia.utils import SeqFile, load_seqinfo, load_txs, load_ann
from igia import GVAR
from igia.linkage import Linkage
from igia.element import identify_element
from igia.transcript import identify_transcript
from igia.skeleton import check_paraclu
import argparse
import sys
import os
import logging
import time
import signal

_logger = logging.getLogger(__name__)


def parse_args(args):
    """Parse command line parameters
    Args:
      args ([str]): command line parameters as list of strings
    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(
        description="Integrative Gene Isoform Assembler")
    parser.add_argument(
        '--version',
        action='version',
        version='igia {ver}'.format(ver=__version__))
    parser.add_argument(
        '-v',
        '--verbose',
        dest="loglevel",
        help="set loglevel to INFO",
        action='store_const',
        const=logging.INFO)
    parser.add_argument(
        '-vv',
        '--very-verbose',
        dest="loglevel",
        help="set loglevel to DEBUG",
        action='store_const',
        const=logging.DEBUG)

    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-o", "--output", type=str, dest="out_dir", metavar="out_dir", required=True)
    base_group.add_argument("--ngs", nargs="*", type=str, dest="ngs_file", metavar="NGS.bam", required=True)
    base_group.add_argument("--tgs", nargs="*", type=str, dest="tgs_file", metavar="TGS.bam", required=True)
    base_group.add_argument("--chrom", type=str, dest="chrom", metavar="chrom", required=True)
    base_group.add_argument("--start", type=int, dest="start", metavar="start", required=True)
    base_group.add_argument("--end", type=int, dest="end", metavar="end", required=True)

    ext_group = parser.add_argument_group("External")
    ext_group.add_argument("--tss", type=str, dest="tss", metavar="tss.csv", default="")
    ext_group.add_argument("--tes", type=str, dest="tes", metavar="tes.csv", default="")
    ext_group.add_argument("--ann", type=str, dest="ann", metavar="NGS_ann.bed12", default="",
                           help="Bed12 format NGS based annotation, needs --size for chromsize")
    ext_group.add_argument("--cfm-ann", type=str, dest="cfm_ann", metavar="comfirmed_ann.bed12", default="",
                           help="Bed12 format comfirmed annotation, needs --size for chromsize")
    ext_group.add_argument("-s", "--size", type=str, dest="size", metavar="chrom.size", default="")
    ext_group.add_argument("-g", "--genome", type=str, dest="f_genome", metavar="genome.fa", default=None)

    opt_group = parser.add_argument_group("Options")
    opt_group.add_argument("-r", type=str, dest="rule", metavar="rule", default="single_end",
                           help="NGS library type", choices=["1++,1--,2+-,2-+", "1+-,1-+,2++,2--", "single_end"])
    opt_group.add_argument("--pir", type=float, dest="pir", metavar="pir", default=0.5,
                           help="PIR cutoff for intron retention")
    opt_group.add_argument("--dtxs", type=int, dest="dtxs", metavar="dtxs", default=500,
                           help="Distance cutoff between two different TSS/TES")
    opt_group.add_argument("--time-out", type=int, dest="time_out", metavar="time_out", default=30*60,
                           help="TimeOut")
    opt_group.add_argument("--paraclu-path", type=str, dest="paraclu_path", metavar="paraclu_path", default=None,
                           help="paraclu_path")
    return parser.parse_args(args)


def setup_logging(loglevel):
    """Setup basic logging
    Args:
      loglevel (int): minimum loglevel for emitting messages
    """
    logformat = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(level=loglevel, stream=sys.stdout,
                        format=logformat, datefmt="%Y-%m-%d %H:%M:%S")


class OutputHandle(object):
    def __init__(self, outdir):
        self.outdir = outdir

    def __enter__(self):
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        self.f_intron = open(os.path.join(self.outdir, "intron.bed6"), "w")
        self.f_internal_exon = open(os.path.join(self.outdir, "internal_exon.bed6"), "w")
        self.f_tss_exon = open(os.path.join(self.outdir, "tss_exon.bed6"), "w")
        self.f_tes_exon = open(os.path.join(self.outdir, "tes_exon.bed6"), "w")
        self.f_isoF = open(os.path.join(self.outdir, "isoF.bed12"), "w")
        self.f_isoA = open(os.path.join(self.outdir, "isoA.bed12"), "w")
        self.f_isoR = open(os.path.join(self.outdir, "isoR.bed12"), "w")
        self.f_isoM = open(os.path.join(self.outdir, "isoM.bed12"), "w")
        self.f_isoC = open(os.path.join(self.outdir, "isoC.bed12"), "w")
        self.f_isoP = open(os.path.join(self.outdir, "isoP.bed12"), "w")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def element_handles(self):
        return self.f_intron, self.f_internal_exon, self.f_tss_exon, self.f_tes_exon

    def isoform_handles(self):
        return self.f_isoF, self.f_isoA, self.f_isoR, self.f_isoM, self.f_isoC, self.f_isoP

    def close(self):
        self.f_intron.close()
        self.f_internal_exon.close()
        self.f_tss_exon.close()
        self.f_tes_exon.close()
        self.f_isoF.close()
        self.f_isoA.close()
        self.f_isoR.close()
        self.f_isoM.close()
        self.f_isoC.close()
        self.f_isoP.close()


class TimeOutError(RuntimeError):
    pass

def time_out_handler(signum, frame):
    raise TimeOutError

def main(args):
    """Main entry point allowing external calls
    Args:
      args ([str]): command line parameter list
    """
    args = parse_args(args)
    check_paraclu(args)
    setup_logging(args.loglevel)
    _logger.debug("Starting IGIA ...")

    ngs_obj_list = [SeqFile(x, "NGS") for x in args.ngs_file]
    tgs_obj_list = [SeqFile(x, "TGS") for x in args.tgs_file]
    ext_tss_list = load_txs(args.tss)
    ext_tes_list = load_txs(args.tes)

    out_dir = args.out_dir
    ann = load_ann(args.ann, args.size, out_dir, "ANN")

    # Update Global variables.
    GVAR.RULE = args.rule
    GVAR.TXS_DIFF = args.dtxs
    GVAR.SPLICED_INTRON_PIR_CUTOFF = args.pir
    f_genome = args.f_genome

    load_seqinfo(ngs_obj_list)

    linkage = Linkage()
    linkage.add_chr_linkage(args.chrom, [(args.start, args.end)])

    cluster_indx = 0
    with OutputHandle(out_dir) as f_out:
        for chrom, start, end in linkage.iterlinkage():
            try:
                signal.signal(signal.SIGALRM, time_out_handler)
                signal.alarm(args.time_out)
                _logger.debug("Start identifying elements in {0}:{1}-{2}".format(chrom, start, end))
                gene_cluster_list = identify_element(
                    chrom, start, end, ngs_obj_list, tgs_obj_list, ext_tss_list, ext_tes_list, ann, f_genome, args.paraclu_path)
                _logger.debug("Finish identifying elements in {0}:{1}-{2}".format(chrom, start, end))

                for gene_cluster in gene_cluster_list:  # list of gene cluster without any common exon
                    if not gene_cluster.has_element():
                        continue
                    print("gene")
                    print(gene_cluster)
                    print("intron_list")
                    print(gene_cluster.intron_list)
                    print("internal_exon_list")
                    print(gene_cluster.internal_exon_list)
                    print("tss_site_list")
                    print(gene_cluster.tss_site_list)
                    print("tes_site_list")
                    print(gene_cluster.tes_site_list)
                    print("tss_exon_list")
                    print(gene_cluster.tss_exon_list)
                    print("tes_exon_list")
                    print(gene_cluster.tes_exon_list)

                    cluster_indx += 1
                    cluster_name = "c_{0}".format(cluster_indx)
                    gene_cluster.write_element2bed6(*f_out.element_handles(), cluster_name)

                    _logger.debug("Start identifying transcript for {0}".format(gene_cluster))
                    trans = identify_transcript(gene_cluster, ann)
                    print("segment")
                    print(trans.segment_list)
                    print("segary")
                    print([iso.segary for iso in trans.iso_list])
                    print("isoF")
                    print(trans.isoF)
                    print("isoR")
                    print(trans.isoR)
                    print("isoA")
                    print(trans.isoA)

                    trans.write2bed12(cluster_name, *f_out.isoform_handles())
                    _logger.debug("Finish identifying transcript for {0}".format(gene_cluster))
                signal.alarm(0)
            except TimeOutError:
                print("TimeOut: {0}\t{1}\t{2}\n".format(chrom, start, end))
                with open(os.path.join(args.out_dir, "igia_debug_timeout.log"), "a") as f:
                    f.write("{0}\t{1}\t{2}\n".format(chrom, start, end))

    _logger.info("End")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
