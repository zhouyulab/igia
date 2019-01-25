from mpi4py import MPI
from igia import __version__
from igia import GVAR
from igia.utils import SeqFile, load_seqinfo, load_txs, bed2bam
from igia.linkage import Linkage, find_linkage_worker
from igia.element import identify_element
from igia.transcript import identify_transcript
from igia.skeleton import parse_args, setup_logging, OutputHandle, check_paraclu
import sys
import os
import logging
import numpy as np
import subprocess
import time, random

_logger = logging.getLogger(__name__)


class LBTask(object):
    """Load balancing tasks distribution system"""
    def __init__(self, comm, finish_msg="done", root=0, master_tag=0, worker_tag=1, sleep=0.1):
        assert isinstance(master_tag, int)
        assert isinstance(worker_tag, int)
        self.comm = comm
        self.master_tag = master_tag
        self.worker_tag = worker_tag
        self.finish_msg = finish_msg
        self.sleep = sleep
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()
        self.root = root
        if self.rank == root:
            self.type = "master"
        else:
            self.type = "worker"

    def do(self, *args, **kwargs):
        raise NotImplementedError()

    def debug(self, text):
        #with open("node_{0}.log".format(self.rank), "a") as f:
        #    f.write(text)
        pass

class LBMaster(LBTask):
    """Master in charge of distributing tasks"""
    def __init__(self, comm, queue_len, finish_msg="done", root=0, master_tag=0, worker_tag=1, sleep=0.1):
        super().__init__(comm, finish_msg, root, master_tag, worker_tag, sleep)
        self.queue_len = queue_len
        self.workers = list(range(self.size))
        self.workers.remove(self.rank)
        self.task_num = np.zeros_like(self.workers)
        self.working_num = 0

    def select_worker(self):
        """Select the worker with least tasks to distribute task"""
        min_task_indx = self.task_num.argmin()
        if self.task_num[min_task_indx] == 0:
            return self.workers[min_task_indx]
        else:
            return None

    def receive_feedback(self, worker):
        """Update the tasks number of each worker"""
        assert self.task_num[self.workers.index(worker)] >= 1
        assert self.working_num >= 1
        self.task_num[self.workers.index(worker)] -= 1
        self.working_num -= 1
        self.debug('LBMaster: received result from worker {0}\n'.format(worker))
        self.debug('LBMaster: there are {0} tasks ongoing\n'.format(self.working_num))

    def record_task(self, worker):
        """Record the tasks number of each worker"""
        self.task_num[self.workers.index(worker)] += 1
        self.working_num += 1

    def send_task(self, data, worker):
        self.comm.send(data, dest=worker, tag=self.master_tag)
        _logger.debug('LBMaster: Sending data "{0}" to worker {1}'.format(data, worker))
        self.debug('LBMaster: Sending data "{0}" to worker {1}\n'.format(data, worker))
        self.record_task(worker)

    def recv_data_and_wait(self):
        """Receive data and """
        try:
            wait = False
            while not self.comm.iprobe(tag=self.worker_tag): # confirm data passage complete
                wait = True # waited once, maybe no queue of results
                time.sleep(self.sleep)

            worker, res = self.comm.recv(tag=self.worker_tag)
            self.receive_feedback(worker)
            self.debug("LBMaster: received data {0} from worker {1}\n".format(res, worker))
            return res, wait
        except Exception as ex:
            self.debug("LBMaster {0}: {1}\n".format(type(ex), ex))
            _logger.error("LBMaster {0}: {1}".format(type(ex), ex))
            raise ex

    def do(self, data_list):
        """Distributing tasks to workers"""
        result_list = list()

        task_indx = 0
        while data_list and (task_indx < len(data_list)):
            # Distribute tasks
            selected_worker = self.select_worker()
            while (selected_worker is not None) and (task_indx < len(data_list)):
                self.debug("LBMaster: send data {0}\n".format(data_list[task_indx]))
                self.send_task(data_list[task_indx], selected_worker)
                task_indx += 1
                selected_worker = self.select_worker()

            # Receive feedback
            while self.working_num:
                res, wait = self.recv_data_and_wait()
                if res != self.finish_msg:
                    result_list.append(res)
                
                if wait:
                    break

        # Collect the final results from workers
        while self.working_num:
            res, _ = self.recv_data_and_wait()
            if res != self.finish_msg:
                result_list.append(res)

        # All tasks have been finished, sending finish massage to all workers
        self.debug("LBMaster: send finish message\n")
        for worker in self.workers:
            self.send_task(self.finish_msg, worker)

        # Collect the finish message from workers
        while self.working_num:
            res, _ = self.recv_data_and_wait()
            assert res == self.finish_msg
        return result_list


class LBWorker(LBTask):
    """Worker in dealing with tasks"""
    def recv_data(self):
        try:
            while not self.comm.iprobe(tag=self.master_tag): # confirm data passage complete
                time.sleep(self.sleep)
            data = self.comm.recv(tag=self.master_tag)
            self.debug("LBWorker: received data {0}\n".format(data))
            return data
        except Exception as ex:
            self.debug("LBWorker {0}: {1}\n".format(type(ex), ex))
            raise ex


class LBIdentifyIsoWorker(LBWorker):
    """Worker for identifying isoform"""
    def do(self, out_dir, ngs_obj_list, tgs_obj_list, ext_tss_list, ext_tes_list, ann, f_genome, paraclu_path=None):
        """Identify isoforms, until receive finish message"""
        gene_cluster_indx = 0
        with OutputHandle(os.path.join(out_dir, "node_{0}".format(self.rank))) as f_out:
            data = self.recv_data()
            while data != self.finish_msg:
                chrom, start, end = data
                gene_cluster_list = identify_element(
                    chrom, start, end, ngs_obj_list, tgs_obj_list,
                    ext_tss_list, ext_tes_list, ann, f_genome, paraclu_path)

                for gene_cluster in gene_cluster_list:
                    if not gene_cluster.has_element():
                        continue
                    gene_cluster_indx += 1
                    gene_cluster_name = "n_{0}_c_{1}".format(self.rank + 1, gene_cluster_indx)
                    gene_cluster.write_element2bed6(*f_out.element_handles(), gene_cluster_name)
                    trans = identify_transcript(gene_cluster, ann)
                    trans.write2bed12(gene_cluster_name, *f_out.isoform_handles())
                self.debug("LBIdentifyIsoWorker: finish identification\n")
                self.comm.send((self.rank, self.finish_msg), dest=self.root, tag=self.worker_tag)
                data = self.recv_data()
        self.comm.send((self.rank, self.finish_msg), dest=self.root, tag=self.worker_tag)


class LBScanLinkageWorker(LBWorker):
    """Worker for scan linkage"""
    def do(self, bam_list):
        """Scan linkages, until receive finish message"""
        res_list = list()
        data = self.recv_data()
        while data != self.finish_msg:
            if data == "err":
                continue
            chrom_size, seq_obj_indx = data
            linkage = find_linkage_worker(chrom_size, bam_list[seq_obj_indx])
            res_list.append(linkage)
            self.comm.send((self.rank, res_list[-1]), dest=self.root, tag=self.worker_tag)
            self.debug("LBScanLinkageWorker: finish scan\n")
            data = self.recv_data()
        self.comm.send((self.rank, self.finish_msg), dest=self.root, tag=self.worker_tag)


def main(args):
    """Main entry point allowing external calls
    Args:
      args ([str]): command line parameter list
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    rank_size = comm.Get_size()

    if rank == 0:
        args = parse_args(args)
        check_paraclu(args)
        out_dir = args.out_dir
        rand_seq = str(hash(time.time() * 100000 + random.random()))
        tmp_dir = os.path.join(out_dir, rand_seq)
        args.tmp_dir = tmp_dir
        setup_logging(args.loglevel)
        _logger.info("Starting IGIA in MPI mode ...")
    else:
        args = None

    args = comm.bcast(args, root=0)

    out_dir = args.out_dir
    tmp_dir = args.tmp_dir

    ngs_obj_list = [SeqFile(x, "NGS") for x in args.ngs_file]
    tgs_obj_list = [SeqFile(x, "TGS") for x in args.tgs_file]
    ext_tss_list = load_txs(args.tss)
    ext_tes_list = load_txs(args.tes)
    f_ann = args.ann
    size = args.size
    f_genome = args.f_genome
    paraclu_path = args.paraclu_path

    ann = None
    if f_ann:
        if rank == 0:
            _logger.info("Loading annotation transcripts ...")
            if not size:
                raise ValueError("Error: missing Chrom Size file")

            ann_bam = bed2bam(f_ann, size, tmp_dir)
        else:
            ann_bam = None
        ann_bam = comm.bcast(ann_bam, root=0)
        ann = SeqFile(ann_bam, "ANN")

    # Update Global variables
    GVAR.RULE = args.rule
    GVAR.TXS_DIFF = args.dtxs
    GVAR.SPLICED_INTRON_PIR_CUTOFF = args.pir

    load_seqinfo(ngs_obj_list)

    _logger.info("Start building linkage ... ")
    bam_list = ngs_obj_list + tgs_obj_list
    if ann is not None:
        bam_list += [ann]

    # Scatter linkage scan tasks
    if rank == 0:
        chrom_size_list = bam_list[0].chromsize()
        scan_linkage_infos = [
            (chrom_size, seq_obj_indx) 
            for chrom_size in chrom_size_list 
            for seq_obj_indx in range(len(bam_list))]
        master = LBMaster(comm, GVAR.MAX_QUEUE_LEN, master_tag=0, worker_tag=1, sleep=GVAR.SLEEP_TIME)
        linkages = master.do(scan_linkage_infos)
        linkage = Linkage()
        for sub_linkage in linkages:
            linkage.add_linkage(sub_linkage)
    else:
        worker = LBScanLinkageWorker(comm, master_tag=0, worker_tag=1)
        worker.do(bam_list)
        linkage = None

    _logger.info("Finish building linkage ... ")

    # with open("node_{0}.log".format(rank), "a") as f:
    #     f.write("finish building linkage\n")

    if rank == 0:
        data = "do"
    else:
        data = None
    data = comm.bcast(data, root=0)
    assert data == "do"

    # Identify transcripts
    if rank == 0:
        _logger.info("Start identifying transcripts")
        linkage_region_list = list(linkage.iterlinkage())
        master = LBMaster(comm, GVAR.MAX_QUEUE_LEN, master_tag=2, worker_tag=3, sleep=GVAR.SLEEP_TIME)
        master.do(linkage_region_list)
    else:
        worker = LBIdentifyIsoWorker(comm, master_tag=2, worker_tag=3)
        worker.do(tmp_dir, ngs_obj_list, tgs_obj_list, ext_tss_list, ext_tes_list, ann, f_genome, paraclu_path)

    # Merge results
    if rank == 0:
        _logger.info("Finish identifying transcripts")
        outfiles = ["intron.bed6", "internal_exon.bed6", "tss_exon.bed6", "tes_exon.bed6",
                    "isoF.bed12", "isoA.bed12", "isoR.bed12", "isoM.bed12", "isoC.bed12", "isoP.bed12"]
        res_dirs = [
            os.path.join(tmp_dir, "node_{0}".format(x)) for x in range(1, rank_size)]
        for name in outfiles:
            out_files = [os.path.join(res_dir, name) for res_dir in res_dirs]
            code = "cat " + " ".join(out_files) + " > " + os.path.join(out_dir, name)
            subprocess.call(code, shell=True)
        for res_dir in res_dirs:
            code = "rm -r {0}".format(res_dir)
            subprocess.call(code, shell=True)

    _logger.info("End")


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
