import time
import multiprocessing
import deeptools.utilities
from deeptools import bamHandler, countReadsPerBin
from deeptoolsintervals import GTF
import pyBigWig
import numpy as np


class CountReadsPerBinWithIntron(countReadsPerBin.CountReadsPerBin):
    def count_reads_in_region_with_intron(self, chrom, start, end, bed_regions_list=None):
        """
        Rewrite deeptools.CountReadsPerBin.count_reads_in_region

        Args:
            chrom (str): Chrom
            start (int): Start position
            end (int): End position
            bed_regions_list (list): List of bed region

        Returns:
            tuple: subnum_reads_per_bin, file_name
        """

        if start > end:
            raise NameError("start %d bigger that end %d" % (start, end))

        if self.stepSize is None:
            raise ValueError("stepSize is not set!")
        # array to keep the read counts for the regions
        subnum_reads_per_bin = []

        start_time = time.time()

        bam_handlers = []
        for fname in self.bamFilesList:
            try:
                bam_handlers.append(bamHandler.openBam(fname))
            except:
                bam_handlers.append(pyBigWig.open(fname))

        blackList = None
        if self.blackListFileName is not None:
            blackList = GTF(self.blackListFileName)

        # A list of lists of tuples
        transcriptsToConsider = []
        if bed_regions_list is not None:
            transcriptsToConsider = [x[1] for x in bed_regions_list]
        else:
            if self.stepSize == self.binLength:
                transcriptsToConsider.append([(start, end, self.binLength)])
            else:
                for i in range(start, end, self.stepSize):
                    if i + self.binLength > end:
                        break
                    if blackList is not None and blackList.findOverlaps(chrom, i, i + self.binLength):
                        continue
                    transcriptsToConsider.append([(i, i + self.binLength)])

        if self.save_data:
            _file = open(deeptools.utilities.getTempFileName(suffix='.bed'), 'w+t')
            _file_name = _file.name
        else:
            _file_name = ''

        for bam in bam_handlers:
            for trans in transcriptsToConsider:
                tcov = self.get_coverage_of_region_with_intron(bam, chrom, trans)
                if bed_regions_list is not None:
                    subnum_reads_per_bin.append(np.sum(tcov))
                else:
                    subnum_reads_per_bin.extend(tcov)

        subnum_reads_per_bin = np.concatenate([subnum_reads_per_bin]).reshape(-1, len(self.bamFilesList), order='F')

        if self.save_data:
            idx = 0
            for i, trans in enumerate(transcriptsToConsider):
                if len(trans[0]) != 3:
                    starts = ",".join([str(x[0]) for x in trans])
                    ends = ",".join([str(x[1]) for x in trans])
                    _file.write("\t".join([chrom, starts, ends]) + "\t")
                    _file.write("\t".join(["{}".format(x) for x in subnum_reads_per_bin[i, :]]) + "\n")
                else:
                    for exon in trans:
                        for startPos in range(exon[0], exon[1], exon[2]):
                            if idx >= subnum_reads_per_bin.shape[0]:
                                # At the end of chromosomes (or due to blacklisted regions), there are bins smaller than the bin size
                                # Counts there are added to the bin before them, but range() will still try to include them.
                                break
                            _file.write("{0}\t{1}\t{2}\t".format(chrom, startPos, startPos + exon[2]))
                            _file.write("\t".join(["{}".format(x) for x in subnum_reads_per_bin[idx, :]]) + "\n")
                            idx += 1
            _file.close()

        if self.verbose:
            endTime = time.time()
            rows = subnum_reads_per_bin.shape[0]
            print("%s countReadsInRegions_worker: processing %d "
                  "(%.1f per sec) @ %s:%s-%s" %
                  (multiprocessing.current_process().name,
                   rows, rows / (endTime - start_time), chrom, start, end))

        return subnum_reads_per_bin, _file_name

    def get_coverage_of_region_with_intron(self, bamHandle, chrom, regions,
                                           fragmentFromRead_func=None):
        """
        Rewrite deeptools.CountReadsPerBin.get_coverage_of_region

        Args:
            bamHandle (AlignmentFile): Bam object
            chrom (str): Chrom
            regions (list): List of block
            fragmentFromRead_func (function): Function to get fragment from read

        Returns:
            float: coverages
        """
        if not fragmentFromRead_func:
            fragmentFromRead_func = self.get_fragment_from_read
        nbins = len(regions)
        if len(regions[0]) == 3:
            nbins = 0
            for reg in regions:
                nbins += (reg[1] - reg[0]) // reg[2]
        coverages = np.zeros(nbins, dtype='float64')

        if self.defaultFragmentLength == 'read length':
            extension = 0
        else:
            extension = self.maxPairedFragmentLength

        blackList = None
        if self.blackListFileName is not None:
            blackList = GTF(self.blackListFileName)

        vector_start = 0
        for idx, reg in enumerate(regions):
            if len(reg) == 3:
                tileSize = int(reg[2])
                nRegBins = (reg[1] - reg[0]) // tileSize
            else:
                nRegBins = 1
                tileSize = int(reg[1] - reg[0])

            # Blacklisted regions have a coverage of 0
            if blackList and blackList.findOverlaps(chrom, reg[0], reg[1]):
                continue
            regStart = int(max(0, reg[0] - extension))
            regEnd = reg[1] + int(extension)

            # If alignments are extended and there's a blacklist, ensure that no
            # reads originating in a blacklist are fetched
            if blackList and reg[0] > 0 and extension > 0:
                o = blackList.findOverlaps(chrom, regStart, reg[0])
                if o is not None and len(o) > 0:
                    regStart = o[-1][1]
                o = blackList.findOverlaps(chrom, reg[1], regEnd)
                if o is not None and len(o) > 0:
                    regEnd = o[0][0]

            start_time = time.time()
            # caching seems faster.
            c = 0
            if chrom in bamHandle.references:
                reads = [r for r in bamHandle.fetch(chrom, regStart, regEnd)
                         if r.flag & 4 == 0]
            else:
                raise NameError("chromosome {} not found in bam file".format(chrom))

            prev_start_pos = None  # to store the start positions
            # of previous processed read pair
            for read in reads:
                if self.minMappingQuality and read.mapq < self.minMappingQuality:
                    continue

                # filter reads based on SAM flag
                if self.samFlag_include and read.flag & self.samFlag_include != self.samFlag_include:
                    continue
                if self.samFlag_exclude and read.flag & self.samFlag_exclude != 0:
                    continue

                # Fragment lengths
                tLen = deeptools.utilities.getTLen(read)
                if self.minFragmentLength > 0 and tLen < self.minFragmentLength:
                    continue
                if self.maxFragmentLength > 0 and tLen > self.maxFragmentLength:
                    continue

                # get rid of duplicate reads that have same position on each of the
                # pairs
                if self.ignoreDuplicates and prev_start_pos \
                        and prev_start_pos == (read.reference_start, read.pnext, read.is_reverse):
                    continue

                # since reads can be split (e.g. RNA-seq reads) each part of the
                # read that maps is called a position block.
                try:
                    position_blocks = fragmentFromRead_func(read)
                except TypeError:
                    # the get_fragment_from_read functions returns None in some cases.
                    # Those cases are to be skipped, hence the continue line.
                    continue

                # Rewrite !!!
                sIdx = vector_start + max((read.reference_start - reg[0]) // tileSize, 0)
                eIdx = vector_start + min(np.ceil(float(read.reference_end - reg[0]) / tileSize).astype('int'),
                                          nRegBins)
                coverages[sIdx:eIdx] += 1

                prev_start_pos = (read.reference_start, read.pnext, read.is_reverse)
                c += 1

            if self.verbose:
                endTime = time.time()
                print("%s,  processing %s (%.1f per sec) reads @ %s:%s-%s" % (
                    multiprocessing.current_process().name, c, c / (endTime - start_time), chrom, reg[0], reg[1]))

            vector_start += nRegBins

        # change zeros to NAN
        if self.zerosToNans:
            coverages[coverages == 0] = np.nan

        return coverages

