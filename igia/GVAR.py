'''
Global variables
'''
SLEEP_TIME = 0.01
N_WORKER = 20
MAX_RAM_PER_ival = 4 * 1024 * 1024 * 1024
NO_PREDICT = False
CONTINUE = False
MAX_QUEUE_LEN = 10

RULE = ""  # library-type: SE or PE RNA strand
ALPHA = 0.01  # for confidence interval
LINKAGE_SEARCH_WINDOW_SIZE = 100
MIN_INTRON_LEN = 10
MAX_INTRON_LEN = 100000
MINEXON_LEN = 3
MAXEXON_LEN = 10000
NGS_FILE_INFO_list = None
NGS_LIMIT_READ_COVER_PER_BASE = None
GENOME_SIZE = None
CHROM_SIZE = None


TXS_DIFF = 500  # TSS/TES distance cutoff for merging in TGS
MAX_GAP_INEXON = 50  # maximum gaps inside exon
TXS_BOUNDARY_MEDIAN_FRACTION = 0.2  # Signal cutoff at TSS/TES boundary
NGS_STRAND_FRACTION = 0.3  # maximum ratio of antisense/(sense+antisense)
TGS_STRAND_FRACTION = 0.5  # maximum ratio of antisense/(sense+antisense)
MIN_DIFF_OUT_TXS = 100  # minimum distance to external TSS/TES to extend as new TXS

# TGS
TGS_LIMIT_READ_COVER = 3
TGS_OVERLAP_RATIO = 0.4  # minimum overlap ratio as TGS cluster for strand decision

# Junction
MAX_MISMATCH = 2  # max mismatch in a intron junction
JUNCTION_READ_ENDEXON_LEN = 8  # minimum length of first/last exonic segment
MIN_JUNCTION_READ_NGS = 2  # Minimum read# in confirming a junction/edge in NGS
MIN_JUNCTION_BOUNDARY_READ_NGS = 5  # Minimum read# in confirming a junction boundary/node in NGS
JUNCTION_BOUNDARY_RATIO = 0.06  # Minimum fraction of being as out-going junction
MIN_JUNCTION_BOUNDARY_READ_TGS = 2  # Minimum read numbers in confirming a junction boundary in TGS
SPLICED_INTRON_PIR_CUTOFF = 0.3  # PIR cutoff for intron
NEIGHBOREXON_LEN = 50  # Neighbor exon search ival flanking intron
MIN_JUNCTION_READ_TGS = 4
MIN_NGS_COV_FOR_ALL_JUNCTION_BY_TGS = 30  # If a gene cluster's max cov less than this value, do not filter.
MAPPING_SLIDING_ERROR = 20

TES_RADIUS_TGS = 50  # TES variation range in TGS
TSS_RADIUS_TGS = 100  # TSS variation range in TGS

# Transcript
SEGMENT_TGS_OVERLAP_FRACTION = 0.5  # minimum coverage of TGS read on segment to be exon
MINOR_AS_PHI_CUTOFF = 0.3
MAJOR_AS_PHT_CUTOFF = 0.7
MINOR_AS_CNT_RATIO = 0.5
MAX_AS_NUM_FOR_ISOP = 50
MAX_ISO_NUM_FOR_ISOP = 200

"""If two overlapped gene linkages are on two strands and have same intron fraction > cutoff, then treat them as one"""
SHARE_INTRON_CUTOFF = 0.1


