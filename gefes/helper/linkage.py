# Futures #
from __future__ import division

# Built-in modules #
from collections import defaultdict

# Internal modules #

# Third party modules #
import pysam

# Orientation constants #
INWARD = 0
OUTWARD = 1
INLINE = 2
INOUTWARD = 3  # unknown orientation

###############################################################################
def read_is_inward(read, regionlength):
    """Determine if read is pointing inward"""
    return (read.pos >= regionlength and not read.is_reverse) or \
           (read.pos <= regionlength and read.is_reverse)

###############################################################################
def mate_is_inward(read, regionlength):
    """Determine if mate is pointing inward"""
    return (read.mpos >= regionlength and not read.mate_is_reverse) or \
           (read.mpos <= regionlength and read.mate_is_reverse)

###############################################################################
def pair_is_inward(read, regionlength):
    """Determine if pair is pointing inward"""
    return read_is_inward(read, regionlength) and mate_is_inward(read, regionlength)

###############################################################################
def pair_is_outward(read, regionlength):
    """Determine if pair is pointing outward"""
    return not read_is_inward(read, regionlength) and not mate_is_inward(read, regionlength)

###############################################################################
def is_in_overlapping_region(readpos, contiglength, readlength, regionlength):
    """Determine if the read is in a location where the searchable region on both tips is ."""
    return readpos - readlength <= regionlength and readpos >= contiglength - regionlength

###############################################################################
def get_orientation(read, regionlength, readlength, contiglength, contigmlength):
    """Get the orientation of a pair."""
    # Check if one of the reads is in an overlapping region or not within a
    # region. In that case link can only follow two orientations
    # inward_or_outward or inline.
    if is_in_overlapping_region(read.pos, contiglength, readlength, regionlength) \
      or is_in_overlapping_region(read.mpos, contigmlength, readlength, regionlength) \
      or not is_within_region(read.pos, contiglength, readlength, regionlength) \
      or not is_within_region(read.mpos, contigmlength, readlength, regionlength):
        if read.is_reverse == read.mate_is_reverse:
            return INLINE
        else:
            return INOUTWARD
    else:
        return get_orientation_tips(read, regionlength)

###############################################################################
def get_orientation_tips(read, regionlength):
    """Determine orientation of pair if the reads are located on the tips of
    the contig."""
    if pair_is_inward(read, regionlength):
        return INWARD
    elif pair_is_outward(read, regionlength):
        return OUTWARD
    else:
        return INLINE

###############################################################################
def is_within_region(readpos, contiglength, readlength, regionlength):
    """Checks if a read is within the given region."""
    return readpos - readlength <= regionlength \
      or readpos >= contiglength - regionlength

###############################################################################
def is_link(read):
    """Checks if the pair is linking to contigs."""
    return read.is_paired \
      and read.tid != read.mrnm

###############################################################################
def parse_linkage_info_bam(bamfile, readlength, min_contig_length, regionlength, fullsearch):
    """Parses a bamfile to retrieve the linkage information.
    Returns a three-dimensional dictionary of linkage information between
    contigs and a dictionary of the read count for each contigs. The linkage
    dictionary stores links only once e.g. if [contig1][contig2] exists,
    [contig2][contig1] does not. The innermost dictionary of linkage has 4 keys
    representing counts for INWARD, OUTWARD, INLINE and INOUTWARD orientation.
    The parameter readlength is the untrimmed readlength. Only contigs longer
    than given min_contig_length are considered. The search space for links is
    only at the tips of the contig in case fullsearch is set to False. The
    parameter regionlength represents the number of bases that should be
    searched over at the tips of the contigs. In case fullsearch is set to True
    the whole contig is searched, but the orientation of a pair cannot be
    determined from one pair that is far away from the tips so all the linkage
    further than regionlength away from the tips will be considered INOUTWARD.
    With very short contigs and large regionlength, the regionlengths might
    overlap, if one of the reads in the pair is in the area INOUTWARD it will
    also be considered INOUTWARD."""
    linkdict = defaultdict(lambda: defaultdict(lambda: [0, 0, 0, 0]))
    read_count_dict = defaultdict(lambda: 0)
    with pysam.Samfile(bamfile, 'rb') as bamh:
        reflens = bamh.lengths
        for read in bamh:
            # check if read is aligned
            try: ref = bamh.getrname(read.tid)
            except ValueError: continue
            # check if contig meets length threshold
            # check if the read falls within specified region
            if reflens[read.tid] >= min_contig_length \
              and (fullsearch
              or is_within_region(read.pos, reflens[read.tid], readlength, regionlength)):
                read_count_dict[ref] += 1
                # look at first read to prevent counting links twice
                # linked contig should be within threshold as well
                # linked mate should be within specified region
                if read.is_read1 \
                  and is_link(read) \
                  and reflens[read.mrnm] >= min_contig_length \
                  and (fullsearch
                  or is_within_region(read.mpos, reflens[read.mrnm], readlength, regionlength)):
                    # check if mate is aligned
                    try: refm = bamh.getrname(read.mrnm)
                    except ValueError: continue
                    # Add one link
                    ori = get_orientation(read, regionlength, readlength, reflens[read.tid], reflens[read.mrnm])
                    if ref not in linkdict[refm]:
                        linkdict[ref][refm][ori] += 1
                    else:
                        linkdict[refm][ref][ori] += 1
    # Remove default generators
    linkdict.default_factory = None
    for v in linkdict.itervalues():
        v.default_factory = None
    read_count_dict.default_factory = None
    return linkdict, read_count_dict