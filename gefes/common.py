# Built-in modules #
import re

# Internal modules #

# First party modules #

# Third party modules #

###############################################################################
def join_paired_filepaths(paths):
    """Every _R1_ goes with every _R2_"""
    assert len(paths) % 2 == 0
    pairs = []
    for path in paths:
        if re.search("_R1_", path):
            pairs.append((path, re.sub("_R1_", "_R2_", path)))
    return pairs