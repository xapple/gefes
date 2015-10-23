# Built-in modules #
import re
from collections import OrderedDict

# Internal modules #

# First party modules #

# Third party modules #

###############################################################################
class Level(object):
    """Describes a taxonomic level."""
    def __init__(self, name, identifier, i):
        self.name = name
        self.identifier = identifier
        self.i = i

levels = [Level('Life',    'l', 1),
          Level('Domain',  'd', 2),
          Level('Kingdom', 'k', 3),
          Level('Phylum',  'p', 4),
          Level('Class',   'c', 5),
          Level('Order',   'o', 6),
          Level('Family',  'f', 7),
          Level('Genus',   'g', 8),
          Level('Species', 's', 9),
          Level('Strain',  't', 10)]

###############################################################################
class Taxon(object):
    """Describes a taxonomic term at a given level."""
    def __nonzero__(self, term, level): return bool(term)
    def __init__(self, term, level):
        self.term = term
        self.level = level

###############################################################################
class Assignment(object):
    """Useful for parsing strings such as the following:
    > d__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Actinomycetaceae.g__Mobiluncus.s__curtisii.t__subsp__holmesii_ATCC_35242"""

    def __str__(self):
        return self.lowest_taxon.term + ' (' + self.lowest_taxon.level.identifier + ')'

    def __init__(self, string, seperator='.'):
        # Save attributes #
        self.string    = string
        self.seperator = seperator
        self.re_sep    = re.escape(seperator)
        # Parse #
        self.taxa = OrderedDict([self.parse_level(l) for l in levels])
        # Lowest taxonomy information #
        self.lowest_taxon = [t for t in self.taxa if t][-1]

    def parse_level(self, level):
        term = re.findall(level.identifier + "__(.+)" + self.re_sep, self.string)[0]
        term = re.match(self.string, level.identifier)
        return Taxon(term, level)

###############################################################################
class Assignments(object):
    """Useful for handling a collection of assignments."""
    def __init__(self, path):
        pass
