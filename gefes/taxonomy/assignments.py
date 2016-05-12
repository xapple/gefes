# Built-in modules #
import re
from collections import OrderedDict

# Internal modules #

# First party modules #

# Third party modules #

###############################################################################
class Level(object):
    """Describes a taxonomic level."""

    def __str__(self):
        return self.name.lower() + ' level'

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
    def __repr__(self): return '<Taxon: "%s" at %s>' % (self.term, self.level)
    def __eq__(self, other): return other == self.term
    def __nonzero__(self): return bool(self.term)
    def __init__(self, term, level):
        self.term = term
        self.level = level

###############################################################################
class Assignment(object):
    """Useful for parsing strings such as the following:
    > d__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Actinomycetaceae.g__Mobiluncus.s__curtisii.t__subsp__holmesii_ATCC_35242"""

    def __str__(self):
        return self.lowest_taxon.term + ' (' + self.lowest_taxon.level.name + ' level)'

    def __init__(self, string, separator='.'):
        # Save attributes #
        self.string    = string.strip('\n')
        self.separator = separator
        self.re_sep    = re.escape(separator)
        # Parse #
        self.taxa = [self.parse_level(l) for l in levels]
        # Lowest taxonomy information #
        self.lowest_taxon = [t for t in self.taxa if t and t.term != '?'][-1]

    def parse_level(self, level):
        found = re.findall(level.identifier +"__([^%s]+)" % self.re_sep, self.string)
        if not found: term = None
        else:         term = found[0]
        return Taxon(term, level)

###############################################################################
class Assignments(object):
    """Useful for handling a collection of assignments."""
    def __init__(self, path):
        pass
