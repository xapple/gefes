# -*- coding: utf-8 -*-

# Built-in modules #

# Internal modules #

# First party modules #

# Third party modules #

# Functions #
def bool_to_unicode(b):
    """Different possibilities for True: ☑️✔︎✓✅👍
       Different possibilities for False: ✕✖︎✗✘✖️❌⛔️❎👎
    """
    if not isinstance(b, bool): b = bool(b)
    if b is True:  return "✅"
    if b is False: return "❎"

###############################################################################
class ProjectStatus(object):
    """Prints information about a project. Will tell you up to which
    point the data has been processed."""

    def __init__(self, project):
        self.project = project
        self.proj    = project
        self.samples = project.samples

    steps = ['raw', 'first_qc', 'cleaned', 'second_qc', 'initial_taxa', 'mono_assembly',
             'co_assembly', 'mono_mapping', 'merged_assembly', 'mappings',
             'binning', 'merged_binning', 'check_m']

    # Default displaying of boolean outcomes #
    bts = lambda b: bool_to_unicode(b) or str(b)

    @property
    def header(self):
        message = "Project '%s' with %i samples"
        message = message % (self.proj.name, len(self.samples))
        return message

    def print_long(self):
        print self.header
        print '-' * 40
        print self.status(details=True)

    def print_short(self):
        print self.header
        print '-' * 40
        print self.status(details=False)

    def status(self, details=False):
        message = ""
        for step in self.steps:
            title, detail, outcome = getattr(self, step)
            message += '-' * 40
            message += title + ':  ' + bool_to_unicode(outcome) + '\n'
            if details: message += detail + '\n'
        return message

    #---------------------------------- Steps --------------------------------#
    @property
    def raw(self):
        title    = "The raw files for each of the samples"
        func     = lambda s: bool(s.pair)
        items    = self.samples
        outcome = all(func(s) for s in items)
        detail   = '\n'.join(str(s) + ': ' + self.bts(func(s)) for s in items)
        return title, detail, outcome

    @property
    def first_qc(self):
        title    = "The first quality control on the raw data"
        func     = lambda s: bool(s.pair.fwd.fastqc)
        items    = self.samples
        outcome  = all(func(s) for s in items)
        detail   = '\n'.join(str(s) + ': ' + self.bts(func(s)) for s in items)
        return title, detail, outcome

    @property
    def cleaned(self):
        title    = "The cleaning of the raw data"
        func     = lambda s: bool(s.quality_checker.results)
        items    = self.samples
        outcome  = all(func(s) for s in items)
        detail   = '\n'.join(str(s) + ': ' + self.bts(func(s)) for s in items)
        return title, detail, outcome

    @property
    def second_qc(self):
        title    = "The second quality control on the cleaned data"
        func     = lambda s: bool(s.clean.fwd.fastqc)
        items    = self.samples
        outcome  = all(func(s) for s in items)
        detail   = '\n'.join(str(s) + ': ' + self.bts(func(s)) for s in items)
        return title, detail, outcome

    @property
    def initial_taxa(self):
        title    = "The initial taxonomic evaluation with Kraken"
        func     = lambda s: bool(s.kraken)
        items    = self.samples
        outcome  = all(func(s) for s in items)
        detail   = '\n'.join(str(s) + ': ' + self.bts(func(s)) for s in items)
        return title, detail, outcome

    @property
    def mono_assembly(self):
        title    = "The mono assembly with just the sample as input"
        func     = lambda s: bool(s.assembly)
        items    = self.samples
        outcome  = all(func(s) for s in items)
        detail   = '\n'.join(str(s) + ': ' + self.bts(func(s)) for s in items)
        return title, detail, outcome

    @property
    def co_assembly(self):
        title    = "The different co-assemblies just with all samples as input"
        func     = lambda a: bool(a)
        items    = [(n,a) for n,a in self.proj.assemblies.items()]
        outcome  = all(func(a) for n,a in items)
        detail   = '\n'.join(n + ': ' + self.bts(func(a)) for n,a in items)
        return title, detail, outcome

    @property
    def mono_mapping(self):
        title    = "The mapping of each sample to their own mono-assembly"
        func     = lambda s: bool(s.mono_mapper.p.coverage)
        items    = self.samples
        outcome  = all(func(s) for s in items)
        detail   = '\n'.join(str(s) + ': ' + self.bts(func(s)) for s in items)
        return title, detail, outcome

    @property
    def merged_assembly(self):
        title    = "The merged assembly with multiple kmer sizes"
        func     = lambda p: bool(p.merged.results)
        items    = [self.proj]
        outcome  = all(func(x) for x in items)
        detail   = '\n'.join(str(x) + ': ' + self.bts(func(x)) for x in items)
        return title, detail, outcome

    @property
    def mappings(self):
        title    = "The mappings of each sample to different assemblies"
        func     = lambda m: bool(m.p.coverage)
        items    = [(s,a,m) for a,m in s.mappers.items() for s in self.samples]
        outcome  = all(func(x) for x in items)
        detail   = '\n'.join("Map %s to %s:"%(s,a) + ': ' + self.bts(func(m)) for s,a,m in items)
        return title, detail, outcome

    @property
    def binning(self):
        title    = "The binning of all the contigs"
        func     = lambda a: bool(a.results.binner.p.clustering)
        items    = [(n,a) for n,a in self.proj.assemblies.items()]
        outcome  = all(func(a) for n,a in items)
        detail   = '\n'.join(n + ': ' + self.bts(func(a)) for n,a in items)
        return title, detail, outcome

    @property
    def merged_binning(self):
        title    = "The binning of the contigs in the merged assembly"
        func     = lambda p: bool(p.merged.results.binner.p.clustering)
        items    = [self.proj]
        outcome  = all(func(x) for x in items)
        detail   = '\n'.join(str(x) + ': ' + self.bts(func(x)) for x in items)
        return title, detail, outcome

    @property
    def check_m(self):
        title    = "The CheckM run on every merged-assembly bin"
        func     = lambda b: bool(b.evaluation.p.stdout)
        items    = [b for b in proj.merged.results.binner.bins]
        outcome  = all(func(s) for s in items)
        detail   = '\n'.join(str(b) + ': ' + self.bts(func(b)) for b in items)
        return title, detail, outcome
