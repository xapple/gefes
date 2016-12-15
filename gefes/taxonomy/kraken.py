# Built-in modules #
import os
from collections import OrderedDict

# Internal modules #

# First party modules #
from plumbing.autopaths import AutoPaths
from plumbing.cache import property_cached
from plumbing.slurm import num_processors

# Third party modules #
import sh, pandas

# Constants #
home = os.environ['HOME'] + '/'
standard_db = home + 'databases/kraken/standard'

###############################################################################
class Kraken(object):
    """Use Kraken at to predict taxonomy on the raw reads.
    Expects version v0.10.5-beta.
    """

    short_name   = 'kraken'
    long_name    = 'Kraken v0.10.5-beta'
    executable   = 'kraken'
    url          = 'https://github.com/DerrickWood/kraken'
    dependencies = ['jellyfish']

    all_paths = """
    /raw_output.txt
    /log.txt
    /summary.tsv
    """

    def __nonzero__(self): return bool(self.p.summary)
    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.source)

    def __init__(self, source, base_dir=None):
        # Basic #
        self.source   = source
        self.base_dir = base_dir
        # Default case #
        if base_dir is None: self.base_dir = self.source.fwd.prefix_path + '.kraken/'
        # Auto paths #
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self, keep_raw=False, cpus=None):
        # Variable threads #
        if cpus is None: cpus = num_processors
        # Check version #
        assert "0.10.5-beta" in sh.kraken('--version')
        # Run the main classification #
        sh.kraken('--preload',
                  '--threads', str(cpus),
                  '--db',      standard_db,
                  '--output',  self.p.raw_output,
                  '--paired',  self.source.fwd, self.source.rev,
                  _err=str(self.p.log))
        # Run the report #
        report = sh.Command('kraken-report')
        report('--db', standard_db, self.p.raw_output, _out=str(self.p.summary))
        # Remove intermediary output #
        if not keep_raw: self.p.raw_output.remove()

    @property_cached
    def results(self):
        results = KrakenResults(self)
        if not results: raise Exception("You can't access results from Kraken before running the tool.")
        return results

###############################################################################
class KrakenResults(object):

    all_paths = """
    /output/lorem
    """

    def __nonzero__(self): return bool(self.kraken.p.summary)
    def __init__(self, kraken):
        self.kraken = kraken

    @property_cached
    def composition(self):
        """1) Percentage of reads covered by the clade rooted at this taxon
           2) Number of reads covered by the clade rooted at this taxon
           3) Number of reads assigned directly to this taxon
           4) A rank code, indicating (U)nclassified, (D)omain, (K)ingdom,
              (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies.
              All other ranks are simply '-'.
           5) NCBI taxonomy ID
           6) indented scientific name
        """
        columns = ['percentage', 'total_reads', 'count_reads', 'rank', 'ncbi_id', 'name']
        result  =  pandas.io.parsers.read_csv(self.kraken.p.summary,
                                              sep='\t', encoding='utf-8',
                                              header=None, names=columns)
        result['name'] = result['name'].apply(lambda x: x.lstrip())
        return result

    @property_cached
    def at_domain_level(self):
        x = self.composition
        return {'Unclassified': x[x['rank']=='U']['percentage'],
                'Bacteria':     x[x['name']=='Bacteria']['percentage'],
                'Archaea':      x[x['name']=='Archaea']['percentage'],
                'Viruses':      x[x['name']=='Viruses']['percentage']}

    @property_cached
    def at_phylum_level(self):
        result = self.composition
        result = result[(result['rank']=='P') & (result['percentage'] != 0)]
        result = result.sort(columns='percentage', ascending=False)
        result = OrderedDict(zip(result['name'], result['percentage']))
        return result

    @property_cached
    def at_species_level(self):
        result = self.composition
        result = result[(result['rank']=='S') & (result['percentage'] != 0)]
        result = result.sort(columns='percentage', ascending=False)
        result = OrderedDict(zip(result['name'], result['percentage']))
        return result
