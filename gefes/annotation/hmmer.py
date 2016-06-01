# Futures #
from __future__ import division

# Built-in modules #
import warnings, multiprocessing

# Internal modules #
from seqsearch.databases.pfam    import pfam
from seqsearch.databases.tigrfam import tigrfam

# First party modules #
from fasta import FASTA
from plumbing.autopaths import FilePath
from plumbing.cache     import property_cached

# Third party modules #
import sh, pandas

# Warnings #
warnings.filterwarnings("ignore", "Bio.SearchIO")
warnings.filterwarnings("ignore", "BiopythonWarning")
from Bio import SearchIO

###############################################################################
class HmmQuery(object):
    """An `hmmsearch` job."""

    short_name = 'hmmsearch'
    long_name  = 'HMMER 3.1b2 (February 2015)'
    executable = 'hmmsearch'
    url        = 'http://hmmer.org/'
    license    = 'GPLv3'
    dependencies = []

    def __nonzero__(self): return bool(self.out_path)
    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.query)

    def __init__(self, query_path,                    # The input sequences
                 db_path      = pfam.hmm_db,          # The database to search
                 seq_type     = 'prot' or 'nucl',     # The seq type of the query_path file
                 e_value      = 1e-10,                # The search threshold
                 params       = None,                 # Add extra params for the command line
                 out_path     = None,                 # Where the results will be dropped
                 executable   = None,                 # If you want a specific binary give the path
                 cpus         = None):                # The number of threads to use
        # Save attributes #
        self.query      = FASTA(query_path)
        self.db         = FilePath(db_path)
        self.params     = params if params else {}
        self.e_value    = e_value
        self.seq_type   = seq_type
        self.executable = FilePath(executable)
        # Cores to use #
        if cpus is None: self.cpus = min(multiprocessing.cpu_count(), 32)
        else:            self.cpus = cpus
        # Auto detect database short name #
        if db_path == 'pfam':    self.db = pfam.hmm_db
        if db_path == 'tigrfam': self.db = tigrfam.hmm_db
        # Output #
        if out_path is None:         self.out_path = FilePath(self.query.prefix_path + '.hmmout')
        elif out_path.endswith('/'): self.out_path = FilePath(out_path + self.query.prefix + '.hmmout')
        else:                        self.out_path = FilePath(out_path)

    @property
    def command(self):
        # Executable #
        if self.executable: cmd = [self.executable.path]
        else:               cmd = ["hmmsearch"]
        # Essentials #
        cmd += ('-o',        '/dev/null',   # direct output to file `f`, not stdout
                '--tblout',  self.out_path, # parsable table of per-sequence hits
                '--seed',    1,             # set RNG seed to `n`
                '-E',        self.e_value,  # report sequences <= this E-value threshold in output
                '--notextw',                # unlimited ASCII text output line width
                '--acc',                    # prefer accessions over names in output
                self.db,
                self.query)
        # Options #
        for k,v in self.params.items(): cmd += [k, v]
        # Return #
        return map(str, cmd)

    def run(self, cpus=None):
        """Simply run the HMM search locally."""
        # Number of threads #
        if cpus is None: cpus = self.cpus
        # Checks #
        assert self.query.exists
        assert self.db.exists
        # Check if query is not empty #
        if self.query.count_bytes == 0:
            message = "Hmm search on a file with no sequences. File at '%s'"
            warnings.warn(message % self.query, RuntimeWarning)
            return False
        # Do it #
        sh.Command(self.command[0])(['--cpu', str(cpus)] + self.command[1:])

    def test_biopython(self):
        """Unfortunately biopython is unable to parse these file."""
        formats = ['hmmer3-tab', 'blast-xml', 'phmmer3-domtab', 'exonerate-text', 'hmmer3-text', 'blast-text', 'hmmer2-text', 'exonerate-vulgar', 'exonerate-cigar', 'blast-tab', 'hmmsearch3-domtab', 'blat-psl', 'hmmscan3-domtab' ]
        for fmt in formats:
            try:
                print SearchIO.read(self.out_path, fmt)
                print "yes: ", fmt
            except ValueError:
                print "no: ", fmt

    @property_cached
    def hits(self):
        """Here the target_name is the contig and the query_accession is the pfam"""
        # Check #
        if not self.out_path:
            raise Exception("You can't access results from HMMER before running the algorithm.")
        # Column names #
        fields = ['target_name', 'target_accession', 'query_name', 'query_accession',
                  'e_value', 'score', 'bias', 'domain_e_value', 'domain_score',
                  'domain_bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc']
        # Parse #
        load = lambda path: pandas.io.parsers.read_csv(path,
                                                       comment          = "#",
                                                       delim_whitespace = True,
                                                       engine           = 'python',
                                                       header           = None,
                                                       index_col        = False,
                                                       names            = fields)
        # Result #
        return load(self.out_path)

    @property_cached
    def distinct_pfams(self):
        """The distinct PFAMS that were found in this search, in a set"""
        pfams = self.hits[self.hits['e_value'] <= self.e_value]
        pfams = set(row['query_accession'].split('.')[0] for i, row in pfams.iterrows())
        return pfams