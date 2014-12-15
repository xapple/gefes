# Built-in modules #
import os

# Internal modules #
from humic.taxonomy import TaxaClassifier
from plumbing.tmpstuff import TmpFile
from humic.helper.trees import convert_megan_tree

# Third party modules #
import sh

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Megan(TaxaClassifier):
    """A megan5 analysis. One per sample.
    http://ab.inf.uni-tuebingen.de/software/megan5/"""

    all_paths = TaxaClassifier.all_paths + """
    /create.txt
    /create.out
    /create.err
    /export.txt
    /export.out
    /export.err
    /taxa.txt
    """

    def __init__(self, sample):
        # Save parent #
        self.parent, self.sample = sample, sample
        # Auto paths #
        self.base_dir = self.parent.p.megan_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Super #
        TaxaClassifier.__init__(self)

    def create_rma(self):
        script = """
            import blastfile='%s' fastafile='%s' meganfile='%s' maxmatches=100 minscore=50.0 toppercent=10.0 minsupport=5 mincomplexity=0.44 useseed=false usekegg=false useidentityfilter=false blastformat=BlastXML
            update
            quit"""
        script = script % (self.sample.searches.nt.p.hits, self.sample.cleaned_fasta, self.sample.searches.nt.p.rma)
        script = '\n'.join(l.lstrip(' ') for l in script.split('\n') if l)
        self.run_script(script, self.p.create_out, self.p.create_err, self.p.create_txt)

    def view(self): self.run_script("open file='%s'"% (self.sample.searches.nt.p.rma))

    def export_taxonomy(self):
        script = """
            open file='%s'
            update
            collapse rank=Species
            update
            select nodes=all
            update
            export what=DSV format=taxonpath_count separator=tab file='%s'
            quit"""
        script = script % (self.sample.searches.nt.p.rma, self.p.taxa)
        script = '\n'.join(l.lstrip(' ') for l in script.split('\n') if l)
        self.run_script(script, self.p.export_out, self.p.export_err, self.p.export_txt)

    def run_script(self, script, out='/dev/null', err='/dev/null', save_script=None):
        # Prepare script #
        if save_script: script_file = TmpFile(save_script, script)
        else: script_file = TmpFile(content=script)
        # Do it #
        megan = sh.Command(home + 'share/megan5/MEGAN')
        megan('+g', '-c', script_file, _out=out, _err=err)
        # Clean up #
        if not save_script: script_file.remove()
        # Check results #
        log = open(out).read()
        if "No X11 DISPLAY variable was set" in log: raise Exception("X11 connection missing.")
        if "Can't connect to X11 window server" in log: raise Exception("X11 connection broken.")

    def make_graph(self):
        # Make new tree #
        self.convert_megan_tree()
        # Super #
        TaxaClassifier.make_graph(self)

    def convert_tree(self):
        convert_megan_tree(self.p.taxa, self.taxa_file)

    @property
    def clade_counts(self):
        # Parse #
        def all_counts():
            with open(self.p.taxa, 'r') as handle:
                for line in handle:
                    line = line.strip('\n').split('\t')
                    count = int(line.pop())
                    path = [name for name in line.pop().strip('"').split(';') if name]
                    level = len(path)
                    name = path[-1].strip('"').replace(',','-')
                    yield count, name, level
        # Return #
        return dict(((n+" (lvl %i)"%l,c) for c,n,l in all_counts()))