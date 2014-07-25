# Built-in modules #
import os

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.common.slurm import nr_threads

# Third party mods #
import sh

# Constants #
kraken_db = "/glob/moritz/data/kraken_data/"
sift = sh.Command(os.environ['PHYLOSIFT_PATH'])
kraken_report = sh.Command(os.environ['KRAKEN_REPORT'])

###############################################################################
class Phylotyper(object):
    """Phylotyping bins with amphora of phylosift."""

    all_paths = """
    /phylosift/
    /kraken/
    /kraken/output_kraken.txt
    /kraken/report_kraken.txt
    /kraken/krona_kraken.html
    """

    def __init__(self, parent):
        self.parent = parent
        self.bini = parent
        # Auto paths #
        self.base_dir = self.parent.p.phylotyping
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        self.phylosift()
        self.kraken()

    def phylosift(self):
        sift("all", "-f", "--threads", nr_threads, "--output=" + self.p.phylosift, self.parent.p.contigs)

    def kraken(self):
        sh.kraken("--db", kraken_db, "--threads", nr_threads, "--fasta-input", "--output", self.p.output, self.parent.p.contigs)
        if os.stat(self.p.output)[6]!=0 : kraken_report("--db", kraken_db, self.p.output, _out = self.p.report)
        sh.ktImportTaxonomy("-t", 5, "-m", 3, "-o", self.p.krona, self.p.report)
