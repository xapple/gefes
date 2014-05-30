
#built-in modules#

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.common.slurm import nr_threads
from pandas import DataFrame

# Third party mods #
import sh
import os

class Phylotyper(object):
    """phylotyping bins with amphora of phylosift"""

    all_paths = """
    /phylosift/
    /kraken/
    /kraken/output_kraken.txt
    /kraken/report_kraken.txt
    /kraken/krona_kraken.html
    """

    kraken_db = "/glob/moritz/data/kraken_data/"
    sift = sh.Command(os.environ['PHYLOSIFT_PATH'])
    kraken_report = sh.Command("kraken-report")

    def __init__(self,parent):
        self.parent = parent
        self.bini = parent
        # Auto paths #
        self.base_dir = self.parent.p.phylotyping
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        self.phylosift()
        self.kraken()

    def phylosift(self):
        self.sift("all", "-f", "--output=" + self.p.phylosift, self.parent.p.contigs)
                  
    def kraken(self):
        sh.kraken("--db", self.kraken_db, "--threads", nr_threads, "--fasta-input", "--output", self.p.output, self.parent.p.contigs)
        if os.stat(self.p.output)[6]!=0 : self.kraken_report("--db", self.kraken_db, self.p.output, _out = self.p.report)
        sh.ktImportTaxonomy("-t", 5, "-m", 3, "-o", self.p.krona, self.p.report)
