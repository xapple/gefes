# Built-in modules #
import os

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.common.slurm import nr_threads

# Third party mods #
import sh
from numpy import array
from pandas import DataFrame
import pandas

# Constants #
kraken_db = "/glob/moritz/data/kraken_large/"
sift = sh.Command(os.environ['PHYLOSIFT_PATH'])
kraken_report = sh.Command(os.environ['KRAKEN_REPORT'])

###############################################################################
class Phylotyper(object):
    """Phylotyping bins with amphora of phylosift."""

    all_paths = """
    /phylosift/
    /kraken/
    /kraken/output_kraken.txt
    /kraken/output_wo_unassigned_kraken.txt
    /kraken/report_kraken.txt
    /kraken/krona_kraken.html
    /kraken/otu_kraken.csv
    """

    def __init__(self, parent):
        self.parent = parent
        self.bini = parent
        # Auto paths #
        self.base_dir = self.parent.p.phylotyping
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        self.phylosift()
#        self.kraken()

    def phylosift(self):
        sift("all", "-f", "--threads", nr_threads/2, "--output=" + self.p.phylosift, self.parent.p.contigs)

    def kraken(self):
        if hasattr(self.parent,"pair"):
            #for a Pool object
            sh.kraken("--preload", "--db", kraken_db, "--threads", nr_threads , "--gzip-compressed", "--fastq-input", "--paired",  parent.fwd_path, parent.rev_path, _out =  str(self.p.output))
        else:
            #for a bin object with pulled reads
            sh.kraken("--preload", "--db", kraken_db, "--threads", nr_threads,  "--fastq-input", "--paired",  self.bini.p.reads_1, self.bini.p.reads_2, "--output",  self.p.output)
            # legacy option for contigs directly --- sh.kraken("--preload", "--db", kraken_db, "--threads", nr_threads, "--fasta-input", "--output", self.p.output, self.parent.p.contigs)
        
        if os.stat(self.p.output)[6]!=0 :
            sh.tail("-n","+2", self.p.output, _out = str(self.p.output_wo_unassigned))
            kraken_report("--db", kraken_db, self.p.output_wo_unassigned, _out = str(self.p.report))
        sh.ktImportTaxonomy("-t", 5, "-m", 3, "-o", self.p.krona, self.p.report)
        self.cumtree2phylacount()
        
    def cumtree2phylacount(self, depth=7):
        tree = pandas.read_table(self.p.report, sep="\t", names=array(["prop","cum_count","count","type","taxa_id","taxa"]),index_col=False)
        tree["depth"] = tree["taxa"].apply(lambda s: (len(s) - len(s.lstrip()))/2)
        tree = tree[tree["depth"] < depth]
        counts = {}
        for i,row in enumerate(tree.iterrows()):
            if i +1 <  len(tree) and tree.iloc[i + 1]['depth'] > row[1]['depth'] :
                count = row[1]['count']
            else :
                count = row[1]['cum_count']
            taxa = row[1]['taxa'].lstrip()
            counts[taxa + " (level " + str(row[1]["depth"]) + ")"] = count
        counts =  { k: v for k,v in counts.iteritems() if v > 0 }
        temp = DataFrame.from_dict(counts,'index')
        temp = temp.sort(column = 0, ascending = False)
        temp.to_csv(self.p.otu, header = [self.bini.name])
        return counts
