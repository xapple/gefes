
#built-in modules#
import os

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.common.slurm import nr_threads
from gefes.common import flatten
from gefes.fasta.single import FASTA

# Third party mods #
import sh
from numpy import array
import numpy
from matplotlib import cm
from matplotlib import pyplot
from pandas import DataFrame
from sklearn import mixture

class GenomeBuilder(object):
    """reassembling a genome from a bins and the reads"""

    all_paths = """
    /scaffolds.fasta
    /filtered.scaffold.fasta
    /contigs
    /graph_contigs.pdf
    /reads.fastq
    /reads.1.fastq
    /reads.2.fastq
    /bin_vs_reassembled.png
    """

    def __init__(self,parent):
        self.parent = parent
        self.bini = parent
        self.fwds=[str(p.fwd.path) for p in self.bini.parent.parent.parent]
        self.revs=[str(p.rev.path) for p in self.bini.parent.parent.parent]
        # Auto paths #
        self.base_dir = self.parent.p.rebuilt
        self.p = AutoPaths(self.base_dir, self.all_paths)

    
    def pull_reads(self):
        sh.bowtie2_build(self.parent.p.contigs,self.p.contigs)
        sh.bowtie2("-p", nr_threads, "-x",self.p.contigs,"-1", ",".join(self.fwds), "-2", ",".join(self.revs), "--al-conc", self.p.fastq, "-S", "/dev/null")

    def assemble_genome(self):
        spades_script = sh.Command("spades.py")
        spades_script("-o", self.base_dir,  "-1", self.p.reads_1, "-2", self.p.reads_2,  "-t", nr_threads, "--careful")

    def filter_assembly(self,cutoff=0.3):
        original = FASTA(self.p.scaffolds)
        data = DataFrame.from_dict({s.id : array(s.id.split("_"))[[3,5]] for s in original}, orient='index')
        data = data.astype(numpy.float)
        previous_min = array([c.length for c in self.parent.contigs]).min()
        data = data[data[0] > previous_min]
        classes = self.clust(numpy.log10(data))
        self.contig_plot(data,classes)
        lens = data[0].groupby(classes).apply(sum)
        data[classes == lens.idxmax()]
        keepers = data[classes == lens.idxmax()].index
        with FASTA(self.p.filtered) as filtered:
            filtered.add_seq([s for s in original if s.id in keepers])
            
    def clust(self,data):
        gmm=mixture.GMM(n_components=2)
        gmm.fit(data[0])
        clust_x = gmm.predict(data[0])
        gmm.fit(data[1])
        clust_y = gmm.predict(data[1])
        return clust_x + 2*clust_y

    def contig_plot(self,data,classes):
        fig = pyplot.figure()
        axes = fig.add_subplot(111)
        axes.set_xscale('log')
        axes.set_yscale('log')
        pyplot.scatter(data[0],data[1],c=classes,cmap=cm.rainbow)
        fig.savefig(self.p.graph)
            
    def graphs(self):
        return 0
        #        nucmer -l 1000  scaffolds.fasta ../../thorsellia_type.fasta
#        mummerplot --layout --png --large -R scaffolds.fasta -Q ../../thorsellia_type.fasta  out.delta
#        gnuplot out.gp > out.png

# echo  $'#!/bin/bash\nstretcher -asequence people/olle/type_strain_correction/scaffold_builder/scaffolds/scaff.fasta -bsequence people/olle/type_strain_correction/scaffold_builder/scaffolds/thorsellia_type.fasta -outfile test.needle' | sbatch -A b2013086 -p core -n 1 -J alignment  -t 1-01:00:00 -e stretcher.err
