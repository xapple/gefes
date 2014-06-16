# Futures #
from __future__ import division

# Built-in modules #
import os, socket

# Internal modules #
from gefes.common import flatten
from gefes.common.autopaths import AutoPaths
from gefes.common.cache import property_cached
from gefes.graphs import assembly_plots
from gefes.helper.contig import Contig
from gefes.helper.metapathways import Metapathways
from gefes.fasta.single import FASTA
from gefes.common.slurm import nr_threads

# Third party modules #
import sh
import pandas

# Constant #
hostname = socket.gethostname()

###############################################################################
class Assembly(object):
    """Will run the co-assembly of several pools by calling the Ray assembler.
    https://github.com/sebhtml/ray"""

    short_name = 'ray'
    executable = 'Ray23'

    all_paths = """
    /graphs/
    /ray_output/frame.csv
    /ray_output/
    /ray_output/Contigs.fasta
    /ray_output/report.txt
    /metapathways/
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, aggregate):
        # Save parent #
        self.parent, self.aggregate = aggregate, aggregate
        # Auto paths #
        self.base_dir = self.parent.p.assembly_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Convenience objects #
        self.contigs_fasta = FASTA(self.p.Contigs)
        # Graphs #
        self.graphs = [getattr(assembly_plots, cls_name)(self) for cls_name in assembly_plots.__all__]
        # Metapath #
        self.metapat = Metapathways(self)

    @property_cached
    def contigs(self):
        """A list of all the contigs produced as custom objects"""
        return [Contig(self, s) for s in self.contigs_fasta]

    def assemble(self):
        # Ray needs a non-existing directory #
        out_dir = self.p.output_dir
        out_dir.remove()
        # Make the pairs of fastq #
        pairs = flatten([('-p', p.cleaner.fwd.path, p.cleaner.rev.path) for p in self.parent])
        # Call Ray on the cray #
        if os.environ.get('CSCSERVICE') == 'sisu':
            stats = sh.aprun('-n', nr_threads, self.executable, '-k', 81, '-o', out_dir, *pairs)
        # Call Ray on Kalkyl #
        elif os.environ.get('SNIC_RESOURCE') == 'kalkyl':
            stats = sh.mpiexec('-n', nr_threads, self.executable, '-k', 81, '-o', out_dir, *pairs)
        # Call Ray just locally #
        else:
            command = sh.Command(self.executable)
            stats = command('-k', 81, '-o', out_dir, *pairs)
        # Print the report #
        with open(self.p.report, 'w') as handle: handle.write(str(stats))

    def index(self):
        """Create two indexes. For both bowtie2 and samtools on the contigs fasta file."""
        sh.bowtie2_build(self.contigs_fasta, self.contigs_fasta)
        sh.samtools('faidx', self.contigs_fasta)

    def make_plots(self):
        for graph in self.graphs: graph.plot()

    def write_contiglist(self,contig_list,path,file_name):
        contig_list = [o for o in self.contigs if o.name in contig_list]
        with FASTA(os.path.join(path,file_name)) as ffile:
            for c in contig_list: ffile.add_seq(c.record)

    @property_cached
    def frame(self):
        columns = ['length'] + ['gc_content'] + [s.id_name for s in self.aggregate] + ['freq_' + t for t in Contig.tetra_cats]
        rows = [c.name for c in self.contigs]
        data = [[c.length] + [c.gc_content] + [s.mapper.coverage[c.name]["cov_mean"] for s in self.aggregate] + c.get_all_tetra_nuc_freqs() for c in self.contigs]
        return pandas.DataFrame(data, columns=columns, index=rows)

    def filtered_frame(self,max_freq=None,min_len=None):
        temp_frame = self.frame
        if(min_len is not None):
            temp_frame = temp_frame[temp_frame.length > min_len]
        if(max_freq is not None):
            good_ones = temp_frame[[c for c in temp_frame if "freq" in c]].apply(lambda x: sum(x>max_freq)==0,1)
            temp_frame = temp_frame[good_ones]
        return temp_frame

    def export_frame(self):
        self.frame.to_csv(self.p.frame)