# Futures #
from __future__ import division

# Built-in modules #
import os, shutil

# Internal modules #
from hiseq.common import AutoPaths
from hiseq.graphs import assembly_plots

# Third party modules #
import sh
from shell_command import shell_output

###############################################################################
class Assembly(object):
    """An assembly analysis."""

    all_paths = """
    /graphs/
    /velvet/contigs.fasta
    /amos/contigs.afg
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, pool):
        # Save parent #
        self.parent, self.pool = pool, pool
        # Auto paths #
        self.base_dir = self.parent.base_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def assemble(self):
        # Clean up #
        shutil.rmtree(self.p.velvet_dir)
        # Params #
        ins = self.pool.report_stats['fwd']['FragmentSize']
        metavelvetg = sh.Command('meta-velvetg')
        contig_paths = []
        # First iteration #
        for kmer in (91,):
            # Message #
            print "-> Kmer set to %i" % kmer
            # Out directory #
            out_dir = self.p.velvet_dir + str(kmer) + '/'
            if not os.path.exists(out_dir): os.mkdir(out_dir)
            # Input #
            fwd_in, rev_in = self.pool.fwd.path, self.pool.rev.path
            # Special case #
            if 'Double' in self.pool.info.get('remarks', ''): fwd_in, rev_in = self.pool.smaller.fwd_path, self.pool.smaller.rev_path
            # Velevet Hashing #
            sh.velveth(out_dir, kmer, '-shortPaired', '-fastq.gz', fwd_in, rev_in,
                       _out=out_dir + 'velveth.out')
            # Velvet Graphing #
            sh.velvetg(out_dir, '-exp_cov', 'auto', '-ins_length', ins, '-unused_reads', 'yes', '-min_contig_lgth', 200,
                      _out=out_dir + 'velvetg.out')
            # Velevet metagenomics #
            metavelvetg(out_dir, '-ins_length', ins,
                       _out=out_dir + 'meta-velvetg.out')
            # Output #
            out_path = out_dir + "meta-velvetg.contigs.fa"
            sh.fasta_prefix_headers(str(kmer) + '_', out_path, '--inplace')
            contig_paths.append(out_path)
            # Next input #
            in_file = out_dir + "UnusedReads.fa"
        # Loop #
        for kmer in (81, 75):
            # Message #
            print "-> Kmer set to %i" % kmer
            # Out directory #
            out_dir = self.p.velvet_dir + str(kmer) + '/'
            if not os.path.exists(out_dir): os.mkdir(out_dir)
            # Velevet Hashing #
            sh.velveth(out_dir, kmer, '-shortPaired', '-fasta', in_file,
                       _out=out_dir + 'velveth.out')
            # Velvet Graphing #
            sh.velvetg(out_dir, '-exp_cov', 'auto', '-ins_length', ins, '-unused_reads', 'yes', '-min_contig_lgth', 200,
                      _out=out_dir + 'velvetg.out')
            # Velevet metagenomics #
            metavelvetg(out_dir, '-ins_length', ins,
                       _out=out_dir + 'meta-velvetg.out')
            # Output #
            out_path = out_dir + "meta-velvetg.contigs.fa"
            sh.fasta_prefix_headers(str(kmer) + '_', out_path, '--inplace')
            contig_paths.append(out_path)
            # Next input #
            in_file = out_dir + "UnusedReads.fa"
        # Merge #
        print "-> Merging outputs"
        shell_output("cat %s > %s" % (' '.join(contig_paths), self.p.contigs_fasta))

    def scaffold(self):
        # Clean #
        shutil.rmtree(self.p.amos_dir)
        # AMOS conversion #
        print "-> Calling amos"
        to_amos = sh.Command('/bubo/home/h3/lucass/share/amos/bin/toAmos')
        to_amos('-s', self.p.contigs_fasta, '-o', self.p.contigs_afg)
        # Scaffolding #
        print "-> Calling minimus"
        minimus = sh.Command('/bubo/home/h3/lucass/share/amos/bin/minimus2')
        minimus(self.p.contigs_afg[0:-4], '-D', 'MINID=98', 'OVERLAP 200', _out=self.p.amos_dir + 'minimus.out')

    def make_plots(self):
        for cls_name in assembly_plots.__all__:
            cls = getattr(assembly_plots, cls_name)
            cls(self).plot()
