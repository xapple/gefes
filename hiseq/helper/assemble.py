# Futures #
from __future__ import division

# Built-in modules #
import os, shutil

# Internal modules #

# Third party modules #
import sh
from shell_command import shell_output

###############################################################################
def assemble(pool):
    # Params #
    ins = pool.report_stats['fwd']['FragmentSize']
    metavelvetg = sh.Command('meta-velvetg')
    contig_paths = []
    # First iteration #
    for kmer in (91,):
        # Message #
        print "-> Kmer set to %i" % kmer
        # Out directory #
        out_dir = pool.p.velvet_dir + str(kmer) + '/'
        if not os.path.exists(out_dir): os.mkdir(out_dir)
        # Input #
        fwd_in, rev_in = pool.fwd.path, pool.rev.path
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
        out_dir = pool.p.velvet_dir + str(kmer) + '/'
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
    shell_output("cat %s > %s" % (' '.join(contig_paths), pool.p.contigs_fasta))

###############################################################################
def scaffold(pool):
    # Clean #
    shutil.rmtree(pool.p.amos_dir)
    # AMOS conversion #
    print "-> Calling amos"
    to_amos = sh.Command('/bubo/home/h3/lucass/share/amos/bin/toAmos')
    to_amos('-s', pool.p.contigs_fasta, '-o', pool.p.contigs_afg)
    # Scaffolding #
    print "-> Calling minimus"
    minimus = sh.Command('/bubo/home/h3/lucass/share/amos/bin/minimus2')
    minimus(pool.p.contigs_afg[0:-4], '-D', 'MINID=98', 'OVERLAP 200', _out=pool.p.amos_dir + 'minimus.out')