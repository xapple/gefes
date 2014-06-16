"""A module to handle metapathways analyses.

https://wiki.bils.se/wiki/Metapathways

Defines an object "Metapathways" with one method "run()"
Takes care of settings things up right before calling the
original Metapathways.py script.

Requires that two environment variables be set:
METAPATHDIR and PATHWAYTOOLSDIR
"""

# Built-in modules #
import os, shutil

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.fasta.single import FASTA

# Third party modules #
from shell_command import shell_output

# Constants #
home = os.environ['HOME'] + '/'
metapath_dir = os.environ['METAPATHDIR'] + '/'
pathwaytools_dir = os.environ['PATHWAYTOOLSDIR'] + '/'

# Config #
default_config = """
PYTHON_EXECUTABLE       %s.pyenv/shims/python
PERL_EXECUTABLE         /usr/bin/perl
METAPATHWAYS_PATH       /proj/b2011138/nobackup/LUCAS/metapathways
PATHOLOGIC_EXECUTABLE   /proj/b2011138/nobackup/LUCAS/pathwaytools/aic-export/pathway-tools/ptools/17.0/pathway-tools
REFDBS                  blastDB/
FORMATDB_EXECUTABLE     executables/linux/bit64/makeblastdb
BLASTP_EXECUTABLE       executables/linux/bit64/blastp
BLASTN_EXECUTABLE       executables/linux/bit64/blastn
EXECUTABLES_DIR         executables/linux/bit64/
LASTDB_EXECUTABLE       executables/linux/bit64/lastdb
LAST_EXECUTABLE         executables/linux/bit64/lastal
ORF_PREDICTION          executables/linux/bit64/prodigal
SCAN_tRNA               executables/linux/bit64/trnascan-1.4
GBK_TO_FNA_FAA_GFF      libs/python_scripts/MetaPathways_parse_genbank.py
GFF_TO_FNA_FAA_GFF      libs/python_scripts/MetaPathways_input_gff.py
PREPROCESS_FASTA        libs/python_scripts/MetaPathways_filter_input.py
GFF_TO_FASTA            libs/python_scripts/MetaPathways_create_amino_sequences.py
COMPUTE_REFSCORE        libs/python_scripts/MetaPathways_refscore.py
PARSE_BLAST             libs/python_scripts/MetaPathways_parse_blast.py
ANNOTATE                libs/python_scripts/MetaPathways_annotate_fast.py
GENBANK_FILE            libs/python_scripts/MetaPathways_create_genbank_ptinput_sequin.py
CREATE_REPORT_FILES     libs/python_scripts/MetaPathways_create_reports_fast.py
STATS_rRNA              libs/python_scripts/MetaPathways_rRNA_stats_calculator.py
MLTREEMAP_IMAGEMAKER    mltreemap/mltreemap_imagemaker/mltreemap_imagemaker.pl
MLTREEMAP_CALCULATION   mltreemap/mltreemap_calculation/mltreemap.pl
""" % home

default_param = """##V.1   do not remove this line
INPUT:format fasta
quality_control:min_length  0
quality_control:delete_replicates  no
orf_prediction:algorithm    prodigal
orf_prediction:min_length   60
annotation:algorithm blast
annotation:dbs   COG_2013-02-05
#COG_2013-02-05,refseq_protein,metacyc-v5-2011-10-21
annotation:min_bsr  0.4
annotation:max_evalue 0.000001
annotation:min_score 20
annotation:min_length 60
annotation:max_hits 5
rRNA:refdbs LSURef_111_tax_silva
rRNA:max_evalue 0.000001
rRNA:min_identity 20
rRNA:min_bitscore 50
ptools_settings:taxonomic_pruning no
grid_engine:batch_size 200
grid_engine:max_concurrent_batches 400
grid_engine:walltime 10:00:00
grid_engine:RAM 8gb
grid_engine:user username
grid_engine:server grid.address.com
metapaths_steps:PREPROCESS_FASTA    yes
metapaths_steps:ORF_PREDICTION  yes
metapaths_steps:GFF_TO_AMINO    yes
metapaths_steps:FILTERED_FASTA  yes
metapaths_steps:COMPUTE_REFSCORE    yes
metapaths_steps:BLAST_REFDB  yes
metapaths_steps:PARSE_BLAST yes
metapaths_steps:SCAN_rRNA   yes
metapaths_steps:STATS_rRNA  yes
metapaths_steps:SCAN_tRNA   yes
metapaths_steps:ANNOTATE    yes
metapaths_steps:PATHOLOGIC_INPUT    yes
metapaths_steps:GENBANK_FILE    yes
metapaths_steps:CREATE_SEQUIN_FILE  yes
metapaths_steps:CREATE_REPORT_FILES yes
metapaths_steps:MLTREEMAP_CALCULATION   skip
metapaths_steps:MLTREEMAP_IMAGEMAKER    skip
metapaths_steps:PATHOLOGIC  yes
"""

###############################################################################
class Metapathways(object):
    """A metapathways analysis."""

    all_paths = """
    /output/
    /param.txt
    /config.txt
    /run.out
    /
    """

    def __repr__(self): return '<%s object of %s>' % (self.__class__.__name__, self.parent)

    def __init__(self, sample):
        # Save parent #
        self.parent, self.sample = sample, sample
        # Auto paths #
        self.base_dir = self.sample.p.metapathways_dir
        self.p = AutoPaths(self.base_dir, self.all_paths)
        # Other #
        self.fasta_path = self.sample.p.Contigs
        self.script_path = metapath_dir + "MetaPathways.py"
        self.pgdb_path = pathwaytools_dir + 'ptools-local/pgdbs/user/' + 'temp.cyc'
        # Command #
        self.params = ['-i', self.fasta_path, '-o', self.p.output_dir,
                       '-c', self.p.config, '-p', self.p.param, '-v', '-r', 'overwrite']
        self.bash_str = self.script_path + ' ' + ' '.join(self.params) + " &>> " + self.p.out
        # Objects #
        self.preprocessed_dir = self.p.output_dir + '/preprocessed/'
        self.preprocessed_fasta = FASTA(self.preprocessed_dir + 'Contigs_preproc.fasta')

    def run(self):
        # Clean #
        if os.path.exists(self.p.output_dir): shutil.rmtree(self.p.output_dir)
        # Make preferences files #
        with open(self.p.config, 'w') as handle: handle.write(default_config)
        with open(self.p.param, 'w') as handle: handle.write(default_param)
        # Make symlink for the input #
        # if os.path.lexists(self.fasta_path): os.remove(self.fasta_path)
        # Remove the PGDB if it exists #
        if os.path.exists(self.pgdb_path): shutil.rmtree(self.pgdb_path)
        # Make a stupid link #
        if os.path.lexists("blastDB"): os.remove("blastDB")
        os.symlink(metapath_dir + "blastDB/", "blastDB")
        # The command #
        with open(self.p.out, 'w') as handle: handle.write(self.bash_str + '\n\n')
        shell_output(self.bash_str)
        # Clean up #
        if os.path.exists("blastDB"): os.remove("blastDB")
        #if os.path.exists(self.fasta_path): os.remove(self.fasta_path)
        # Check for errors #
        with open(self.p.out, 'r') as handle: log = handle.read()
        if 'Error!' in log:
            raise Exception('Metapathway tools reported an error in file "%s".') % self.p.out
        # Check that the preprocessing didn't affect the data #
        assert self.sample.cleaned_fasta.count == self.preprocessed_fasta.count
        assert len(self.sample.cleaned_fasta.first) == len(self.preprocessed_fasta.first) + 1
