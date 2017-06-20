# Built-in modules #

# Internal modules #
import gefes
from gefes.groups.aggregates import Aggregate

# First party modules #
from plumbing.cache     import property_cached
from plumbing.autopaths import DirectoryPath

# Third party modules #

###############################################################################
class Bundle(Aggregate):
    """Regroup result files and reports from one or several projects for delivery."""

    all_paths = """
    /projects/
    /samples/
    /metadata/samples.xlsx
    """

    def __init__(self, name, samples, out_dir=None):
        # Directory #
        if out_dir is None: out_dir = gefes.bundles_dir
        # Super #
        super(self.__class__,self).__init__(name, samples, out_dir)
        # Figure out the projects within #
        proj_names = sorted(list(set([s.project_short_name for s in samples])))
        self.projects = [gefes.projects[p] for p in proj_names]

    def run(self):
        self.base_dir.remove()
        self.base_dir.create()
        for p in self.projects:
            # Directories #
            proj_dir = DirectoryPath(self.p.projects_dir + p.name)
            # Reports directory #
            reports_dir = DirectoryPath(proj_dir + 'reports')
            reports_dir.create(safe=True)
            # Reports for samples #
            smpl_rprts_dir = DirectoryPath(reports_dir + 'samples')
            smpl_rprts_dir.create(safe=True)
            for s in p: s.report.output_path.link_to(smpl_rprts_dir + s.short_name + '.pdf')
            # Reports for assemblies #
            ably_rprts_dir = DirectoryPath(reports_dir + 'assemblies')
            ably_rprts_dir.create(safe=True)
            p.merged.results.report.output_path.link_to(ably_rprts_dir + p.merged.short_name + '.pdf')
            # Reports for bins #
            bin_rprts_dir = DirectoryPath(reports_dir + 'bins')
            bin_rprts_dir.create(safe=True)
            for b in p.merged.results.binner.results.bins:
                if not b.good: continue
                b.report.output_path.link_to(bin_rprts_dir + b.name + '.pdf')
            # Extras #
            extras_dir = DirectoryPath(proj_dir + 'extras')
            extras_dir.create(safe=True)
            # PFAM table #
            p.merged.results.hit_profile.p.norm_samples_x_pfams.link_to(extras_dir)
            # Coverage table #
            p.merged.results.binner.coverage_matrix_tsv.link_to(extras_dir)
            # Reports for project (not super useful) #
            p.report.output_path.link_to(extras_dir + 'project_report.pdf')
            # Data files #
            data_dir = DirectoryPath(proj_dir + 'data')
            data_dir.create(safe=True)
            # Contigs #
            contigs_dir = DirectoryPath(data_dir + 'contigs')
            contigs_dir.create(safe=True)
            # Contig fasta #
            p.merged.results.contigs_fasta.link_to(contigs_dir + 'all_contigs.fasta')
            # TSV which contigs goes in which bins #
            lines = (k + '\t' + v + '\n' for k,v in p.merged.results.binner.results.contig_id_to_bin_id.items())
            (contigs_dir + 'contig_id_to_bin_id.tsv').writelines(lines)
            lines = (k + '\t' + ','.join(v) + '\n' for k,v in p.merged.results.binner.results.bin_id_to_contig_ids.items())
            (contigs_dir + 'bin_id_to_contig_ids.tsv').writelines(lines)
            # Bins #
            bins_data_dir = DirectoryPath(data_dir + 'bins')
            bins_data_dir.create(safe=True)
            for b in p.merged.results.binner.results.bins:
                if not b.good: continue
                current_dir = DirectoryPath(bins_data_dir + b.name)
                current_dir.create(safe=True)
                b.fasta.link_to(current_dir + 'contigs.fasta')
            # Bowtie output #
            mapping_dir = DirectoryPath(data_dir + 'mapping')
            mapping_dir.create(safe=True)
            for s in p:
                s.mapper_merged.results.p.map_smds_bam.link_to(mapping_dir + s.name + '.bam')
                (s.mapper_merged.results.p.map_smds_bam + '.bai').link_to(mapping_dir + s.name + '.bam.bai')

    @property_cached
    def results(self):
        results = BundleResults(self)
        message = "You can't access results from a bundle before making the bundle."
        if not results: raise Exception(message)
        return results

###############################################################################
class BundleResults(object):

    def __nonzero__(self): return bool(self.p.x)

    def __init__(self, parent):
        self.parent   = parent
        self.base_dir = parent.base_dir
        self.p        = parent.p
