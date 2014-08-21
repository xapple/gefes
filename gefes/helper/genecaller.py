# Built-in modules #

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.fasta.single import FASTA
from gefes.common.slurm import nr_threads


# Third party modules #
import sh, shutil

###############################################################################
class GeneCaller(object):
    """Genecalling pipelines:
    * glimmer3 inspired by the g3-iterated script from glimmer
    * prodigal
    """

    all_paths = """
    /glimmer/
    /prodigal/
    /prokka/
    """

    def __repr__(self): return '<%s object of bin "%s" with %i predictions >' % (self.__class__.__name__,self.parent.name, len(self))
    def __iter__(self): return iter(self.callers)
    def __len__(self): return len(self.callers)
    def __getitem__(self, index): return self.callers[index]

    def __init__(self,parent):
        self.parent = parent
        # Auto paths #
        self.base_dir = self.parent.p.genes
        self.p = AutoPaths(self.base_dir, self.all_paths)
        self.callers = [ Glimmer(self), Prodigal(self) , Prokka(self)]

    def run(self):
        for cal in self: cal.run()

###################################################################
class Glimmer(object):

    glimmer_params={'updist' : 24}

    all_paths= """
    /genes.fasta
    """

    def __init__(self,parent):
        self.parent = parent
        self.bini = self.parent.parent
        # Auto paths #
        
        self.base_dir = self.parent.p.glimmer
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        seqs=FASTA(self.bini.p.contigs)
        motif_count = sh.Command('get-motif-counts.awk')
        sh.touch('temp.genes.fasta')
        # generate single fasta entry files to compute long-orfs for gene prediction
        for s in seqs:
            temp=FASTA('./temp.fasta')
            s.description=s.name
            temp.write(s)
            try:
                sh.long_orfs('-l','-g', 30 ,'-n', '-t', 1.15, './temp.fasta', './temp.longorfs')
                sh.extract('-t', './temp.fasta', './temp.longorfs', _out="temp2.fasta")
                sh.sed( '-n', 'wtemp.fasta', 'temp.genes.fasta', 'temp2.fasta')
                sh.mv('temp.fasta', 'temp.genes.fasta')
            except sh.ErrorReturnCode:
                sh.rm('temp.fasta')
        seqs.close()
        sh.rm('temp2.fasta')
        sh.rm('temp.longorfs')
        # generate ICM from long-orfs
        with open('temp.genes.fasta') as inpo:
            sh.build_icm('-r','temp.icm',_in=inpo)
        # use ICM to predict genes
        sh.glimmer3('-l','-o50', '-g110', '-t30', self.bini.p.contigs, 'temp.icm', 'temp')
        # extract and correct upstream corrdinated to extract from multi-fasta file
        with open('temp.predict') as f_in:
                content = f_in.readlines()
                with open('temp.gene.coord','w') as g_out:
                        with open('temp.upstream.coord','w') as up_out:
                            for line in content:
                                    if '>' in line:
                                        iD=line[1:-1].split(" ")[0]
                                    else:
                                        line = [w for w in line.split(" ") if w!='']
                                        line_gene = "\t".join([line[0],iD,"\t".join(line[1:])])
                                        g_out.write(line_gene)
                                        start = int(line[1])
                                        end = int(line[2])
                                        if start > end:
                                            start = str(int(start)+1)
                                            end = str(int(start)+self.glimmer_params['updist'])
                                        else:
                                            start = str(int(start)-1)
                                            end = str(int(start)-self.glimmer_params['updist'])
                                        line[1] = start
                                        line[2] = end
                                        line= "\t".join([line[0],iD,"\t".join(line[1:3])])+"\n"
                                        up_out.write(line)
        # extract upstream regions for binding motif stats
        sh.multi_extract(self.bini.p.contigs, 'temp.upstream.coord', _out='temp.upstream')
        motif_count(sh.elph('temp.upstream', 'LEN=6'), _out='temp.motif')
        # TODO : need still to modify start codon use stat thing they use
        # second glimmer run
        sh.glimmer3('-l', '-o50', '-g110', '-t30','-b','temp.motif', self.bini.p.contigs, 'temp.icm', 'final')
        # extract it from the data again
        with open('final.predict') as f_in:
            content = f_in.readlines()
            with open('gene.coord','w') as g_out:
                    for line in content:
                        if '>' in line:
                                iD=line[1:-1].split(" ")[0]
                        else:
                            line=[w for w in line.split(" ") if w!='']
                            line_gene= "\t".join([line[0],iD,"\t".join(line[1:])])
                            g_out.write(line_gene)
        sh.multi_extract('-d', '-t',self.bini.p.contigs, 'gene.coord', _out= 'temp.genes')
        # rewrite it properly so that each gene seq as a unique name of type orf_name;contig
        in_file=FASTA('./temp.genes')
        with FASTA(self.p.genes) as out_file:
            for seq in in_file:
                    seq.description=seq.description.replace("  ",";",1)
                    seq.name=seq.description.split("  ")[0]
                    seq.id=seq.name
                    out_file.add_seq(seq)
        # clean-up
        sh.rm('gene.coord')
        sh.rm('temp.motif')
        sh.rm('temp.genes')
        sh.rm('temp.genes.fasta')
        sh.rm('temp.icm')
        sh.rm('temp.predict')
        sh.rm('temp.gene.coord')
        sh.rm('temp.upstream.coord')
        sh.rm('temp.upstream')
        sh.rm('final.predict')
        sh.rm("final.detail")
        sh.rm("temp.detail")

###################################################################
class Prodigal(object):

    all_paths= """
    /genes.fasta
    """

    def __init__(self,parent):
        self.parent = parent
        self.bini = self.parent.parent
        # Auto paths #
        self.base_dir = self.parent.p.prodigal
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        sh.prodigal('-q', '-c', '-i', self.bini.p.contigs, '-d', self.p.genes, _out = "/dev/null")


class Prokka(object):

    all_paths= """
    /genes.fasta
    /proteins.fasta
    /output/
    """

    def __init__(self,parent):
        self.parent = parent
        self.bini = self.parent.parent
        # Auto paths #
        self.base_dir = self.parent.p.prokka
        self.p = AutoPaths(self.base_dir, self.all_paths)

    def run(self):
        sh.prokka( "--outdir", self.p.output, "--prefix",  self.bini.name, "--compliant", "--locustag", self.bini.name, "--force", "--cpus", nr_threads, self.bini.p.contigs)    
        shutil.copy(self.p.output + "/" + self.bini.name + ".fna", self.p.genes)
        shutil.copy(self.p.output + "/" + self.bini.name + ".faa", self.p.proteins)
