
#built-in modules#

# Internal modules #
from gefes.common.autopaths import AutoPaths
from gefes.helper.contig import Contig
from gefes.helper.genecaller import GeneCaller
from gefes.helper.blast_dbs import blast_dbs
from gefes.helper.blast_dbs import blast_header
from parallelblast import BLASTquery
from gefes.common.slurm import nr_threads
from pandas import DataFrame
from pandas import Series
from numpy import array
from numpy import argsort

# Third party mods #
import sh
import os

class BinAnnotater(object):
    """A bin is a object containing a multifasta of contigs that are ready to be annotated"""

    blast_params={
        'single_copy_cogs' :  {'algorithm': 'tblastx' , '-m': '8', '-a': nr_threads},
        }

    blast_cutoff={
        'single_copy_cogs' : 0.01
        }
    
    def __init__(self,parent):
        self.parent = parent
        self.bini = parent

    def single_copy_cog_blast(self):
        scc_dbs = [dbs for dbs in blast_dbs.iteritems() if 'SCG' in dbs[1]['sets'] and  'cog' in dbs[1]['sets']]
        db_max_size = max([os.stat(lise[1]['path']).st_size for lise in scc_dbs])
        
        output = []
        single = []
        weird = [] # count of perfectly similiar high quality hits
        idx = -1

        for caller in self.bini.caller:
            ipath = caller.p.genes
            path = caller.p._base_dir + "/SingleCopyCogs/"
            opath = path + "raw/"
            ppath = path +"processed/"
            print "annotating for " + str(caller)

        


            if not os.path.exists(opath): os.makedirs(opath)
            if not os.path.exists(ppath): os.makedirs(ppath)

            search_space_size=db_max_size + os.stat(ipath).st_size
            search_space_size = search_space_size*search_space_size
            scc_blast_params = self.blast_params['single_copy_cogs']
            scc_blast_params['-Y'] = search_space_size

            idx = idx + 1

            output.append({})
            single.append({})
            weird.append({})
            for db in scc_dbs:
                    b_query = BLASTquery(ipath,db[1]['path'],scc_blast_params, opath + db[0])
                    b_query.run()
                    if os.stat(b_query.out_path).st_size > 0:
                            raw = DataFrame.from_csv(b_query.out_path,sep='\t',header=-1)
                            orfs = list(set(raw.index))
                            processed=DataFrame(columns=raw.columns)
                    
                            for orf in orfs:
                                    max_es = raw.loc[orf][10] == raw.loc[orf][10].min()
                                    if isinstance(max_es,Series):
                                            best_hits = raw.loc[orf][max_es]
                                            if len(max_es) > 1:
                                                max_len = best_hits[3] == best_hits[3].max()
                                                best_hits = best_hits[max_len]
                                            processed = processed.append(best_hits.iloc[0])
                                    else:
                                            processed = processed.append(raw.loc[orf])

                                            
                            output[idx][db[0]] = processed[processed[10] < self.blast_cutoff['single_copy_cogs']]
                            output[idx][db[0]].to_csv(ppath + db[0],sep="\t",header=blast_header)

                            if len(output[idx][db[0]]) == 0 : single[idx][db[0]] = 'abscent'
                            elif len(output[idx][db[0]]) == 1 : single[idx][db[0]] = 'single'
                            else :
                                data = output[idx][db[0]][[8,9]]

                                # remove self-chimeras
                                tups = [row[1:4] for row in output[idx][db[0]][[8,9,10]].itertuples()]
                                not_dups = [t not in tups[0:i] for i,t in enumerate(tups)]
                                data = data.iloc[not_dups]
                                weird[idx][db[0]] = sum([not b for b in not_dups])
                                 
                                coords = array([sorted(list(z)[1:]) for z in data.itertuples()])
                                rnk = argsort([z[0] for z in coords])
                                f_coords = coords[rnk].flatten()
                                bads = 0
                                for i in range(1,len(f_coords)-1):
                                    if f_coords[i] < f_coords[i-1]:
                                        bads = bads +1
                                if bads > 0 : single[idx][db[0]] = 'multiple'
                                else : single[idx][db[0]] = 'single'
                    else:
                            single[idx][db[0]] = 'abscent'
                            
                    DataFrame(single[idx],index=["presence"]).transpose().to_csv(path + "single_copy_cogs.tsv",sep="\t")
                    DataFrame(weird[idx],index=["presence"]).transpose().to_csv(path + "auto_chimeras.tsv",sep="\t")

        self.single = single
        self.weird = weird
                    

 
