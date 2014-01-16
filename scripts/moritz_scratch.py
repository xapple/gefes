import gefes
from gefes.fasta.single import FASTA
import csv

maps=[ p.mapper.coverage for p in gefes.projects['humic']]
contig_lists=[[m for m in sample if sample[m]['cov_mean']>1] for sample in maps]

pat="/home/moritz/people/sari/assemblies/"
names=[p.label for p in gefes.projects["humic"]]
fpath= [pat+n+".fasta" for n in names]

    
def write_contig_list(contig_list,path,file):
    contig_list=[o for o in gefes.projects["acI"].assembly.contigs if o.name in contig_list]
    ffile=FASTA(path+file)
    print "saving to "+path+file
    ffile.create()
    for c in contig_list: ffile.add_seq(c.record)
    ffile.close()


for i in range(0,3):
    write_contig_list(contig_lists[i],fpath[i])

with open('/home/moritz/people/sari/assemblies/ray/concoct-output_clustering_gt1000.csv', 'rb') as csvfile:
    bins = csv.reader(csvfile, delimiter=',', quotechar='|')
    binning={row[0] : row[1] for row in bins}

bins={s:[contig for contig in binning if binning[contig] == s] for s in set(binning.viewvalues())}

with open('/home/moritz/people/sari/assemblies/ray/Contigs.fasta.dir/phylotype.result', 'rb') as csvfile:
    phylo = csv.reader(csvfile, delimiter='\t', quotechar='|')
    phylotyping={row[0] : {'marker': row[1], 'phylo': [row[r] for r in range(2,len(row)-1)]}  for row in phylo if row[0]!='Query'}

temp= [ mark.split("_") for mark in  phylotyping.keys()]
markermao={contig:[] for contig in set([cs[0] for cs in temp ]) }

for k in markermao:
    t=['_'.join(m) for m in temp if m[0]==k]
    markermao[k]=[{p:phylotyping[p]} for p in phylotyping if p in t]

binned_markers={}
for bin_id in bins:
    binned_markers[bin_id]=[markermao[m] for m in markermao if m in bins[bin_id]]
    binned_markers[bin_id]=[item for sublist in binned_markers[bin_id] for item in sublist]

def write_marker_list(list_of_markers,path,name):
    ffile = open(path+name, 'w')
    buff = csv.writer(ffile)
    data=[]
    print 'saving markers to '+path+name
    for m in list_of_markers:
        line=[]
        line.append(m.keys()[0])
        m=m[m.keys()[0]]
        line.append(m['marker'])
        line.extend(m['phylo'])
        data.append(line)
    buff.writerows(data)
    ffile.close()
    
    
for bin_id in binned_markers:
    write_marker_list(binned_markers[bin_id],'/home/moritz/people/sari/assemblies/','markers_of_bin_'+bin_id+'.csv')

for bin_id in bins:
    write_contig_list(bins[bin_id],'/home/moritz/people/sari/assemblies/','contigs_of_bin_'+bin_id+'.csv')

    
