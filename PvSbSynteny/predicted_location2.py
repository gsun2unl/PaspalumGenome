refn = 'sorghum3'

class bed_object:
    def __init__(self,afile):
        self.genes = {}
        self.positions = []
        fh = open(afile)
        count = 0
        for x in fh:
            y = x.strip().split('\t')
            mygene = y[3]
            self.genes[mygene] = y[:]
            self.positions.append(mygene)
        self.positions.sort(key=lambda a:int(self.genes[a][1]))
        self.positions.sort(key=lambda a:self.genes[a][0])
        for xind,x in enumerate(self.positions):
            self.genes[x].append(xind)
            
class last_object:
    def __init__(self,afile):
        self.mylastfile = afile
    def gene_hits(self,query_gene,target_genes):
        import subprocess as sp
        import os
        tempfile = open("temp_hits2.lastz",'w')
        proc = sp.Popen(['grep',query_gene,self.mylastfile],stdout=tempfile)
        proc.wait()
        tempfile.close()
        tempfile = open("temp_hits2.lastz")
        hit_dict = {}
        for x in tempfile:
            if not x: continue
            y = x.strip().split('\t')
            if len(y) < 3: continue
#            print y
            #if y[7] != '+': continue
            myscore = float(y[-1])
            mygene = y[0] #target gene, query gene is sorghum
            if '.' in mygene[-4:]:
                mygene = '.'.join(mygene.split('.')[:-1])
            if not mygene in target_genes: continue
            hit_dict[mygene] = myscore
        return hit_dict

class synteny_object:
    def __init__(self,afile):
        self.lookup = {}
        fh = open(afile)
        aline = fh.readline()
        self.names_list = aline.strip().split(',')[1:]
        for x in fh:
            y = x.strip().split(',')
            mygene = y[0]
            self.lookup[mygene] = {}
            for other_gene,other_genome in zip(y[1:],self.names_list):
                if other_gene == 'No Gene':
                    self.lookup[mygene][other_genome] = False
                else:
                    self.lookup[mygene][other_genome] = other_gene
 
def position_lookup(mystart,mystop,myind,synteny_file,genome,tgene_list):
    for x in range(mystart,mystop,myind):
        if x < 0: continue
        if x >= len(tgene_list): continue
        mygene = tgene_list[x]
        if not synteny_file.lookup[mygene][genome]: continue
        return synteny_file.lookup[mygene][genome]
    return ''

def anchor_search(synteny_file,genome,position,genome_data_sets,bounds=20):
#    if position < bounds: position = bounds
    tgene_list = genome_data_sets[refn]['Bed'].positions
    lower_bound = position_lookup(position+1,position+bounds,1,synteny_file,genome,tgene_list)
    upper_bound = position_lookup(position-1,position-bounds,-1,synteny_file,genome,tgene_list)
    target_genes_set = set([])
    if lower_bound:
        lb_data = genome_data_sets[genome]['Bed'].genes[lower_bound]
        geno_max = len(genome_data_sets[genome]['Bed'].positions)
        lpos = lb_data[-1]
        lchr = lb_data[0]
        for x in range(lpos-bounds,lpos+bounds):
            if x >= geno_max: continue
            candidate_gene = genome_data_sets[genome]['Bed'].positions[x]
            if genome_data_sets[genome]['Bed'].genes[candidate_gene][0] != lchr: continue
            target_genes_set.add(candidate_gene)
    if upper_bound:
        lb_data = genome_data_sets[genome]['Bed'].genes[upper_bound]
        geno_max = len(genome_data_sets[genome]['Bed'].positions)
        lpos = lb_data[-1]
        lchr = lb_data[0]
        for x in range(lpos-bounds,lpos+bounds):
            if x >= geno_max: continue
            candidate_gene = genome_data_sets[genome]['Bed'].positions[x]
            if genome_data_sets[genome]['Bed'].genes[candidate_gene][0] != lchr: continue
            target_genes_set.add(candidate_gene)
    return target_genes_set

import sys
from pyfasta import Fasta

naive_pairs = sys.argv[1]
data_files = sys.argv[2]

fh = open(data_files)
fh.readline()
gfiles = {}
for x in fh:
    y = x.strip().split(',')
    mygenome = y[0]
    gfiles[mygenome] = {'Bed':bed_object(y[1]),'CDS':Fasta(y[2]),'Genome':Fasta(y[3]),'lastz':last_object(y[4])}
#    except: print mygenome
my_syn_obj = synteny_object(naive_pairs)
plist = [refn]
plist.extend(my_syn_obj.names_list)
plist.append("GEvo Link")
print ",".join(plist)

#print list(gfiles)
#1/0

for position,gene in enumerate(gfiles[refn]['Bed'].positions):
    plist = [gene]
    for genome in my_syn_obj.names_list:
        target_genes = anchor_search(my_syn_obj,genome,position,gfiles)
        hit_dict = gfiles[genome]['lastz'].gene_hits(gene,target_genes)
        if len(hit_dict) == 0: 
            plist.append('No Gene')
            continue
        candidates = list(hit_dict)
        if len(hit_dict) == 1:
            plist.append(candidates[0])
            continue
        candidates.sort(key=lambda a:hit_dict[a])
        plist.append(candidates[-1])
    proofing_link = []
    num_seqs = 0
    for x in plist:
        if x == 'No Gene': continue
        num_seqs += 1
        proofing_link.append("accn{0}={1}".format(num_seqs,x))
    proofing_link.append("num_seqs={0};autogo=1".format(num_seqs))
    plist.append("http://genomevolution.org/CoGe/GEvo.pl?" + ";".join(proofing_link))
    print ",".join(plist)
    
