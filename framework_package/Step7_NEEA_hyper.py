import sys,os
import numpy as np
from scipy import stats
from scipy.stats import hypergeom
import statsmodels.stats.multitest
from collections import defaultdict
import scipy.sparse
import json
import networkx as nx

def check_file(expres):
    checkset = set(["", "NA", "Na", "na", "nan", "null"])
    for c in checkset:
        loc = np.where(expres == c)
        if loc[0].size:
            expres[loc] = "0"
            print(f"There is {c} in the 'gene expression matrix' file and it will be assigned to 0.")
    return expres

### ez -> ens
data_e = {}
with open("./example_input/geneid_ens2ez.txt",'r')as f:
    _g = f.readline()
    for line in f:
        g_e, ent = line.strip('\n').split(',')     
        g_e = [(g_e)]
        if ent == '':
            continue
        else:   
            if ent in data_e.keys() :
                # ent = ent.split()    
                data_e[ent].extend(g_e)
            else :
                data_e[ent] = g_e
# Pathway edges (ens)
def calculate_pathway_edges_ens(pathway_file):
    # Record genes
    pathway_gene = defaultdict(list)
    with open(pathway_file, 'r') as input_term:
        for line_term in input_term:
            id, name, genes = line_term.strip().split('\t')
            for g in genes.split(','):
                if g in data_e :
                    for k in data_e[g] :
                        pathway_gene[name].append(k)
    # Record edges
    pathway_edge = defaultdict(list)
    for i,j in pathway_gene.items():
            for x in j:
                for y in j:
                    if(x != y):
                        pathway_edge[i].append(f'{x}-{y}') #including A-B, B-A
    return pathway_edge

# Pathway edges (ez)
def calculate_pathway_edges_ez(pathway_file):
    # Record genes
    pathway_gene = defaultdict(list)
    with open(pathway_file, 'r') as input_term:
        for line_term in input_term:
            id, name, genes = line_term.strip().split('\t')
            for g in genes.split(','):
                pathway_gene[name].append(g)
    # Record edges
    pathway_edge = defaultdict(list)
    for i,j in pathway_gene.items():
            for x in j:
                for y in j:
                    if(x != y):
                        pathway_edge[i].append(f'{x}-{y}') #including A-B, B-A
    return pathway_edge
                
# Calculate network edges (N,n)

def calculate_network_edges(load_path,pvalue):
    
    zscore = stats.norm.isf(float(pvalue)/2)
    
    node_set = set()
    edge_set = set()
    
    if load_path.endswith('.npz') :
        csc = scipy.sparse.load_npz(load_path)
    
        json_path = load_path.replace('.npz','_index.json')
        with open(json_path, 'r') as f:
            gene_map = json.load(f)
    
        gene_map_reversed = {v: k for k, v in gene_map.items()}
    
        G = nx.from_scipy_sparse_array(csc, create_using=nx.Graph) 
        edges = [
            (gene_map_reversed[u], gene_map_reversed[v], d['weight'])
            for u, v, d in G.edges(data=True) if u < v ]  
        
        for nline in edges :
            g1, g2, w = nline[0],nline[1],nline[2]
            if np.abs(float(w)) >= float(zscore) :
                node_set.add(g1)
                node_set.add(g2)
                edge_set.add(g1+'-'+g2)
                edge_set.add(g2+'-'+g1)
    
    else :
        with open(load_path,'r') as f :
            _ = f.readline()
            for line in f :
                g1, g2, w = line.strip('\n').split('\t')
                # g1=g1.split('_')[0]
                # g2=g2.split('_')[0]
                if np.abs(float(w)) >= float(zscore) :
                    node_set.add(g1)
                    node_set.add(g2)
                    edge_set.add(g1+'-'+g2)
                    edge_set.add(g2+'-'+g1)
                
    num_nodes = len(node_set)
    num_edges = len(edge_set)/2 # A-B, B-A only need to calculate in one time
    return num_edges, edge_set

# Calculate overlapping edges (M,m)
def calculate_overlapping_edges(network_edge, pathway_edge):
    network_overlapping = {}
    for i,j in pathway_edge.items():
        if len(set(j) & network_edge) > 0: overlapping_edges = len(set(j) & network_edge)/2 # A-B, B-A only need to calculate in one time
        else: overlapping_edges = 0
        network_overlapping[i] = overlapping_edges
    return network_overlapping
    
# Hypergeometric p-value
def Record(sample, save, N, network_overlapping, n, subnetwork_overlapping):
    m1 , m3,  m4 , m5 , m6 , m7  = [] , [] , [] , [] , [] , [] 
    with open(save, 'w') as w_line :
        w_line.write('Pathway name\tP-value\tFDR q-value\tN\tM\tn\tm\n')
        for i,j in zip(network_overlapping.items(),subnetwork_overlapping.items()):
            M = i[1]
            m = j[1]
            hyper = hypergeom.sf(m-1,N,M,n)
            m1.append(i[0])
            m3.append(hyper)
            m4.append(N)
            m5.append(M)
            m6.append(n)
            m7.append(m)
        if len(m3) == 0:
            print("No pathways!")
            return False
        m8 = statsmodels.stats.multitest.multipletests(m3,method='fdr_bh',alpha=0.05)[1]
        for result in zip(m1 , m3 , m8 , m4 , m5 , m6 , m7) :
            w_line.write('\t'.join(str(s) for s in result)+'\n')

# if __name__ == '__main__':
    
def neea(file_e, pvalue, file_p, rate, save_path) : 
    
    #open patient file
    patlist = []
    if file_p != None :                    
        with open(file_p,mode='r') as rline :
            for nline in rline :
                tem = nline.strip('\n').split('\t')
                patlist.append(tem[0])
        del rline , nline , tem
        if not patlist :
            print("patient file doesn't have patients")
            sys.exit()
    else :
        with open(file_e,mode='r') as rline :
            pat = rline.readline().strip('\n').split('\t')[1:]
            for i in pat :
                patlist.append(i)

    for p in patlist[0:]:
        
        # Pathway
        pathway_file = "./example_input/KEGG_api_347_hsa_pathway_20220518_v102.txt"
        
        
        # Whole network
        # N
        ori_load_path=f'{save_path}/{p}_zscore.npz'
        N, edge_set = calculate_network_edges(ori_load_path,pvalue)
        
        if str(next(iter(edge_set))).startswith('ENS') :
            pathway_edge = calculate_pathway_edges_ens(pathway_file)
        elif not str(next(iter(edge_set))).startswith('ENS') :
            pathway_edge = calculate_pathway_edges_ez(pathway_file)
        
        # M
        network_overlapping = calculate_overlapping_edges(edge_set, pathway_edge)

        # sub-network 
        # n
        sub_load_path = f'{save_path}/{p}_{pvalue}_rwr{rate}_subnetwork.txt'
        n, sub_edge_set = calculate_network_edges(sub_load_path,pvalue)
        
        # m
        subnetwork_overlapping = calculate_overlapping_edges(sub_edge_set, pathway_edge)

        if not os.path.isdir(f'{save_path}/neea_results'):
            os.makedirs(f'{save_path}/neea_results')
        
        save = f'{save_path}/neea_results/{p}_{pvalue}_rwr{rate}_edge_based_hyper.txt'
        Record(p, save, N, network_overlapping, n, subnetwork_overlapping)
        print(f'{p} NEEA Done!')

                