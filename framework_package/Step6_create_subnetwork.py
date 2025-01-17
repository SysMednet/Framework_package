import os, sys
import numpy as np
from scipy import stats
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

def create_subnetwork(file_e, pvalue, file_p, rate, save_path) :

    knee_point = {}
    with open(f'{save_path}/{pvalue}_rwr{rate}_knee_point.txt','r') as f1 :
        for line in f1 :
            i, j = line.strip().split('\t')
            knee_point[i]=j
    
    zscore = stats.norm.isf(float(pvalue)/2)
    
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
        
        genes = []
        with open(f'{save_path}/{p}_{pvalue}_rwr{rate}_score.txt','r') as f2 :
            for line in f2 :
                if float(line.strip().split('\t')[1]) >= float(knee_point[f'{p}']) :
                    genes.append(line.strip().split('\t')[0])
               
        with  open(f'{save_path}/{p}_{pvalue}_rwr{rate}_subnetwork.txt','w') as w :
        
            csc = scipy.sparse.load_npz(f'{save_path}/{p}_zscore.npz')
    
            json_file = f"{p}_zscore_index.json"
            json_path = os.path.join(f'{save_path}', json_file)
            if os.path.exists(json_path):
                with open(json_path, 'r') as f:
                    gene_map = json.load(f)
            else:
                print(f"{json_file} not found!")
    
            gene_map_reversed = {v: k for k, v in gene_map.items()}
    
            G = nx.from_scipy_sparse_array(csc, create_using=nx.Graph) 
            edges = [
                (gene_map_reversed[u], gene_map_reversed[v], d['weight'])
                for u, v, d in G.edges(data=True) if u < v ]  
            
            for nline in edges :
                i, j, k = nline[0],nline[1],nline[2]
                if np.abs(float(k)) >= float(zscore) :
                    if nline[0] in genes and nline[1] in genes :
                        w.write(nline[0]+'\t'+nline[1]+'\t'+str(nline[2])+'\n')
        
        
        # with open(f'{save_path}/{p}_zscore.txt','r') as f3 , open(f'{save_path}/{p}_{pvalue}_rwr{rate}_subnetwork.txt','w') as w :
        #     _ = f3.readline()
        #     for line in f3 :
        #         i, j, k = line.strip('\n').split('\t')
        #         if np.abs(float(k)) >= float(zscore) :
        #             if line.strip().split('\t')[0] in genes and line.strip().split('\t')[1] in genes :
        #                 w.write(line)
    # print('Finish')                    