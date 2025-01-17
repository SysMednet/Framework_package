import numpy as np
import sys
import os
import networkx as nx
import pandas as pd
from scipy import stats
import scipy.sparse
import json


def check_file(expres):
    checkset = set(["", "NA", "Na", "na", "nan", "null"])
    for c in checkset:
        loc = np.where(expres == c)
        if loc[0].size:
            expres[loc] = "0"
            print(f"There is {c} in the 'gene expression matrix' file and it will be assigned to 0.")
    return expres

def edge_to_adjacency_matrix(edge_list, num_nodes):
    adjacency_matrix = np.zeros((num_nodes, num_nodes), dtype=bool)
    for edge in edge_list:
        i, j = edge
        adjacency_matrix[i, j] = True
        adjacency_matrix[j, i] = True
    return adjacency_matrix

def create_gene_weight_dictionary(nodes,genelist,genum):
    dictionary = {}
    for g in nodes :
        g = genum[g]
        if g in genelist :
            dictionary[g]=1
        else :
            dictionary[g]=0
    return dictionary

def calculate_pagerank(G, dictionary,rate):
    rwr_score = nx.pagerank(G, personalization=dictionary, alpha = (1-float(rate)))
    return rwr_score


def rwr(file_e, pvalue, file_p, file_i, rate, save_path) :

    gem_file = file_e
    gem = pd.read_csv(gem_file,sep='\t',index_col=0)
    
    genum={}
    genumr={}
    for n,m in zip(range(len(gem.index)),gem.index) :
        # m = m.split('_')[0]
        genum[m]=n
        genumr[n]=m
    num_nodes = len(gem.index)
        
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
    
    
    zscore = stats.norm.isf(float(pvalue)/2)
    
    for p in patlist[0:]:
        
        file = f"{save_path}/{p}_zscore.npz"
        if not os.path.exists(file):
            print(f"{file} not found")
            sys.exit()

        with open(f'{save_path}/{p}_{pvalue}_rwr{rate}_score.txt','w') as w:
            
            edge_list = []
            nodes = set()
            
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
                if np.abs(nline[2]) >= float(zscore) :
                    nodes.add(nline[0])
                    nodes.add(nline[1])
                    I = genum[nline[0]]
                    J = genum[nline[1]]
                    edge_list.append((I,J))
            
            # with open(f'{save_path}/{p}_zscore.txt', 'r') as f:
            #     _ = f.readline()
            #     for line in f:
            #         i, j, k = line.strip('\n').split('\t')
            #         # i = i.split('_')[0]
            #         # j = j.split('_')[0]
            #         if np.abs(float(k)) >= float(zscore) :
            #             nodes.add(str(i))
            #             nodes.add(str(j))
            #             I = genum[str(i)]
            #             J = genum[str(j)]
            #             edge_list.append((I, J))
            
            adj_matrix = edge_to_adjacency_matrix(edge_list, num_nodes)
            G = nx.from_numpy_array(adj_matrix)
            # print('nodes='+str(len(G.nodes)))
            # print('edges='+str(len(G.edges)))
            
            genelist = []
            with open(f'{file_i}','r') as f :
                for line in f :
                    if line.strip() in genum :
                        gene = genum[line.strip()]
                        genelist.append(gene)
            
            dictionary = create_gene_weight_dictionary(nodes,genelist,genum)
            
            rwr_score = calculate_pagerank(G, dictionary, rate)
            
            for s in rwr_score :
                w.write(str(genumr[s])+'\t'+str(rwr_score[s])+'\n')
    
