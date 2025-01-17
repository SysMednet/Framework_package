import numpy as np
import pandas as pd
import sys
import os
import scipy.sparse
import json
import networkx as nx

def check_file(expres):
    checkset = set(["", "NA", "Na", "na", "nan", "null"])
    if expres.dtype != object:
        expres = expres.astype(str)
    for c in checkset:
        loc = np.where(expres == c)
        if loc[0].size:
            expres[loc] = "0"
            print(f"There is {c} in the 'gene expression matrix' file and it will be assigned to 0.")
    return expres

def sweet3(file_e, file_p, save, output_network) :

    #open patient file
    patlist = []
    gene = []
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
            for nline in rline :
                g , *v = nline.strip('\n').split('\t')
                # if g in geneset :
                gene.append(g)
    
    geneset = set()
    pair = []
    
    # file = f"{save}/{patlist[0]}.npz"
    # if not os.path.exists(file):
    #     print(f"{file} not found")
    #     sys.exit()
    
    # # load npz file
    # csc = scipy.sparse.load_npz(f'{save}/{patlist[0]}.npz')

    # json_file = f"{patlist[0]}_index.json"
    # json_path = os.path.join(f'{save}', json_file)
    # if os.path.exists(json_path):
    #     with open(json_path, 'r') as f:
    #         gene_map = json.load(f)
    # else:
    #     print(f"{json_file} not found!")

    # gene_map_reversed = {v: k for k, v in gene_map.items()}

    # G = nx.from_scipy_sparse_array(csc, create_using=nx.Graph) 
    # edges = [
    #     (gene_map_reversed[u], gene_map_reversed[v], d['weight'])
    #     for u, v, d in G.edges(data=True) if u < v ]  
    
    # for nline in edges :
    #     geneset.add(nline[0]+'\t'+nline[1])
    #     pair.append(nline[2])
        
    # with open(f"{save}/{patlist[0]}.txt", mode='r') as rline:
    #     _ = rline.readline()
    #     for nline in rline:
    #         if nline != '\n':
    #             val = nline.strip('\n').split('\t')
    #             geneset.add(val[0]+'\t'+val[1])
    #             pair.append(val[2])
    
    for p in patlist[0:]:
        file = f"{save}/{p}.npz"
        if not os.path.exists(file):
            print(f"{file} not found")
            sys.exit()
        
        # load npz file
        csc = scipy.sparse.load_npz(f'{save}/{p}.npz')

        json_file = f"{p}_index.json"
        json_path = os.path.join(f'{save}', json_file)
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
            # if (nline[0]+'\t'+nline[1]) not in geneset:
            #     print(f"Warning! In the sample {p}, there are gene pair(s) that cannot be found in the other samples.")
            pair.append(nline[2])
            
        # file = f"{save}/{p}.txt"
        # if not os.path.exists(file):
        #     print(f"{file} not found")
        #     sys.exit()
        # with open(f"{save}/{p}.txt", mode='r') as rline:
        #     _ = rline.readline()
        #     for nline in rline:
        #         if nline != '\n':
        #             val = nline.strip('\n').split('\t')
        #             if (val[0]+'\t'+val[1]) not in geneset:
        #                 print(f"Warning! In the sample {p}, there are gene pair(s) that cannot be found in the other samples.")
        #             pair.append(val[2])
    
    pair = np.array(pair)
    pair = check_file(pair)
    pair = pair.astype(float)
    vmean, vstd = np.mean(pair), np.std(pair)
    
    with open(f'{save}/mean_std.txt', mode='w') as wline:
        wline.write(f"mean\t{vmean}\nstd\t{vstd}\n")
    
    for p in patlist:
        
        csc = scipy.sparse.load_npz(f'{save}/{p}.npz')
        
        csc = csc.tocoo()  
        values = csc.data
        values = np.abs((values - vmean) / vstd)  # Z-score 
        csc.data = values
        
        csc_upper = scipy.sparse.triu(csc, k=1)

        json_file = f"{p}_index.json"
        json_path = os.path.join(f'{save}', json_file)
        with open(json_path, 'r') as f:
            gene_map = json.load(f)

        gene_map_reversed = {v: k for k, v in gene_map.items()}
    
        matrix_output_path = f"{save}/{p}_zscore.npz"
        scipy.sparse.save_npz(matrix_output_path, csc_upper)
    
        index_output_path = f"{save}/{p}_zscore_index.json"
        with open(index_output_path, 'w') as f:
            json.dump(gene_map, f)
        
        # Output network txt file
        if output_network == 'yes' :
            # write txt file
            edge_pairs = []
            with open(f"{save}/{p}_zscore.txt", mode='w') as wline:
                wline.write('gene1\tgene2\tz_score\n')
                for nline in edges :
                    z = str((float(nline[2])-vmean)/vstd)
                    wline.write(f'{nline[0]}\t{nline[1]}\t{z}\n')
                    edge_pairs.append((nline[0], nline[1], z))
         
        print(f'{p} Done!\n')
        
        
        
        
        # file = f"{save}/{p}.txt"
        # if not os.path.exists(file):
        #     print(f"{file} not found")
        #     sys.exit()
        # with open(f"{file}", mode='r') as rline, open(f"{save}/{p}_zscore.txt", mode='w') as wline:
        #     _ = rline.readline()
        #     wline.write('gene1\tgene2\tz_score\n')
        #     for nline in rline:
        #         if nline != '\n':
        #             val = nline.strip('\n').split('\t')
        #             z = str((float(val[2])-vmean)/vstd)
        #             wline.write(f'{val[0]}\t{val[1]}\t{z}\n')
    
    # print("Finish")