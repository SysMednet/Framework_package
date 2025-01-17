import numpy as np
import sys, os
import pandas as pd
import scipy.sparse
import json
from scipy.sparse import csr_matrix

def check_file(expres) :
    checkset = set(["","NA","Na","na","nan","null"])
    for c in checkset :
        loc = np.where(expres==c)
        if loc[0].size :
            expres[loc] = "0"
            print(f"{c} in expres and transfer to 0")
    return expres

def sweet2(file_e, file_p, save_path) :

    #open patient file
    patset = set()
    if file_p != None :                    
        with open(file_p,mode='r') as rline :
            for nline in rline :
                tem = nline.strip('\n').split('\t')
                patset.add(tem[0])
        del rline , nline , tem
        if not patset :
            print("patient file doesn't have patients")
            sys.exit()
    else :
        with open(file_e,mode='r') as rline :
            pat = rline.readline().strip('\n').split('\t')[1:]
            for i in pat :
                patset.add(i)
        
    #open expression file
    gene , value = [] , []
    with open(file_e,mode='r') as rline :
        pat = rline.readline().strip('\n').split('\t')[1:]
        for nline in rline :
            g , *v = nline.strip('\n').split('\t')
            # if g in geneset :
            value += v
            gene.append(g)
    del rline , nline , g , v
    patlen , genelen = len(pat) , len(gene)
    if (not patlen) :
        print("expression file are empty")
        sys.exit()    
    
    #open weight file   
    weight = {}
    with open(f'{save_path}/'+file_e.split('/')[-1].strip('.txt')+'_weight.txt',mode='r') as rline :
        _ = rline.readline()
        for nline in rline :
            p , w , *_ = nline.strip('\n').split('\t')
            weight.update({p:float(w)})
    del rline , nline , p , w , _
    if not weight :
        print("weight file doesn't have pateint")
        sys.exit()
    
    #check the pateints and genes in expression file
    patloc = [ l for l,p in enumerate(pat) if p in patset]
    if (not genelen) or (len(patloc)!=len(patset)) :
        print("expression file doesn't map to pateint or gene file")
        sys.exit()
    if len(set(pat)&weight.keys()) != patlen :
        print("expression file doesn't map to weight file")
        sys.exit()
    del patset 
    print(f"patient : {len(patloc)}\ngene : {genelen}")
    
    value = np.array(value).reshape(genelen,patlen)
    value = check_file(value)
    # value = np.array(value,dtype=float)
    value = value.astype(float)
    
    agg = np.corrcoef(value)
    
    for l in patloc :
        p = pat[l]
        value_s = np.c_[value,value[:,l]]
        value_s = np.corrcoef(value_s)
        value_s = weight[p] * (value_s - agg) + agg
        
        # with open(f"{save_path}/{p}.txt",mode='w') as wline :
        #     wline.write("gene1\tgene2\tedge_weight\n")
        #     for l , g1 , v1 in zip(range(genelen),gene,value_s) :        
        #         wline.write('\n'.join(g1+'\t'+g2+'\t'+str(v2) for g2, v2 in zip(gene[(l+1):], v1[(l+1):])))                        
        #         wline.write('\n')
                
        # edge_pairs = [ (g1, g2, v2) for i, (g1, v1) in enumerate(zip(gene, value_s)) for g2, v2 in zip(gene[(i + 1):], v1[(i + 1):])]
        # 
        # edge_pairs = (
        #     (g1, g2, v2)
        #     for i, (g1, v1) in enumerate(zip(gene, value_s))
        #     for g2, v2 in zip(gene[(i + 1):], v1[(i + 1):])
        # )
        
        # pair_file = pd.DataFrame(edge_pairs, columns=['gene1', 'gene2', 'edge_weight'])
        
        # all_genes = sorted(set(pair_file['gene1'].astype(str)).union(set(pair_file['gene2'].astype(str))))
        
        all_genes = sorted(set(gene))
        gene_map = {g: idx for idx, g in enumerate(all_genes)}
    
        # 創建稀疏矩陣
        # rows = pair_file['gene1'].astype(str).map(gene_map).values
        # cols = pair_file['gene2'].astype(str).map(gene_map).values
        # values = pair_file['edge_weight'].values
        # values = np.array(values, dtype=float)
        
        rows, cols = np.triu_indices_from(value_s, k=1)
        values = value_s[rows, cols]
        
        matrix = scipy.sparse.csr_matrix((values, (rows, cols)), shape=(len(all_genes), len(all_genes)))
        # matrix = scipy.sparse.triu(matrix, k=1)  # 取上三角 (排除對角線和下三角)
    
        # 儲存矩陣和索引字典
        matrix_output_path = f"{save_path}/{p}.npz"
        scipy.sparse.save_npz(matrix_output_path, matrix)
    
        # 儲存索引字典
        index_output_path = f"{save_path}/{p}_index.json"
        with open(index_output_path, 'w') as f:
            json.dump(gene_map, f)
            
    del value , agg , patloc , value_s , p , weight , l
    