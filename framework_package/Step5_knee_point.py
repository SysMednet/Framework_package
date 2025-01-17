import os, sys
import numpy as np
from kneed import KneeLocator

def check_file(expres):
    checkset = set(["", "NA", "Na", "na", "nan", "null"])
    for c in checkset:
        loc = np.where(expres == c)
        if loc[0].size:
            expres[loc] = "0"
            print(f"There is {c} in the 'gene expression matrix' file and it will be assigned to 0.")
    return expres

def knee_point(file_e, pvalue, file_p, rate, save_path) :

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
    
    file = f"{save_path}/{patlist[0]}_{pvalue}_rwr{rate}_score.txt"
    if not os.path.exists(file):
        print(f"{file} not found")
        sys.exit()
    
    xlist,ylist=[],[]
    with open(f'{save_path}/{pvalue}_rwr{rate}_knee_point.txt','w') as w :    
        for p in patlist[0:]:
            
            rwr_scores = {}
            scores = []
            with open(f'{save_path}/{p}_{pvalue}_rwr{rate}_score.txt','r')as f :
                for line in f :
                    g, val = line.strip().split('\t')
                    rwr_scores[g]=float(val)
                    scores.append(float(val))
    
            scores.sort(reverse=True)
                    
            gene_nums = np.arange(1,len(scores)+1,1) 
    
            node_scores_list = [float(scores[0])]
            for i in range(1, len(scores)) :
                node_scores_list.append(float(node_scores_list[i-1]) + float(scores[i]))
            
            node_scores_list_sort = list(node_scores_list)
            x = list(n for n in range(1,len(node_scores_list_sort)+1))
            # x = list(gene_nums)
            y = node_scores_list_sort
            kneedle = KneeLocator(x, y, S=1, curve='concave', direction='increasing')
            xlist.append(kneedle.knee)
            ylist.append(kneedle.knee_y)
            # find cutoff
            knee_point = kneedle.knee_y
            index = node_scores_list_sort.index(knee_point)
            value = scores[index]
            w.write(f"{p}\t{value}\n")
    