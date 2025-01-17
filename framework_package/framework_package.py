import os 
import Step1_SWEET_sample_weight_calculating_HNSCC as SWEET1
import Step2_SWEET_edge_weight_calculating_HNSCC as SWEET2
import Step3_SWEET_calculating_mean_std_zscore_HNSCC as SWEET3
import Step4_RWR_algorithm as RWR
import Step5_knee_point as KNEE
import Step6_create_subnetwork as SUBNETWORK
import Step7_NEEA_hyper as NEEA

def framework(GEM, 
              k=0.1, 
              output_edgescore_network='no',
              Samples=None, 
              Interest_genes=None,
              save_path=None, 
              pvalue=0.05,                
              rate=0.3) :
    
    #check GEM
    if not os.path.isfile(GEM): #確認檔案存在
        print('GEM file does not exist')
        return
    
    #check interest gene list file
    if Interest_genes == None: #確認檔案存在
        print('interest gene list file does not exist')
        return
    
    # check output folder
    if save_path != None:
        if not os.path.isdir(save_path): #確認資料夾存在
            print('Output folder does not exist')
            return
    else:
        print('Please input the output directory')
        return
    
        
    print('Constructing SWEET networks....')
    SWEET1.sweet1(GEM, k, save_path)
    SWEET2.sweet2(GEM, Samples, save_path)
    SWEET3.sweet3(GEM, Samples, save_path, output_edgescore_network)
    print('Complete\n')
    
    
    print('Calculate RWR algorithm....')
    RWR.rwr(GEM, pvalue, Samples, Interest_genes, rate, save_path)
    KNEE.knee_point(GEM, pvalue, Samples, rate, save_path)
    print('Complete\n')
    
    
    print('Constructing subnetworks....')
    SUBNETWORK.create_subnetwork(GEM, pvalue, Samples, rate, save_path)
    print('Complete\n')
    
    
    print('Calculating NEEA....')
    NEEA.neea(GEM, pvalue, Samples, rate, save_path)
    print('Complete\n')
    
    
    
    
    
    
    
    
    
    
    
    