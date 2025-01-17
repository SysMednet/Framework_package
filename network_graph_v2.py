import numpy as np
import os
import shutil
import pkg_resources
import pandas as pd
import argparse


def json_generate_for_pathway(pathway_name):
    # print(pathway_name)
    genes = pathways.get(pathway_name, [])
    op = '[\n'
    elements = []  

    # Add gene nodes
    for gene in genes:
        if gene in subnetwork_genes:
            elements.append(
                '\t' * 9 + '{"group":"nodes", "data": { "id": "' + gene + '","shape": "ellipse","color":"#FF9797"}}'
            )

    # Add edges between genes
    for i, gene in enumerate(genes):
        for j in range(i + 1, len(genes)):
            if genes[i] + '-' + genes[j] in subnetwork_edges:
                elements.append(
                    '\t' * 9 + '{"group":"edges", "data": { "id": "e' + str(i) + '_' + str(j) +
                    '", "source": "' + genes[i] + '", "target": "' + genes[j] + '" }}'
                )

    op += ',\n'.join(elements) + ']\n'
    # op += "]\n"
    # print(op)
    return op

def gene_set_network(DEG_list_result_df, criterion, output_dir):
    if criterion=='fdrq':
        legend_figure = './for_graphs/graph_files/legend_fdrq.png'
        enrich_DEGlist_df = DEG_list_result_df[DEG_list_result_df['FDR q-value']<0.05]
        enrich_DEGlist_df['p_log'] = [-np.log10(x) for x in enrich_DEGlist_df['FDR q-value']]
        rank_by_tag = '(selected by FDR <i>q</i> value)'
    else:
        legend_figure = './for_graphs/graph_files/legend_rawp.png'
        enrich_DEGlist_df = DEG_list_result_df[DEG_list_result_df['P-value']<0.05]
        enrich_DEGlist_df['p_log'] = [-np.log10(x) for x in enrich_DEGlist_df['P-value']]
        rank_by_tag = '(selected by raw <i>p</i> value)'
    

    total_select=[len(enrich_DEGlist_df)]
    function_list = ['one']

    html_selection = ''
    for index, row in enrich_DEGlist_df.iterrows():
        pathway_name = row['Pathway name']
        
        html_selection += f'<option value='+str(pathway_name.replace(',', '').replace(' ', '_').replace('-','_'))+' selected>'+str(pathway_name)+'</option>\n'

    
    # print(html_selection)
    #copy graph need files
    graph_files_database = './for_graphs/graph_files/'
    if not os.path.isdir(os.path.join('./for_graphs','graph_files')):
        shutil.copytree(graph_files_database, os.path.join('./for_graphs','graph_files'))
    
    #write html file
    fileo_html = open(os.path.join(f'./for_graphs/{sample}_{pvalue}_rwr{rate}_subnetwork_pathway_graph.html'),'w')
    ###html backbone
    filef1_1_html_path = './for_graphs/graph_code_fragment/html_f1_1.txt'
    filef1_1_html = open(filef1_1_html_path)
    html_f1_1 = filef1_1_html.read()
    filef1_1_html.close()
    
    filef1_2_html_path = './for_graphs/graph_code_fragment/html_f1_2.txt'
    filef1_2_html = open(filef1_2_html_path)
    html_f1_2 = filef1_2_html.read()
    filef1_2_html.close()
    
    filef2_html_path = './for_graphs/graph_code_fragment/html_f2.txt'
    filef2_html = open(filef2_html_path)
    html_f2 = filef2_html.read()
    filef2_html.close()
    
    filef3_html_path = './for_graphs/graph_code_fragment/html_f3.txt'
    filef3_html = open(filef3_html_path)
    html_f3 = filef3_html.read()
    filef3_html.close()
    
    fileo_html.write(html_f1_1)
    # fileo_html.write(legend_figure)
    fileo_html.write(f'    <script src="graph_files/script{sample}_{pvalue}_{rate}.js"></script>\n')
    fileo_html.write(html_f1_2)
    fileo_html.write(html_f2)
    fileo_html.write(html_selection)
    fileo_html.write('\t</select>\n\tpathway\n')
    fileo_html.write('\t&emsp;&emsp;\n\t<button onclick="reset()">Refresh graph</button>\n')
    fileo_html.write('\t<br />\n\t'+rank_by_tag+'\n')
    fileo_html.write(html_f3)
    fileo_html.close()

    #write js file
    fileo_js = open(os.path.join('./for_graphs',f'graph_files/script{sample}_{pvalue}_{rate}.js'),'w')
    ###js backbone
    filef1_js_path = './for_graphs/graph_code_fragment/js_f1.txt'
    filef1_js = open(filef1_js_path)
    js_f1 = filef1_js.read()
    filef1_js.close()
    
    filef2_js_path = './for_graphs/graph_code_fragment/js_f2.txt'
    filef2_js = open(filef2_js_path)
    js_f2 = filef2_js.read()
    filef2_js.close()

    filet1_cyto_path = './for_graphs/graph_code_fragment/cyto_t1.txt'
    filet1_cyto = open(filet1_cyto_path)
    cyto_t1 = filet1_cyto.read()
    filet1_cyto.close()
    
    filet2_cyto_path = './for_graphs/graph_code_fragment/cyto_t2.txt'
    filet2_cyto = open(filet2_cyto_path)
    cyto_t2 = filet2_cyto.read()
    filet2_cyto.close()
    
    fileo_js.write(js_f1)

    c=0
    for index, row in enrich_DEGlist_df.iterrows():
        pathway_name = row['Pathway name']
        fileo_js.write('\t' * 4 + pathway_name.replace(',', '').replace(' ', '_').replace('-','_') + cyto_t1)
        json_data = json_generate_for_pathway(pathway_name)
        fileo_js.write(json_data)
        fileo_js.write(cyto_t2)
        c+=1
        if c<len(enrich_DEGlist_df):
            fileo_js.write(',')
        fileo_js.write('\n')
    fileo_js.write(js_f2)
    # fileo_js.write('}\n')
    fileo_js.close()
    
#%%
parser = argparse.ArgumentParser(description="Manual")
parser.add_argument("-c", type=str, default="0.05" , help='network pvalue cutoff value')# network pvalue cutoff value
parser.add_argument("-cri", type=str, default="pv" , help='hyper p-value')# enrich hyper p-value
parser.add_argument("-sample", type=str , default="13-IT" , help="sample name")# sample name
parser.add_argument("-r", type=str , default="0.3" , help="restart rate")# 1-restart rate
parser.add_argument("-s", type=str , default="./example_output" , help="save path")# save path

args = parser.parse_args()
pvalue = args.c
cri = args.cri
sample = args.sample
rate = args.r
save_path = (args.s).rstrip('/')

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

subnetwork_genes = set()
subnetwork_edges = set()
with open(f'{save_path}/{sample}_{pvalue}_rwr{rate}_subnetwork.txt','r') as f :
    for line in f:
        g1, g2, w = line.strip().split('\t')
        # g1 = g1.split('_')[1]
        # g2 = g2.split('_')[1]
        subnetwork_genes.add(g1)
        subnetwork_genes.add(g2)
        subnetwork_edges.add(g1+'-'+g2)
        subnetwork_edges.add(g2+'-'+g1)


pathways = {}
input_dir = save_path.replace('output','input')
with open(f'{input_dir}/KEGG_api_347_hsa_pathway_20220518_v102.txt','r') as f :
    for line in f :
        id, name, genes = line.strip().split('\t')
        for g in genes.split(',') :
            if not str(next(iter(subnetwork_genes))).startswith('ENS') :
                if name not in pathways :
                    pathways[name] = [g]
                else :
                    pathways[name].extend([g])
            else :
                if g in data_e :
                    for k in data_e[g] :
                        if name not in pathways :
                            pathways[name] = [k]
                        else :
                            pathways[name].extend([k])



DEG_list_result_df = pd.read_csv(f'{save_path}/neea_results/{sample}_{pvalue}_rwr{rate}_edge_based_hyper.txt',sep='\t')
criterion = f'{cri}'
output_dir = f'{save_path}'
gene_set_network(DEG_list_result_df, criterion, output_dir)


