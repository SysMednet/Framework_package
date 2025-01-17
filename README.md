# Framework_package
This framework package is a single-sample network-based framework to analyze RNA sequencing data from a cohort of HNSCC patients. It is designed to uncover biological differences in different tumor samples by analyzing differential gene expression and network structure changes across samples. 

1. Single-Sample Networks (SINs) Construction: Used the sample-specific weighted correlation network (SWEET) method to construct SINs.
For more information, you can visit:
[[Publication]](https://pmc.ncbi.nlm.nih.gov/articles/PMC10025435)
[[SWEET repository on GitHub]](https://github.com/SysMednet/SWEET/tree/main).

2. Subnetwork Extraction: Applied the random walk with restart (RWR) algorithm to extract subnetworks centered on user-specified gene sets.
3. Pathway Enrichment Analysis: Performed network edge-based enrichment analysis (NEEA) to identify enriched biological pathways in subnetworks.

## Framework
<img src="https://github.com/user-attachments/assets/f261edb8-91ce-48e3-b2e7-aed5dffc5172" width="400x900">

## Dependencies and requirement
The framework package code is written in python 3.9.17. Users also need to install following python package listed below:
  - numpy
  - pandas
  - shutil
  - scipy
  - statsmodels
  - matplotlib
  - argparse
  - json
  - networkx

## Installation
<p>To install the framework package, run:</p>
<pre><code>pip install git+https://github.com/HsinYu-Hsu/Framework_package</code></pre>
<p>If the error message "cannot find command git" was shown, please install "git".</p>
<pre><code>conda install git</code></pre>

<a name="variable-table"></a>
## Variable
### Framework function
<pre><code>def framework(GEM, 
              k=0.1, 
              output_edgescore_network='no',
              Samples=None, 
              Interest_genes=None,
              save_path=None, 
              pvalue=0.05,                
              rate=0.3)
</code></pre>

| Variable | Description | 
| ---- | ----- |
| GEM | A path to the "gene expression data" file. |
| k | Balance parameter. Default: 0.1. |
| output_edgescore_network | Decide whether to export the created network as a .txt file ('no', 'yes'). Default: 'no'. If 'no', only output the .npz file. |
| Samples | A path to the "sample list" file. If "None", calculate all samples. |
| Interest_genes | A path to the "interest gene list" file. |
| save_path | A path to the output files. |
| pvalue | The cutoff for constructing networks will be set to different p-values. Default: 0.05. |
| rate | The restart rate for calculating random walk with restart algorithm, which should be set between 0 and 1. Default: 0.3. |

## Gene expression data format
[Back to variable table](#variable-table)
   - Each elements should be delimited by tab (\t).
     - Columns: Samples.
     - Rows: Genes.
   - The gene ID should be entrez ID or ensembl ID.
   - As building a network demands substantial computational resources, it is recommended to avoid including too many genes.

| Gene | 13-IC | 13-N | ... |
| ---- | ----- | ---- | ---- |
| ENSG00000117152 | 1.526069 | 1.321928 | ... |
| ENSG00000179632 | 3.986411 | 3.643856 | ... |
| ENSG00000127314 | 7.681309 | 7.565978 | ... |
| ... | ... | ... | ... |

## Sample list format
[Back to variable table](#variable-table)
   - Each elements should be seperated with \n.
   - The sample name should be include in gene expression data.

| Sample |
| ---- |
| 13-IC |
| 13-N |
| ... |

## Interest gene list format
[Back to variable table](#variable-table)
  - Each elements should be seperated with \n.
  - Gene ID should be the same for gene expression data.

| Gene |
| ---- |
| ENSG00000117152 |
| ENSG00000179632 |
| ENSG00000127314 |
| ... |

## Usage and outputs
### Import and execute framework function with python

  - All example input files can be found in "example_input" folder.
  - All example output files can be found in "example_output" folder.

Every output files from this function will be saved in the folder set by "save_path" variable.
<pre><code>from framework_package import framework  
    
framework('./example_input/gene_expression_1000genes.txt',   
            k=0.1, 
            output_edgescore_network='yes',
            Samples='./example_input/samples.txt', 
            Interest_genes='./example_input/interest_genes.txt',
            save_path='./example_output', 
            pvalue=0.05,                
            rate=0.3)
</code></pre>

### Gene expression data file
Gene expression data file (Example file: "gene_expression_1000genes.txt" in "example_input" folder)

To reduce the consumption of computing time, the gene expression matrix file in the "example_input" folder contains only 1000 genes.

| Gene | 13-IC | 13-N | ... |
| ---- | ----- | ---- | ---- |
| ENSG00000117152 | 1.526069 | 1.321928 | ... |
| ENSG00000179632 | 3.986411 | 3.643856 | ... |
| ENSG00000127314 | 7.681309 | 7.565978 | ... |
| ... | ... | ... | ... |

### Sample weight file
Sample weight result output (Example file: "gene_expression_weight.txt" in "example_output" folder)

Provides weights for each sample based on the SWEET method.

| sample | sample_weight |
| ---- | ---- |
| 13-IT | 4.17586578629123 |
| 13-N | 4.302891217808787 |
| ... | ... |

### Mean and std file
Mean and std result output (Example file: "mean_std.txt" in "example_output" folder)

|  |  |
| ---- | ---- |
| mean | 0.08817455347984886 |
| std | 0.25319309455468386 |

### Sample subnetwork file
Sample subnetwork result output (Example file: "{sample}_{pvalue}_rwr{rate}_subnetwork.txt" in "example_output" folder)

Lists edges within extracted subnetworks and their z-scores.

| gene1 | gene2 | z-score |
| ---- | ---- | ---- |
| ENSG00000070610 | ENSG00000166387 | -2.9484170173090085 |
| ENSG00000070610 | ENSG00000172236 | -2.873478066857216 |
| ENSG00000117152 | ENSG00000166387 | -2.8527161492919006 |
| ... | ... | ... |

### Pathway enrichment analysis result
Pathway enrichment analysis result output (Example file: "{sample}_{pvalue}_rwr{rate}_edge_based_hyper.txt" in "example_output/neea_results" folder)

The result shows the significant associations between each subnetwork and KEGG pathway.

| Pathway name | P-value | FDR q-value | N | M | n | m |
| ---- | ---- | ---- | --- | --- | --- | --- |
| Glycolysis / Gluconeogenesis | 0.1625 | 1 | 12712 | 7 | 318 | 1 |
| Citrate cycle (TCA cycle) | 0.9151 | 1 | 12712 | 97 | 318 | 1 |
| Pentose phosphate pathway | 0.0121 | 1 | 12712 | 7 | 318 | 2 |
| ... | ... | ... | ... | ... | ... | ... |


### Network Graph

<pre><code>
python3 network_graph.py -sample {sample} -c 0.05 -r 0.3 -cri 'pv' -s ./example_output
</code></pre>

  - `-sample` : Sample name.
  - `-c` : The cutoff for constructing networks.
  - `-r` : The restart rate for calculating random walk with restart algorithm.
  - `-cri` : The criterion for filtered pathways. Choose p-value or fdr q-value ('pv', 'fdrq').
  - `-s` : A path to the output files of the previous subnetwork and pathway enrichment results.

Network graph output (Example file: "{sample}_{pvalue}_rwr{rate}_subnetwork_pathway_graph.html" in "for_graphs" folder)

The offline webpage shows the enriched pathway genes on subnetwork by network graph with cytoscape.js.

![image](https://github.com/user-attachments/assets/fb1a6694-0ce5-42ee-8a6a-e62f4db468fc)
