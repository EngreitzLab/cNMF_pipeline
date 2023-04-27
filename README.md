# Variant-to-Gene-to-Program (V2G2P) Approach

## Description

A pipeline to connect GWAS variants to genes to disease-associate gene programs. This pipeline uses snakemake and cNMF from Kotliar et al. The V2G2P approach could be applied to any GWAS studies with the right cell type(s). 

## Usage
### Step 1: Clone this github repository

### Step 2: Install conda environment
Install Snakemake and conda environment using conda:

```bash
conda create -f conda_env/cnmf_env.yml
conda create -f conda_env/cnmf_analysis_R.yml
```
cnmf_env contains snakemake. If you do not have snakemake installed already, you can activate the environment via `conda activate cnmf_env`, then run the pipeline.


### Step 3: Gather all input data
#### Necessary inputs:
Config file slots:
| field         | meaning                   |
|---------------|---------------------------|
| total_workers | Number of processes to run in parallel |
| seed | A number to set seed for reproducibility |
| num_runs | Number of NMF run (recommend 100 for the actual data analysis, and 10 for testing the pipeline) |
| run_per_worker | Number of run assign to each process/worker to run in series (make sure that num_runs = total_workers * run_per_worker)|
| k | number of components to factorize |
| thresholds | Threshold for filtering outlier components, 2 for no filtering, recommend 0.2 |
| analysisDir | Directory for outputting results |
| figDir | Directory for outputting figures |
| dataDir | If using 10X scRNA-seq matrix as input, the directory with 10X files (file names must be features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz). This can be left blank if supplying the input matrix as an h5ad file. |
| input_h5ad_mtxDir | Path to input matrix in h5ad format, can leave this blank if supplying the input matrix in other formats |
| scratchDir | Directory for cNMF pipeline to store temporary files |
| barcodeDir | (Necessary for Perturb-seq, optional for scRNA-seq without perturbation) Path to cell barcode / identifier tsv file, the sequence of appearance should be the same as in the input cell x gene matrix, and it must contain columns named "Gene" and "Guide" if doing Perturb-seq |
| motif_meme | MEME file for matching sequence to motifs |
| fasta_file | Fasta file to find sequences based on coordinates |
| promoter_fimo_formatted | (Optional) Path to promoter motif matches from FIMO output |
| enhancer_fimo_formatted | (Optional) Path to ABC enhancer motif matches from FIMO output | 
| ABC_enhancers | Coordinates for ABC enhancers for the approapriate cell type |
| PoPS_raw_featureDir (not used in this version) | Directory for storing feature txt files (one row per gene in ENSGID) for PoPS input |
| PoPS_gene_annotation_path (not used in this version) | Gene annotation file for PoPS pipeline |
| magma_dir (not used in this version) | Path to MAGMA results for PoPS pipeline input |
| magma_prefix (not used in this version) | Desired magma_prefix for pointing to files in magma_dir |
| PoPS_control_features (not used in this version) | A list of features that serves as controls for the PoPS pipeline |
| PoPS_features_metadata_path (not used in this version) | PoPS external feature's annotation |
| pipelineDir (not used in this version) | The directory to this snakemake pipeline |
| K_spectra_threshold_table (Optional) | Table used for specifying spectra cut-off for each K tested in the pipeline. default: 0.2 for all K |
| GWAS_traits | GWAS traits to conduct V2G2P |
| sampleName | Name for this run |
| cNMF_gene_selection | "all_genes" for using all expressed genes as input, "top3000VariableGenes" for using the top 3,000 most variable genes as an input, the number of top variable genes can be adjusted by the user, e.g. "top2000VariableGenes" for using the top 2,000 most variable genes |
| Perturb-seq | True for an Perturb-seq experiment. The pipeline will generate additional results based on comparison between perturbation and control. False for any scRNA-seq experiment. |
| num_Genes_per_MAST_runGroup | (Optional, required for Perturb-seq datasets) Number of perturbations to be tested against control in one job submission. This is to speed up the code. Recommend 20. |
| numCtrls_for_MAST | (Optional, required for Perturb-seq datasets) Number of (randomly selected) controls to use in statstical test. Recommend 5000 or less for speed. |
| num_cells | Expected number of cells in the input dataset. It doesn't have to be exact. This is for optimizing memory requests and job scheduling. |


## Outputs
The output files can be found in the folders specified in analysisDir and figDir fields in the config file.

### Analysis files (in analysisDir)
\*\_ProgramSummary\_k\_\*\_dt\_\*.xlsx
```
primerid	Pr(>Chisq)	coef	ci.hi	ci.lo	fdr	perturbation	zlm.model.name	fdr.across.ptb
topic_1		0.931958981499932	-0.00308004078521135	0.0570387818054911		-0.0631988633759138	0.997587143723284 A1BG batch.correction 0.980597200884276
topic_10	0.997587143723284	0.00216003169008006	0.0715983146968009		-0.0672782513166408	0.997587143723284 A1BG batch.correction	0.999373860636769
topic_11	0.996904491792907	-6.48785292078902e-06	0.0407059273499904		-0.040718903055832	0.997587143723284 A1BG batch.correction	0.99916248171282
topic_12	0.938744741036869	0.00810003101565848	0.0544634189392725		-0.0382633569079555	0.997587143723284 A1BG batch.correction	0.982593736237456  

Column description:
primerid    program name
Pr(>Chisq)  p-value
coef    natural log fold change
ci.hi   confidence interval upper bound
ci.lo   confidence interval lower bound
fdr	perturbation    FDR across K programs per perturbation
zlm.model.name  Model name (default is ~ProgramExpression + SampleIndex)
fdr.across.ptb  FDR across K programs x all perturbations

```

Perturb-seq outputs:
MAST_DEtopics.txt


### Figures (in figDir)

### Perturb-seq data outputs


### Running the pipeline:

