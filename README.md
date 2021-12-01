# cNMF pipeline


### Config file:
| field         | meaning                   |
|---------------|---------------------------|
| total_workers | number of processes to run in parallel |
| seed | a number to set seed for reproducibility |
| num_runs | number of NMF run (recommend 100 for analyzing the data, and 10 for testing the pipeline) |
| run_per_worker | number of run assign to each process/worker to run in series |
| k | number of topics to factorize |
| thresholds | spectra threshold for filtering outlier spectras, 2 is no filtering |
| analysisDir | directory for outputting results |
| figDir | directory for outputting figures |
| dataDir | directory with 10X output as the input of the pipeline (features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz |
| scratchDir | directory for cNMF pipeline to store temporary files |
| cNMF_gene_selection | "all_genes" for using all expressed genes as input, "top3000VariableGenes" for using the top 3,000 most variable genes as an input, the number of top variable genes can be adjusted by the user, e.g. "top2000VariableGenes" for using the top 2,000 most variable genes |
| motif_meme | MEME file for matching sequence to motifs |
| fasta_file | fasta file to find sequences based on coordinates |
| fimo_formatted | ABC enhancer motif matches (if available) | 
| ABC_enhancers | coordinates for ABC enhancers for the correct cell type |
| PoPS_raw_featureDir | directory for storing feature txt files (one row per gene in ENSGID) for PoPS input |
| PoPS_gene_annotation_path | gene annotation file for PoPS pipeline |
| magma_dir | path to MAGMA results for PoPS pipeline input |
| magma_prefix | desired magma_prefix for pointing to files in magma_dir |
| PoPS_control_features | a list of features that serves as controls for the PoPS pipeline |
| PoPS_features_metadata_path | PoPS external feature's annotation |
| pipelineDir | the directory to this snakemake pipeline |
| sampleName | the name of this run |


### Conda Environment:
Download the .yml files from conda_envs/ and install via
```bash
conda create -f cnmf_env.yml
conda create -f cnmf_analysis_R.yml
```
cnmf_env contains snakemake. If you do not have snakemake installed already, you can activate the environment via `conda activate cnmf_env`, then run the pipeline.

### Running the pipeline:

