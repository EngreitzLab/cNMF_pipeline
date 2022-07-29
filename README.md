# cNMF pipeline


### Config file:
| field         | meaning                   |
|---------------|---------------------------|
| total_workers | number of processes to run in parallel |
| seed | a number to set seed for reproducibility |
| num_runs | number of NMF run (recommend 100 for analyzing the data, and 10 for testing the pipeline) |
| run_per_worker | number of run assign to each process/worker to run in series (num_runs = total_workers * run_per_worker)|
| k | number of topics to factorize |
| thresholds | spectra threshold for filtering outlier spectras, 2 means no filtering, recommend 0.2 |
| analysisDir | directory for outputting results |
| figDir | directory for outputting figures |
| dataDir | directory with 10X output as the input of the pipeline (features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz) |
| input_h5ad_mtxDir | directory for input matrix in h5ad format, can leave this blank if supplying the input matrix in other formats |
| scratchDir | directory for cNMF pipeline to store temporary files |
| barcodeDir | path to cell barcode / identifier tsv file, the sequence of appearance should be the same as in the input cell x gene matrix, and it must contain columns named "Gene" and "Guide" if doing Perturb-seq |
| motif_meme | MEME file for matching sequence to motifs |
| fasta_file | fasta file to find sequences based on coordinates |
| promoter_fimo_formatted | path to promoter motif matches |
| enhancer_fimo_formatted | path to ABC enhancer motif matches | 
| ABC_enhancers (not used in this version) | coordinates for ABC enhancers for the correct cell type |
| PoPS_raw_featureDir (not used in this version) | directory for storing feature txt files (one row per gene in ENSGID) for PoPS input |
| PoPS_gene_annotation_path (not used in this version) | gene annotation file for PoPS pipeline |
| magma_dir (not used in this version) | path to MAGMA results for PoPS pipeline input |
| magma_prefix (not used in this version) | desired magma_prefix for pointing to files in magma_dir |
| PoPS_control_features (not used in this version) | a list of features that serves as controls for the PoPS pipeline |
| PoPS_features_metadata_path (not used in this version) | PoPS external feature's annotation |
| pipelineDir (not used in this version) | the directory to this snakemake pipeline |
| K_spectra_threshold_table (Optional) | Table used for specifying spectra cut-off for each K tested in the pipeline. default: 0.2 for all K |
| sampleName | Name for this run |
| subsample_type (not used in this version) | Type of subsampling (by gene or by guide) |
| cNMF_gene_selection | "all_genes" for using all expressed genes as input, "top3000VariableGenes" for using the top 3,000 most variable genes as an input, the number of top variable genes can be adjusted by the user, e.g. "top2000VariableGenes" for using the top 2,000 most variable genes |
| Perturb-seq | True for an Perturb-seq experiment. The pipeline will generate additional results based on comparison between perturbation and control. False for any scRNA-seq experiment. |


### Conda Environment:
Download the .yml files from conda_envs/ and install via
```bash
conda create -f cnmf_env.yml
conda create -f cnmf_analysis_R.yml
```
cnmf_env contains snakemake. If you do not have snakemake installed already, you can activate the environment via `conda activate cnmf_env`, then run the pipeline.

### Running the pipeline:
Example script in log.sh
Recommend doing a dry-run to double check that the workflow is working properly.
