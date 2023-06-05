### Run cNMF on Perturb-seq data and analyze the results

import subprocess
import numpy as np
import math
import pandas as pd
import os

threshold_txt = [ s.replace(".", "_") for s in config["thresholds"] ]
# num_workers = math.ceil(config["num_runs"] / config["run_per_worker"])
run_index_list = [n for n in range(config["run_per_worker"])]


def get_rule_prepare_cNMF_memory(wildcards):
	num_cells = config["num_cells"]
	if num_cells > 1e6:
		return "300"
	else:
		return "196"

def get_rule_prepare_cNMF_partition(wildcards):
	num_cells = config["num_cells"]
	if num_cells > 1e6:
		return "bigmem"
	else:
		return "owners,normal"


# def get_output_files(wildcards):
# 	output_files = []

# 	h5ad_mtx = expand(os.path.join(config["analysisDir"], "data/{sample}.h5ad"), sample=config["sampleName"]),
# 	gene_name_txt = expand(os.path.join(config["analysisDir"], "data/{sample}.h5ad.all.genes.txt"), sample=config["sampleName"]),
# 	# seurat_object_withUMAP = expand(os.path.join(config["analysisDir"], "data/{sample}.withUMAP_SeuratObject.RDS"), sample=config["sampleName"]),
# 	usages = expand(os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.usages.k_{k}.dt_{threshold}.consensus.txt"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]]),
# 	findK_plot_cNMF = expand(os.path.join(config["figDir"],"{folder}/{sample}/acrossK/{sample}.k_selection.png"), folder = config["cNMF_gene_selection"], sample=config["sampleName"]),
	
# 	aggregated_output = expand(os.path.join(config["analysisDir"], "{folder}/{sample}/acrossK/aggregated.outputs.findK.RData"), folder = config["cNMF_gene_selection"], sample=config["sampleName"]),
		
# 	cNMF_Results = expand(os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/cNMF_results.k_{k}.dt_{threshold}.RData"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]]),
# 	topic_zscore_top50 = expand(os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_top50GeneInTopics.zscore.pdf"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]]),
# 	topic_zscore_correlation = expand(os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_topic.Pearson.correlation.pdf"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]]),
# 	motif_plot = expand(os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_zscore.{ep_type}.motif.{test_type_and_threshold}.pdf"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]], ep_type = ["enhancer", "promoter"], test_type_and_threshold = ["count.ttest.enrichment_motif.thr.qval0.1"]),
# 	# motif_plot = expand(os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_zscore.{ep_type}.motif.enrichment{test_type_and_threshold}.pdf"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]], ep_type = ["enhancer", "promoter"], test_type_and_threshold = ["", "_motif.thr.10e-6", ".by.count.ttest", ".by.count.ttest_motif.thr.10e-6", ".by.count.ttest.labelPos", ".by.count.ttest_motif.thr.10e-6.labelPos"]),
# 	fgsea_result = expand(os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/fgsea/fgsea_all_pathways_df_{ranking_type}_k_{k}.dt_{threshold}.RData"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]], ranking_type=["raw.score", "z.score"]),
# 	# factor_expression_UMAP = expand(os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_Factor.Expression.UMAP.pdf"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]]),
# 	# GO_raw_plot = expand(os.path.join(config["figDir"], "{folder}/{sample}/acrossK/GO.enrichment.on.raw.score.ranking.threshold0.1.pdf"), folder = config["cNMF_gene_selection"], sample=config["sampleName"]),
	

# 	## Across K
# 	clustering_plot = expand(os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.clustering.k_{k}.dt_{threshold}.png"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]]),
# 	topic_clustering_plot = expand(os.path.join(config["figDir"], "{folder}/{sample}/acrossK/cluster.topic.zscore.by.Pearson.corr.pdf"), folder = config["cNMF_gene_selection"], sample=config["sampleName"])
 	
#  	# output_files.extend([h5ad_mtx, gene_name_txt, usages, findK_plot_cNMF, aggregated_output, cNMF_Results, topic_zscore_top50, topic_zscore_correlation, motif_plot, fgsea_result])
# 	# print(h5ad_mtx)
# 	# print(type(h5ad_mtx))
# 	# print([i for i in list(h5ad_mtx)])
# 	output_files.extend(list(h5ad_mtx))
# 	output_files.extend(list(gene_name_txt))

# 	output_files.extend(list(usages))
# 	output_files.extend(list(findK_plot_cNMF))
# 	output_files.extend(list(aggregated_output))
# 	output_files.extend(list(cNMF_Results))
# 	output_files.extend(list(topic_zscore_top50))
# 	output_files.extend(list(topic_zscore_correlation))
# 	output_files.extend(list(motif_plot))
# 	output_files.extend(list(fgsea_result))

#  	if config["Perturb-seq"] == "True":
#  		batch_correlation_pdf = expand(os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_batch.correlation.heatmap.pdf"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]]),
# 		percent_batch_topics_plot = expand(os.path.join(config["figDir"], "{folder}/{sample}/acrossK/percent.batch.topics.pdf"), folder = config["cNMF_gene_selection"], sample=config["sampleName"]),
# 		batch_topic_list = expand(os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/batch.topics.txt"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]]),

# 		## motif enrichment
# 		motif_plot_perturbseq = expand(os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_zscore.{ep_type}.motif.{test_type_and_threshold}.pdf"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]], ep_type = ["enhancer", "promoter"], test_type_and_threshold = ["count.ttest.enrichment_motif.thr.qval0.1"]),
		
# 		## Statistical Test
# 		# cNMF_Analysis = expand(os.path.join(config["analysisDir"], "{sample}/{folder}/K{k}/threshold_{threshold}/cNMFAnalysis.k_{k}.dt_{threshold}.RData"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]]),
# 		wilcoxon_txt_results = expand(os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/all.test.k_{k}.dt_{threshold}.minGuidePerPtb_{min_guide_per_ptb}.minCellPerGuide_{min_cell_per_guide}.txt"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]], min_guide_per_ptb = [1], min_cell_per_guide = [2]),
# 		MAST_txt_result = expand(os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/{sample}_MAST_DEtopics.txt"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]]),

# 		## aggregated across K results for Perturb-seq
# 		aggregated_output_perturbseq = expand(os.path.jn(config["analysisDir"], "{folder}/{sample}/acrossK/aggregated.outputs.findK.perturb-seq.RData"), folder = config["cNMF_gene_selection"], sample=config["sampleName"]),


# 		## PoPS
# 		feature_x_gene_importance_score_RDS = expand(os.path.join(config["scratchDir"], 
# 			"{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}_coefs.marginals.feature.outer.prod.RDS"), folder=config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]], magma_prefix=config["magma_prefix"]),
# 		feature_x_gene_component_importance_score_coefs_plot = expand(os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{magma_prefix}_cNMF{k}_dt_{threshold}_feature_x_gene.component.importance.score.coefs.pdf"), folder = config["cNMF_gene_selection"], sample=config["sampleName"], k=config["k"], threshold=[n.replace(".","_") for n in config["thresholds"]], magma_prefix=config["magma_prefix"])

# 		output_files.extend(list(batch_correlation_pdf))
# 		output_files.extend(list(percent_batch_topics_plot))
# 		output_files.extend(list(batch_topic_list))
# 		output_files.extend(list(motif_plot_perturbseq))
# 		output_files.extend(list(wilcoxon_txt_results))
# 		output_files.extend(list(MAST_txt_result))
# 		output_files.extend(list(aggregated_output_perturbseq))
# 		output_files.extend(list(feature_x_gene_importance_score_RDS))
# 		# output_files.extend(list(feature_x_gene_component_importance_score_coefs_plot))

# 	all_output = []
# 	for output_list in output_files:
# 		for item in output_list:
# 			all_output.extend([item])

# 	# print(all_output)
# 	# print(output_files)
# 	return all_output


# ## todo: rethink the structure of UMAP part: what are the possible input file types?
# rule create_Seurat_Object:
# 	input:
# 		mtx = os.path.join(config["dataDir"], "matrix.mtx.gz"),
# 		features = os.path.join(config["dataDir"], "features.tsv.gz"),
# 		barcodes = os.path.join(config["dataDir"], "barcodes.tsv.gz")
# 	output:
# 		seurat_object = os.path.join(config["analysisDir"], "data/{sample}.SeuratObject.RDS")
# 	params:
# 		time = "2:00:00",
# 		mem_gb = "200",
# 		datadir = config["dataDir"],
# 		outdir = os.path.join(config["analysisDir"], "data"),
# 		partition = "owners,normal"
# 	shell:
# 		"bash -c ' source $HOME/.bashrc; \
# 		conda activate cnmf_analysis_R; \
# 		Rscript workflow/scripts/create_seurat_object.R \
# 		--outdir {params.outdir}/ \
# 		--datadir {params.datadir} \
# 		--sampleName {wildcards.sample} \
# 		' "


## convert Seurat Object to h5ad file
rule Seurat_Object_to_h5ad:
	input:
		seurat_object = os.path.join(config["analysisDir"], "data/{sample}.SeuratObject.RDS")
	output:
		h5ad_mtx = os.path.join(config["analysisDir"], "data/{sample}.h5ad"),
		gene_name_txt = os.path.join(config["analysisDir"], "data/{sample}.h5ad.all.genes.txt")
	params:
		time = "2:00:00",
		mem_gb = "64",
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/seurat_to_h5ad.R \
		--inputSeuratObject {input.seurat_object} \
		--output_h5ad {output.h5ad_mtx} \
		--output_gene_name_txt {output.gene_name_txt} ' "


rule raw_h5ad_to_filtered_h5ad:
	input:
		raw_h5ad_file = os.path.join(config["analysisDir"], "data/{sample}.raw.h5ad")
	output:	
		h5ad_mtx = os.path.join(config["analysisDir"], "data/{sample}.h5ad"),
		gene_name_txt = os.path.join(config["analysisDir"], "data/{sample}.h5ad.all.genes.txt")
	params:
		time = "2:00:00",
		mem_gb = "64"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		python workflow/scripts/raw_h5ad_to_h5ad.py \
		--input_h5ad {input.seurat_object} \
		--output_h5ad {output.h5ad_mtx} \
		--output_gene_name_txt {output.gene_name_txt} ' "
	

# ## need a rule to convert txt to h5ad or RDS to h5ad ## outdated but keeping it here for reference
# rule txt_to_h5ad:
# 	input:
# 		mtx_txt_file = os.path.join(config["input_h5ad_mtxDir"], "{sample}.txt")
# 	output:
# 		h5ad_mtx = os.path.join(config["analysisDir"], "{sample}.h5ad"),
# 		gene_name_txt = os.path.join(config["analysisDir"], "{sample}.h5ad.all.genes.txt")
# 	shell:
# 		" bash -c ' source $HOME/.bashrc; \
# 		conda activate cnmf_env; \
# 		python workflow/scripts/txt_to_h5ad.py \
# 		--txtPath {input.mtx_txt_file} \
# 		--output_h5ad_mtx {output.h5ad_mtx} \
# 		--output_gene_name_txt {output.gene_name_txt} ' "


def get_num_cells(wildcards):
	return config["num_cells"]



# rule prepare_varGenes_cNMF:
# 	input:
# 		# h5ad = os.path.join(config["analysisDir"], "{folder}/{sample}/{sample}.h5ad")
# 		h5ad_mtx = os.path.join(config["analysisDir"], "data/{sample}.h5ad")
# 	output:
# 		tpm_h5ad = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
# 		tpm_stats = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
# 		nmf_yaml = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
# 		nmf_params = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
# 		norm_counts = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
# 		overdispersed_genes = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker{workerIndex}/{sample}/{sample}.overdispersed_genes.txt")
# 	params:
# 		time = "3:00:00",
# 		mem_gb = get_rule_prepare_cNMF_memory, #"64"
# 		seed = config["seed"],
# 		run_per_worker = config["run_per_worker"],
# 		total_workers = config["total_workers"],
# 		outdir = os.path.join(config["scratchDir"], "top{num_genes}VariableGenes/K{k}/worker{workerIndex}/"),
# 		partition = get_rule_prepare_cNMF_partition #"owners,normal"
# 	# resources: 
# 	# 	mem_mb=128*1000,
# 	# 	time = "3:00:00"
# 	# threads: config["total_workers"]
# 	shell:
# 		" bash -c ' source $HOME/.bashrc; \
# 		conda activate cnmf_env; \
# 		mkdir -p {params.outdir}/{wildcards.sample}; \
# 		python workflow/scripts/cNMF/cnmf.py prepare \
# 		--output-dir {params.outdir} \
# 		--name {wildcards.sample} \
# 		-c {input.h5ad_mtx} \
# 		-k {wildcards.k} \
# 		--n-iter {params.run_per_worker} \
# 		--total-workers {params.run_per_worker} \
# 		--seed {params.seed} \
# 		--numgenes {wildcards.num_genes} ' "


rule prepare_varGenes_cNMF_oneRun:
	input:
		# h5ad = os.path.join(config["analysisDir"], "{folder}/{sample}/{sample}.h5ad")
		h5ad_mtx = os.path.join(config["analysisDir"], "data/{sample}.h5ad")
	output:
		tpm_h5ad = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
		tpm_stats = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
		nmf_yaml = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
		nmf_params = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
		norm_counts = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
		overdispersed_genes = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/{sample}.overdispersed_genes.txt")
	params:
		time = "3:00:00",
		mem_gb = get_rule_prepare_cNMF_memory, #"64"
		seed = config["seed"],
		run_per_worker = config["run_per_worker"],
		total_workers = config["total_workers"],
		outdir = os.path.join(config["scratchDir"], "top{num_genes}VariableGenes/K{k}/worker0/"),
		partition = get_rule_prepare_cNMF_partition #"owners,normal"
	# resources: 
	# 	mem_mb=128*1000,
	# 	time = "3:00:00"
	# threads: config["total_workers"]
	shell:
		" bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		mkdir -p {params.outdir}/{wildcards.sample}; \
		python workflow/scripts/cNMF/cnmf.py prepare \
		--output-dir {params.outdir} \
		--name {wildcards.sample} \
		-c {input.h5ad_mtx} \
		-k {wildcards.k} \
		--n-iter {params.run_per_worker} \
		--total-workers {params.run_per_worker} \
		--seed {params.seed} \
		--numgenes {wildcards.num_genes} ' "



rule distribute_prepare_varGenes_cNMF:
	input:
		tpm_h5ad = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
		tpm_stats = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
		nmf_yaml = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
		nmf_params = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
		norm_counts = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
		overdispersed_genes = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/{sample}.overdispersed_genes.txt")
	output:
		tpm_h5ad = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
		tpm_stats = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
		nmf_yaml = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
		nmf_params = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
		norm_counts = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
		overdispersed_genes = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker{workerIndex}/{sample}/{sample}.overdispersed_genes.txt")
	params:
		time = "1:00:00",
		mem_gb = "12",
		seed = config["seed"],
		run_per_worker = config["run_per_worker"],
		total_workers = config["total_workers"],
		fromdir = os.path.join(config["scratchDir"], "top{num_genes}VariableGenes/K{k}/worker0"),
		todir = os.path.join(config["scratchDir"], "top{num_genes}VariableGenes/K{k}/worker{workerIndex}"),
		partition = "owners,normal"
	shell:
		" bash -c 'source $HOME/.bashrc; \
		conda activate cnmf_env; \
		mkdir -p {params.todir}/{wildcards.sample}; \
		cp -r {params.fromdir}/{wildcards.sample}/cnmf_tmp {params.todir}/{wildcards.sample}/; \
		cp {params.fromdir}/{wildcards.sample}/{wildcards.sample}.overdispersed_genes.txt {params.todir}/{wildcards.sample}/ ' "



# rule distribute_prepare_varGenes_cNMF_secondMethod:
# 	input:
# 		tpm_h5ad = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
# 		tpm_stats = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
# 		nmf_yaml = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
# 		nmf_params = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
# 		norm_counts = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
# 		overdispersed_genes = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes/K{k}/worker0/{sample}/{sample}.overdispersed_genes.txt")
# 	output:
# 		tpm_h5ad = expand(os.path.join(config["scratchDir"],"top{{num_genes}}VariableGenes/K{{k}}/worker{workerIndex}/{{sample}}/cnmf_tmp/{{sample}}.tpm.h5ad"), workerIndex = [i+1 for i in range(config["total_workers"]-1)]),
# 		tpm_stats = expand(os.path.join(config["scratchDir"],"top{{num_genes}}VariableGenes/K{{k}}/worker{workerIndex}/{{sample}}/cnmf_tmp/{{sample}}.tpm_stats.df.npz"), workerIndex = [i+1 for i in range(config["total_workers"]-1)]),
# 		nmf_yaml = expand(os.path.join(config["scratchDir"],"top{{num_genes}}VariableGenes/K{{k}}/worker{workerIndex}/{{sample}}/cnmf_tmp/{{sample}}.nmf_idvrun_params.yaml"), workerIndex = [i+1 for i in range(config["total_workers"]-1)]),
# 		nmf_params = expand(os.path.join(config["scratchDir"],"top{{num_genes}}VariableGenes/K{{k}}/worker{workerIndex}/{{sample}}/cnmf_tmp/{{sample}}.nmf_params.df.npz"), workerIndex = [i+1 for i in range(config["total_workers"]-1)]),
# 		norm_counts = expand(os.path.join(config["scratchDir"],"top{{num_genes}}VariableGenes/K{{k}}/worker{workerIndex}/{{sample}}/cnmf_tmp/{{sample}}.norm_counts.h5ad"), workerIndex = [i+1 for i in range(config["total_workers"]-1)]),
# 		overdispersed_genes = expand(os.path.join(config["scratchDir"],"top{{num_genes}}VariableGenes/K{{k}}/worker{workerIndex}/{{sample}}/{{sample}}.overdispersed_genes.txt"), workerIndex = [i+1 for i in range(config["total_workers"]-1)])
# 	params:
# 		time = "3:00:00",
# 		mem_gb = 64
# 		seed = config["seed"],
# 		run_per_worker = config["run_per_worker"],
# 		total_workers = config["total_workers"],
# 		outdir = os.path.join(config["scratchDir"], "top{num_genes}VariableGenes/K{k}/worker0/"),
# 		partition = get_rule_prepare_cNMF_partition #"owners,normal"


rule prepare_geneSet_cNMF:
	input:
		# h5ad_mtx = os.path.join(config["input_h5ad_mtxDir"], "{sample}.h5ad"),
		# genes = os.path.join(config["input_h5ad_mtxDir"],"{sample}.h5ad.{gene_selection_method}.genes.txt")
		h5ad_mtx = os.path.join(config["analysisDir"], "data/{sample}.h5ad"),
		genes = os.path.join(config["analysisDir"], "data/{sample}.h5ad.all.genes.txt")
	output:
		tpm_h5ad = os.path.join(config["scratchDir"],"{gene_selection_method}_genes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
		tpm_stats = os.path.join(config["scratchDir"],"{gene_selection_method}_genes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
		nmf_yaml = os.path.join(config["scratchDir"],"{gene_selection_method}_genes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
		nmf_params = os.path.join(config["scratchDir"],"{gene_selection_method}_genes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
		norm_counts = os.path.join(config["scratchDir"],"{gene_selection_method}_genes/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
		overdispersed_genes = os.path.join(config["scratchDir"],"{gene_selection_method}_genes/K{k}/worker{workerIndex}/{sample}/{sample}.overdispersed_genes.txt")
	params:
		time = "5:00:00",
		mem_gb = get_rule_prepare_cNMF_memory, #"196"
        seed = config["seed"],
		run_per_worker = config["run_per_worker"],
		total_workers = config["total_workers"],
		outdir = os.path.join(config["scratchDir"], "{gene_selection_method}_genes/K{k}/worker{workerIndex}/"),
		partition = get_rule_prepare_cNMF_partition #"owners,normal"
	# resources: 
	# 	mem_mb=128*1000,
	# 	time = "3:00:00"
	# threads: config["total_workers"]
	shell:
		" bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		mkdir -p {params.outdir}/{wildcards.sample}; \
		python workflow/scripts/cNMF/cnmf.py prepare \
		--output-dir {params.outdir} \
		--name {wildcards.sample} \
		-c {input.h5ad_mtx} \
		-k {wildcards.k} \
		--n-iter {params.run_per_worker} \
		--total-workers {params.total_workers} \
		--seed {params.seed} \
		--genes-file {input.genes} ' "


def get_cNMF_memory(wildcards): ## add condition to allot more memory for large matrices
	if config["num_cells"] > 1e6:
		return "64"
	else:
		if "VariableGenes" in wildcards.folder:
			return "16"
		elif int(wildcards.k) >= 100:
			return "64"
		else:
			return "32"

def get_cNMF_time(wildcards):
	if int(wildcards.k) > 100 and wildcards.folder == "all_genes":
		return "48:00:00" ## temporarily changed to 48 to accommodate for owners partition
	elif int(wildcards.k) > 25:
		return "48:00:00"
	else:
		return "24:00:00"

def get_cNMF_memory_slurm(wildcards):
	if "VariableGenes" in wildcards.folder:
		return 16000
	elif int(wildcards.k) > 100:	
		return 64000
	else:
		return 32000

def get_cNMF_partition_slurm(wildcards):
	if ("all_genes" in wildcards.folder) and (int(wildcards.k) > 100):
		return "normal,engreitz"
	else:
		return "owners,normal"


rule run_cNMF: # we don't know which worker will be working on which run, and snakemake gets confused by this
    # https://github.com/dylkot/cNMF/blob/master/Tutorials/analyze_pbmc_example_data.ipynb
    # try GNU parallel and manually submit jobs?
	input:
		tpm_h5ad = os.path.join(config["scratchDir"],"{folder}/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
		tpm_stats = os.path.join(config["scratchDir"],"{folder}/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
		nmf_yaml = os.path.join(config["scratchDir"],"{folder}/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
		nmf_params = os.path.join(config["scratchDir"],"{folder}/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
		norm_counts = os.path.join(config["scratchDir"],"{folder}/K{k}/worker{workerIndex}/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
		overdispersed_genes = os.path.join(config["scratchDir"],"{folder}/K{k}/worker{workerIndex}/{sample}/{sample}.overdispersed_genes.txt")
		# dummy_file = os.path.join(config["scratchDir"], "K{k}/worker{workerIndex}/{sample}/tmp.worker.{workerIndex}.txt")
	output:
		individual_runs = expand(os.path.join(config["scratchDir"],"{{folder}}/K{{k}}/worker{{workerIndex}}/{{sample}}/cnmf_tmp/{{sample}}.spectra.k_{{k}}.iter_{run_index}.df.npz"), run_index = run_index_list)
	# dummy_file = os.path.join(config["scratchDir"], "K{k}/{sample}/worker{workerIndex}/tmp.worker.{workerIndex}.started.txt")
	# individual_runs = expand(os.path.join(config["scratchDir"],"K{{k}}/{{sample}}/cnmf_tmp/no_IL1B.spectra.k_{{k}}.iter_{run_index}.df.npz"), run_index=range(config["num_runs"]))
	params:
		time = get_cNMF_time, ## explore how to make snakemake automatically retry with more time if the time allocated is not enough
		mem_gb = get_cNMF_memory,
		total_workers = config["total_workers"],
		num_runs = config["num_runs"],
		outdir = os.path.join(config["scratchDir"],"{folder}/K{k}/worker{workerIndex}/"),
		partition = "owners,normal"
	# resources: 
	# 	mem_mb=get_cNMF_memory_slurm,
	# 	time = get_cNMF_time,
	# 	partition = get_cNMF_partition_slurm
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		python workflow/scripts/cNMF/cnmf.py factorize \
		--output-dir {params.outdir} \
		--name {wildcards.sample} ' "


rule pool_together_results:
	input:
		individual_runs = expand(os.path.join(config["scratchDir"],"{{folder}}/K{{k}}/worker{workerIndex}/{{sample}}/cnmf_tmp/{{sample}}.spectra.k_{{k}}.iter_{run_index}.df.npz"), workerIndex = [n for n in range(config["total_workers"])], run_index = [n for n in range(config["run_per_worker"])])
	output:
		pooled_individual_runs = expand(os.path.join(config["scratchDir"],"{{folder}}_combined/K{{k}}/{{sample}}/cnmf_tmp/{{sample}}.spectra.k_{{k}}.iter_{pooled_run_index}.df.npz"), pooled_run_index = [n for n in range(config["num_runs"])]) # should use all config["num_runs"], but DAG takes a while to build if all are used
	params:
		time = "3:00:00",
		mem_gb = "64",
		inputdir = os.path.join(config["scratchDir"],"{folder}/K{k}"),
		outdir = os.path.join(config["scratchDir"],"{folder}_combined/K{k}"),
		num_workers = config["total_workers"],
		run_per_worker = config["run_per_worker"],
		partition = "owners,normal"
	# resources: 
	# 	mem_mb=64*1000,
	# 	time = "3:00:00"
	run:
		cmd = "mkdir -p " + os.path.join(params.inputdir, wildcards.sample, "cnmf_tmp")
		shell(cmd)
		for worker in range(params.num_workers):
			for run in range(params.run_per_worker):
				from_file = os.path.join(params.inputdir, "worker" + str(worker), wildcards.sample, "cnmf_tmp/" + wildcards.sample + ".spectra.k_" + wildcards.k + ".iter_" + str(run) + ".df.npz")
				index_here = worker * params.run_per_worker + run  ### give the runs a new index
				to_file = os.path.join(params.outdir, wildcards.sample, "cnmf_tmp/" + wildcards.sample + ".spectra.k_" + wildcards.k + ".iter_" + str(index_here) + ".df.npz")
				cmd = "cp " + from_file + " " + to_file
				shell(cmd)


rule prepare_combine_varGene: # this step causes snakemake hanging at "Building DAG of jobs..."
	input:
		# h5ad_mtx = os.path.join(config["input_h5ad_mtxDir"], "{sample}.h5ad"),
		h5ad_mtx = os.path.join(config["analysisDir"], "data/{sample}.h5ad"),
		pooled_individual_runs = expand(os.path.join(config["scratchDir"],"top{{num_genes}}VariableGenes_combined/K{{k}}/{{sample}}/cnmf_tmp/{{sample}}.spectra.k_{{k}}.iter_{pooled_run_index}.df.npz"), pooled_run_index = [n for n in range(config["num_runs"])])
	output:
		tpm_h5ad = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes_combined/K{k}/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
		tpm_stats = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes_combined/K{k}/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
		nmf_yaml = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes_combined/K{k}/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
		nmf_params = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes_combined/K{k}/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
		norm_counts = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes_combined/K{k}/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
		overdispersed_genes = os.path.join(config["scratchDir"],"top{num_genes}VariableGenes_combined/K{k}/{sample}/{sample}.overdispersed_genes.txt")
	params:
		time = "3:00:00",
		mem_gb = get_rule_prepare_cNMF_memory, #"64",
		seed = config["seed"],
		num_runs = config["num_runs"],
		outdir = os.path.join(config["scratchDir"], "top{num_genes}VariableGenes_combined/K{k}"),
		partition = get_rule_prepare_cNMF_partition
	# resources: 
	# 	mem_mb=64*1000,
	# 	time = "3:00:00"
	shell:
		" bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		mkdir -p {params.outdir}/K{wildcards.k}/; \
		python workflow/scripts/cNMF/cnmf.py prepare \
		--output-dir {params.outdir} \
		--name {wildcards.sample} \
		-c {input.h5ad_mtx} \
		-k {wildcards.k} \
		--n-iter {params.num_runs} \
		--total-workers 1 \
		--seed {params.seed} \
		--numgenes {wildcards.num_genes} ' "


rule prepare_combine_geneSet: # this step causes snakemake hanging at "Building DAG of jobs..."
	input:
		# h5ad_mtx = os.path.join(config["input_h5ad_mtxDir"], "{sample}.h5ad"),
		pooled_individual_runs = expand(os.path.join(config["scratchDir"],"{{gene_selection_method}}_genes_combined/K{{k}}/{{sample}}/cnmf_tmp/{{sample}}.spectra.k_{{k}}.iter_{pooled_run_index}.df.npz"), pooled_run_index = [n for n in range(config["num_runs"])]),
		# genes = os.path.join(config["dataDir"],"{sample}.h5ad.{gene_selection_method}.genes.txt")
		h5ad_mtx = os.path.join(config["analysisDir"], "data/{sample}.h5ad"),
		genes = os.path.join(config["analysisDir"], "data/{sample}.h5ad.all.genes.txt")
	output:
		tpm_h5ad = os.path.join(config["scratchDir"],"{gene_selection_method}_genes_combined/K{k}/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
		tpm_stats = os.path.join(config["scratchDir"],"{gene_selection_method}_genes_combined/K{k}/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
		nmf_yaml = os.path.join(config["scratchDir"],"{gene_selection_method}_genes_combined/K{k}/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
		nmf_params = os.path.join(config["scratchDir"],"{gene_selection_method}_genes_combined/K{k}/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
		norm_counts = os.path.join(config["scratchDir"],"{gene_selection_method}_genes_combined/K{k}/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
		overdispersed_genes = os.path.join(config["scratchDir"],"{gene_selection_method}_genes_combined/K{k}/{sample}/{sample}.overdispersed_genes.txt")
	params:
		time = "3:00:00",
		mem_gb = get_rule_prepare_cNMF_memory, #"160",
		seed = config["seed"],
		num_runs = config["num_runs"],
		outdir = os.path.join(config["scratchDir"], "{gene_selection_method}_genes_combined/K{k}"),
		partition = get_rule_prepare_cNMF_partition #"owners,normal"
	# resources: 
	# 	mem_mb=160*1000,
	# 	time = "3:00:00"
	shell:
		" bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		mkdir -p {params.outdir}/K{wildcards.k}/; \
		python workflow/scripts/cNMF/cnmf.py prepare \
		--output-dir {params.outdir} \
		--name {wildcards.sample} \
		-c {input.h5ad_mtx} \
		-k {wildcards.k} \
		--n-iter {params.num_runs} \
		--total-workers 1 \
		--seed {params.seed} \
		--genes-file {input.genes} ' "


rule combine_runs_varGenes:
	input:
		tpm_h5ad = os.path.join(config["scratchDir"],"{folder}_combined/K{k}/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
		tpm_stats = os.path.join(config["scratchDir"],"{folder}_combined/K{k}/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
		nmf_yaml = os.path.join(config["scratchDir"],"{folder}_combined/K{k}/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
		nmf_params = os.path.join(config["scratchDir"],"{folder}_combined/K{k}/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
		norm_counts = os.path.join(config["scratchDir"],"{folder}_combined/K{k}/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
		overdispersed_genes = os.path.join(config["scratchDir"],"{folder}_combined/K{k}/{sample}/{sample}.overdispersed_genes.txt"),
		pooled_individual_runs = expand(os.path.join(config["scratchDir"],"{{folder}}_combined/K{{k}}/{{sample}}/cnmf_tmp/{{sample}}.spectra.k_{{k}}.iter_{pooled_run_index}.df.npz"), pooled_run_index = [n for n in range(config["num_runs"])])
	# individual_runs = dynamic(os.path.join(config["scratchDir"],"K{{k}}/{{sample}}/cnmf_tmp/{{sample}}.spectra.k_{{k}}.iter_{pooled_run_index}.df.npz"))
	# # individual_runs = expand(os.path.join(config["scratchDir"],"K{{k}}/{{sample}}/cnmf_tmp/no_IL1B.spectra.k_{{k}}.iter_{run_index}.df.npz"), run_index=range(config["num_runs"]))
	output:
		merged_result = os.path.join(config["scratchDir"],"{folder}_combined/K{k}/{sample}/cnmf_tmp/{sample}.spectra.k_{k}.merged.df.npz")
	params:
		time = "3:00:00",
		mem_gb = "64",
		outdir = os.path.join(config["scratchDir"],"{folder}_combined/K{k}/"),
		partition = "owners,normal"
	# resources: 
	# 	mem_mb=64*1000,
	# 	time = "3:00:00"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		python workflow/scripts/cNMF/cnmf.py combine \
		--output-dir {params.outdir} \
		--name {wildcards.sample} ' "


rule aggregate_combined_runs:
	input:
		merged_result = os.path.join(config["scratchDir"],"{folder}_combined/K{k}/{sample}/cnmf_tmp/{sample}.spectra.k_{k}.merged.df.npz")
	output:
		merged_copied_result = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.spectra.k_{k}.merged.df.npz")
	# resources: time_min = 60
	params:
		time = "3:00:00",
		mem_gb = "64",
		outdir = config["scratchDir"],
		partition = "owners,normal"
	# resources: 
	# 	mem_mb=64*1000,
	# 	time = "3:00:00"
	run:
		cmd = "mkdir -p " + os.path.join(params.outdir, wildcards.folder + "_acrossK", wildcards.sample, "cnmf_tmp")
		shell(cmd)
		cmd = "cp " + input.merged_result + " " + output.merged_copied_result
		shell(cmd)



rule prepare_findK_varGene:
	input:
		# h5ad_mtx = os.path.join(config["input_h5ad_mtxDir"], "{sample}.h5ad")
		h5ad_mtx = os.path.join(config["analysisDir"], "data/{sample}.h5ad")
		# merged_copied_result = expand(os.path.join(config["scratchDir"],"top{{num_genes}}VariableGenes_acrossK/{{sample}}/cnmf_tmp/{{sample}}.spectra.k_{k}.merged.df.npz"), k=config["k"])
	output:
		tpm_h5ad = os.path.join(config["analysisDir"],"top{num_genes}VariableGenes_acrossK/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
		tpm_stats = os.path.join(config["analysisDir"],"top{num_genes}VariableGenes_acrossK/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
		nmf_yaml = os.path.join(config["analysisDir"],"top{num_genes}VariableGenes_acrossK/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
		nmf_params = os.path.join(config["analysisDir"],"top{num_genes}VariableGenes_acrossK/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
		norm_counts = os.path.join(config["analysisDir"],"top{num_genes}VariableGenes_acrossK/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
		overdispersed_genes = os.path.join(config["analysisDir"],"top{num_genes}VariableGenes_acrossK/{sample}/{sample}.overdispersed_genes.txt")
	params:
		time = "3:00:00",
		mem_gb = get_rule_prepare_cNMF_memory,
		seed = config["seed"],
		num_runs = config["num_runs"],
		outdir = os.path.join(config["analysisDir"], "top{num_genes}VariableGenes_acrossK"),
		klist = " ".join(str(k) for k in config["k"]),
		partition = get_rule_prepare_cNMF_partition
	threads: config["total_workers"]
	# resources: 
	# 	mem_mb=64*1000,
	# 	time = "3:00:00"
	shell:
		" bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		mkdir -p {params.outdir}; \
		python workflow/scripts/cNMF/cnmf.py prepare \
		--output-dir {params.outdir} \
		--name {wildcards.sample} \
		-c {input.h5ad_mtx} \
		-k {params.klist} \
		--n-iter {params.num_runs} \
		--total-workers 1 \
		--seed {params.seed} \
		--numgenes {wildcards.num_genes} ' "


rule prepare_findK_geneSet:
	input:
		# h5ad_mtx = os.path.join(config["input_h5ad_mtxDir"], "{sample}.h5ad"),
		# genes = os.path.join(config["dataDir"],"{sample}.h5ad.{gene_selection_method}.genes.txt")
		h5ad_mtx = os.path.join(config["analysisDir"], "data/{sample}.h5ad"),
		genes = os.path.join(config["analysisDir"], "data/{sample}.h5ad.all.genes.txt")
		# merged_copied_result = expand(os.path.join(config["scratchDir"],"top{{num_genes}}VariableGenes_acrossK/{{sample}}/cnmf_tmp/{{sample}}.spectra.k_{k}.merged.df.npz"), k=config["k"])
	output:
		tpm_h5ad = os.path.join(config["analysisDir"],"{gene_selection_method}_genes_acrossK/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
		tpm_stats = os.path.join(config["analysisDir"],"{gene_selection_method}_genes_acrossK/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
		nmf_yaml = os.path.join(config["analysisDir"],"{gene_selection_method}_genes_acrossK/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
		nmf_params = os.path.join(config["analysisDir"],"{gene_selection_method}_genes_acrossK/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
		norm_counts = os.path.join(config["analysisDir"],"{gene_selection_method}_genes_acrossK/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
		overdispersed_genes = os.path.join(config["analysisDir"],"{gene_selection_method}_genes_acrossK/{sample}/{sample}.overdispersed_genes.txt")
	params:
		time = "3:00:00",
		mem_gb = get_rule_prepare_cNMF_memory, #"200",
		seed = config["seed"],
		num_runs = config["num_runs"],
		outdir = os.path.join(config["analysisDir"], "{gene_selection_method}_genes_acrossK"),
		klist = " ".join(str(k) for k in config["k"]),
		partition = get_rule_prepare_cNMF_partition
	threads: config["total_workers"]
	# resources: 
	# 	mem_mb=128*1000,
	# 	time = "3:00:00"
	shell:
		" bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		mkdir -p {params.outdir}; \
		python workflow/scripts/cNMF/cnmf.py prepare \
		--output-dir {params.outdir} \
		--name {wildcards.sample} \
		-c {input.h5ad_mtx} \
		-k {params.klist} \
		--n-iter {params.num_runs} \
		--total-workers 1 \
		--seed {params.seed} \
		--genes-file {input.genes} ' "


def get_findK_cNMF_partition(wildcards):
	if "all_genes" in wildcards.folder or config["num_cells"] > 1e6:
		return "bigmem,owners" ## 2kG.library 800G, 12hrs
	else:
		return 'owners,normal'

def get_findK_cNMF_mem_gb(wildcards):
	if config["num_cells"] > 1e6:
		return "800" 
	elif "all_genes" in wildcards.folder:
		return "500" ## 256 for 210,000 cells
	else:
		return "96"

def get_findK_cNMF_mem_gb_slurm(wildcards):
	if "all_genes" in wildcards.folder:
		return 256000
	else:
		return 96000


rule findK_cNMF:
	input:
		tpm_h5ad = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
		tpm_stats = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
		nmf_yaml = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
		nmf_params = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
		norm_counts = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
		overdispersed_genes = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.overdispersed_genes.txt"),
		merged_copied_result = expand(os.path.join(config["analysisDir"],"{{folder}}_acrossK/{{sample}}/cnmf_tmp/{{sample}}.spectra.k_{k}.merged.df.npz"), k=[str(k) for k in config["k"]])
	output:
		plot = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.k_selection.png"),
		plot_new_location = os.path.join(config["figDir"],"{folder}/{sample}/acrossK/{sample}.k_selection.png")
	# resources: mem_mb=32000
	params:
		time = "24:00:00",
		mem_gb = get_findK_cNMF_mem_gb,
		outdir = os.path.join(config["analysisDir"], "{folder}_acrossK"),
		partition = get_findK_cNMF_partition
	# resources: 
	# 	mem_mb=get_findK_cNMF_mem_gb_slurm,
	# 	time = "24:00:00"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		python workflow/scripts/cNMF/cnmf.modified.py k_selection_plot --output-dir {params.outdir} --name {wildcards.sample}; \
		cp {output.plot} {output.plot_new_location} ' "


def get_concensus_factors_time(wildcards):
	if config["num_cells"] > 1e6:
		return "24:00:00"
	else:
		if int(wildcards.k) >= 70:
			return "36:00:00"
		elif (int(wildcards.k) >= 40 and wildcards.folder == "all_genes") or config["num_cells"] > 1e6:
			return "9:00:00"
		else:
			return "3:00:00"

def get_concensus_factors_time_slurm(wildcards):
	if int(wildcards.k) >= 100:
		return 36*3600
	elif int(wildcards.k) >= 40 and wildcards.folder == "all_genes":
		return 9*3600
	else:
		return 3*3600

def get_concensus_factors_memory(wildcards):
	if config["num_cells"] > 1e6:
		if int(wildcards.k) >= 150:
			return "1200"
		elif int(wildcards.k) >= 50: 
			return "800"
		else:
			return "500"
	else:
		return "255"

def get_concensus_factors_partition(wildcards):
	if config["num_cells"] > 1e6:
		if int(wildcards.k) >= 50:
			return "bigmem"	
		else:
			return "normal,owners"
	else:
		return "normal,owners"


rule get_concensus_factors:
	input:
		tpm_h5ad = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
		tpm_stats = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
		nmf_yaml = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
		nmf_params = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
		norm_counts = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
		overdispersed_genes = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.overdispersed_genes.txt"),
		individual_k = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.spectra.k_{k}.merged.df.npz")
		# dummy_file = os.path.join(config["scratchDir"],"{sample}/cnmf_tmp/{sample}.dummy.k_{k}.dt_{threshold}.txt")
	output:
		tmp_spectra_tpm = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.gene_spectra_tpm.k_{k}.dt_{threshold}.df.npz"),
		tmp_spectra_score = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.gene_spectra_score.k_{k}.dt_{threshold}.df.npz"),
		tmp_spectra_consensus = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.spectra.k_{k}.dt_{threshold}.consensus.df.npz"),
		tmp_stats = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.stats.k_{k}.dt_{threshold}.df.npz"),
		tmp_usages = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.usages.k_{k}.dt_{threshold}.consensus.df.npz"),
		# clustering_plot = os.path.join(config["scratchDir"],"{folder}_acrossK/{sample}/{sample}.clustering.k_{k}.dt_{threshold}.png"),
		spectra_tpm = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.gene_spectra_tpm.k_{k}.dt_{threshold}.txt"),
		spectra_score = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.gene_spectra_score.k_{k}.dt_{threshold}.txt"),
		spectra_consensus = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.spectra.k_{k}.dt_{threshold}.consensus.txt"),
		usages = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.usages.k_{k}.dt_{threshold}.consensus.txt")
	# resources: mem_mb=96000
	params:
		time = get_concensus_factors_time,
		mem_gb = get_concensus_factors_memory, # 96 for top 3000 genes
		outdir = os.path.join(config["analysisDir"], "{folder}_acrossK"),
		partition = get_concensus_factors_partition
	# resources: 
	# 	mem_mb=196000,
	# 	time = get_concensus_factors_time_slurm
	run:
		threshold_here = wildcards.threshold.replace("_",".")
		shell("bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		python workflow/scripts/cNMF/cnmf.modified.py consensus \
		--output-dir {params.outdir} \
		--name {wildcards.sample} \
		--components {wildcards.k} \
		--local-density-threshold {threshold_here} ' ") # --show-clustering 


rule get_concensus_factors_plot:
	input:
		tpm_h5ad = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
		tpm_stats = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.tpm_stats.df.npz"),
		nmf_yaml = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.nmf_idvrun_params.yaml"),
		nmf_params = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.nmf_params.df.npz"),
		norm_counts = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
		overdispersed_genes = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.overdispersed_genes.txt"),
		individual_k = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/cnmf_tmp/{sample}.spectra.k_{k}.merged.df.npz")
		# dummy_file = os.path.join(config["scratchDir"],"{sample}/cnmf_tmp/{sample}.dummy.k_{k}.dt_{threshold}.txt")
	output:
		clustering_plot = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.clustering.k_{k}.dt_{threshold}.png")
	# resources: 
	# 	mem_mb=128000,
	# 	time = get_concensus_factors_time_slurm
	params:
		time = get_concensus_factors_time,
		mem_gb = "256", # 96 for top 3000 genes
		outdir = os.path.join(config["analysisDir"], "{folder}_acrossK"),
		figdir = os.path.join(config["figDir"]),
		partition = "owners,normal"
	run:
		threshold_here = wildcards.threshold.replace("_",".")
		shell("bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		python workflow/scripts/cNMF/cnmf.modified.py consensus \
		--output-dir {params.outdir} \
		--name {wildcards.sample} \
		--components {wildcards.k} \
		--local-density-threshold {threshold_here} \
		--show-clustering; \
		cp {params.outdir}/{wildcards.sample}/{wildcards.sample}.clustering.k_{wildcards.k}.dt_{wildcards.threshold}.png {params.figdir}/{wildcards.folder}/{wildcards.sample}/K{wildcards.k}/ ' ")


def get_cNMF_filter_threshold_double(wildcards):
        return wildcards.threshold.replace("_",".")


# def get_topicModelAnalysis_time_slurm(wildcards):
# 	if int(wildcards.k) > 60:
# 		return "12:00:00" # 12* 60*60
# 	else:
# 		return "6:00:00" # 6 * 60*60

# def get_topicModelAnalysis_memory_slurm(wildcards):
# 	if int(wildcards.k) > 60:
# 		return 200*1000
# 	else:
# 		return 128*1000

# def get_topicModelAnlaysis_partition_slurm(wildcards):
# 	if int(wildcards.k) > 60:
# 		return "bigmem"
# 	else:
# 		return "owners,normal"



## same here: do we need a separate step to create the Seurat Object?
rule calc_UMAP:
	input:
		input_seurat_object = os.path.join(config["analysisDir"], "data/{sample}.SeuratObject.RDS")
	output:
		seurat_object_withUMAP = os.path.join(config["analysisDir"], "data/{sample}.withUMAP_SeuratObject.RDS")
	params:
		time = "6:00:00",
		mem_gb = "200",
		analysisdir = config["analysisDir"],
		datadir = config["dataDir"],
		outdir = os.path.join(config["analysisDir"], "data"),
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/calcUMAP.only.R \
			--outdir {params.outdir}/ \
			--inputSeuratObject {input.input_seurat_object} \
			--sampleName {wildcards.sample} \
			--maxMt 50 \
			--maxCount 25000 \
			--minUniqueGenes 0 \
			--UMAP.resolution 0.6' " ## can make this a wildcard

rule plot_UMAP:
	input:
		seurat_object_withUMAP = os.path.join(config["analysisDir"], "data/{sample}.withUMAP_SeuratObject.RDS"),
		cNMF_Results = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/cNMF_results.k_{k}.dt_{threshold}.RData")
	output:
		factor_expression_UMAP = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_Factor.Expression.UMAP.pdf")
	params:
		time = "6:00:00",
		mem_gb = "200",
		datadir = config["dataDir"],
		outdir = os.path.join(config["analysisDir"], "{folder}_acrossK/{sample}"),
		figdir = os.path.join(config["figDir"], "{folder}"), 
		analysisdir = os.path.join(config["analysisDir"], "{folder}"), # K{k}/threshold_{threshold}
		threshold = get_cNMF_filter_threshold_double,
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/cNMF_UMAP_plot.R \
		--sampleName {wildcards.sample} \
		--inputSeuratObject {input.seurat_object_withUMAP} \
		--figdir {params.figdir}/ \
		--outdir {params.analysisdir} \
		--K.val {wildcards.k} \
		--density.thr {params.threshold} \
		--recompute F ' "


# rule FIMO: # refer to /oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2104_remove_lincRNA/
# 	input:
# 	output:
# 	params:
# 	shell:


rule variance_explained:
	input:
		tpm_counts = os.path.join(config["analysisDir"], "{folder}_acrossK/{sample}/cnmf_tmp/{sample}.tpm.h5ad"),
		spectra_tpm = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.gene_spectra_tpm.k_{k}.dt_{threshold}.txt"),
		spectra_zscore = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.gene_spectra_score.k_{k}.dt_{threshold}.txt"),
		spectra_consensus = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.spectra.k_{k}.dt_{threshold}.consensus.txt")
	output:
		Var_k_txt = os.path.join(config["analysisDir"],"{folder}/{sample}/K{k}/threshold_{threshold}/metrics.varianceExplained.df.txt"),
		Var_k_summary_txt = os.path.join(config["analysisDir"],"{folder}/{sample}/K{k}/threshold_{threshold}/summary.varianceExplained.df.txt")
	params:
		time = "6:00:00",
		mem_gb = "160",
		path_to_topics = os.path.join(config["analysisDir"], "{folder}_acrossK/"),
		# tpm_counts_path = os.path.join(config["analysisDir"], "{folder}_acrossK/{sample}/cnmf_tmp/"),
		X_normalized_path = os.path.join(config["analysisDir"], "{folder}_acrossK/{sample}/cnmf_tmp/{sample}.norm_counts.h5ad"),
		outdir = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/"),
		threshold = get_cNMF_filter_threshold_double,
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		python workflow/scripts/variance_explained_v2.py \
		--path_to_topics {params.path_to_topics} \
		--topic_sampleName {wildcards.sample} \
		--X_normalized {params.X_normalized_path} \
		--outdir {params.outdir} \
		--k {wildcards.k} \
		--density_threshold {params.threshold} ' "
# --tpm_counts_path {input.tpm_counts} \


rule analysis:
	input:
		# mtx_RDS = config["inputRDSmtx"], # os.path.join(config["dataDir"],"aggregated.{sample}.mtx.cell_x_gene.RDS"),
		#todo: add input files
		spectra_tpm = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.gene_spectra_tpm.k_{k}.dt_{threshold}.txt"),
		spectra_zscore = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.gene_spectra_score.k_{k}.dt_{threshold}.txt"),
		spectra_consensus = os.path.join(config["analysisDir"],"{folder}_acrossK/{sample}/{sample}.spectra.k_{k}.dt_{threshold}.consensus.txt")
	output:
		#todo: add output files
		cNMF_Results = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/cNMF_results.k_{k}.dt_{threshold}.RData"),
		# cNMF_Analysis = os.path.join(config["analysisDir"], "{sample}/{folder}/K{k}/threshold_{threshold}/cNMFAnalysis.k_{k}.dt_{threshold}.RData")
		cNMF_ENSG_topic_zscore = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/topic.zscore.ensembl_k_{k}.dt_{threshold}.txt"),
		cNMF_ENSG_topic_zscore_scaled = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/topic.zscore.ensembl.scaled_k_{k}.dt_{threshold}.txt"),
		cNMF_ENSG_topic_raw = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/topic.tpm.ensembl_k_{k}.dt_{threshold}.txt"),
		cNMF_ENSG_topic_raw_scaled = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/topic.tpm.ensembl.scaled_k_{k}.dt_{threshold}.txt"),
		cNMF_median_spectra_zscore = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/median.spectra.zscore.df_k_{k}.dt_{threshold}.txt"),
		cNMF_median_spectra_zscore_ensembl_scaled = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/median.spectra.zscore.ensembl.scaled_k_{k}.dt_{threshold}.txt")
	params:
		time = "3:00:00",
		mem_gb = "64",
		outdir = os.path.join(config["analysisDir"], "{folder}_acrossK/{sample}"),
		figdir = os.path.join(config["figDir"], "{folder}"), 
		analysisdir = os.path.join(config["analysisDir"], "{folder}"), # K{k}/threshold_{threshold}
		# barcode = os.path.join(config["barcodeDir"], "{sample}.barcodes.tsv"),
		threshold = get_cNMF_filter_threshold_double,
		barcode_names = config["barcodeDir"],
		partition = "owners,normal"
		# subsample_type = config["subsample_type"]
	# resources:
	# 	mem_mb=get_topicModelAnalysis_memory_slurm,
	# 	time = get_topicModelAnalysis_time_slurm, ## 6 hours
	# 	partition = get_topicModelAnlaysis_partition_slurm  
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		conda env info; \
		Rscript workflow/scripts/cNMF_analysis.R \
		--topic.model.result.dir {params.outdir}/ \
		--sampleName {wildcards.sample} \
		--barcode.names {params.barcode_names} \
		--figdir {params.figdir}/ \
		--outdir {params.analysisdir}/ \
		--K.val {wildcards.k} \
		--density.thr {params.threshold} \
		--recompute F \
		--motif.enhancer.background /oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/data/fimo_out_ABC_TeloHAEC_Ctrl_thresh1.0E-4/fimo.formatted.tsv \
		--motif.promoter.background /oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2104_remove_lincRNA/data/fimo_out_all_promoters_thresh1.0E-4/fimo.tsv \
		' "


rule topic_plot: 
	input:
		cNMF_Results = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/cNMF_results.k_{k}.dt_{threshold}.RData")
	output:
		topic_zscore_top50 = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_top50GeneInTopics.zscore.pdf"),
		topic_raw_score_top50 = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_top50GeneInTopics.rawWeight.pdf"),
		topic_raw_score_top10 = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_top10GeneInTopics.rawWeight.pdf"),
		topic_zcore_top10 = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_top10GeneInTopics.zscore.pdf"),
		topic_zscore_correlation = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_topic.Pearson.correlation.pdf")
	params:
		time = "3:00:00",
		mem_gb = "64",
		outdir = os.path.join(config["analysisDir"], "{folder}_acrossK/{sample}"),
		figdir = os.path.join(config["figDir"], "{folder}"), 
		analysisdir = os.path.join(config["analysisDir"], "{folder}"), # K{k}/threshold_{threshold}
		threshold = get_cNMF_filter_threshold_double,
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/cNMF_analysis_topic_plot.R \
		--sampleName {wildcards.sample} \
		--figdir {params.figdir}/ \
		--outdir {params.analysisdir}/ \
		--K.val {wildcards.k} \
		--density.thr {params.threshold} \
		--recompute F ' "


rule batch_topic_correlation:
	input:
		cNMF_Results = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/cNMF_results.k_{k}.dt_{threshold}.RData"),
		barcode_names = os.path.join(config["barcodeDir"])		
	output:
		batch_correlation_mtx_RDS = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/batch.correlation.RDS"),
		batch_correlation_pdf = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_batch.correlation.heatmap.pdf"),
		batch_topic_list = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/batch.topics.txt")
	params:
		time = "3:00:00",
		mem_gb = "64",
		figdir = os.path.join(config["figDir"], "{folder}"), 
		analysisdir = os.path.join(config["analysisDir"], "{folder}"), # K{k}/threshold_{threshold}
		threshold = get_cNMF_filter_threshold_double,
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/batch.topic.correlation.R \
		--figdir {params.figdir} \
		--outdir {params.analysisdir} \
		--sampleName {wildcards.sample} \
		--K.val {wildcards.k} \
		--density.thr {params.threshold} \
		--barcode.names {input.barcode_names} \
		--recompute F ' "

def get_fimo_results(wildcards):
	if os.path.isfile(config["fimo_formatted"]):
		return(config["fimo_formatted"])
	else:
		return(os.path.join(config["analysisDir"], "{folder}/{sample}/fimo/fimo_out/fimo.txt"))


rule motif_enrichment_analysis:
	input:
		cNMF_Results = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/cNMF_results.k_{k}.dt_{threshold}.RData"),
		fimo_formatted = get_fimo_results # os.path.join(config["analysisDir"], "{folder}/{sample}/fimo/fimo_out/fimo.formatted.tsv")
	output: 
		# motif_enrichment = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/cNMFAnalysis.factorMotifEnrichment.k_{k}.dt_{threshold}.RData")
		# motif_table_qval = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/{ep_type}.topic.top.300.zscore.gene_motif.count.ttest.enrichment_motif.thr.qval0.1_k_{k}.dt_{threshold}.txt")
		motif_table = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/{ep_type}.topic.top.300.zscore.gene_motif.count.ttest.enrichment_motif.thr.{motif_match_thr}_k_{k}.dt_{threshold}.txt")
	params:
		time = "6:00:00",
		mem_gb = "64",
		figdir = os.path.join(config["figDir"], "{folder}"), 
		analysisdir = os.path.join(config["analysisDir"], "{folder}"), # K{k}/threshold_{threshold}
		threshold = get_cNMF_filter_threshold_double,
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/motif_enrichment.R \
		--sampleName {wildcards.sample} \
		--figdir {params.figdir}/ \
		--outdir {params.analysisdir} \
		--K.val {wildcards.k} \
		--density.thr {params.threshold} \
		--recompute F \
		--ep.type {wildcards.ep_type} \
		--motif.match.thr.str {wildcards.motif_match_thr} \
		--motif.enhancer.background {input.fimo_formatted} \
		--motif.promoter.background /oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2104_remove_lincRNA/data/fimo_out_all_promoters_thresh1.0E-4/fimo.tsv \
		' "
		# --motif.enhancer.background /oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/data/fimo_out_ABC_TeloHAEC_Ctrl_thresh1.0E-4/fimo.formatted.tsv \
		

## old
# rule motif_enrichment_analysis_plot:
# 	input:
# 		# motif_enrichment = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/cNMFAnalysis.factorMotifEnrichment.k_{k}.dt_{threshold}.RData")
# 		motif_table_qval = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/{ep_type}.topic.top.300.zscore.gene_motif.count.ttest.enrichment_motif.thr.qval0.1_k_{k}.dt_{threshold}.txt")
# 	output:
# 		# motif_plot = expand(os.path.join(config["figDir"], "{{folder}}/{{sample}}/K{{k}}/{{sample}}_K{{k}}_dt_{{threshold}}_zscore.{ep_type}.motif.enrichment{test_type_and_threshold}.pdf"), ep_type = ["enhancer", "promoter"], test_type_and_threshold = ["", "_motif.thr.10e-6", ".by.count.ttest", ".by.count.ttest_motif.thr.10e-6", ".by.count.ttest.labelPos", ".by.count.ttest_motif.thr.10e-6.labelPos"])
# 		motif_plot = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_zscore.{ep_type}.motif.{test_type_and_threshold}.pdf")
# 	params:
# 		time = "3:00:00",
# 		mem_gb = "64",
# 		figdir = os.path.join(config["figDir"], "{folder}"), 
# 		analysisdir = os.path.join(config["analysisDir"], "{folder}"), # K{k}/threshold_{threshold}
# 		threshold = get_cNMF_filter_threshold_double
# 	shell:
# 		"bash -c ' source $HOME/.bashrc; \
# 		conda activate cnmf_analysis_R; \
# 		Rscript workflow/scripts/cNMF_analysis_motif.enrichment_plot.R \
# 		--sampleName {wildcards.sample} \
# 		--ep.type {wildcards.ep_type} \
# 		--figdir {params.figdir}/ \
# 		--outdir {params.analysisdir} \
# 		--K.val {wildcards.k} \
# 		--density.thr {params.threshold} \
# 		--recompute F ' "
## end of old

rule motif_enrichment_analysis_plot:
	input:
		# motif_enrichment = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/cNMFAnalysis.factorMotifEnrichment.k_{k}.dt_{threshold}.RData")
		motif_table = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/{ep_type}.topic.top.300.zscore.gene_motif.count.ttest.enrichment_motif.thr.{motif_match_thr}_k_{k}.dt_{threshold}.txt")
	output:
		# motif_plot = expand(os.path.join(config["figDir"], "{{folder}}/{{sample}}/K{{k}}/{{sample}}_K{{k}}_dt_{{threshold}}_zscore.{ep_type}.motif.enrichment{test_type_and_threshold}.pdf"), ep_type = ["enhancer", "promoter"], test_type_and_threshold = ["", "_motif.thr.10e-6", ".by.count.ttest", ".by.count.ttest_motif.thr.10e-6", ".by.count.ttest.labelPos", ".by.count.ttest_motif.thr.10e-6.labelPos"])
		motif_plot = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_zscore.{ep_type}.motif.{test_type}_motif.thr.{motif_match_thr}.pdf")
	params:
		time = "3:00:00",
		mem_gb = "64",
		figdir = os.path.join(config["figDir"], "{folder}"), 
		analysisdir = os.path.join(config["analysisDir"], "{folder}"), # K{k}/threshold_{threshold}
		threshold = get_cNMF_filter_threshold_double,
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/cNMF_analysis_motif.enrichment_plot.R \
		--sampleName {wildcards.sample} \
		--ep.type {wildcards.ep_type} \
		--figdir {params.figdir}/ \
		--outdir {params.analysisdir} \
		--K.val {wildcards.k} \
		--density.thr {params.threshold} \
		--motif.match.thr.str {wildcards.motif_match_thr} \
		--recompute F ' "


# rule fgsea:
# 	input:
# 		cNMF_Results = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/cNMF_results.k_{k}.dt_{threshold}.RData")
# 	output:
# 		fgsea_result = expand(os.path.join(config["analysisDir"], "{{folder}}/{{sample}}/K{{k}}/threshold_{{threshold}}/fgsea/fgsea_all_pathways_df_{ranking_type}_k_{{k}}.dt_{{threshold}}.RData"), ranking_type=["raw.score", "z.score"])
# 	params:
# 		time = "6:00:00",
# 		mem_gb = "64",
# 		outdir = os.path.join(config["analysisDir"], "{folder}_acrossK/{sample}"),
# 		figdir = os.path.join(config["figDir"], "{folder}"), 
# 		analysisdir = os.path.join(config["analysisDir"], "{folder}"), # K{k}/threshold_{threshold}
# 		threshold = get_cNMF_filter_threshold_double
# 	shell:
# 		"bash -c ' source $HOME/.bashrc; \
# 		conda activate cnmf_analysis_R; \
# 		Rscript workflow/scripts/cNMF_analysis_fgsea.R \
# 		--topic.model.result.dir {params.outdir} \
# 		--sampleName {wildcards.sample} \
# 		--figdir {params.figdir}/ \
# 		--outdir {params.analysisdir} \
# 		--K.val {wildcards.k} \
# 		--density.thr {params.threshold} \
# 		--recompute F ' "
	


rule clusterProfiler_GSEA:
	input:
		cNMF_Results = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/cNMF_results.k_{k}.dt_{threshold}.RData")
	output:
		clusterProfiler_result = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/clusterProfiler_GeneRankingType{ranking_type}_EnrichmentType{GSEA_type}.txt")
		# clusterProfiler_top300Gene_GO_result = expand(os.path.join(config["analysisDir"], "{{folder}}/{{sample}}/K{{k}}/threshold_{{threshold}}/clusterProfiler_top300Genes{ranking_type_here}_GOEnrichment.txt"), ranking_type_here = ["zscore", "raw"]),
		# clusterProfiler_allGeneWeight_gsea_result = expand(os.path.join(config["analysisDir"], "{{folder}}/{{sample}}/K{{k}}/threshold_{{threshold}}/clusterProfiler_allGene{ranking_type_here}_ByWeight_GSEA.txt"), ranking_type_here = ["zscore", "raw"]),
		# clusterProfiler_top300Gene_gsea_result = expand(os.path.join(config["analysisDir"], "{{folder}}/{{sample}}/K{{k}}/threshold_{{threshold}}/clusterProfiler_top300Genes{ranking_type_here}_GSEA.txt"), ranking_type_here = ["zscore", "raw"])
	params:
		time = "7:00:00",
		mem_gb = "64", ## 64 is fine too
		outdir = os.path.join(config["analysisDir"], "{folder}_acrossK/{sample}"),
		figdir = os.path.join(config["figDir"], "{folder}"), 
		analysisdir = os.path.join(config["analysisDir"], "{folder}"), # K{k}/threshold_{threshold}
		threshold = get_cNMF_filter_threshold_double,
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/cNMF_analysis_gsea_clusterProfiler.R \
		--topic.model.result.dir {params.outdir} \
		--sampleName {wildcards.sample} \
		--figdir {params.figdir}/ \
		--outdir {params.analysisdir} \
		--K.val {wildcards.k} \
		--density.thr {params.threshold} \
		--ranking.type {wildcards.ranking_type} \
		--GSEA.type {wildcards.GSEA_type} \
		--recompute F ' "


rule clusterProfiler_GSEA_plot:
	input:
		clusterProfiler_result = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/clusterProfiler_GeneRankingType{ranking_type}_EnrichmentType{GSEA_type}.txt")
	output:
		clusterProfiler_plot = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{sample}_K{k}_dt_{threshold}_top10EnrichedPathways_GeneRankingType{ranking_type}_EnrichmentType{GSEA_type}.pdf")
	params:
		time = "1:00:00",
		mem_gb = "16",
		figdir = os.path.join(config["figDir"], "{folder}"), 
		analysisdir = os.path.join(config["analysisDir"], "{folder}"), # K{k}/threshold_{threshold}
		threshold = get_cNMF_filter_threshold_double,
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/plot_gsea_clusterProfiler.R \
		--sampleName {wildcards.sample} \
		--figdir {params.figdir}/ \
		--outdir {params.analysisdir} \
		--K.val {wildcards.k} \
		--density.thr {params.threshold} \
		--ranking.type {wildcards.ranking_type} \
		--GSEA.type {wildcards.GSEA_type} '"
		

# rule clusterTopics:
# 	input: #add input
# 	output: #add output
# 	params:
# 		outdir = config["scratchDir"],
# 		figdir = config["figDir"],
# 		analysisdir = config["analysisDir"],
# 		subset_lib_name = config["subset_lib_name"]
# 	shell:
# 		"Rscript workflow/scripts/cluster.topics.R \
# 		--outdir {params.outdir}/clsuter.topics/ \
# 		--figdir {params.figdir}/cluster.topics/ \
# 		--topic.model.result.dir.list {params.analysisdir}/{sample},{params.analysisdir}/{sample}_subsample_Gene_0_95/,{params.analysisdir}/{sample}_subsample_Guide_0_95/ \
# 		--topic.model.result.dir.list.name {params.subset_lib_name} \
# 		--K.val {k} \
# 		--sampleName {sample} \
# 		--density.thr 0.2 \
# 		--cluster.topic.expression T "


## Topic Summary Table ##here
rule TopicSummary:
	input:
		batch_correlation_mtx_RDS = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/batch.correlation.RDS"),
		cNMF_Results = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/cNMF_results.k_{k}.dt_{threshold}.RData"),
		clusterProfiler_median_spectra_zscore_result = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/clusterProfiler_GeneRankingTypemedian_spectra_zscore_EnrichmentTypeGOEnrichment.txt"),
		clusterProfiler_result = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/clusterProfiler_GeneRankingTypezscore_EnrichmentTypeGOEnrichment.txt")
	output:
		ProgramSummary = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/{sample}_ProgramSummary_k_{k}.dt_{threshold}.txt")
	params:
		time = "1:00:00",
		mem_gb = "16",
		analysisdir = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/"), # K{k}/threshold_{threshold}
		threshold = get_cNMF_filter_threshold_double,
		perturbseq = "T" if config["Perturb-seq"] == "True" else "F",
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
			conda activate cnmf_analysis_R; \
			Rscript workflow/scripts/create_program_summary_table.R \
			--sampleName {wildcards.sample} \
			--outdir {params.analysisdir} \
			--K.val {wildcards.k} \
			--density.thr {params.threshold} \
			--perturbSeq {params.perturbSeq} '"



rule aggregate_over_K:
	input:
		cNMF_Results = expand(os.path.join(config["analysisDir"], "{{folder}}/{{sample}}/K{k}/threshold_{threshold}/cNMF_results.k_{k}.dt_{threshold}.RData"), k=[str(k) for k in config["k"]], threshold=[n.replace(".","_") for n in config["thresholds"]]),
		motif_table = expand(os.path.join(config["analysisDir"], "{{folder}}/{{sample}}/K{k}/threshold_{threshold}/{ep_type}.topic.top.300.zscore.gene_motif.count.ttest.enrichment_motif.thr.{motif_match_thr}_k_{k}.dt_{threshold}.txt"), k=[str(k) for k in config["k"]], threshold=[n.replace(".","_") for n in config["thresholds"]], ep_type = ["enhancer", "promoter"], motif_match_thr = ["qval0.1", "pval1e-4", "pval1e-6"]),
		clusterProfiler_result = expand(os.path.join(config["analysisDir"], "{{folder}}/{{sample}}/K{k}/threshold_{threshold}/clusterProfiler_GeneRankingType{ranking_type}_EnrichmentType{GSEA_type}.txt"), k=[str(k) for k in config["k"]], threshold=[n.replace(".","_") for n in config["thresholds"]], ranking_type = ["zscore"], GSEA_type = ["GOEnrichment", "ByWeightGSEA", "GSEA"]), ## skipping raw x ByWeightGSEA
		variance_explained_txt = expand(os.path.join(config["analysisDir"],"{{folder}}/{{sample}}/K{k}/threshold_{threshold}/metrics.varianceExplained.df.txt"), k=[str(k) for k in config["k"]], threshold=[n.replace(".","_") for n in config["thresholds"]]),
		variance_explained_summary_txt = expand(os.path.join(config["analysisDir"],"{{folder}}/{{sample}}/K{k}/threshold_{threshold}/summary.varianceExplained.df.txt"), k=[str(k) for k in config["k"]], threshold=[n.replace(".","_") for n in config["thresholds"]])
		# clusterProfiler_result = expand(os.path.join(config["analysisDir"], "{{folder}}/{{sample}}/K{k}/threshold_{threshold}/clusterProfiler_GeneRankingType{ranking_type}_EnrichmentType{GSEA_type}.txt"), k=[str(k) for k in config["k"]], threshold=[n.replace(".","_") for n in config["thresholds"]], ranking_type = ["zscore", "raw"], GSEA_type = ["GOEnrichment", "ByWeightGSEA", "GSEA"])
		# clusterProfiler_top300Gene_GO_result = expand(os.path.join(config["analysisDir"], "{{folder}}/{{sample}}/K{k}/threshold_{threshold}/clusterProfiler_top300Genes{ranking_type_here}_GOEnrichment.txt"), k=[str(k) for k in config["k"]], threshold=[n.replace(".","_") for n in config["thresholds"]], ranking_type_here = ["zscore", "raw"]),
		# clusterProfiler_allGeneWeight_gsea_result = expand(os.path.join(config["analysisDir"], "{{folder}}/{{sample}}/K{k}/threshold_{threshold}/clusterProfiler_allGene{ranking_type_here}_ByWeight_GSEA.txt"), k=[str(k) for k in config["k"]], threshold=[n.replace(".","_") for n in config["thresholds"]], ranking_type_here = ["zscore", "raw"]),
		# clusterProfiler_top300Gene_gsea_result = expand(os.path.join(config["analysisDir"], "{{folder}}/{{sample}}/K{k}/threshold_{threshold}/clusterProfiler_top300Genes{ranking_type_here}_GSEA.txt"), k=[str(k) for k in config["k"]], threshold=[n.replace(".","_") for n in config["thresholds"]], ranking_type_here = ["zscore", "raw"])
		# fgsea_result = expand(os.path.join(config["analysisDir"], "{{folder}}/{{sample}}/K{k}/threshold_{threshold}/fgsea/fgsea_all_pathways_df_{ranking_type}_k_{k}.dt_{threshold}.RData"), ranking_type=["raw.score", "z.score"], k=[str(k) for k in config["k"]], threshold=[n.replace(".","_") for n in config["thresholds"]]),
		# motif_table_qval = expand(os.path.join(config["analysisDir"], "{{folder}}/{{sample}}/K{k}/threshold_{threshold}/{ep_type}.topic.top.300.zscore.gene_motif.count.ttest.enrichment_motif.thr.qval0.1_k_{k}.dt_{threshold}.txt"), ep_type = ["enhancer", "promoter"], k=[str(k) for k in config["k"]], threshold=[n.replace(".","_") for n in config["thresholds"]])
		# motif_enrichment = expand(os.path.join(config["analysisDir"], "{{folder}}/{{sample}}/K{k}/threshold_{threshold}/cNMFAnalysis.factorMotifEnrichment.k_{k}.dt_{threshold}.RData"), k=[str(k) for k in config["k"]], threshold=[n.replace(".","_") for n in config["thresholds"]])
		# cNMF_Analysis = expand(os.path.join(config["analysisDir"], "{{sample}}/{{folder}}/K{k}/threshold_0_2/cNMFAnalysis.k_{k}.dt_0_2.RData"), k=config["k"])
	output:
		aggregated_output = os.path.join(config["analysisDir"], "{folder}/{sample}/acrossK/aggregated.outputs.findK.RData")
	params:
		time = "3:00:00",
		mem_gb = "64",
		figdir = os.path.join(config["figDir"], "{folder}"),
		analysisdir = os.path.join(config["analysisDir"], "{folder}"),
		datadir = config["dataDir"],
		klist_comma = ",".join(str(k) for k in config["k"]),
		K_spectra_threshold_table = config["K_spectra_threshold_table"],
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/aggregate_across_K.R \
		--figdir {params.figdir} \
		--outdir {params.analysisdir} \
		--sampleName {wildcards.sample} \
		--K.list {params.klist_comma} \
		--K.table {params.K_spectra_threshold_table} ' " ## how to create this automatically?


rule findK_plot:
	input:
		toplot = os.path.join(config["analysisDir"], "{folder}/{sample}/acrossK/aggregated.outputs.findK.RData")
	output:
		GSEA_plots = os.path.join(config["figDir"], "{folder}/{sample}/acrossK/All_GSEA.pdf"),
		TFMotifEnrichment_plots = os.path.join(config["figDir"], "{folder}/{sample}/acrossK/All_TFMotifEnrichment.pdf"),
		topic_clustering_plot = os.path.join(config["figDir"], "{folder}/{sample}/acrossK/cluster.topic.zscore.by.Pearson.corr.pdf"),
		variance_explained_plot = os.path.join(config["figDir"], "{folder}/{sample}/acrossK/variance.explained.by.model.pdf")
	params:
		time = "3:00:00",
		mem_gb = "64",
		figdir = os.path.join(config["figDir"], "{folder}/{sample}/acrossK/"),
		analysisdir = os.path.join(config["analysisDir"], "{folder}/{sample}/acrossK/"),
		GO_threshold = 0.1,
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/cNMF_findK_plots.R \
		--figdir {params.figdir} \
		--outdir {params.analysisdir} \
		--sampleName {wildcards.sample} \
		--p.adj.threshold 0.1 \
		--aggregated.data {input.toplot} \
		' "


##########################################################################################
### run PoPS and analysis

checkpoint munge_features_without_cNMF:
	input:
		magma_genes_out=os.path.join(config["magma_dir"], "{magma_prefix}.genes.out"),
		magma_genes_raw=os.path.join(config["magma_dir"], "{magma_prefix}.genes.raw")
	output:
		munged_features = directory(os.path.join(config["scratchDir"], "pops/features/pops_features_munged_{magma_prefix}"))
	params:
		time = "2:00:00",
		mem_gb = "16",
		raw_features_dir = os.path.join(config["scratchDir"], "pops/features/pops_features_raw"),
		external_features = config["PoPS_raw_featureDir"],
		gene_annot_path = config["PoPS_gene_annotation_path"],
		munged_features_dir = os.path.join(config["scratchDir"], "pops/features/pops_features_munged_{magma_prefix}"),
		magma_prefix = config["magma_prefix"],
		pipelineDir = config["pipelineDir"],
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		mkdir -p {output.munged_features} {params.raw_features_dir}; \
		cp -r {params.external_features}/* {params.raw_features_dir}/; \
		cd {params.munged_features_dir}; \
		python {params.pipelineDir}/workflow/scripts/pops/munge_feature_directory.py \
       		--gene_annot_path {params.gene_annot_path} \
       		--feature_dir {params.raw_features_dir} \
       		--save_prefix {wildcards.magma_prefix} \
       		--nan_policy zero  ' "


checkpoint munge_features_with_cNMF: ## need to test pops in cnmf_env 
	input:
		magma_genes_out=os.path.join(config["magma_dir"], "{magma_prefix}.genes.out"),
		magma_genes_raw=os.path.join(config["magma_dir"], "{magma_prefix}.genes.raw"),
		cNMF_ENSG_topic_zscore_scaled = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/topic.zscore.ensembl.scaled_k_{k}.dt_{threshold}.txt")
	output:
		munged_features = directory(os.path.join(config["scratchDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/features/pops_features_munged_{magma_prefix}_cNMF"))
	params:
		time = "2:00:00",
		mem_gb = "16",
		raw_features_dir = os.path.join(config["scratchDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/features/pops_features_raw"),
		external_features = config["PoPS_raw_featureDir"],
		gene_annot_path = config["PoPS_gene_annotation_path"],
		munged_features_dir = os.path.join(config["scratchDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/features/pops_features_munged_{magma_prefix}_cNMF"),
		pipelineDir = config["pipelineDir"],
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		rm -r {params.raw_features_dir}; \
		mkdir -p {output.munged_features} {params.raw_features_dir}; \
 		cp -r {params.external_features}/* {params.raw_features_dir}/; \
		cp {input.cNMF_ENSG_topic_zscore_scaled} {params.raw_features_dir}/; \
		cd {params.munged_features_dir}; \
		python {params.pipelineDir}/workflow/scripts/pops/munge_feature_directory.py \
       		--gene_annot_path {params.gene_annot_path} \
       		--feature_dir {params.raw_features_dir} \
       		--save_prefix {wildcards.magma_prefix}_cNMF{wildcards.k} \
       		--nan_policy zero  ' "



## helper to get number of chunck pops mung_feature_directory.py created
def get_pops_feature_chunk_number_with_cNMF(wildcards):
    checkpoint_output = checkpoints.munge_features_with_cNMF.get(**wildcards).output[0]
    # print(os.path.join(checkpoint_output, wildcards.magma_prefix + "_cNMF" + wildcards.k + ".cols.{i}.txt"))
    chunks = glob_wildcards(os.path.join(checkpoint_output, wildcards.magma_prefix + "_cNMF" + wildcards.k + ".cols.{i}.txt")).i
    # print(max([int(i) for i in chunks])+1)
    return max([int(i) for i in chunks])+1

def get_pops_feature_chunk_number_without_cNMF(wildcards):
    checkpoint_output = checkpoints.munge_features_without_cNMF.get(**wildcards).output[0]
    chunks = glob_wildcards(os.path.join(checkpoint_output, wildcards.magma_prefix + ".cols.{i}.txt")).i
    # print(max([int(i) for i in chunks])+1)
    return max([int(i) for i in chunks])+1


## helper to get munged features files as inputs to run_PoPS_cNMF
def get_pops_feature_checkpoint_output_with_cNMF(wildcards): # if "cNMF" in wildcards.prefix
	checkpoint_output = checkpoints.munge_features_with_cNMF.get(**wildcards).output[0]
	return expand(os.path.join(config["scratchDir"], "{{wildcards.folder}}/{{wildcards.sample}}/K{{wildcards.k}}/threshold_{{wildcards.threshold}}/pops/features/pops_features_munged_{{magma_prefix}}_cNMF/{{wildcards.magma_prefix}}_cNMF{{wildcards.k}}.cols.{i}.txt"), i=glob_wildcards(os.path.join(checkpoint_output, "pops_features.cols.{i}.txt")).i)

def get_pops_feature_checkpoint_output_without_cNMF(wildcards): # if "cNMF" in wildcards.prefix
	checkpoint_output = checkpoints.munge_features_without_cNMF.get(**wildcards).output[0]
	return expand(os.path.join(config["scratchDir"], "pops/features/pops_features_munged_{{wildcards.magma_prefix}}/{{wildcards.magma_prefix}}.cols.{i}.txt"), i=glob_wildcards(os.path.join(checkpoint_output, "{{wildcards.magma_prefix}}.cols.{i}.txt")).i)


rule run_PoPS_with_cNMF: ## num_feature_chunks: needs specification
	input:
		get_pops_feature_checkpoint_output_with_cNMF
	output:
		coefs = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.coefs"),
		marginals = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.marginals"),
		preds = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.preds")
	params:
		time = "2:00:00",
		mem_gb = "32",
		gene_annot_path = config["PoPS_gene_annotation_path"],
		munged_features_dir = os.path.join(config["scratchDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/features/pops_features_munged_{magma_prefix}_cNMF"),
		magma_dir = config["magma_dir"],
		PoPS_control_features = config["PoPS_control_features"],
		outdir = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops"),
		num_munged_feature_chunks = get_pops_feature_chunk_number_with_cNMF,
		pipelineDir = config["pipelineDir"],
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		python {params.pipelineDir}/workflow/scripts/pops/pops.py \
			--gene_annot_path {params.gene_annot_path} \
			--feature_mat_prefix {params.munged_features_dir}/{wildcards.magma_prefix}_cNMF{wildcards.k} \
			--num_feature_chunks {params.num_munged_feature_chunks} \
			--control_features_path {params.PoPS_control_features} \
			--magma_prefix {params.magma_dir}/{wildcards.magma_prefix} \
			--out_prefix {params.outdir}/{wildcards.magma_prefix}_cNMF{wildcards.k} \
			--verbose ' "


rule run_PoPS_without_cNMF: ## num_feature_chunks: needs specification
	input:
		get_pops_feature_checkpoint_output_without_cNMF
	output:
		coefs = os.path.join(config["analysisDir"], "pops/{magma_prefix}.coefs"),
		marginals = os.path.join(config["analysisDir"], "pops/{magma_prefix}.marginals"),
		preds = os.path.join(config["analysisDir"], "pops/{magma_prefix}.preds")
	params:
		time = "2:00:00",
		mem_gb = "32",
		gene_annot_path = config["PoPS_gene_annotation_path"],
		munged_features_dir = os.path.join(config["scratchDir"], "pops/features/pops_features_munged_{magma_prefix}"),
		magma_dir = config["magma_dir"],
		PoPS_control_features = config["PoPS_control_features"],
		outdir = os.path.join(config["analysisDir"], "pops"),
		num_munged_feature_chunks = get_pops_feature_chunk_number_without_cNMF,
		pipelineDir = config["pipelineDir"],
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		python {params.pipelineDir}/workflow/scripts/pops/pops.py \
			--gene_annot_path {params.gene_annot_path} \
			--feature_mat_prefix {params.munged_features_dir}/{wildcards.magma_prefix} \
			--num_feature_chunks {params.num_munged_feature_chunks} \
			--control_features_path {params.PoPS_control_features} \
			--magma_prefix {params.magma_dir}/{wildcards.magma_prefix} \
			--out_prefix {params.outdir}/{wildcards.magma_prefix} \
			--verbose ' "


def get_features_txt_stream(wildcards):
	raw_features_dir = os.path.join(config["scratchDir"], "pops/features/pops_features_raw")
	import glob
	return glob.glob(os.path.join(raw_features_dir, "*.txt"))


rule aggregate_features_txt_to_RDS: ## aggregate all external features
	input:
		get_features_txt_stream
	output:
		all_features_RDS = os.path.join(config["analysisDir"], "pops/features/processed_features/full_external_features.RDS"),
		all_features_txt=os.path.join(config["analysisDir"], "pops/features/processed_features/full_external_features.txt")
	params:
		time = "3:00:00",
		mem_gb = "128",
		outdir = os.path.join(config["analysisDir"], "pops/features/processed_features"),
		raw_features_dir = os.path.join(config["scratchDir"], "pops/features/pops_features_raw"),
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		mkdir -p {params.outdir}; \
		Rscript workflow/scripts/PoPS_aggregate_features.R \
		--feature.dir {params.raw_features_dir} \
		--output {params.outdir}/ ' "


rule aggregate_features_txt_to_RDS_with_cNMF:
	input:
		all_features_RDS = os.path.join(config["analysisDir"], "pops/features/processed_features/full_external_features.RDS"),
		cNMF_ENSG_topic_zscore_scaled = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/topic.zscore.ensembl.scaled_k_{k}.dt_{threshold}.txt")
	output:
		all_features_with_cNMF_RDS=os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/full_features_{magma_prefix}_cNMF{k}.RDS")
	params:
		time = "3:00:00",
		mem_gb = "128",
		outdir = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops"),
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		mkdir -p {params.outdir}; \
		Rscript workflow/scripts/PoPS_aggregate_features_with_cNMF.R \
		--feature.RDS {input.all_features_RDS} \
		--cNMF.features {input.cNMF_ENSG_topic_zscore_scaled} \
		--output {params.outdir}/ \
		--prefix {wildcards.magma_prefix}_cNMF{wildcards.k} ' "


rule analyze_PoPS_results: ## need RDS file with all features
	input:
		coefs_with_cNMF = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.coefs"),
		marginals_with_cNMF = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.marginals"),
		preds_with_cNMF = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.preds"),
		coefs_without_cNMF = os.path.join(config["analysisDir"], "pops/{magma_prefix}.coefs"),
		marginals_without_cNMF = os.path.join(config["analysisDir"], "pops/{magma_prefix}.marginals"),
		preds_without_cNMF = os.path.join(config["analysisDir"], "pops/{magma_prefix}.preds"),
		all_features_with_cNMF_RDS=os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/full_features_{magma_prefix}_cNMF{k}.RDS"),
		cNMF_ENSG_topic_zscore_scaled = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/topic.zscore.ensembl.scaled_k_{k}.dt_{threshold}.txt")
		# all_features_txt=os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/full_features_{magma_prefix}_cNMF{k}.txt")
	output:
		feature_x_gene_importance_score_RDS = os.path.join(config["scratchDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}_coefs.marginals.feature.outer.prod.RDS"),
		combined_preds = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.combined.preds"),
		coefs_defining_top_topic_RDS = os.path.join(config["scratchDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}_coefs.defining.top.topic.RDS"),
		preds_importance_score_key_columns = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}_PoPS_preds.importance.score.key.columns.txt")
	params:
		time = "3:00:00",
		mem_gb = "128",
		outdir = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops"),
		scratch_outdir = os.path.join(config["scratchDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops"),
		external_features_metadata = config["PoPS_features_metadata_path"],
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/PoPS.data.processing.R \
		--output {params.outdir}/ \
		--scratch.output {params.scratch_outdir}/ \
		--prefix {wildcards.magma_prefix}_cNMF{wildcards.k} \
		--external.features.metadata {params.external_features_metadata} \
		--coefs_with_cNMF {input.coefs_with_cNMF} \
		--preds_with_cNMF {input.preds_with_cNMF} \
		--marginals_with_cNMF {input.marginals_with_cNMF} \
		--coefs_without_cNMF {input.coefs_without_cNMF} \
		--preds_without_cNMF {input.preds_without_cNMF} \
		--marginals_without_cNMF {input.marginals_with_cNMF} \
		--cNMF.features {input.cNMF_ENSG_topic_zscore_scaled} \
		--all.features {input.all_features_with_cNMF_RDS} \
		--recompute F ' "


rule PoPS_plots:
	input:
		coefs_with_cNMF = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.coefs"),
		marginals_with_cNMF = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.marginals"),
		preds_with_cNMF = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.preds"),
		coefs_without_cNMF = os.path.join(config["analysisDir"], "pops/{magma_prefix}.coefs"),
		marginals_without_cNMF = os.path.join(config["analysisDir"], "pops/{magma_prefix}.marginals"),
		preds_without_cNMF = os.path.join(config["analysisDir"], "pops/{magma_prefix}.preds"),
		all_features_with_cNMF_RDS=os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/full_features_{magma_prefix}_cNMF{k}.RDS"),
		cNMF_ENSG_topic_zscore_scaled = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/topic.zscore.ensembl.scaled_k_{k}.dt_{threshold}.txt"),
		feature_x_gene_importance_score_RDS = os.path.join(config["scratchDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}_coefs.marginals.feature.outer.prod.RDS"),
		combined_preds = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.combined.preds"),
		coefs_defining_top_topic_RDS = os.path.join(config["scratchDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}_coefs.defining.top.topic.RDS"),
		preds_importance_score_key_columns = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}_PoPS_preds.importance.score.key.columns.txt")
	output: ## double check: threshold
		feature_x_gene_component_importance_score_coefs_plot = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{magma_prefix}_cNMF{k}_dt_{threshold}_feature_x_gene.component.importance.score.coefs.pdf"),
		topic_coef_beta_plot = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{magma_prefix}_cNMF{k}_dt_{threshold}_Topic.coef.beta.pdf"),
		all_coef_beta_plot = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{magma_prefix}_cNMF{k}_dt_{threshold}_all.coef.beta.pdf"),
		PoPS_score_list_plot = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{magma_prefix}_cNMF{k}_dt_{threshold}_PoPS_score_list.pdf"),
		with_and_without_cNMF_comparison_plot = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{magma_prefix}_cNMF{k}_dt_{threshold}_before.vs.after.cNMF.pdf")
	params:
		time = "3:00:00",
		mem_gb = "64",
		outdir = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops"),
		figdir = os.path.join(config["figDir"], "{folder}/{sample}/K{k}"),
		scratch_outdir = os.path.join(config["scratchDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops"),
		external_features_metadata = config["PoPS_features_metadata_path"],
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/PoPS.plots.R \
		--sampleName {wildcards.sample} \
		--output {params.outdir}/ \
		--figure {params.figdir}/ \
		--scratch.output {params.scratch_outdir}/ \
		--prefix {wildcards.magma_prefix}_cNMF{wildcards.k} \
		--k.val {wildcards.k} \
		--coefs_with_cNMF {input.coefs_with_cNMF} \
		--preds_with_cNMF {input.preds_with_cNMF} \
		--marginals_with_cNMF {input.marginals_with_cNMF} \
		--coefs_without_cNMF {input.coefs_without_cNMF} \
		--preds_without_cNMF {input.preds_without_cNMF} \
		--marginals_without_cNMF {input.marginals_with_cNMF} \
		--cNMF.features {input.cNMF_ENSG_topic_zscore_scaled} \
		--all.features {input.all_features_with_cNMF_RDS} \
		--external.features.metadata {params.external_features_metadata} \
		--combined.preds {input.combined_preds} \
		--coefs.defining.top.topic.RDS {input.coefs_defining_top_topic_RDS} \
		--preds.importance.score.key.columns {input.preds_importance_score_key_columns} ' "


######################################################################################### END OF PoPS


## IGVF formatting:
rule IGVF_formatting_model_programGenes:
	input: cNMF_Results = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/cNMF_results.k_{k}.dt_{threshold}.RData"),
	output: model_yaml = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/IGVF_format/{sample}.k_{k}.dt_{threshold}.modelYAML.yaml")
	params:
		time = "1:00:00",
		mem_gb = "8",
		analysisdir = os.path.join(config["analysisDir"], "{folder}"), # K{k}/threshold_{threshold}
		partition = "owners,normal",
		threshold = get_cNMF_filter_threshold_double
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/output_IGVF_format.R \
			--sampleName {wildcards.sample} \
			--outdir {params.analysisdir}/ \
			--K.val {wildcards.k} \
			--density.thr {params.threshold} \
		' "

rule IGVF_formatting_model_cellxgene:
	input: cNMF_usage_output = os.path.join(config["analysisDir"], "{folder}_acrossK/{sample}/{sample}.usages.k_{k}.dt_{threshold}.consensus.txt")
	output: cellxgene_h5ad = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/IGVF_format/{sample}.k_{k}.dt_{threshold}.cellxgene.h5ad")
	params:
		time = "1:00:00",
		mem_gb = "8",
		analysisdir = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/IGVF_format/"), # K{k}/threshold_{threshold}
		cNMF_outdir = os.path.join(config["analysisDir"], "{folder}_acrossK"),
		partition = "owners,normal",
		threshold = get_cNMF_filter_threshold_double,
		barcode_dir = config["barcodeDir"]
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		python workflow/scripts/create_cellxgene_h5ad_IGVF_format.py \
			--path_to_topics {params.cNMF_outdir} \
			--topic_sampleName {wildcards.sample} \
			--outdir {params.analysisdir} \
			--k {wildcards.k} \
			--density_threshold {params.threshold} \
			--barcode_dir {params.barcode_dir} \
		' "


# rule concatenate_findK:
# 	input:
# 		large_summary_file = os.path.join(config["analysisDir"], "{sample}/aggregated.outputs.findK.RData")
# 	output:
# 		toplot = os.path.join(config["analysisDir"], "{sample}/findK.toPlot.RData")
# 	shell:
# 		"Rscript workflow/scripts/concatenate.findK.input.R --input.path {input.large_summary_file} --output.path {output.toplot}"
#
# ## RuleException: CalledProcessError in line 165 of /oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/Perturb-seq_CAD/workflow/rules/Perturb-seq_analysis.smk: Command 'set -euo pipefail;  Rscript --vanilla -e 'rmarkdown::render("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/Perturb-seq_CAD/.snakemake/scripts/tmp0ezqneja.findK.plots.Rmd", output_file="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2105_snakemake_output/analysis/no_IL1B/findK.html", quiet=TRUE, knit_root_dir = "/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/Perturb-seq_CAD", params = list(rmd="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/Perturb-seq_CAD/.snakemake/scripts/tmp0ezqneja.findK.plots.Rmd"))'' returned non-zero exit status 1.
# ## File "/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/Perturb-seq_CAD/workflow/rules/Perturb-seq_analysis.smk", line 165, in __rule_findK_plot


	# shell:
	# 	"""
	# 		Rscript -e " .libPaths("/home/groups/engreitz/Software/R_3.6.1") \\
	# 		library("rmarkdown") \\
	# 		rmarkdown::render('../../workflow/scripts/findK.plots.Rmd')"
	# 	"""
