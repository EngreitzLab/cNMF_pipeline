### Run cNMF on Perturb-seq data and analyze the results

import subprocess
import numpy as np
import math


### run PoPS and analysis
checkpoint munge_features: # add script to allow for curating features from mutiple directories based on a list in config file
	input:
		magma_genes_out=os.path.join(config["magma_dir"], "{magma_prefix}.genes.out"),
		magma_genes_raw=os.path.join(config["magma_dir"], "{magma_prefix}.genes.raw")
	output:
		munged_features = directory(os.path.join(config["scratchDir"], "pops/features/pops_features_munged_{magma_prefix}_{sample}"))
	params:
		time = "2:00:00",
		mem_gb = "16",
		raw_features_dir = os.path.join(config["scratchDir"], "pops/features/pops_features_raw_{magma_prefix}_{sample}"),
		external_features = config["PoPS_raw_featureDir"],
		gene_annot_path = config["PoPS_gene_annotation_path"],
		munged_features_dir = os.path.join(config["scratchDir"], "pops/features/pops_features_munged_{magma_prefix}_{sample}"),
		magma_prefix = config["magma_prefix"],
		pipelineDir = config["pipelineDir"]
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		mkdir -p {output.munged_features} {params.raw_features_dir}; \
		cp -r {params.external_features}/* {params.raw_features_dir}/; \
		cd {params.munged_features_dir}; \
		python {params.pipelineDir}/workflow/scripts/pops/munge_feature_directory.py \
       		--gene_annot_path {params.gene_annot_path} \
       		--feature_dir {params.raw_features_dir} \
       		--save_prefix {wildcards.magma_prefix}_{wildcards.sample} \
       		--nan_policy zero  ' "



## helper to get number of chunck pops mung_feature_directory.py created
def get_pops_feature_chunk_number(wildcards):
    checkpoint_output = checkpoints.munge_features.get(**wildcards).output[0]
    print(checkpoint_output)
    chunks = glob_wildcards(os.path.join(checkpoint_output, wildcards.magma_prefix + "_" + wildcards.sample + ".cols.{i}.txt")).i
    print(chunks)
    # print(max([int(i) for i in chunks])+1)
    return max([int(i) for i in chunks])+1


## helper to get munged features files as inputs to run_PoPS_cNMF
def get_pops_feature_checkpoint_output(wildcards): # if "cNMF" in wildcards.prefix
	checkpoint_output = checkpoints.munge_features.get(**wildcards).output[0]
	return expand(os.path.join(config["scratchDir"], "pops/features/pops_features_munged_{{wildcards.magma_prefix}}_{{wildcards.sample}}/{{wildcards.magma_prefix}}_{{wildcards.sample}}.cols.{i}.txt"), i=glob_wildcards(os.path.join(checkpoint_output, "{{wildcards.magma_prefix}}_{{wildcards.sample}}.cols.{i}.txt")).i)


rule run_PoPS: ## num_feature_chunks: needs specification
	input:
		get_pops_feature_checkpoint_output
	output:
		coefs = os.path.join(config["analysisDir"], "pops/{magma_prefix}_{sample}.coefs"),
		marginals = os.path.join(config["analysisDir"], "pops/{magma_prefix}_{sample}.marginals"),
		preds = os.path.join(config["analysisDir"], "pops/{magma_prefix}_{sample}.preds")
	params:
		time = "2:00:00",
		mem_gb = "32",
		gene_annot_path = config["PoPS_gene_annotation_path"],
		munged_features_dir = os.path.join(config["scratchDir"], "pops/features/pops_features_munged_{magma_prefix}_{sample}"),
		magma_dir = config["magma_dir"],
		PoPS_control_features = config["PoPS_control_features"],
		outdir = os.path.join(config["analysisDir"], "pops"),
		num_munged_feature_chunks = get_pops_feature_chunk_number,
		pipelineDir = config["pipelineDir"]
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		python {params.pipelineDir}/workflow/scripts/pops/pops.py \
			--gene_annot_path {params.gene_annot_path} \
			--feature_mat_prefix {params.munged_features_dir}/{wildcards.magma_prefix} \
			--num_feature_chunks {params.num_munged_feature_chunks} \
			--control_features_path {params.PoPS_control_features} \
			--magma_prefix {params.magma_dir}/{wildcards.magma_prefix} \
			--out_prefix {params.outdir}/{wildcards.magma_prefix}_{wildcards.sample} \
			--verbose ' "


def get_features_txt_stream(wildcards):
	raw_features_dir = os.path.join(config["scratchDir"], "pops/features/pops_features_raw_{magma_prefix}_{wildcards.sample}")
	import glob
	return glob.glob(os.path.join(raw_features_dir, "*.txt"))


rule aggregate_features_txt_to_RDS: ## aggregate all external features
	input:
		get_features_txt_stream
	output:
		all_features_RDS = os.path.join(config["analysisDir"], "pops/features/processed_features/full_external_features_{magma_prefix}_{sample}.RDS"),
		all_features_txt=os.path.join(config["analysisDir"], "pops/features/processed_features/full_external_features_{magma_prefix}_{sample}.txt")
	params:
		time = "3:00:00",
		mem_gb = "128",
		outdir = os.path.join(config["analysisDir"], "pops/features/processed_features"),
		raw_features_dir = os.path.join(config["scratchDir"], "pops/features/pops_features_raw_{magma_prefix}_{sample}")
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		mkdir -p {params.outdir}; \
		Rscript workflow/scripts/PoPS_aggregate_features.R \
		--feature.dir {params.raw_features_dir} \
		--output {params.outdir}/ \
		--prefix {wildcards.magma_prefix} \
		--sampleName {wildcards.sample} ' "


rule analyze_PoPS_results: ## need RDS file with all features
	input:
		coefs = os.path.join(config["analysisDir"], "pops/{magma_prefix}_{sample}.coefs"),
		marginals = os.path.join(config["analysisDir"], "pops/{magma_prefix}_{sample}.marginals"),
		preds = os.path.join(config["analysisDir"], "pops/{magma_prefix}_{sample}.preds"),
		all_features_RDS = os.path.join(config["analysisDir"], "pops/features/processed_features/full_external_features_{magma_prefix}_{sample}.RDS"),
		all_features_txt=os.path.join(config["analysisDir"], "pops/features/processed_features/full_external_features_{magma_prefix}_{sample}.txt")
	output:
		feature_x_gene_importance_score_RDS = os.path.join(config["scratchDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}_coefs.marginals.feature.outer.prod.RDS"),
		combined_preds = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.combined.preds"),
		coefs_defining_top_topic_RDS = os.path.join(config["scratchDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}_coefs.defining.top.topic.RDS"),
		preds_importance_score_key_columns = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}_PoPS_preds.importance.score.key.columns.txt")
	params:
		time = "3:00:00",
		mem_gb = "128",
		outdir = os.path.join(config["analysisDir"], "pops"),
		scratch_outdir = os.path.join(config["scratchDir"], "pops"),
		external_features_metadata = config["PoPS_features_metadata_path"]
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/PoPS.data.processing.R \
		--output {params.outdir}/ \
		--scratch.output {params.scratch_outdir}/ \
		--prefix {wildcards.magma_prefix} \
		--sampleName {wildcards.sample} \
		--external.features.metadata {params.external_features_metadata} \
		--coefs {input.coefs} \
		--preds {input.preds} \
		--marginals {input.marginals} \
		--all.features {input.all_features_with_cNMF_RDS} \
		--recompute F ' "


rule PoPS_plots:
	input:
		coefs = os.path.join(config["analysisDir"], "pops/{magma_prefix}_{sample}.coefs"),
		marginals = os.path.join(config["analysisDir"], "pops/{magma_prefix}_{sample}.marginals"),
		preds = os.path.join(config["analysisDir"], "pops/{magma_prefix}_{sample}.preds"),
		all_features_RDS = os.path.join(config["analysisDir"], "pops/features/processed_features/full_external_features_{magma_prefix}_{sample}.RDS"),
		all_features_txt=os.path.join(config["analysisDir"], "pops/features/processed_features/full_external_features_{magma_prefix}_{sample}.txt")
		# feature_x_gene_importance_score_RDS = os.path.join(config["scratchDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}_coefs.marginals.feature.outer.prod.RDS"),
		# combined_preds = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.combined.preds"),
		# coefs_defining_top_topic_RDS = os.path.join(config["scratchDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}_coefs.defining.top.topic.RDS"),
		# preds_importance_score_key_columns = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}_PoPS_preds.importance.score.key.columns.txt")
	output: ## double check: threshold
		# feature_x_gene_component_importance_score_coefs_plot = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{magma_prefix}_cNMF{k}_dt_{threshold}_feature_x_gene.component.importance.score.coefs.pdf"),
		# topic_coef_beta_plot = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{magma_prefix}_cNMF{k}_dt_{threshold}_Topic.coef.beta.pdf"),
		all_coef_beta_plot = os.path.join(config["figDir"], "{magma_prefix}_{sample}_all.coef.beta.pdf"),
		PoPS_score_list_plot = os.path.join(config["figDir"], "{magma_prefix}_{sample}_PoPS_score_list.pdf")
		# with_and_without_cNMF_comparison_plot = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/{magma_prefix}_cNMF{k}_dt_{threshold}_before.vs.after.cNMF.pdf")
	params:
		time = "3:00:00",
		mem_gb = "64",
		outdir = os.path.join(config["analysisDir"], "pops"),
		figdir = os.path.join(config["figDir"]),
		scratch_outdir = os.path.join(config["scratchDir"], "pops"),
		external_features_metadata = config["PoPS_features_metadata_path"]
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/PoPS.plots.R \
		--output {params.outdir}/ \
		--figdir {params.figdir}/ \
		--sampleName {wildcards.sample} \
		--scratch.output {params.scratch_outdir}/ \
		--prefix {wildcards.magma_prefix}_cNMF{wildcards.k} \
		--coefs {input.coefs} \
		--preds {input.preds} \
		--marginals {input.marginals} \
		--all.features {input.all_features_with_cNMF_RDS} ' "

		# --external.features.metadata {params.external_features_metadata} \
		# --combined.preds {input.combined_preds} \
		# --coefs.defining.top.topic.RDS {input.coefs_defining_top_topic_RDS} \
		# --preds.importance.score.key.columns {input.preds_importance_score_key_columns} ' "




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
