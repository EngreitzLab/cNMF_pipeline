## Prioritize Topics

def get_cNMF_filter_threshold_double(wildcards):
        return wildcards.threshold.replace("_",".")

rule cNMF_PoPS_result_table:
	input:
		# coefs_with_cNMF = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.coefs"),
		# marginals_with_cNMF = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.marginals"),
		# preds_with_cNMF = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.preds"),
		coefs_with_cNMF = os.path.join("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/{magma_prefix}_cNMF60.coefs"),
		marginals_with_cNMF = os.path.join("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/{magma_prefix}_cNMF60.marginals"),
		preds_with_cNMF = os.path.join("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/{magma_prefix}_cNMF60.preds"),
		batch_topic_list = os.path.join(config["analysisDir"], "{sample}/{folder}/{sample}/K{k}/threshold_{threshold}/batch.topics.txt")
	output:
		input_table_for_compute_enrichment = os.path.join(config["analysisDir"], "{sample}/{folder}/{sample}/K{k}/threshold_{threshold}/cNMF.PoPS.input.to.large.table.pipeline_k_{k}.dt_{threshold}.minGuidePerPtb_1.minCellPerGuide_2.{magma_prefix}.txt")
	params:
		time = "6:00:00",
		mem_gb = "64",
		figdir = os.path.join(config["figDir"], "{sample}/{folder}/"), # need to distinguish k and dt
		analysisdir = os.path.join(config["analysisDir"], "{sample}/{folder}/"), # K{k}/threshold_{threshold}
		threshold = get_cNMF_filter_threshold_double,
		popsdir = "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/create_cNMF_PoPS_input_table.R \
		--sampleName {wildcards.sample} \
		--outdir {params.analysisdir} \
		--pops.dir {params.popsdir} \
		--prefix {wildcards.magma_prefix}_cNMF60 \
		--K.val {wildcards.k} \
		--density.thr {params.threshold} ' " # popsdir =  "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/"



rule compute_enrichment:
	input:
		input_table_for_compute_enrichment = os.path.join(config["analysisDir"], "{sample}/{folder}/{sample}/K{k}/threshold_{threshold}/cNMF.PoPS.input.to.large.table.pipeline_k_{k}.dt_{threshold}.minGuidePerPtb_1.minCellPerGuide_2.{magma_prefix}.txt"),
		input_GWAS_table = os.path.join(config["GWASTableDir"], "CAD_GWAS_gene_incl_ubq_genes.{version}.csv")
	output:
		statistical_test_result = os.path.join(config["analysisDir"], "{sample}/{folder}/{sample}/K{k}/threshold_{threshold}/statistical.test.to.prioritize.topics.{magma_prefix}.{version}.txt")
	params:
		time = "6:00:00",
		mem_gb = "64",
		figdir = os.path.join(config["figDir"], "{sample}/{folder}/"), # need to distinguish k and dt
		analysisdir = os.path.join(config["analysisDir"], "{sample}/{folder}/"), # K{k}/threshold_{threshold}
		threshold = get_cNMF_filter_threshold_double
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/compute_enrichment.R \
		--input.GWAS.table {input.input_GWAS_table} \
		--version {wildcards.version} \
		--sampleName {wildcards.sample} \
		--figdir {params.figdir}/ \
		--outdir {params.analysisdir} \
		--K.val {wildcards.k} \
		--prefix {wildcards.magma_prefix} \
		--density.thr {params.threshold} \
		--cNMF.PoPS.input.table {input.input_table_for_compute_enrichment} ' "


rule combine_statistical_tests:
	input:
		statistical_test_result = expand(os.path.join(config["analysisDir"], "{{sample}}/{{folder}}/{{sample}}/K{k}/threshold_{{threshold}}/statistical.test.to.prioritize.topics.{{magma_prefix}}.{{version}}.txt"), k=config["k"])
	output:
		combined_statistical_test_result = os.path.join(config["analysisDir"], "{sample}/{folder}/{sample}/acrossK/unique.K.topic.pairs_statistical.test.to.prioritize.topic.{threshold}.{magma_prefix}.{version}.txt")
	params:
		time = "6:00:00",
		mem_gb = "64",
		figdir = os.path.join(config["figDir"], "{sample}/{folder}/{sample}/acrossK/"), # need to distinguish k and dt
		analysisdir = os.path.join(config["analysisDir"], "{sample}/{folder}/"), # K{k}/threshold_{threshold}
		threshold = get_cNMF_filter_threshold_double,
		klist_comma = ",".join(str(k) for k in config["k"])
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		echo \" {params.klist_comma} \"; \
		Rscript workflow/scripts/combine_statistical_test_results.R \
		--version {wildcards.version} \
		--sampleName {wildcards.sample} \
		--figdir {params.figdir}/ \
		--outdir {params.analysisdir} \
		--prefix {wildcards.magma_prefix} \
		--density.thr {params.threshold} ' "



#
