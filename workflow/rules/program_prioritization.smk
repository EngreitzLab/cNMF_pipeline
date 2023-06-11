## Prioritize Topics

def get_cNMF_filter_threshold_double(wildcards):
        return wildcards.threshold.replace("_",".")

rule prepare_compute_enrichment:
	input:
		# # coefs_with_cNMF = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.coefs"),
		# # marginals_with_cNMF = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.marginals"),
		# # preds_with_cNMF = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/pops/{magma_prefix}_cNMF{k}.preds"),
		# coefs_with_cNMF = os.path.join("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/{magma_prefix}_cNMF60.coefs"),
		# marginals_with_cNMF = os.path.join("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/{magma_prefix}_cNMF60.marginals"),
		# preds_with_cNMF = os.path.join("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/{magma_prefix}_cNMF60.preds")
		os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/{sample}_MAST_DEtopics.txt") if config["Perturb-seq"] == "True" else os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/cNMF_results.k_{k}.dt_{threshold}.RData")
	output:
		# input_table_for_compute_enrichment = os.path.join(config["analysisDir"], "{sample}/{folder}/{sample}/K{k}/threshold_{threshold}/cNMF.PoPS.input.to.large.table.pipeline_k_{k}.dt_{threshold}.magma_{magma_prefix}.txt")
		input_table_for_compute_enrichment = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/prepare_compute_enrichment.txt")
	params:
		time = "1:00:00",
		mem_gb = "64",
		figdir = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/"), # need to distinguish k and dt
		analysisdir = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/"), # K{k}/threshold_{threshold}
		threshold = get_cNMF_filter_threshold_double,
		popsdir = "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/",
		perturbseq = "T" if config["Perturb-seq"] == "True" else "F",
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/program_prioritization/create_input_table.R \
		--sampleName {wildcards.sample} \
		--outdir {params.analysisdir} \
		--K.val {wildcards.k} \
		--density.thr {params.threshold} \
		--perturbSeq {params.perturbseq}' " # popsdir =  "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/"


		# --pops.dir {params.popsdir} \
		# --prefix {wildcards.magma_prefix}_cNMF60 \


rule compute_enrichment:
	input:
		# input_table_for_compute_enrichment = os.path.join(config["analysisDir"], "{sample}/{folder}/{sample}/K{k}/threshold_{threshold}/cNMF.PoPS.input.to.large.table.pipeline_k_{k}.dt_{threshold}.magma_{magma_prefix}.txt"),
		input_table_for_compute_enrichment = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/prepare_compute_enrichment.txt"),
		input_GWAS_table = os.path.join(config["GWAS_gene_priorization_table_path"], "overlap/{GWAS_trait}/{GWAS_trait}_GWAS_gene_incl_ubq_genes.txt"),
		coding_variant_table = os.path.join(config["GWAS_gene_priorization_table_path"], "{GWAS_trait}/{GWAS_trait}variant.list.1.coordinate.txt"),
		batch_topic_list = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/batch.topics.txt")
	output:
		program_prioritization_statistical_test_result = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/program_prioritization_{regulator_analysis_type}/{GWAS_trait}/{GWAS_trait}.program_prioritization.txt"),
		program_prioritization_plot = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/program_prioritization_{regulator_analysis_type}/{GWAS_trait}/{sample}_K{k}_dt_{threshold}_{GWAS_trait}_ProgramGene_EnrichmentBarPlot.pdf")
	params:
		time = "6:00:00",
		mem_gb = "64",
		figdir = os.path.join(config["figDir"], "{folder}/{sample}/K{k}/program_prioritization_{regulator_analysis_type}/{GWAS_trait}/"), # need to distinguish k and dt
		analysisdir = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/program_prioritization_{regulator_analysis_type}/{GWAS_trait}/"),
		outdirsample = os.path.join(config["analysisDir"], "{folder}/{sample}/K{k}/threshold_{threshold}/"),
		celltype = config["celltype"],
		threshold = get_cNMF_filter_threshold_double,
		perturbseq = "T" if config["Perturb-seq"] == "True" else "F",
		partition = "owners,normal"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_analysis_R; \
		Rscript workflow/scripts/program_prioritization/compute_enrichment.R \
		--input.GWAS.table {input.input_GWAS_table} \
		--coding.variant.df {input.coding_variant_table} \
		--sampleName {wildcards.sample} \
		--outdirsample {params.outdirsample}/ \
		--celltype {params.celltype} \
		--figdir {params.figdir}/ \
		--outdir {params.analysisdir} \
		--K.val {wildcards.k} \
		--trait.name {wildcards.GWAS_trait} \
		--density.thr {params.threshold} \
		--cNMF.table {input.input_table_for_compute_enrichment} \
		--regulator.analysis.type {wildcards.regulator_analysis_type} \
		--perturbSeq {params.perturbseq}' "


# rule combine_statistical_tests:
# 	input:
# 		statistical_test_result = expand(os.path.join(config["analysisDir"], "{{sample}}/{{folder}}/{{sample}}/K{k}/threshold_{{threshold}}/statistical.test.to.prioritize.topics.{{magma_prefix}}.{{version}}.txt"), k=config["k"])
# 	output:
# 		combined_statistical_test_result = os.path.join(config["analysisDir"], "{sample}/{folder}/{sample}/acrossK/unique.K.topic.pairs_statistical.test.to.prioritize.topic.{threshold}.{magma_prefix}.{version}.txt")
# 	params:
# 		time = "6:00:00",
# 		mem_gb = "64",
# 		figdir = os.path.join(config["figDir"], "{sample}/{folder}/{sample}/acrossK/"), # need to distinguish k and dt
# 		analysisdir = os.path.join(config["analysisDir"], "{sample}/{folder}/"), # K{k}/threshold_{threshold}
# 		threshold = get_cNMF_filter_threshold_double,
# 		klist_comma = ",".join(str(k) for k in config["k"])
# 	shell:
# 		"bash -c ' source $HOME/.bashrc; \
# 		conda activate cnmf_analysis_R; \
# 		echo \" {params.klist_comma} \"; \
# 		Rscript workflow/scripts/combine_statistical_test_results.R \
# 		--version {wildcards.version} \
# 		--sampleName {wildcards.sample} \
# 		--figdir {params.figdir}/ \
# 		--outdir {params.analysisdir} \
# 		--prefix {wildcards.magma_prefix} \
# 		--density.thr {params.threshold} ' "



#
