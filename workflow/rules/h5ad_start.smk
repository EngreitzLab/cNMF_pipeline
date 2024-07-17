



def get_raw_h5ad_file(wildcards):
	if os.path.isfile(config["input_h5ad_mtxDir"]):
		return(config["input_h5ad_mtxDir"])
	else:
		return(os.path.join(config["analysisDir"], "data/raw/{sample}.h5ad"))


rule raw_h5ad_to_filtered_h5ad:
	input:
		raw_h5ad_mtx = get_raw_h5ad_file
	output:	
		h5ad_mtx = os.path.join(config["analysisDir"], "data/{sample}.h5ad"),
		gene_name_txt = os.path.join(config["analysisDir"], "data/{sample}.h5ad.all.genes.txt")
	params:
		time = "2:00:00",
		mem_gb = "64"
	shell:
		"bash -c ' source $HOME/.bashrc; \
		conda activate cnmf_env; \
		python workflow/scripts/filter_to_h5ad.py \
		--inputPath {input.raw_h5ad_mtx} \
		--output_h5ad {output.h5ad_mtx} \
		--output_gene_name_txt {output.gene_name_txt} ' "