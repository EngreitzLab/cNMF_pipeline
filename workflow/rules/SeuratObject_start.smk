## convert Seurat Object to h5ad file
rule Seurat_Object_to_h5ad:
	input:
		seurat_object = os.path.join(config["analysisDir"], "data/" + config["sampleName"] + ".SeuratObject.RDS")
	output:
		h5ad_mtx = os.path.join(config["analysisDir"], "data/" + config["sampleName"] + ".h5ad"),
		gene_name_txt = os.path.join(config["analysisDir"], "data/" + config["sampleName"] + ".h5ad.all.genes.txt")
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