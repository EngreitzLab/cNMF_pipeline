
rule get_ABC_enhancer_fasta:
	input:
		fasta = config["fasta_file"],
		coord = config["ABC_enhancers"]
	output:
		fasta = os.path.join(config["analysisDir"], "{folder}/{sample}/fimo/fasta_to_fimo.fa")
	params:
		time = "3:00:00",
		mem_gb = "16",
		partition = "owners,normal"
	shell:
		"bash -c ' source ~/.bashrc; \
		conda activate cnmf_env; \
		bash workflow/scripts/fimo_motif_match.sh {input.coord} {input.fasta} {output.fasta} ' "


rule FIMO_ABC_enhancers:
	input:
		fasta = os.path.join(config["analysisDir"], "{folder}/{sample}/fimo/fasta_to_fimo.fa"),
		motif_meme = os.path.join(config["motif_meme"])
	output:
		fimo_result = os.path.join(config["analysisDir"], "{folder}/{sample}/fimo/fimo_out/fimo.txt")
	params:
		time = "12:00:00",
		mem_gb = "128",
		partition = "owners,normal",
		fimo_outdir = os.path.join(config["analysisDir"], "{folder}/{sample}/fimo/fimo_out")
	shell: ## need FIMO wrapper
		"bash -c ' source ~/.bashrc; \
		conda activate cnmf_env; \
		fimo -oc {params.fimo_outdir} --verbosity 1 --thresh 1.0E-4 {input.motif_meme} {input.fasta}' "

	# shell: ## need FIMO wrapper
	# 	"bash -c ' source ~/.bashrc; \
	# 	conda activate cnmf_env; \
	# 	fimo -oc {params.fimo_outdir} --verbosity 1 --thresh 1.0E-4 {input.motif_meme} {input.fasta}; \
	# 	echo {output.fimo_result} | tr \"|\" \"\\t\" > {output.fimo_formatted} ' "

