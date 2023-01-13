process collate_results {
	publishDir params.output_dir, mode: params.publish_mode

	input:
	path(results)
	path(collate_script)
	path(gtdb_markers)

	output:
	path("collated/*.rds"), emit: collated, optional: true

	script:
	"""
	mkdir -p collated/

	mkdir -p kraken2/
	(mv *kraken2_report.txt kraken2/) || :

	mkdir -p pathseq/
	(mv *pathseq.scores pathseq/) || :

	mkdir -p libsize/
	(mv *libsize.txt libsize/) || :

	mkdir -p liblayout/
	(mv *is_paired.txt liblayout/) || :

	mkdir -p flagstats/
	(mv *flagstats.txt flagstats/) || :

	mkdir -p raw_counts/
	(mv *.txt raw_counts/) || :

	Rscript ${collate_script} \
		--libdir \$(dirname \$(readlink ${collate_script})) \
		--gtdb_markers ${gtdb_markers} \
		--kraken2_res_path kraken2/ \
		--PathSeq_res_path pathseq/ \
		--libsize_res_path libsize/ \
		--lib_layout_res_path liblayout/ \
		--N_raw_counts_path raw_counts/ \
		--out_folder collated/
	"""
}