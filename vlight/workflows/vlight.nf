#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { kraken2 } from  "../modules/profilers/kraken2"
include { pathseq } from "../modules/profilers/pathseq"
include { fq2fa } from "../../nevermore/modules/converters/fq2fa"


if (!params.publish_mode) {
	params.publish_mode = "symlink"
}

if (!params.output_dir) {
	params.output_dir = "vlight_out"
}

output_dir = "${params.output_dir}"

if (!params.pathseq_min_clipped_read_length) {
	params.pathseq_min_clipped_read_length = 31
}

def run_kraken2 = (!params.skip_kraken2 || params.run_kraken2);
def run_pathseq = (!params.skip_pathseq || params.run_pathseq);


workflow bam_analysis {
	take:
		bam_ch

	main:
		out_ch = Channel.empty()
    	if (run_pathseq) {
			pathseq(bam_ch, params.pathseq_database)
			out_ch = out_ch.concat(pathseq.out.scores)
    	}

		out_ch = out_ch
			.map { sample, files -> return files }

	emit:
		results = out_ch
}


workflow fastq_analysis {
	take:
		fastq_ch

	main:
		out_ch = Channel.empty()

		if (run_kraken2) {
			kraken2(fastq_ch, params.kraken_database)
			out_ch = out_ch.concat(kraken2.out.kraken2_out)
		}

		out_ch = out_ch
			.map { sample, files -> return files }

	emit:
		results = out_ch
}
