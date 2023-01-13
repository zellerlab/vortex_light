#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { fastq_input; bam_input } from "./nevermore/workflows/input"
include { vlight_main } from "./vlight/workflows/vlight"


if (params.input_dir && params.remote_input_dir) {
	log.info """
		Cannot process both --input_dir and --remote_input_dir. Please check input parameters.
	""".stripIndent()
	exit 1
} else if (!params.input_dir && !params.remote_input_dir) {
	log.info """
		Neither --input_dir nor --remote_input_dir set.
	""".stripIndent()
	exit 1
}


params.do_bam2fq_conversion = true // vlight bam input requires this

def input_dir = (params.input_dir) ? params.input_dir : params.remote_input_dir


workflow {

	if (params.bam_input_pattern) {
		bam_input(
			Channel.fromPath(bam_input_pattern)
		)
		fastq_ch = bam_input.out.bamfiles
	} else {
		fastq_input(
			Channel.fromPath(input_dir + "/*", type: "dir")
		)
		fastq_ch = fastq_input.out.fastqs
	}

	vlight_main(fastq_ch)

}
