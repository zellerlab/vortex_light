#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { kraken2 } from  "../modules/profilers/kraken2"
include { pathseq } from "../modules/profilers/pathseq"

include { fq2fa } from "../../nevermore/modules/converters/fq2fa"
include { fastqc } from "../../nevermore/modules/qc/fastqc"
include { multiqc } from "../../nevermore/modules/qc/multiqc"
include { fq2bam } from "../../nevermore/modules/converters/fq2bam"
include { nevermore_simple_preprocessing } from "../../nevermore/workflows/nevermore"
include { remove_host_kraken2; remove_host_kraken2_individual } from "../../nevermore/modules/decon/kraken2"
include { flagstats } from "../../nevermore/modules/stats"
include { collate_results } from "../modules/collate"


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

def asset_dir = "${projectDir}/nevermore/assets"

def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)
def get_basecounts = (!params.skip_basecounts || params.run_basecounts);

def run_bam_analysis = run_pathseq
def run_fastq_analysis = run_kraken2



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


workflow vlight_main {
	take:
		fastq_ch
	main:
		results_ch = Channel.empty()

		if (do_preprocessing) {

			nevermore_simple_preprocessing(fastq_ch)

			preprocessed_ch = nevermore_simple_preprocessing.out.main_reads_out
			results_ch = results_ch
				.concat(nevermore_simple_preprocessing.out.raw_counts)
				.map { sample, files -> files }

			if (params.remove_host) {

				if (params.remove_host == "individual") {

					remove_host_kraken2_individual(nevermore_simple_preprocessing.out.main_reads_out, params.remove_host_kraken2_db)

					preprocessed_ch = remove_host_kraken2_individual.out.reads


				} else if (params.remove_host == "pair") {

					remove_host_kraken2(nevermore_simple_preprocessing.out.main_reads_out, params.remove_host_kraken2_db)

					preprocessed_ch = remove_host_kraken2.out.reads

				}

			}

		} else {

			preprocessed_ch = fastq_ch

		}

		//
		/*	perform post-qc fastqc analysis and generate multiqc report on merged single-read and paired-end sets */

		fastqc(preprocessed_ch, "qc")

		multiqc(
			fastqc.out.stats
				.map { sample, report -> return report }.collect(),
			"${asset_dir}/multiqc.config",
			"qc"
		)

		//

		if (get_basecounts || run_bam_analysis) {

			fq2bam(preprocessed_ch)

			if (get_basecounts) {

				flagstats(fq2bam.out.reads)

				flagstat_results_ch = flagstats.out.flagstats
					.concat(flagstats.out.counts)
					.concat(flagstats.out.is_paired)
					.map { sample, files -> files }
				results_ch = results_ch.concat(flagstat_results_ch)

			}

			if (run_bam_analysis) {

				bam_analysis(fq2bam.out.reads)
				results_ch = results_ch.concat(bam_analysis.out.results)

			}

		}

		if (run_fastq_analysis) {

			fastq_analysis(preprocessed_ch)
			results_ch = results_ch.concat(fastq_analysis.out.results)

		}

		if (!params.skip_collate) {
			collate_results(
				results_ch.collect(),
				"${projectDir}/scripts/ExtractProfiledCounts_210823.R",
				params.GTDB_markers
			)
		}

	emit:
		results = results_ch
}
