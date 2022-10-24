#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { bam2fq } from "./nevermore/modules/converters/bam2fq"
include { fq2bam } from "./nevermore/modules/converters/fq2bam"
include { prepare_fastqs } from "./nevermore/modules/converters/prepare_fastqs"
include { nevermore_simple_preprocessing } from "./nevermore/workflows/nevermore"
include { classify_sample } from "./nevermore/modules/functions"
include { remove_host_kraken2; remove_host_kraken2_individual } from "./nevermore/modules/decon/kraken2"
include { flagstats } from "./nevermore/modules/stats"
include { bam_analysis; fastq_analysis } from "./vlight/workflows/vlight"


def run_kraken2 = false
def run_pathseq = (!params.skip_pathseq || params.run_pathseq)

def get_basecounts = (!params.skip_basecounts || params.run_basecounts);
def convert_fastq2bam = (run_pathseq || get_basecounts);

def do_preprocessing = (!params.skip_preprocessing || params.run_preprocessing)

def run_bam_analysis = run_pathseq
def run_fastq_analysis = run_kraken2


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



workflow {

	fastq_ch = Channel
		//.fromPath(params.input_dir + "/" + "**.{fastq,fq,fastq.gz,fq.gz}")
		.fromPath(params.input_dir + "/" + "**[._]{fastq.gz,fq.gz}")
		.map { file ->
				def sample = file.name.replaceAll(/.?(fastq|fq)(.gz)?$/, "")
				sample = sample.replaceAll(/_R?[12]$/, "")
				return tuple(sample, file)
		}
		.groupTuple(sort: true)
        .map { classify_sample(it[0], it[1]) }

	bam_ch = Channel
		.fromPath(params.input_dir + "/" + "**.bam")
		.map { file ->
			def sample = file.name.replaceAll(/.bam$/, "")
			return tuple(sample, file)
		}
		.groupTuple(sort: true)
		.map { classify_sample(it[0], it[1]) }

	bam2fq(bam_ch)

	bfastq_ch = bam2fq.out.reads
		.map { classify_sample(it[0].id, it[1]) }

	prepare_fastqs(fastq_ch)

	results_ch = Channel.empty()

	if (do_preprocessing) {

		raw_fastq_ch = prepare_fastqs.out.reads.concat(bfastq_ch)

		nevermore_simple_preprocessing(raw_fastq_ch)

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

		preprocessed_ch = prepare_fastqs.out.reads
			.concat(bfastq_ch)

	}


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

}
