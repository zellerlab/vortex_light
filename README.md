# vortex_knight

### Installing locally and running from local installation

1. Clone the repo from GitHub.

```
git clone https://github.com/zellerlab/vortex_light.git
```

2. Create a conda environment with NextFlow, e.g. by using the provided `environment.yml`.

```
cd vortex_light
conda env create -f environment.yml
conda activate vortex_light
```

3. Make a copy of the `config/run.config` file and adjust it to your environment.

4. Run the pipeline 

``` 
nextflow run /path/to/vortex_light/main.nf --input_dir /path/to/input_files --output_dir /path/to/output_dir -c /path/to/run.config
```

*Note: Nextflow itself requires at least `5GB` of memory.*

### Running from GitHub

This requires a local nextflow installation. If you don't have one, see Steps 1/2 above.

1. Make a local copy of the [run configuration file](https://raw.githubusercontent.com/zellerlab/vortex_light/main/nextflow/run.config) and adjust to your environment.

2. Run the pipeline

```
nextflow run zellerlab/vortex_light --input_dir /path/to/input_files --output_dir /path/to/output_dir -c /path/to/run.config
```

*Note: Nextflow itself requires at least `5GB` of memory.*

### Input parameters

* `--input_dir` should be a folder with bam files or with gzipped fastq files. For fastq files, individual samples should be separated into individual folders.
* `--output_dir` is `vknight_out` in the local directory by default.
* `--skip_<analysis>`, `--run_<analysis>` skips, resp. explicitly requires execution of the specified analysis (`pathseq`, `base_counts` (read counts post pre-processing), `kraken2`)
* `--publishMode` allows to switch between various modes of how results files are placed in the `output_dir` (cf. NextFlow documentation)

#### Notes
* `kraken2` can only run when the parameter `kraken_database` is set.
* `pathseq` can only run when the parameter `pathseq_database` is set.


5. Outputs

The output folder contains:

* a subdirectory per sample (named `<sample>`) with
  * the kraken2 report `<sample>.kraken2_report.txt`
  * the library size `<sample>.libsize.txt`
  * pathseq output
    - `<sample>.pathseq.bam`
    - `<sample>.pathseq.bam.sgi`
    - `<sample>.pathseq.score_metrics`
    - `<sample>.pathseq.scores`

**Note that by default, all files in the output folder are symlinks into the work dir! Before you delete the work dir, ensure you have dereferenced copies. Alternatively, change the --publishMode parameter to `copy` or `link` (if the target file system supports hard links).**
