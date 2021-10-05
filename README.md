# Pipeline for making 3'UTR predictions in mouse using F3UTER

This `snakemake` pipeline uses the F3UTER trained model (without DNA structural properties) to classify input regions from the mouse genome into potential 3'UTRs and non-3'UTRs. F3UTER was trained on a set of known 3’UTRs (positive examples) and non-3’UTRs (negative examples) from Ensembl human genome annotation (v94), and features derived from the human data. F3UTER preforms remarkably well in classifying 3’UTRs in M.musculus as well, achieving an AUROC of 0.96 and AUPRC of 0.8 in classifying known 3'UTRs in the mouse genome. The trained F3UTER model is provided in `/data/rf_model_no_dsp.rda`. This pipeline generates omic features for a set of input candidate regions, compiles the different omic features into a feature matrix and applies F3UTER to make predictions.

# Getting Started

## Input

- Regions of interest: can be a standard BED file, it must contain at least the following columns - `seqnames`, `start`, `end`, `strand` and `id`. Please note that contig names in BED file should be UCSC style: `chr1`, `chr2` .... `chrM`.
- Aligned RNA-seq reads in bigwig format to calculate the expression features. Multiple RNA-seq replicates can be provided. In case of multiple bigwigs, mean expression will be calculated. Please note that contig names in the bigwig should be UCSC style: `chr1`, `chr2` .... `chrM`.
- chromosome lengths are required by the code. A file containing chromosome lengths for mm10 is provided in `/data` and is automatically used by the pipeline.

External data resources:
- PhastCons scores from UCSC (bigwig format): https://hgdownload.soe.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons60wayEuarchontoGlire.bw
- Repeats data from repeatmasker.org: https://www.repeatmasker.org/genomes/mm10/RepeatMasker-rm406-dfam2.0/mm10.fa.out.gz

Process the repeatmasker file using the commands below. Please see `/test_data` for a demo file.
```bash
gunzip mm10.fa.out.gz
awk -F " " '{print $5" "$6" "$7" "$11}' mm10.fa.out | sed 1,3d > mm10.repeatMasker.mod.fa.out
```

## Output

The main output files generated by the pipeline are:

- `F3UTER_features/<sample_name>_nt_freq.txt` - mono- and di-nucleotide frequencies (n=20)
- `F3UTER_features/<sample_name>_phastcons.txt` - sequence conservation (phastCons) scores (n=1)
- `F3UTER_features/<sample_name>_polyA_signal.txt` - poly(A) signal occurrence (n=1; binary outcome)
- `F3UTER_features/<sample_name>_repeats.txt` - transposon occurrence (n=1)
- `F3UTER_features/<sample_name>_exp_feat.txt` - expression features (n=2)
- `F3UTER_prediction/<sample_name>_ml_table.rds` - compiled feature matrix for F3UTER (R object)
- `F3UTER_prediction/<sample_name>_predictions_meta.txt` - predictions from F3UTER and meta-data added from BED file


## Depedencies

- [miniconda](https://conda.io/miniconda.html)
- [snakemake](http://snakemake.readthedocs.io/en/latest/) - can be installed via conda (`snakemake>=5.3`)
- The rest of the dependencies (R packages) are installed via conda.

## Installation

Clone the pipeline:

```bash
git clone --recursive https://github.com/sid-sethi/Generate-F3UTER-predictions-mouse.git
```

## Usage

Edit `config.yml` to set up the working directory and input files. `Snakemake` command should be issued from within the pipeline directory.

```bash
cd Generate-F3UTER-predictions-mouse
snakemake --use-conda -j <num_cores> all
```
If you provide more than one core, independent snakemake rules will be processed simultaneously. This pipeline uses 5 cores at most. It is a good idea to do a dry run (using -n parameter) to view what would be done by the pipeline before executing the pipeline.

```bash
snakemake --use-conda -n all
```
Snakemake can be run to only install the required conda environments without running the full workflow. Subsequent runs with --use-conda will make use of the local environments without requiring internet access. This is suitable for running the pipeline offline.

```bash
snakemake --use-conda --conda-create-envs-only
```

## Licence

Copyright 2020 Astex Therapeutics Ltd.

This repository is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the [LICENSE](LICENSE) file (GNU General Public License) for more details.
