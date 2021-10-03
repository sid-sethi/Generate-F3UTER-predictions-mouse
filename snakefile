import os
from os import path

from snakemake.utils import min_version
min_version("5.3")

# ----------------------------------------------------------------

configfile: "config.yml"
workdir: config["workdir"]

WORKDIR = config["workdir"]
SNAKEDIR = path.dirname(workflow.snakefile)

er_file = config["ers"]
sample = config["sample_name"]

chr = list(range(1,19))
chr.append("X")
chr.append("Y")
chr.append("M")

# ----------------------------------------------------------------

rule all:
    input:
        "F3UTER_features/" + sample + "_nt_freq.txt",
        "F3UTER_features/" + sample + "_polyA_signal.txt",
        "F3UTER_features/" + sample + "_repeats.txt",
        "F3UTER_features/" + sample + "_phastcons.txt",
        "F3UTER_features/" + sample + "_exp_feat.txt",
        "F3UTER_prediction/" + sample + "_ml_table.rds",
        "F3UTER_prediction/" + sample + "_predictions_meta.txt"



#########################################################################


rule cal_nt_freq:
    input:
        ers = er_file

    output:
        out = "F3UTER_features/" + sample + "_nt_freq.txt"

    params:
        resPath = "F3UTER_features",
        prefix = sample,
        script = SNAKEDIR + "/scripts/nt_freq.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.ers} {params.resPath} {params.prefix}
        """



rule cal_pas:
    input:
        ers = er_file

    output:
        out = "F3UTER_features/" + sample + "_polyA_signal.txt"

    params:
        resPath = "F3UTER_features",
        prefix = sample,
        script = SNAKEDIR + "/scripts/polyA_signal.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.ers} {params.resPath} {params.prefix}
        """



rule cal_repeats:
    input:
        ers = er_file,
        repeats_data = config["repeats"]

    output:
        out = "F3UTER_features/" + sample + "_repeats.txt"

    params:
        resPath = "F3UTER_features",
        prefix = sample,
        script = SNAKEDIR + "/scripts/repeats.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.ers} {params.resPath} {params.prefix} {input.repeats_data}
        """



rule cal_phastcons:
    input:
        ers = er_file,
        bw = config["phastcons"]

    output:
        out = "F3UTER_features/" + sample + "_phastcons.txt"

    params:
        resPath = "F3UTER_features",
        prefix = sample,
        script = SNAKEDIR + "/scripts/phastCons.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.ers} {params.resPath} {params.prefix} {input.bw}
        """




rule generate_coverage:
    input:
        chrom_info = SNAKEDIR + "/data/mm10.chrom.sizes"

    output:
        out = expand("Coverage_files/" + sample + "_chr{seqnames}_mean_cov.rda", seqnames=chr)

    params:
        bw_dir = config["bigwig_dir"],
        resPath = "Coverage_files",
        prefix = sample,
        script = SNAKEDIR + "/scripts/generate_coverage.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {params.bw_dir} {params.resPath} {params.prefix} {input.chrom_info}
        """



rule cal_exp_feat:
    input:
        ers = er_file,
        coverage_files = rules.generate_coverage.output.out

    output:
        out = "F3UTER_features/" + sample + "_exp_feat.txt"

    params:
        resPath = "F3UTER_features",
        prefix = sample,
        script = SNAKEDIR + "/scripts/exp_features.R",
        coverage_dir = "Coverage_files"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.ers} {params.resPath} {params.prefix} {params.coverage_dir}
        """




rule make_ml_table:
    input:
        nt = rules.cal_nt_freq.output.out,
        pas = rules.cal_pas.output.out,
        cons = rules.cal_phastcons.output.out,
        repeats = rules.cal_repeats.output.out,
        exp = rules.cal_exp_feat.output.out

    output:
        ml_table = "F3UTER_prediction/" + sample + "_ml_table.rds"

    params:
        feat_dir = "F3UTER_features",
        resPath = "F3UTER_prediction",
        prefix = sample,
        script = SNAKEDIR + "/scripts/make_ML_table.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {params.feat_dir} {params.resPath} {params.prefix}
        """




rule make_f3uter_predictions:
    input:
        ml_table = rules.make_ml_table.output.ml_table,
        prediction_model = SNAKEDIR + "/data/rf_model_no_dsp.rda"

    output:
        preds = "F3UTER_prediction/" + sample + "_predictions.txt"

    params:
        resPath = "F3UTER_prediction",
        prefix = sample,
        script = SNAKEDIR + "/scripts/f3uter_predict.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.ml_table} {params.resPath} {params.prefix} {input.prediction_model}
        """



rule add_metadata_predictions:
    input:
        preds = rules.make_f3uter_predictions.output.preds,
        ers = er_file

    output:
        meta = "F3UTER_prediction/" + sample + "_predictions_meta.txt",
        bed = "F3UTER_prediction/" + sample + "_predictions.bed"

    params:
        resPath = "F3UTER_prediction",
        prefix = sample,
        script = SNAKEDIR + "/scripts/add_metadata_predictions.R"

    conda: "envs/r_v3.yml"

    shell:
        """
        Rscript {params.script} {input.preds} {input.ers} {params.resPath} {params.prefix}
        """
