configfile: "config/config.yaml"

POP1_NAME = config["populations"]["pop1_name"]
POP2_NAME = config["populations"]["pop2_name"]

POP1_FILE = config["populations"]["pop1_file"]
POP2_FILE = config["populations"]["pop2_file"]

REGION = config["region"]

REF_GENOME = config["ref_genome"]

#################################
# Parameters
#################################

MIN_MAP_Q = config["filters"]["minMapQ"]
MIN_Q = config["filters"]["minQ"]
MIN_IND = config["filters"]["minInd"]
UNIQUE_ONLY = config["filters"]["uniqueOnly"]
ONLY_PROPER_PAIRS = config["filters"]["only_proper_pairs"]
MIN_DEPTH = config["filters"]["minDepth"]
MAX_DEPTH = config["filters"]["maxDepth"]

FILTER_TAG = (
    f"mapq{MIN_MAP_Q}"
    f"_q{MIN_Q}"
    f"_minInd{MIN_IND}"
    f"_uniq{UNIQUE_ONLY}"
    f"_proper{ONLY_PROPER_PAIRS}"
    f"_minD{MIN_DEPTH}"
    f"_maxD{MAX_DEPTH}"
)

DEPTH_FLAGS = ""

if MIN_DEPTH > 0:
    DEPTH_FLAGS += f" -setMinDepth {MIN_DEPTH}"

if MAX_DEPTH > 0:
    DEPTH_FLAGS += f" -setMaxDepth {MAX_DEPTH}"

WIN  = config["window"]["size"]
STEP = config["window"]["step"]


#################################
# Final target
#################################

rule all:
    input:
        f"results/{FILTER_TAG}/fst_{POP1_NAME}_{POP2_NAME}_{REGION}.pdf",
        f"results/{FILTER_TAG}/fst_global_{POP1_NAME}_{POP2_NAME}_{REGION}",
        f"results/{FILTER_TAG}/fst_{POP1_NAME}_{POP2_NAME}_{REGION}_win_{WIN//1000}kb_step_{STEP//1000}kb.sliding_window"

#################################
# Final plot
#################################

rule plot_sliding_fst:
    input:
        f"results/{FILTER_TAG}/fst_global_{POP1_NAME}_{POP2_NAME}_{REGION}",
        f"results/{FILTER_TAG}/fst_{POP1_NAME}_{POP2_NAME}_{REGION}_win_{WIN//1000}kb_step_{STEP//1000}kb.sliding_window"
    output:
        f"results/{FILTER_TAG}/fst_{POP1_NAME}_{POP2_NAME}_{REGION}.pdf"
    conda:
        "envs/r.yaml"
    shell:
        """
        Rscript scripts/fst_sex_diff_chrom6.r \
            {input[0]} \
            {input[1]} \
            {output} \
            {POP1_NAME} \
            {POP2_NAME} \
            {REGION} \
            {WIN} \
            {STEP}    
        """

#################################
# Sliding window FST
#################################

rule fst_sliding:
    input:
        idx = f"intermediate/{FILTER_TAG}/{POP1_NAME}_{POP2_NAME}_{REGION}.fst.idx",
        fst_file = f"intermediate/{FILTER_TAG}/{POP1_NAME}_{POP2_NAME}_{REGION}.fst.gz"
    output:
        f"results/{FILTER_TAG}/fst_{POP1_NAME}_{POP2_NAME}_{REGION}_win_{WIN//1000}kb_step_{STEP//1000}kb.sliding_window"
    threads:
        config["threads"]
    conda:
        "envs/angsd.yaml"
    shell:
        """
            realSFS fst stats2 \
            {input.idx} \
            -fold 1 \
            -win {WIN} \
            -step {STEP} \
            -P {threads} \
            > {output}
        """

#################################
# Global FST
#################################        

rule global_fst:
    input:
        idx = f"intermediate/{FILTER_TAG}/{POP1_NAME}_{POP2_NAME}_{REGION}.fst.idx",
        fst_file = f"intermediate/{FILTER_TAG}/{POP1_NAME}_{POP2_NAME}_{REGION}.fst.gz"
    output:
        f"results/{FILTER_TAG}/fst_global_{POP1_NAME}_{POP2_NAME}_{REGION}"
    threads:
        config["threads"]
    conda:
        "envs/angsd.yaml"
    shell:
        """
            realSFS fst stats \
            {input.idx} \
            -fold 1 \
            -P {threads} \
            > {output}
        """


#########################################################
# Intermediate per-site FST components in indexed format
#########################################################

rule fst_index:
    input:
        pop1_idx=f"intermediate/{FILTER_TAG}/{POP1_NAME}_{REGION}.saf.idx",
        pop1_saf=f"intermediate/{FILTER_TAG}/{POP1_NAME}_{REGION}.saf.gz",
        pop1_pos=f"intermediate/{FILTER_TAG}/{POP1_NAME}_{REGION}.saf.pos.gz",
        pop2_idx=f"intermediate/{FILTER_TAG}/{POP2_NAME}_{REGION}.saf.idx",
        pop2_saf=f"intermediate/{FILTER_TAG}/{POP2_NAME}_{REGION}.saf.gz",
        pop2_pos=f"intermediate/{FILTER_TAG}/{POP2_NAME}_{REGION}.saf.pos.gz",
        sfs_file=f"intermediate/{FILTER_TAG}/{POP1_NAME}_{POP2_NAME}_{REGION}.ml"
    output:
        idx = f"intermediate/{FILTER_TAG}/{POP1_NAME}_{POP2_NAME}_{REGION}.fst.idx",
        fst_file = f"intermediate/{FILTER_TAG}/{POP1_NAME}_{POP2_NAME}_{REGION}.fst.gz"
    threads:
        config["threads"]
    conda:
        "envs/angsd.yaml"
    shell:
        """
            realSFS fst index \
            {input.pop1_idx} \
            {input.pop2_idx} \
            -sfs {input.sfs_file} \
            -fold 1 \
            -P {threads} \
            -fstout intermediate/{FILTER_TAG}/{POP1_NAME}_{POP2_NAME}_{REGION}
        """

#########################################################
# 2D SFS estimation
#########################################################

rule sfs_2d:
    input:
        pop1_idx=f"intermediate/{FILTER_TAG}/{POP1_NAME}_{REGION}.saf.idx",
        pop1_saf=f"intermediate/{FILTER_TAG}/{POP1_NAME}_{REGION}.saf.gz",
        pop1_pos=f"intermediate/{FILTER_TAG}/{POP1_NAME}_{REGION}.saf.pos.gz",
        pop2_idx=f"intermediate/{FILTER_TAG}/{POP2_NAME}_{REGION}.saf.idx",
        pop2_saf=f"intermediate/{FILTER_TAG}/{POP2_NAME}_{REGION}.saf.gz",
        pop2_pos=f"intermediate/{FILTER_TAG}/{POP2_NAME}_{REGION}.saf.pos.gz"
    output:
        sfs_file=f"intermediate/{FILTER_TAG}/{POP1_NAME}_{POP2_NAME}_{REGION}.ml"
    threads:
        config["threads"]
    conda:
        "envs/angsd.yaml"
    shell:
        """
            realSFS \
            {input.pop1_idx} \
            {input.pop2_idx} \
            -fold 1 \
            -P {threads} \
            > {output.sfs_file}
        """

#########################################################
# Site allele freqs for pop 2
#########################################################

rule safs_pop2:
    input:
        ref = REF_GENOME,
        pop2_samples = POP2_FILE
    output:
        pop2_idx=f"intermediate/{FILTER_TAG}/{POP2_NAME}_{REGION}.saf.idx",
        pop2_saf=f"intermediate/{FILTER_TAG}/{POP2_NAME}_{REGION}.saf.gz",
        pop2_pos=f"intermediate/{FILTER_TAG}/{POP2_NAME}_{REGION}.saf.pos.gz"
    log:
        f"log/{FILTER_TAG}/safs_pop2_{POP2_NAME}_{REGION}.log"
    threads:
        config["threads"]
    conda:
        "envs/angsd.yaml"
    shell:
        """
            angsd -b {input.pop2_samples} -r {REGION} -doSaf 1 \
            -out intermediate/{FILTER_TAG}/{POP2_NAME}_{REGION} \
            -anc {input.ref} \
            -GL 2 \
            -minMapQ {MIN_MAP_Q} \
            -minQ {MIN_Q} \
            -minInd {MIN_IND} \
            -uniqueOnly {UNIQUE_ONLY} \
            -only_proper_pairs {ONLY_PROPER_PAIRS} \
            {DEPTH_FLAGS} \
            -P {threads} \
            2> {log}
        """

#########################################################
# Site allele freqs for pop 1
#########################################################

rule safs_pop1:
    input:
        ref = REF_GENOME,
        pop1_samples = POP1_FILE
    output:
        pop1_idx=f"intermediate/{FILTER_TAG}/{POP1_NAME}_{REGION}.saf.idx",
        pop1_saf=f"intermediate/{FILTER_TAG}/{POP1_NAME}_{REGION}.saf.gz",
        pop1_pos=f"intermediate/{FILTER_TAG}/{POP1_NAME}_{REGION}.saf.pos.gz"
    log:
        f"log/{FILTER_TAG}/safs_pop1_{POP1_NAME}_{REGION}.log"
    threads:
        config["threads"]
    conda:
        "envs/angsd.yaml"
    shell:
        """
            angsd -b {input.pop1_samples} -r {REGION} -doSaf 1 \
            -out intermediate/{FILTER_TAG}/{POP1_NAME}_{REGION} \
            -anc {input.ref} \
            -GL 2 \
            -minMapQ {MIN_MAP_Q} \
            -minQ {MIN_Q} \
            -minInd {MIN_IND} \
            -uniqueOnly {UNIQUE_ONLY} \
            -only_proper_pairs {ONLY_PROPER_PAIRS} \
            {DEPTH_FLAGS} \
            -P {threads} \
            2> {log}
        """










