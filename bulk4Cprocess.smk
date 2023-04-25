####################################
#          4C_pipline              #
#@author Z Liu                     #
#@Ver 0.0.1                        #
#@date 2023/4/9                    #
####################################

import os

SAMPLES = [i.split(sep='_')[0] for i in os.listdir("./Rawdata")]

configfile: "bulk4Cprocess/config.yaml"


rule all:
    input:
    	expand("processed/pairs/{sample}.pairs.gz",sample=SAMPLES),
        expand("processed/cools/{sample}.5k.cool",sample=SAMPLES),


rule bwa_map:
    input:
        DNA1="Rawdata/{sample}/{sample}_R1.fq.gz",
        DNA2="Rawdata/{sample}/{sample}_R2.fq.gz",
        BWA_REF_GENOME=config["refs"][config["ref_genome"]]["BWA_REF_GENOME"],
    output:
        bamfile = "processed/mapping/{sample}.aln.bam"
    threads: 8
    resources:
        nodes = 8
    params:
        extra=r"-R '@RG\tID:{sample}\tPL:ILLUMINA\tSM:{sample}'",
    shell:"""
        set +u
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate py3
        set -u

        bwa-mem2 mem -5SP -t{threads} {params.extra} {input.BWA_REF_GENOME} {input.DNA1} {input.DNA2} | samtools sort -@{threads} -n -o {output.bamfile} -

        set +u
        conda deactivate
        set -u
        """

rule bam2pairs:
    input:
    	bamfile = rules.bwa_map.output.bamfile,
    	chromsizes = config["refs"][config["ref_genome"]]["chr_len"],
    output:
    	dedupPairs = "processed/pairs/{sample}.pairs.gz",
    threads: 10
    resources:
        nodes = 10
    shell:"""
    	set +u
        source activate
        conda activate hic
        set -u

        pairtools parse --chroms-path {input.chromsizes} {input.bamfile} | \
        pairtools sort --nproc 10 --tmpdir=./ | pairtools dedup | pairtools split --output-pairs {output.dedupPairs}

        set +u
        conda deactivate
        set -u
    """

rule pairs2cool:
    input:
        pairs = "processed/pairs/{sample}.pairs.gz",
        chr_len = config["refs"][config["ref_genome"]]["chr_len"],
    output:
        balancedCool = "processed/cools/{sample}.5k.cool",
    params:
        resolution = 5000,
    threads: 4
    shell:"""
        set +u
        source activate
        conda activate hic2
        set -u

        cooler cload pairs -c1 2 -c2 4 -p1 3 -p2 5 {input.chr_len}:{params.resolution} {input.pairs} {output.balancedCool}

        set +u
        conda deactivate
        set -u
    """

    
