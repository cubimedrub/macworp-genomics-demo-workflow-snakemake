import pandas as pd
import os

# Load config and sample sheet
configfile: "config.yaml"
samples = pd.read_csv(config["sample_sheet"])
print(samples)

# Reference genome URL from config (with default fallback)
refGenomePathOrUrl = config.get(
    "refGenomePathOrUrl",
    "https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
)

# File paths
refGenome = "resources/genome.fa.gz"
genomeFasta = "resources/genome.fa"
starIndexDir = "resources/star_index"

# Prepare sample input mapping
SAMPLES = {
    row["sampleID"]: {
        "fastq1": row["fastq1"],
        "fastq2": row["fastq2"]
    }
    for _, row in samples.iterrows()
}

rule all:
    input:
        expand("results/aligned/{sample}.bam", sample=SAMPLES.keys())

rule download_reference:
    output:
        ref=refGenome
    container:
        "docker://quay.io/medbioinf/base-tools:1.0.0"
    shell:
        """
        mkdir -p resources
        curl -v -k -L -o {output.ref} {refGenomePathOrUrl}
        """

rule unzip_reference:
    input:
        ref=refGenome
    output:
        fasta=genomeFasta
    container:
        "docker://quay.io/medbioinf/base-tools:1.0.0"
    shell:
        """
        gunzip -c {input.ref} > {output.fasta}
        """

rule generate_star_index:
    input:
        fasta=genomeFasta
    output:
        index_dir=directory(starIndexDir)
    threads: 8
    container:
        "docker://quay.io/biocontainers/star:2.6.1d--0"
    shell:
        """
        mkdir -p {output.index_dir}
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output.index_dir} \
            --genomeFastaFiles {input.fasta}
        """

rule star_align:
    input:
        fastq1=lambda wildcards: SAMPLES[wildcards.sample]["fastq1"],
        fastq2=lambda wildcards: SAMPLES[wildcards.sample]["fastq2"],
        genome_dir=starIndexDir
    output:
        bam="results/aligned/{sample}.bam",
        log="logs/{sample}_STAR.log"
    container:
        "docker://quay.io/biocontainers/star:2.6.1d--0"
    params:
        extra=config.get("star_extra", "--outSAMtype BAM SortedByCoordinate"),
    shell:
        """
        mkdir -p results/aligned logs
        STAR \
            --genomeDir {input.genome_dir} \
            --readFilesIn {input.fastq1} {input.fastq2} \
            --readFilesCommand zcat \
            {params.extra} \
            --outFileNamePrefix results/aligned/{wildcards.sample}_ \
            > {output.log} 2>&1
        mv results/aligned/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output.bam}
        """
