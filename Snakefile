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
    shell:
        """
        mkdir -p resources
        wget -O {output.ref} {refGenomePathOrUrl}
        """

rule unzip_reference:
    input:
        ref=refGenome
    output:
        fasta=genomeFasta
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
    params:
        docker=config.get("docker", True),
        sjdbOverhang=100
    shell:
        """
        {{"docker run --rm -v $(pwd):/data quay.io/biocontainers/star:2.6.1d--0 STAR" if params.docker else "STAR"}} \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output.index_dir} \
            --genomeFastaFiles {input.fasta} \
            --sjdbOverhang {params.sjdbOverhang}
        """

rule star_align:
    input:
        fastq1=lambda wildcards: SAMPLES[wildcards.sample]["fastq1"],
        fastq2=lambda wildcards: SAMPLES[wildcards.sample]["fastq2"],
        genome_dir=starIndexDir
    output:
        bam="results/aligned/{sample}.bam",
        log="logs/{sample}_STAR.log"
    params:
        extra=config.get("star_extra", "--outSAMtype BAM SortedByCoordinate"),
        docker=config.get("docker", True)
    shell:
        """
        mkdir -p results/aligned logs
        {{"docker run --rm -v $(pwd):/data quay.io/biocontainers/star:2.6.1d--0 STAR" if params.docker else "STAR"}} \
            --genomeDir {input.genome_dir} \
            --readFilesIn {input.fastq1} {input.fastq2} \
            --readFilesCommand zcat \
            {params.extra} \
            --outFileNamePrefix results/aligned/{wildcards.sample}_ \
            > {output.log} 2>&1
        mv results/aligned/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output.bam}
        """
