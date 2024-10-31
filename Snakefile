#STAR and snakemake need to be installed
#navigate to Snakefile location and start pipeline with 
# "snakemake --cores 12" 
#to align reads of two fastq samples to hg38 using STAR aligner


import pandas as pd

# Read the sample sheet
configfile: "config.yaml"
samples = pd.read_csv(config["sample_sheet"])
print(samples)
# Define input files for each sample
SAMPLES = {
    row["sampleID"]: {
        "fastq1": row["fastq1"],
        "fastq2": row["fastq2"]
    }
    for _, row in samples.iterrows()
}


rule all:
    input:
        # Collect all BAM files and the final MultiQC report
        expand("results/aligned/{sample}.bam", sample=SAMPLES.keys())

rule star_align:
    input:
        fastq1=lambda wildcards: SAMPLES[wildcards.sample]["fastq1"],
        fastq2=lambda wildcards: SAMPLES[wildcards.sample]["fastq2"],
        genome_dir=config["star_index"]
    output:
        bam="results/aligned/{sample}.bam",
        log="logs/{sample}_STAR.log"
    params:
        # Adjust parameters as needed
        extra=config.get("star_extra", "--outSAMtype BAM SortedByCoordinate")
    shell:
        """
        STAR --genomeDir {input.genome_dir} \
             --readFilesIn {input.fastq1} {input.fastq2} \
             --readFilesCommand zcat \
             {params.extra} \
             --outFileNamePrefix results/aligned/{wildcards.sample}_ \
             > {output.log} 2>&1
        mv results/aligned/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output.bam}
        """

