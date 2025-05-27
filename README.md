# macworp-genomics-demo-workflow-snakemake
## Workflow
Tasks included in this workflow:
1) Align reads to a reference genome
    * Input: reference genome index, fastq files, sample sheet
    * Tool: STAR (needs to be installed)
    * Output: Aligned reads from fastq files as .bam files
## How to run
Quick start: 
1) Navigate to folder containing "Snakefile"
2) Adjust parameters in config.yaml
3) Start with `snakemake --sdm apptainer --cores all` (using default parameters and example input files) 
    * Attention: snakemake and STAR need to be installed and a reference genome index built with STAR

