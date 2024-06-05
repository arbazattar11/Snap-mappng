# SNP Mapping Pipeline

This repository contains a script for mapping single nucleotide polymorphisms (SNPs) from raw sequencing data using commonly used bioinformatics tools. The pipeline includes alignment, post-alignment processing, variant calling, and SNP filtering.

## Overview

The analysis pipeline performs the following steps:
1. Index the reference genome for alignment.
2. Align reads to the reference genome using BWA-MEM.
3. Convert SAM to BAM, sort, and index.
4. Mark duplicates using GATK.
5. Create a sequence dictionary and index the reference genome.
6. Call variants using GATK HaplotypeCaller.
7. Filter SNPs using bcftools.

## Requirements

### Software
- [BWA](http://bio-bwa.sourceforge.net/)
- [SAMtools](http://www.htslib.org/)
- [GATK](https://gatk.broadinstitute.org/hc/en-us)
- [bcftools](http://samtools.github.io/bcftools/)

### Data
- Raw FASTQ files.
- Reference genome in FASTA format.

## Setup

1. Install the required software tools.
2. Ensure the raw FASTQ files and reference genome are prepared.

## Usage

1. Clone this repository or download the script.

    ```bash
    git clone https://github.com/arbazattar11/Snap-mappng
    cd Snap-mappng
    ```

2. Modify the variables in the script (`snp_mapping.sh`) to match your data and paths:
    - `SAMPLE_NAME`: Name of your sample.
    - `FASTQ_R1`: Path to the forward reads FASTQ file.
    - `FASTQ_R2`: Path to the reverse reads FASTQ file.
    - `REF_GENOME`: Path to the reference genome.
    - `OUTPUT_DIR`: Directory where output files will be saved.
    - `BWA`: Path to the BWA executable.
    - `SAMTOOLS`: Path to the SAMtools executable.
    - `GATK`: Path to the GATK executable.
    - `BCFTOOLS`: Path to the bcftools executable.
    - `THREADS`: Number of threads to use for alignment.

3. Make the script executable:

    ```bash
    chmod +x snp_mapping.sh
    ```

4. Run the script:

    ```bash
    ./snp_mapping.sh
    ```

## Script Details

### snp_mapping.sh

This bash script performs the following steps:

1. **Index Reference Genome**: Prepares the reference genome for alignment with BWA.
2. **Align Reads**: Aligns paired-end reads to the reference genome using BWA-MEM.
3. **Post-Alignment Processing**: Converts SAM to BAM, sorts the BAM file, marks duplicates, and indexes the BAM file.
4. **Reference Genome Preparation**: Creates a sequence dictionary and indexes the reference genome.
5. **Variant Calling**: Calls variants using GATK HaplotypeCaller and generates a GVCF file.
6. **SNP Filtering**: Filters SNPs based on quality and depth using bcftools.

## Output

The script generates several output files and directories:
- `*.sam`: SAM file from BWA-MEM alignment.
- `*.sorted.bam`: Sorted BAM file.
- `*.dedup.bam`: Deduplicated BAM file.
- `*.g.vcf`: GVCF file with called variants.
- `*.filtered.snps.vcf`: Filtered VCF file with SNPs.

## Troubleshooting

If you encounter any issues, ensure that:
- All paths and filenames are correct.
- All required software tools are installed and accessible.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

- [BWA](http://bio-bwa.sourceforge.net/) for read alignment.
- [SAMtools](http://www.htslib.org/) for BAM processing.
- [GATK](https://gatk.broadinstitute.org/hc/en-us) for variant calling.
- [bcftools](http://samtools.github.io/bcftools/) for SNP filtering.

## Contact

For any questions or issues, please contact [Arbaz] at [arbazattar1137@gmail.com].
