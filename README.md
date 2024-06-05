# Detection and analysis of SDs and CNVs 

This repository contains the Snakemake pipeline and scripts used for the analysis of segmental duplications (SDs) and copy number variations (CNVs) in the European amphioxus. The purpose of this study is to explore the distribution, length, and prevalence of these genomic features, providing insights into the genomic architecture and evolutionary processes of this species.

## Requirements

To conduct this analysis, you need:
- A reference genome
- Population data

These files should be placed in the `data/` directory.

## Usage

1. **Prepare the Data:**
   - Ensure the reference genome and population data are available in the `data/` directory.
   - You may need to modify the names of the input genomes in the Snakemake pipeline to match your specific file names.

2. **Run the Pipeline:**
   - Execute the Snakemake pipeline to reproduce the results and plots from the study.

## Reproducibility

The Snakemake pipeline and scripts provided here ensure that all results and plots are fully reproducible. Detailed instructions and all necessary files are included to replicate the analyses and visualizations.

## Repository Contents

- `Snakefile`: The Snakemake pipeline definition.
- `scripts/`: Directory containing all custom scripts used in the analysis.
- `data/`: Directory where the reference genome and population data should be placed.
- `results/`: Directory where output results and plots will be generated.

## Contact

For any questions or issues, contact me at antonio.iuliano@unil.ch
