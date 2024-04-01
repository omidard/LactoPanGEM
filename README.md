# Pangenome-Scale Reconstruction of Lactobacillaceae Metabolism

This repository contains the code and documentation for reconstructing strain-specific Genome-Scale Metabolic Models (GEMs) of the Lactobacillaceae family. The workflow is designed to automate the generation of GEMs, identify and fill gaps in the models, and perform various analyses to understand the metabolic capabilities and differences among strains.

## Workflow Overview

### 1. GEM Generation (`GEMgenerator.py`)
The main workflow for generating strain-specific GEMs. It requires four types of input:
- A Reactome model in JSON format.
- A nucleotide FASTA file containing sequences of genes to be incorporated into Reactome Gene-Protein-Reaction associations (GPRs).
- An amino acid FASTA file containing sequences of genes for Reactome GPRs.
- GenBank files (.gbk) of target strains.

The script is optimized for high-performance computing, utilizing multiprocessing over 96 cores on an Azure virtual machine. Adjustments may be necessary to optimize performance on machines with different core counts. The output is draft strain-specific GEMs.

### 2. Gap Filling (`dgap.py` and `run_dgap.py`)
- `dgap.py` contains functions for gap filling.
- `run_dgap.py` leverages these functions to identify and fill gaps in the GEMs, using a Reactome model and the draft GEMs as inputs. This produces gap-filled GEMs.

### 3. Manual Curation
Draft GEMs undergo manual curation to ensure accuracy and completeness. Details of the curation process are discussed in our accompanying paper.

### 4. Computation Requirements
Due to the computational demands of large-scale analysis, the code is optimized for high-performance environments. It may require modifications to adapt to the computational power available for smaller or larger projects.

### 5. Reporting and Analysis Tools
- `species_all_gaps.py`: Generates a report on gap-filled reactions. Requires a directory of GEMs and a dataframe listing GEM IDs and their species labels.
- `species_knockout_fluxes.py`: Calculates flux distributions for all reactions across GEMs, producing a .csv report. Useful for identifying essential compounds and reactions.
- `basic_inf.py`: Generates a .csv file with basic GEM information, such as the number of reactions, gaps, growth rates, and gene counts.
- `lacto_pangem_fva.py`: Performs Flux Variability Analysis (FVA) on panGEMs, outputting a dataframe with minimum and maximum fluxes for all reactions/strains.
- `niche_reactions.py`: Identifies reactions enriched in each niche, using two .csv files: one containing GEM IDs with species names and isolation sources, and another listing reaction presence (1) or absence (0) across all GEMs.
