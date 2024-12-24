# stereo-xcr-seq

## Description

**Stereo-XCR-seq** is an innovative approach designed to overcome the lack of tools for **in situ single-cell T/BCR (XCR) sequencing**. This efficient strategy retrieves and sequences **TCR** and **BCR** from Stereo-seq cDNA libraries at subcellular resolution. Stereo-XCR-seq provides unbiased full-length XCR reads alongside **spatial transcriptomics**, enabling comprehensive insights into immune repertoire and spatial gene expression.

## Repository Structure

This repository is organized into four main directories to streamline data integration, analysis, and visualization:

1. **`meta_build`**  
   - Scripts for integrating second-generation (NGS) and third-generation (TGS) sequencing data.  
   - Constructs meta information and adata objects for downstream analysis.

2. **`cycclone_immune_repertoire`**  
   - Tools and scripts for analyzing **third-generation sequencing** (TGS) data.

3. **`NGS_immune_repertoire`**  
   - Scripts for processing and analyzing **second-generation sequencing** (NGS) data.

4. **`analysis_script`**  
   - Includes the full analysis pipeline used in the study.  
   - Contains code for generating **figures** used in the publication.
