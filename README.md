# Antibody Affinity Birth Through Somatic Hypermutation

This repository supports the publication **"Antibody affinity birth through somatic hypermutation"**[1]. It contains raw FASTA files, Python scripts, intermediate results, and results used in the analysis described in the paper. The **“raw_fasta_files”** folder contains all raw FASTA sequences analyzed in this project, except for the passenger dataset published by Yeap, Leng-Siew et al. (2015)[2]. The **“analyses”** folder contains Python scripts, intermediate results, and results.

Before running any script, ensure that the data files from the previous stage are decompressed to allow the script to locate its input files.

## Analysis Workflow

### Stage A: Data Preparation

Each 10X Genomics VDJ (BCR) library is accompanied by CSP (cell surface protein) data. By experimental design, each library contains samples from three organs—spleen (SPL), mesenteric lymph nodes (MLN), and Peyer’s patches (PPs)—multiplexed using Hashtag Oligos (HTOs)[3]. HTO barcodes are used to demultiplex the FASTA files for these organs using the **htoDemux** function from the Seurat R package[4], as summarized in the table below:

|     CSP     |     BCR     |   Mouse   |  SPL  |  MLN  |  PPs   |
|:-----------:|:-----------:|:---------:|:-----:|:-----:|:-----:|
| WL-46-CSP   | WL-47-BCR   | B18-383   | HTO1  | HTO2  | HTO3  |
| WL-48-CSP   | WL-49-BCR   | HA-uMT    | HTO4  | HTO5  | HTO6  |
| WL-52-CSP   | WL-53-BCR   | B18-383   | HTO1  | HTO2  | HTO3  |
| WL-54-CSP   | WL-55-BCR   | HA-uMT    | HTO4  | HTO5  | HTO6  |

The passenger dataset, published by Yeap, Leng-Siew et al. (2015)[2], is downloaded from **NCBI BioSample SAMN04014807**, and the R1 and R2 reads are merged to create a single dataset. This merged passenger dataset is added to the 10X and Sanger sequencing datasets for further analysis.

All datasets are annotated using a local version of IgBLAST, **run through the Immcantation pipeline**[5]. The annotated results are stored in the **Stage A data folder**.

### Stage B: Quality Filtering

The datasets from Stage A are inspected for quality based on predefined filtering criteria such as correct locus alignment, heavy and light chain expression, as well as consistent gene annotation and length. Only high-quality sequences are retained. The retained datasets are exported as FASTA files.

### Stage C: Custom Annotation

The FASTA files from Stage B are annotated using a custom IgBLAST database. This database is built from monoclonal antibody heavy and light chain reference sequences from B18-383 and HA mouse models used in this project. The annotated sequences are prepared for additional filtering and downstream analyses.

### Stage D: Additional Quality Filtering

Further quality filtering steps are applied to refine the results from Stage C. These quality filtering steps ensure that only the most reliable and relevant sequences are retained for alignment and data analysis in the next stage.

### Stage E: Data Analysis and Visualization

- **Script E0:** Aligned sequences from Stage D are exported as FASTA alignments.
- **Script E1:** The aligned sequences are processed to generate multiple dataframes, which serve as input data for subsequent scripts.
- **Script E2:** Tables summarizing the number of nucleotide and amino acid mutations per sequence are created.
- **Script E3:** Frequencies of nucleotide and amino acid mutations at each position in different datasets are calculated and stored in tables.
- **Script E4:** Tables for donut plots visualizing the proportions of unmutated and mutated sequences in the datasets are generated.
- **Script E5:** Sequence logos are prepared and visualized using Logomaker[6].
- **Scripts E6 and E7:** R/S (replacement-to-silent mutation) ratio tables are generated.
- **Scripts E8 and E9:** Privacy index-related tables and figures are created and visualized to analyze the distribution of mutations across datasets.

## References

1. "Antibody affinity birth through somatic hypermutation."
2. Yeap, L. S., et al. (2015). "Sequence-Intrinsic Mechanisms that Target AID Mutational Outcomes on Antibody Genes." *Cell*. DOI: [10.1016/j.cell.2015.10.042](https://doi.org/10.1016/j.cell.2015.10.042).
3. Stoeckius, M., et al. (2018). "Cell Hashing with barcoded antibodies enables multiplexing and doublet detection for single cell genomics." *Genome Biology*. DOI: [10.1186/s13059-018-1603-1](https://doi.org/10.1186/s13059-018-1603-1).
4. Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., Fernandez-Granda, C., & Satija, R. (2024). "Dictionary learning for integrative, multimodal and scalable single-cell analysis." *Nature Biotechnology*, 42(2), 293–304. DOI: [10.1038/s41587-023-01767-y](https://doi.org/10.1038/s41587-023-01767-y).
5. Gupta, N. T., et al. (2015). "Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data." *Bioinformatics*. DOI: [10.1093/bioinformatics/btv359](https://doi.org/10.1093/bioinformatics/btv359).
6. Tareen, A., and Kinney, J. B. (2020). "Logomaker: beautiful sequence logos in Python." *Bioinformatics*. DOI: [10.1093/bioinformatics/btz921](https://doi.org/10.1093/bioinformatics/btz921).
