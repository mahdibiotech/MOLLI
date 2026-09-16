MOLLI

MOLLI is an R/Shiny application for the interactive exploration and integration of regulatory protein-binding data and transcriptomic differential-expression results in liver biology.

The application brings together gene-level regulatory information, CRM/TSS-associated protein-binding profiles, and RNA-seq differential-expression results in a single interactive interface.

Project status: research prototype / work in progress.

Overview

MOLLI is designed to help explore relationships between:

genes;

regulatory regions associated with CRM/TSS coordinates;

proteins detected as bound in these regions;

transcriptomic responses measured across multiple experimental comparisons.

The current version provides several complementary views for interactome exploration, transcriptomic visualization, and integrative analyses.

Importantly, an observed protein-binding signal combined with differential expression is compatible with a potential regulatory relationship but does not by itself demonstrate direct transcriptional regulation.

Main features

1. Global gene selection

A shared gene-selection system is used across the application so that the same genes can be explored in several views.

Selected genes can be:

added manually;

removed interactively;

explored in regulatory and transcriptomic views;

transferred to the detailed gene view.

2. Interactome exploration

The interactome module explores a binary gene × protein binding matrix.

Available visualizations include:

binary binding heatmap;

protein–protein Jaccard similarity heatmap;

protein binding-frequency plot;

protein-combination frequency view;

gene–protein bipartite network.

The current implementation aggregates binding information at the gene level.

3. CRM / TSS exploration

The CRM / TSS view displays regulatory-region information associated with selected genes and proteins.

Available information includes:

gene;

protein;

CRM chromosome;

CRM start and end coordinates;

CRM size;

associated TSS chromosome and position.

The current binding dataset contains one regulatory entry per gene. Therefore, this module should not be interpreted as a transcript- or isoform-level promoter annotation.

4. Transcriptome exploration

The transcriptomic module visualizes differential-expression results across experimental comparisons.

Available views include:

Bubble Plot;

signed log2 fold-change heatmap;

Volcano Plot;

multi-study expression profile;

Top differentially expressed genes.

The application supports filtering by:

gene;

experimental comparison (StudiesCond);

FDR threshold;

significance status.

5. Detailed gene view

The Gene Card brings together information for one selected gene:

CRM/TSS coordinates;

associated proteins;

number of bound proteins;

differential-expression profile across studies;

RNA-seq result table.

This view is intended to provide a compact biological summary for gene-centered exploration.

6. Integrative analyses

MOLLI integrates protein-binding and RNA-seq information for genes shared between both datasets.

Current analyses include:

enrichment of differentially expressed genes among bound genes;

Fisher's exact test;

Benjamini–Hochberg correction for multiple testing in the protein × study enrichment grid;

separate enrichment analyses for up- and down-regulated genes;

comparison of signed log2 fold-change distributions between bound and unbound genes using a Wilcoxon test;

bound versus unbound expression visualization;

integrated gene–protein network.

These analyses are exploratory and should be interpreted together with the biological design of each experiment.

7. Dataset and QC summary

A dedicated tab provides technical summaries of the loaded datasets, including:

number of rows;

number of genes;

number of proteins;

number of experimental comparisons;

duplicate identifiers;

missing values;

overlap between regulatory and transcriptomic datasets;

protein-binding frequencies.

These checks are intended as technical input and consistency controls, not as biological quality-control analyses.

Input files

MOLLI currently expects two input files in the application directory.

Regulatory binding dataset

Expected filename:

pool_MultiBind_ProseqActiveLiverPromoters_Annotated_withColnames.bed

Required metadata columns:

chrCRM
startCRM
endCRM
tssC
tssS
tssE
Name
inter

All additional columns are interpreted as protein-binding columns.

Protein values are expected to represent binary binding information.

RNA-seq differential-expression dataset

Expected filename:

RNAseq_DE_results.tsv

Required columns:

gene
StudiesCond
abslogFC
FDR
State

State is expected to identify the direction of regulation (up or down).

MOLLI reconstructs a signed fold-change value internally:

up   -> +abslogFC
down -> -abslogFC

Data availability

The biological input datasets used during development are not distributed in this repository.

They are excluded through .gitignore.

To run MOLLI with the current configuration, place compatible input files in the same directory as app.R, or modify the file paths defined at the beginning of the script.

Installation

Requirements

R

RStudio is recommended but not required

The application currently uses the following R packages:

install.packages(c(
  "shiny",
  "data.table",
  "dplyr",
  "pheatmap",
  "DT",
  "ggplot2",
  "colourpicker",
  "shinyWidgets",
  "scales"
))

Package and R version information used during development is provided in:

MOLLI_sessionInfo.txt

Running the application

Clone the repository:

git clone https://github.com/mahdibiotech/MOLLI.git
cd MOLLI

Place the required input files in the project directory.

Then start R in the project directory and run:

shiny::runApp()

Alternatively:

shiny::runApp("app.R")

If Rscript is available from the command line, the application can also be started with:

Rscript app.R

Repository structure

MOLLI/
├── app.R
├── MOLLI_sessionInfo.txt
├── .gitignore
└── README.md

The current version is intentionally distributed as a single-file Shiny application.

A future development step is to modularize the code into separate components for data loading, interactome analyses, transcriptomic analyses, integration, statistics, and visualization.

Current analysis workflow

Regulatory binding data
          │
          ├───────────────┐
          │               │
          ▼               │
  Interactome / CRM-TSS   │
                          │
                          ▼
                    Integration
                          ▲
                          │
          ┌───────────────┘
          │
          ▼
RNA-seq differential-expression data
          │
          ▼
 Transcriptome exploration

Current limitations

The current version has several important limitations:

analyses are primarily gene-level;

the present regulatory dataset does not provide transcript/isoform-level promoter resolution;

results depend on the structure and naming conventions of the input files;

gene identifiers must be compatible across regulatory and transcriptomic datasets;

the application currently loads local files directly;

the application remains monolithic and has not yet been converted into Shiny modules;

integrative statistical results are exploratory and do not establish causal or direct regulation;

broader validation on additional datasets remains necessary.

Planned developments

Potential future developments include:

modularization of the Shiny codebase;

direct import of standardized CRESCENT output;

improved provenance and metadata tracking;

transcript- and isoform-level analyses;

additional DET, DTU and transposable-element views;

interactive selection directly from additional plots;

improved network layouts;

downloadable analysis settings and metadata;

automated testing;

deployment on a Shiny-compatible server.

Development context and attribution

This repository contains a redevelopment and extension of an existing application.

The database, the general concept of the application, and part of the visualization framework pre-dated the work represented by this repository. The current development work focused on restructuring and extending the R/Shiny application, improving the backend, integrating transcriptomic results, expanding interactive exploration, adding export and visualization capabilities, and developing integrative analyses between regulatory binding and differential-expression data.

The work was carried out in the context of a bioinformatics project on integrated transcriptional analysis of liver disease.

Reproducibility

For reproducibility, the repository includes the R session information used during development:

MOLLI_sessionInfo.txt

For a more fully reproducible software environment, future versions may use a dependency-management system such as renv.

Citation

If MOLLI is used in academic work, please cite the repository and the associated laboratory/project publications when they become available.

A formal software citation will be added if a versioned public release is archived.

License

No open-source license is currently provided.

Reuse, redistribution, or modification should therefore not be assumed to be permitted until the ownership and licensing conditions of the pre-existing application components have been clarified.

Author

Mahdi
Bioinformatics

GitHub: @mahdibiotech

Disclaimer

MOLLI is a research prototype intended for exploratory analysis.

The software is provided for research purposes and should not be used as a standalone basis for clinical or medical decision-making.
