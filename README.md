# Spatial Proximity of UTR Variants of Uncertain Significance to ClinVar-Classified Pathogenic Variants

Analysis code accompanying the manuscript:

**Spatial Proximity of UTR Variants of Uncertain Significance to ClinVar-Classified Pathogenic Variants and Its Association with Review Status**

## Overview

This repository contains the analysis code used to investigate the spatial proximity of untranslated region (UTR) variants of uncertain significance (VUS) to ClinVar-classified pathogenic or likely pathogenic (P/LP) variants.

The analyses examine whether apparent spatial proximity is robust to alternative null models, gene-level structure, review status, population allele frequency, transcript overlap, conservation, and other sensitivity analyses.

The study does not treat spatial proximity as independent evidence of pathogenicity. Instead, the results are interpreted in the context of ClinVar's spatial sampling and curation structure.

## Repository structure

- code/ — analysis scripts used for the manuscript
- data/ — local input data directory; raw public datasets are not included
- esults/ — local output directory; generated result files are not included
- .gitignore — excludes raw datasets, generated outputs, and local files

## Data sources

The analyses use publicly available variant resources, including:

- ClinVar
- MANE
- gnomAD

Raw source datasets are not redistributed in this repository. Users should obtain the relevant source datasets directly from their official repositories and place local input files in the data/ directory.

## Reproducibility

The scripts in code/ correspond to the main statistical analyses and sensitivity analyses reported in the manuscript.

Because the source datasets are external and may be updated over time, exact reproduction should use the same dataset versions and annotations described in the manuscript.

## License

License information will be added before publication.
