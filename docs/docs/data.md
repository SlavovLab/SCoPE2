---
layout: default
title: Download data
nav_order: 2
permalink: docs/data
---

# SCoPE2 data download

&nbsp;


### Quantifying proteins in single cells at high-throughput with mass spectrometry

[SCoPE2 Preprint](https://www.biorxiv.org/content/10.1101/665307v1){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }[Perspective](https://pubs.acs.org/doi/10.1021/acs.jproteome.8b00257){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2}[Data Repository](ftp://massive.ucsd.edu/MSV000083945){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }[GitHub](https://github.com/SlavovLab/){: .btn .fs-5 .mb-4 .mb-md-0 }

------------



## SCoPE2 data processed to ASCII text matrices


* [Peptides-raw.csv](http://slavovlab.net/scope2/data/Peptides-raw.csv)
  - Peptides by single cells at 1% FDR based on [DART-ID](https://dart-id.slavovlab.net/). The first 2 columns list the corresponding protein identifiers and peptide sequences and each subsequent column corresponds to a single cell.

&nbsp;

* [Proteins-processed.csv](http://slavovlab.net/scope2/data/Proteins-processed.csv)
   - Proteins by single cells at 1% FDR, imputed and batch corrected.

&nbsp;

* [Cells.csv](http://slavovlab.net/scope2/data/Cells.csv)
   - Annotation by single cells. Each column corresponds to a single cell and the rows include relevant metadata, such as, cell type if known, measurements from the isolation of the cell, and derivative quantities, i.e., rRI, CVs, reliability.

&nbsp;


## SCoPE2 RAW data and search results from MaxQuant

* **MassIVE Repository:**
  - [http: MSV000083945](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=de6aace2096845378ab9ef288e43aa75)
  - [ftp: MSV000083945](ftp://massive.ucsd.edu/MSV000083945)
