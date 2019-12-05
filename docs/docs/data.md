---
layout: default
title: Download data
nav_order: 2
permalink: docs/data
description: "download single-cell proteomics data from SCoPE2, a second generation SCoPE-MS" 
---

# SCoPE2 data download

&nbsp;


### Quantifying proteins in single cells at high-throughput with mass spectrometry

[SCoPE2 Preprint](https://www.biorxiv.org/content/10.1101/665307v1){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }[Perspective](https://pubs.acs.org/doi/10.1021/acs.jproteome.8b00257){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2}[A dream of single-cell proteomics](https://www.nature.com/articles/s41592-019-0540-6){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }[Data Repository 1](ftp://massive.ucsd.edu/MSV000083945){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }[Data Repository 1](ftp://massive.ucsd.edu/MSV000084660){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }[GitHub](https://github.com/SlavovLab/SCoPE2/tree/master/code){: .btn .fs-5 .mb-4 .mb-md-0 }

------------



## SCoPE2 data processed to ASCII text matrices


* [Peptides-raw.csv](https://drive.google.com/drive/folders/1bx4sBkKHe1R4VLH_8U8x4A5LFOjNluSO?usp=sharing)
  - `Peptides` **x** `single cells` at 1% FDR.  The first 2 columns list the corresponding protein identifiers and peptide sequences and each subsequent column corresponds to a single cell. Peptide identification is based on spectra analyzed by [MaxQuant](https://www.maxquant.org/)  and is enhanced by using [DART-ID](https://dart-id.slavovlab.net/) to incorporate retention time information. See [Specht et al., 2019](https://www.biorxiv.org/content/10.1101/665307v1) for details.   

&nbsp;

* [Proteins-processed.csv](https://drive.google.com/drive/folders/1bx4sBkKHe1R4VLH_8U8x4A5LFOjNluSO?usp=sharing)
   - `Proteins` **x** `single cells` at 1% FDR, imputed and batch corrected.

&nbsp;

* [Cells.csv](https://drive.google.com/drive/folders/1bx4sBkKHe1R4VLH_8U8x4A5LFOjNluSO?usp=sharing)
   - `Annotation` **x**  `single cells`. Each column corresponds to a single cell and the rows include relevant metadata, such as, cell type if known, measurements from the isolation of the cell, and derivative quantities, i.e., rRI, CVs, reliability.

&nbsp;

* [Joint protein-RNA data](https://drive.google.com/drive/folders/1bx4sBkKHe1R4VLH_8U8x4A5LFOjNluSO?usp=sharing)
   - `Gene` **x**  `single cells`. Boths sets imputed and batch-corrected separately then combined, taking only genes common to both data sets. Uniprot accession numbers used to denote gene. 

&nbsp;

* [mRNA biological replicate two, raw data](https://drive.google.com/drive/folders/1bx4sBkKHe1R4VLH_8U8x4A5LFOjNluSO?usp=sharing)
   - `Transcript` **x**  `single cells`. Raw UMI counts. 

&nbsp;

* [mRNA biological replicate two, raw data](https://drive.google.com/drive/folders/1bx4sBkKHe1R4VLH_8U8x4A5LFOjNluSO?usp=sharing)
   - `Transcript` **x**  `single cells`. Raw UMI counts. 

&nbsp;


## SCoPE2 RAW data and search results from MaxQuant

* **MassIVE Repository 1:**
  - [**http:**  MSV000083945](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=de6aace2096845378ab9ef288e43aa75)
  - [**ftp:** &nbsp; MSV000083945](ftp://massive.ucsd.edu/MSV000083945)

* **MassIVE Repository 2:**
  - [**http:**  MSV000084660](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000084660)
  - [**ftp:** &nbsp; MSV000084660](ftp://massive.ucsd.edu/MSV000084660)