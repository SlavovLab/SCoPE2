---
layout: default
title: nPOP
nav_order: 1
permalink: nPOP
description: "nano-ProteOmic sample Preparation (nPOP)"
parent: Sample preparation

---
{% include social-media-links.html %}
# Droplet sample preparation

## nano ProteOmic sample Preparation (nPOP)

&nbsp;

<span class="text-center"></span>
[bioRxiv Preprint](https://doi.org/10.1101/2020.08.24.264994){: .btn .fs-5 .mr-4 }


**Table of Contents**

1. [Abstract](#abstract)
2. [RAW Data](#raw_data)
3. [Processed Data](#proc_data)


## Abstract
Mass spectrometry methods have enabled quantifying thousands of proteins at the single cell level. These methods open the door to tackling many biological challenges, such as characterizing heterogeneity in the tumor micro-environment and better understanding signaling pathways driving stem cell differentiation. To further advance single-cell MS analysis, we developed an automated nano-ProteOmic sample Preparation (nPOP).



<h2 style="letter-spacing: 2px; font-size: 26px;" id="raw_data" >Raw Data from experiments benchmarking nPOP</h2>

* **MassIVE Repository:**
  - [**http:**  MSV000082841](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=0374fefddfc64cb8b400f77e4c19536e)
  - [**ftp:** &nbsp; MSV000082841](ftp://massive.ucsd.edu/MSV000087152)

  &nbsp;

  &nbsp;

<h2 style="letter-spacing: 2px; font-size: 26px;" id="proc_data" >Processed Data from experiments benchmarking nPOP</h2>

* [Peptides-raw.csv]()
   - `Peptides` **x** `single cells` at 1% FDR.  The first 2 columns list the corresponding protein identifiers and peptide sequences and each subsequent column corresponds to a single cell. Peptide identification is based on spectra analyzed by [MaxQuant](https://www.maxquant.org/)  and is enhanced by using [DART-ID](https://dart-id.slavovlab.net/) to incorporate retention time information. See [Specht et al., 2019](https://www.biorxiv.org/content/10.1101/665307v3) for details.

&nbsp;

* [Proteins-processed.csv](https://drive.google.com/file/d/1sF5STkofF_f2msnYaaYdWabou84Qf2Xr/view?usp=sharing)
   - `Proteins` **x** `single cells` at 1% FDR, imputed and batch corrected.
&nbsp;

* [HeLa-proteins.csv](https://drive.google.com/file/d/1sF5STkofF_f2msnYaaYdWabou84Qf2Xr/view?usp=sharing)
   - `Proteins` **x** `single cells` for HeLa cells at 1% FDR, unimputed and zscored.

&nbsp;

* [U-937-proteins.csv](https://drive.google.com/file/d/1sF5STkofF_f2msnYaaYdWabou84Qf2Xr/view?usp=sharing)
   - `Proteins` **x** `single cells` for U-937 cells at 1% FDR, unimputed and zscored.

&nbsp;

* [Cells.csv](https://drive.google.com/file/d/1sF5STkofF_f2msnYaaYdWabou84Qf2Xr/view?usp=sharing)
   - `Annotation` **x**  `single cells`. Each column corresponds to a single cell and the rows include relevant metadata, such as, cell type if known, measurements from the isolation of the cell, and derivative quantities, i.e., rRI, CVs, reliability.


&nbsp;

* [CellCycle-Proteins.csv](https://drive.google.com/file/d/1sF5STkofF_f2msnYaaYdWabou84Qf2Xr/view?usp=sharing)
   - A list of cell cycle dependent proteins used in analysis. Proteins used taken form [Botstein et al, 2002](http://genome-www.stanford.edu/Human-CellCycle/Hela/)

&nbsp;


<!--
<span class="text-center"></span>
[bioRxiv Preprint](https://doi.org/10.1101/2020.08.24.264994){: .btn .fs-5 .mr-4 }

**Table of Contents**

1. [Abstract](#abstract)
2. [RAW Data](#data)


## Abstract

Mass spectrometry methods have enabled quantifying thousands of proteins at the single cell level. These methods open the door to tackling many biological challenges, such as characterizing heterogeneity in the tumor micro-environment and better understanding signaling pathways driving stem cell differentiation. To further advance single-cell MS analysis, we developed an automated nano-ProteOmic sample Preparation (nPOP). nPOP isolates individual cells in 300 picoliter volumes and performs all subsequent preparation steps in small droplets on a hydrophobic glass slide, which allows to keep sample volumes below 15 nl.


 


&nbsp;


<h2 style="letter-spacing: 2px; font-size: 26px;" id="data" >Data from experiments with increasing isobaric carriers</h2>

* **MassIVE Repository:**
  - [**http:**  MSV000082841](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=bfd7f21d718940fdbaccc0d58ad6b122)
  - [**ftp:** &nbsp; MSV000082841](ftp://massive.ucsd.edu/MSV000082841)

  &nbsp;

  &nbsp;


&nbsp;

&nbsp;  

&nbsp;

## About the project

This project on characterizing the isobaric carrier was conducted in the [Slavov Laboratory](https://slavovlab.net) and [SCP Center](https://center.single-cell.net) at [Northeastern University](https://www.northeastern.edu/), and was authored by [Harrison Specht](http://harrisonspecht.com) and [Nikolai Slavov](https://coe.northeastern.edu/people/slavov-nikolai/). Learn more about [single-cell mass-spectrometry analysis](https://scope2.slavovlab.net/mass-spec/single-cell-proteomics).  


This project was supported by funding from the [NIH Director's Award](https://projectreporter.nih.gov/project_info_description.cfm?aid=9167004&icde=31336575).

-->

&nbsp;  

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;

&nbsp;
