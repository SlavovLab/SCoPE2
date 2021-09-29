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

## nano-ProteOmic sample Preparation (nPOP)

&nbsp;

<span class="text-center"></span>
[bioRxiv Preprint](https://www.biorxiv.org/content/10.1101/2021.04.24.441211v2){: .btn .fs-5 .mr-4 }

<span class="text-center"></span>
[Step by Step Protocol](https://www.protocols.io/view/npop-bwy7pfzn){: .btn .fs-5 .mr-4 }


**Table of Contents**

1. [Introduction to droplet sample preparation for single-cell proteomics](#introduction)
2. [Application to the cell division cycle](#application)
3. [RAW Data](#raw_data)
4. [Processed Data](#proc_data)
4. [Video presentations](#talks)


## Introduction
Mass spectrometry methods have enabled quantifying thousands of proteins at the single cell level. These methods open the door to many biological opportunities, such as characterizing heterogeneity in the tumor micro-environment, signaling pathways driving stem cell differentiation, and intrinsically single-cell processes, such as the cell division cycle. To further advance single-cell MS analysis, we developed an automated nano-ProteOmic sample Preparation (nPOP). nPOP uses piezo acoustic dispensing to isolate individual cells in 300 picoliter volumes and performs all subsequent preparation steps in small droplets on a hydrophobic slide. This allows massively parallel sample preparation, including lysing, digesting, and labeling individual cells in volumes below 20 nanoliters.

## Application
Single-cell protein analysis using nPOP classified cells by cell type and by cell cycle phase. Furthermore, the data allowed us to quantify the covariation between cell cycle protein markers and thousands of proteins. Based on this covariation, we identify cell cycle associated proteins and functions that are shared across cell types and those that differ between cell types.   


<h2 style="letter-spacing: 2px; font-size: 26px;" id="raw_data" >Raw Data from experiments benchmarking nPOP</h2>

* **MassIVE Repository:**
  - [**http:**  MSV000082841](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=0374fefddfc64cb8b400f77e4c19536e)
  - [**ftp:** &nbsp; MSV000082841](ftp://massive.ucsd.edu/MSV000087152)
  - [**http:**  MSV000088167](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=2f82c5f336a441d7a7aee378d84f7a58)
  - [**ftp:** &nbsp; MSV000088167](ftp://massive.ucsd.edu/MSV000088167)

  &nbsp;

  &nbsp;

<h2 style="letter-spacing: 2px; font-size: 26px;" id="proc_data" >Processed Data from experiments benchmarking nPOP</h2>

* [Peptides-raw.csv](https://drive.google.com/file/d/1cVEq5KIHdyhVfDObo31W2GbFH5XgDne8/view?usp=sharing)
   - `Peptides` **x** `single cells` at 1% FDR.  The first columns list the peptide sequences and each subsequent column corresponds to a single cell. Peptide identification is based on spectra analyzed by [MaxQuant](https://www.maxquant.org/)  and is enhanced by using [DART-ID](https://dart-id.slavovlab.net/) to incorporate retention time information. See [Specht et al., 2019](https://www.biorxiv.org/content/10.1101/665307v3) for details.

&nbsp;

* [Proteins-processed.csv](https://drive.google.com/file/d/1LHyHfE0WoyVWbMyhYD1DtnpVxfRz5a79/view?usp=sharing)
   - `Proteins` **x** `single cells` at 1% FDR, imputed and batch corrected.

&nbsp;

* [Cells.csv](https://drive.google.com/file/d/1PKaGrIOizIxx9zmrM7ShZYiTaM49p_9R/view?usp=sharing)
   - `Annotation` **x**  `single cells`. Each column corresponds to a single cell and the rows include relevant metadata, such as, cell type, measurements from the isolation of the cell, and derivative quantities, i.e., rRI, CVs, reliability. This file corresponds to Proteins-processed.csv and Peptides-raw.csv files.

&nbsp;


* [HeLa-proteins.csv](https://drive.google.com/file/d/1BMj5YF_qVu34JXkcBn54GhdB2fS2-_0z/view?usp=sharing)
   - `Proteins` **x** `single cells` for HeLa cells at 1% FDR, unimputed and z-scored.

&nbsp;

* [U-937-proteins.csv](https://drive.google.com/file/d/1BLNher4z0agGGoJjM2VRRGAGS077ONYi/view?usp=sharing)
   - `Proteins` **x** `single cells` for U-937 cells at 1% FDR, unimputed and z-scored.

&nbsp;


* [CellCycle-Proteins.csv](https://drive.google.com/file/d/1BM4ffkpu0vW_9rfSnmkPwA66RTGKRd21/view?usp=sharing)
   - A list of cell cycle dependent proteins used in analysis. Cell cycle protein and phases information was used from [Whitfield et al, 2002](http://genome-www.stanford.edu/Human-CellCycle/Hela/).

&nbsp;


<h2 style="letter-spacing: 2px; font-size: 26px;" id="talks" >Recorded video presentations</h2>

<iframe width="560" height="315" src="https://www.youtube.com/embed/DJ1U_KpMNcY" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>


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
