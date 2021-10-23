---
layout: default
title: Prioritization
nav_order: 7
description: "Prioritized single-cell proteomics for high-throughput and multiplexed single-cell proteomics by mass-spectrometry using pSCoPE"
permalink: prioritization
---

# Prioritized single-cell proteomics

&nbsp;

&nbsp;


<iframe width="560" height="315" src="https://www.youtube.com/embed/SP0x3gAALtg" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>


<!--

[SCoPE2 Protocol Preprint](https://www.biorxiv.org/content/10.1101/2021.03.12.435034v1){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
[SCoPE2 code on GitHub](https://github.com/SlavovLab/SCoPE2/tree/master/code){: .btn .fs-5 .mb-4 .mb-md-0 }



------------

Many biological systems are composed of diverse single cells. This diversity necessitates functional and molecular single-cell analysis. Single-cell protein analysis has long relied on affinity reagents, but emerging mass-spectrometry methods (either label-free or multiplexed) have enabled quantifying over 1,000 proteins per cell while simultaneously increasing the specificity of protein quantification. Isobaric carrier based multiplexed single-cell proteomics is a scalable, reliable, and cost-effective method that can be fully automated and implemented on widely available equipment. It uses inexpensive reagents and is applicable to any sample that can be processed to a single-cell suspension. Here we describe an automated SCoPE2 workflow that allows analyzing about 200 single cells per 24 hours using only standard commercial equipment. We emphasize experimental steps and benchmarks required for achieving quantitative protein analysis.

[Protocol data]({{site.baseurl}}#data){: .btn .fs-3 .mr-2}
[Recorded workshops]({{site.baseurl}}#workshops){: .btn .fs-3 .mr-2}

&nbsp;

## Data

* [RAW Files](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=66e7837857194b67b3050099747833e3)

* [Peptides x Samples](Protocol_data/SingleCell_PeptidesXsamples.txt)
  - `Peptides` **x** `single cells` at 1% FDR.  The first 2 columns list the corresponding protein identifiers and peptide sequences and each subsequent column corresponds to a single cell. Peptide identification is based on spectra analyzed by [MaxQuant](https://www.maxquant.org/)  and is enhanced by using [DART-ID](https://dart-id.slavovlab.net/) to incorporate retention time information. See [Specht et al., 2019](https://www.biorxiv.org/content/10.1101/665307v3) for details.   

&nbsp;

* [Peptides x Samples](Protocol_data/ProteinsXSamples_BulkANDSSC.txt)
   - `Proteins` **x** `single cells` at 1% FDR, imputed and batch corrected.

&nbsp;

* [Cells.txt](Protocol_data/SingleCell_ids.txt)
   - `Annotation` **x**  `single cells`. Each column corresponds to a single cell and the rows include relevant metadata, such as, cell type if known, measurements from the isolation of the cell, and derivative quantities, i.e., rRI, CVs, reliability.


&nbsp;

* Source data for figures
  - Figure 2: [LC-MS/MS setup for SCoPE2 experiments](https://doi.org/10.6084/m9.figshare.15060720.v1)

  - Figure 3: [Evaluating data acquisition and interpretation using DO-MS diagnostic plots](https://doi.org/10.6084/m9.figshare.15060774.v1)

  - Figure 4: [Principal Component Analysis of results from SCoPE2 analysis](https://doi.org/10.6084/m9.figshare.15060789.v1)

  - Extended Data Figure 1: [The ABIRD device may suppress contaminant ions](https://doi.org/10.6084/m9.figshare.15060846.v1)



## Workshops
The presentations below describe sample preparation and data analysis for the latest version of Single Cell ProtEomics by Mass Spectrometry (SCoPE-MS), [SCoPE2](https://www.biorxiv.org/content/10.1101/665307v5). They were presented at the workshop of the second [Single Cell Proteomics Conference](https://single-cell.net/proteomics/scp2019).


* [Design of single-cell proteomics experiments using SCoPE2](https://youtu.be/mz6Yq2XSu-8)
* [Sample preparation for single-cell MS analysis by SCoPE2](https://youtu.be/Eq_s6Jlzfnk)
* [High-throughput single-cell proteomics quantifies the emergence of macrophage heterogeneity](https://youtu.be/NNLh4nE687I)

&nbsp;

<iframe width="560" height="315" src="https://www.youtube.com/embed/Eq_s6Jlzfnk" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

-->

&nbsp;

&nbsp;

&nbsp;
