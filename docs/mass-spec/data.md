---
layout: default
title: Download data
nav_order: 2
permalink: mass-spec/data
description: "download single-cell proteomics data from SCoPE2, a second generation SCoPE-MS"
---
{% include social-media-links.html %}

# Download single-cell protein and RNA data

&nbsp;

[Processed SCoPE2 Data]({{site.baseurl}}#processed-single-cell-protein-data){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
[RAW SCoPE2 Data]({{site.baseurl}}#RAW-single-cell-protein-data){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
[10x Genomics Data]({{site.baseurl}}#single-cell-RNA-data){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }

&nbsp;

<h2 style="letter-spacing: 2px; font-size: 26px;" id="processed-single-cell-protein-data" >SCoPE2 data processed to ASCII text matrices</h2>

* [Peptides-raw.csv](https://drive.google.com/file/d/15DwDzAKFuRDTV31EnU83aRRhuAyz85_v/view?usp=sharing)
  - `Peptides` **x** `single cells` at 1% FDR.  The first 2 columns list the corresponding protein identifiers and peptide sequences and each subsequent column corresponds to a single cell. Peptide identification is based on spectra analyzed by [MaxQuant](https://www.maxquant.org/)  and is enhanced by using [DART-ID](https://dart-id.slavovlab.net/) to incorporate retention time information. See [Specht et al., 2019](https://www.biorxiv.org/content/10.1101/665307v3) for details.   

&nbsp;

* [Proteins-processed.csv](https://drive.google.com/file/d/1sF5STkofF_f2msnYaaYdWabou84Qf2Xr/view?usp=sharing)
   - `Proteins` **x** `single cells` at 1% FDR, imputed and batch corrected.

&nbsp;

* [Cells.csv](https://drive.google.com/file/d/16vf6rjIsk-oK9naAH6BQnCFrlWnYtJsS/view?usp=sharing)
   - `Annotation` **x**  `single cells`. Each column corresponds to a single cell and the rows include relevant metadata, such as, cell type if known, measurements from the isolation of the cell, and derivative quantities, i.e., rRI, CVs, reliability.

&nbsp;

* [sdrf_meta_data.tsv](https://drive.google.com/file/d/1T8BTfNDlYQkBTs8La6YRSCyD1RwNTvqk/view?usp=sharing)
   -  Meta data following the [Sample to Data file format (SDRF) for Proteomics project guidelines](https://github.com/bigbio/proteomics-metadata-standard) for  for all single cells used in analysis constituting all figures.

&nbsp;

* [Joint protein-RNA data](https://drive.google.com/file/d/130FWc-s-Pd-mx3ymg22bI1qH5fiT7Ktv/view?usp=sharing)
   - `Gene` **x**  `single cells`. Both sets imputed and batch-corrected separately then combined, taking only genes common to both data sets. Uniprot accession numbers used to denote gene.

&nbsp;

* [Signal-to-noise data](https://drive.google.com/file/d/16dmI7qNdpJlPOn83dOZFhHfXv0du5Dip/view?usp=sharing)
  - `Peptides` and `Proteins` **x** `single cells` at 1% FDR.  The first 2 columns list the corresponding protein identifiers and peptide sequences and each subsequent column corresponds to a single cell. The quantitation is the Signal-to-noise (S/N) ratio for each single cell's corresponding reporter ion extracted from the RAW file. The single cell identification numbers are [mapped](https://drive.google.com/file/d/1PUfiGhmInYP3JW5Xoul7Tikl9RSHyQcN/view?usp=sharing) to cell type and RAW file. Complete extracted S/N for each RAW file can be found [here](https://drive.google.com/drive/folders/18_BQ15_JQKzbDt1JZo36MaJuOhN3tJCX?usp=sharing).  


&nbsp;
* [DART-ID input](https://drive.google.com/drive/folders/1ohLco5KHX95jyXIZUAZDvrrbip1RzZ_1?usp=sharing)


&nbsp;
* [GSEA: GOrilla output](https://drive.google.com/drive/folders/1DCp_euY0Cj_NWWG5xQsx7CTN3ju5LI_O?usp=sharing)
&nbsp;

* [Minimal data files](https://drive.google.com/drive/folders/10pOMMlxHsFIyPa9X2auq6xKJssqFgo-D?usp=sharing) necessary for generating Peptides-raw.csv and Proteins-processed.csv

&nbsp;
* [Additional data files](https://drive.google.com/drive/folders/1Zhjik_JFjCQNIVjg63-fooJ4K0HZxWjV?usp=sharing) necessary for generating figures from the [SCoPE2 article](https://doi.org/10.1101/665307).

&nbsp;

* [Processed Data](https://drive.google.com/drive/folders/1NJODxiKrnfW2_nTP-_n_UDvIpwcDEz4C?usp=sharing) from the [second version (v2) of the SCoPE2 preprint](https://www.biorxiv.org/content/10.1101/665307v3)

&nbsp;

* [Processed Data](https://drive.google.com/open?id=1cMQ-SIGpHwSfx9wJF2fIa-t8yX329LPM) from the [first version (v1) of the SCoPE2 preprint](https://www.biorxiv.org/content/10.1101/665307v1)

&nbsp;

* [Single cell proteomics data processing](https://uclouvain-cbio.github.io/scp/): The analysis of the data described here has been replicated by Christophe Vanderaa and Laurent Gatto with the scp Bioconductor package: The scp package is used to process and analyze mass spectrometry-based single cell proteomics data and is freely available from their [Github repository](https://github.com/UCLouvain-CBIO/scp/). The scp package and the replication are described in this [video](https://youtu.be/XMxZkw8yorY).



&nbsp;


<h2 style="letter-spacing: 2px; font-size: 26px;" id="RAW-single-cell-protein-data" >SCoPE2 RAW data and search results from MaxQuant</h2>
The repositories below contain RAW mass-spectrometry data files generated by a Q exactive instrument as well as the search results from analyzing the  RAW files by [MaxQuant](https://www.maxquant.org/)  and by [DART-ID](https://dart-id.slavovlab.net/). The files in Repository 1 were generated in the spring of 2019 and described in the first version of [Specht et al., 2019](https://www.biorxiv.org/content/10.1101/665307v2). The files in Repository 2 were generated in the fall of 2019 from biological replicates and are described in the third version of [Specht et al., 2019](https://www.biorxiv.org/content/10.1101/665307v3). The data in Repository 1 were generated with 11-plex TMT while data in Repository 2 were generated with 16-plex TMT pro.

* [MaxQuant results descriptions](https://drive.google.com/open?id=1qXThKpGPx1tBcxvYFvNM0zCSeyILDzE6)

* [Raw file descriptions](https://drive.google.com/open?id=1-RPN6FOk3ULhkmIH7uc3pIgyQjRkdIdJ)

* **MassIVE Repository 1:**
  - [**http:**  MSV000083945](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=de6aace2096845378ab9ef288e43aa75)
  - [**ftp:** &nbsp; MSV000083945](ftp://massive.ucsd.edu/MSV000083945)

* **MassIVE Repository 2:**
  - [**http:**  MSV000084660](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000084660)
  - [**ftp:** &nbsp; MSV000084660](ftp://massive.ucsd.edu/MSV000084660)



&nbsp;


<h2 style="letter-spacing: 2px; font-size: 26px;" id="single-cell-RNA-data" >scRNA-seq 10x Genomics RAW and processed data</h2>

A cellular mixture identical to that used for the single-cell proteomics was assessed with scRNA-seq using 10x Genomics Chromium platform and the Single Cell 3’ Library & Gel Bead Kit (v2). Two biological replicates of the cell suspension (stock concentration: 1200 cells/μl) were loaded into independent lanes of the device. An average number of about 10,000 cells/lane were recovered. Following the library preparation and sample QC by Agilent BioAnalyzer High Sensitivity chip, the two libraries were pooled together, quantified by KAPA Library Quantification kit and sequenced using the Illumina Novaseq 6000 system (Nova S1 100 flow cell) with the following run parameters: Read1: 26 cycles, i7 index: 8 cycles, Read2: 93 cycles.

* **GEO Repository**
  - [**http:**  GSE142392](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142392)

  &nbsp;  

* **Processed data matrices: `Transcript` x  `single cells` in UMI counts.**
  - [mRNA biological replicate one](https://drive.google.com/open?id=1cN6UgSrZfqKdOjwJ0VyEYp6m_Fy9eANR)

  - [mRNA biological replicate two](https://drive.google.com/open?id=1cuoYiqKgzVnUoFnFmrpXWKVfaFwiboeo)

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

&nbsp;

&nbsp;
