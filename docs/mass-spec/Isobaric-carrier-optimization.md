---
layout: default
title: Isobaric carrier optimization
nav_order: 1
permalink: mass-spec/Isobaric-carrier-optimization
description: "Optimizing accuracy and depth of protein quantification in experiments using isobaric carriers"
parent: Optimization

---
{% include social-media-links.html %}

# Isobaric carrier optimization

&nbsp;

<span class="text-center"></span>
[JPR Article](https://doi.org/10.1021/acs.jproteome.0c00675){: .btn .fs-5 .mr-4 }
[bioRxiv Preprint](https://doi.org/10.1101/2020.08.24.264994){: .btn .fs-5 .mr-4 }
[GitHub Code](https://github.com/SlavovLab/isobaric-carrier){: .btn .fs-5 .mr-4 }

**Table of Contents**

1. [Abstract](#abstract)
2. [RAW Data](#data)
3. [DO-MS Reports](#do-ms-reports)
4. [Guidelines](#guidelines)

## Abstract

The isobaric carrier approach, which combines small isobarically-labeled samples with a larger isobarically-labeled carrier sample, is finding diverse applications in ultrasensitive mass-spectrometry analysis of very small samples, such as single cells. To inform the growing use of isobaric carriers, we characterized the trade-offs of using isobaric carriers in controlled experiments with complex human proteomes. The data indicate that isobaric carriers directly enhances peptide sequence identification without simultaneously increasing the number of protein copies sampled from small samples. Furthermore, the results indicate strategies for optimizing the amount of isobaric carrier and analytical parameters, such as ion accumulation time, for different priorities such as improved quantification or increased number of identified proteins. Balancing these trade-offs enables adapting isobaric carrier experiments to different applications, such as quantifying the proteome of limited biopsies or organoids, building single-cell atlases, or modeling protein networks at single-cell resolution. In all cases, the reliability of protein quantification should be estimated and incorporated in all subsequent data analysis. We expect that these guidelines will aid explicit incorporation of the characterized trade-offs in experimental designs and transparent error propagation in the data analysis.

![]({{site.baseurl}}Figures/Single-cell-Proteomics_Applications_iCarrier.png){: width="80%" .center-image}
<!--
To increase the throughput and quantitative accuracy of single-cell protein analysis by [SCoPE-MS](https://doi.org/10.1101/102681), we introduced many technical improvements in both the sample preparation and in the mass-spectrometry analysis. The [synergistic effect](https://www.biorxiv.org/content/biorxiv/early/2019/12/05/665307/T1.medium.gif) is to increase quantitative accuracy by 4-fold and the throughput of data acquisition about 8-fold. Below, we outline controlled experiments that illustrate the benefits of **individual** improvements. To comprehensively compare the mass-spec data at all levels (including chromatographic separation, precursor abundance, ion isolation, spectral purity, and peptide sequence identification), we include the full [Data-driven Optimization of MS (DO-MS)](https://do-ms.slavovlab.net) reports for each set of experiments. 
-->

&nbsp;


<h2 style="letter-spacing: 2px; font-size: 26px;" id="data" >Data from experiments with increasing isobaric carriers</h2>

* **MassIVE Repository:**
  - [**http:**  MSV000086004](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=b432a22b241e4c4881d63f1a97db4a4e)
  - [**ftp:** &nbsp; MSV000086004](ftp://massive.ucsd.edu/MSV000086004)

  &nbsp;

  &nbsp;

## DO-MS Reports
These reports contain plots with detailed information for the chromatographic separation (elution peaks), precursor ion intensities, ion accumulation times, apex sampling, number of MS2 scans per experiment and the confidence of the assigned PSMs, and contaminants. This information is provided as distributions across all ions as previously described. For more information, see [Data-driven Optimization of Mass-Spectrometry (DO-MS)](https://do-ms.slavovlab.net).  
* [DO-MS Report for the low MS2 AGC target](DO-MS_Reports/DO-MS_Report_lowAGC.html)
* [DO-MS Report for the high MS2 AGC target](DO-MS_Reports/DO-MS_Report_highAGC.html)

  &nbsp;


  <h2 style="letter-spacing: 2px; font-size: 26px;" id="guidelines" >Guidelines for optimizing isobaric carrier experiments</h2>

  - **The size of the isobaric carrier should reflect the project priorities.** The [data](https://www.biorxiv.org/content/10.1101/2020.08.24.264994v4) reinforce previous suggestions that isobaric carriers that are about 200-fold larger than the small samples provide most of the needed increase in peptide fragments to enhance sequence identification without adverse affects on quantification. Larger carrier sizes can further benefit peptide sequence identification but may result in reduced copy number sampling and noisier data.   

  - **Evaluate whether the isobaric carrier affects peptide sampling.** We recommend estimating whether the carrier levels and AGC target result in reduced accumulation times and sampling of proteins from the small sample. This can be visualized by plotting the distributions of accumulation times and reporter ion intensities. These distributions can be automatically generated by [Data-driven Optimization of MS (DO-MS)](https://do-ms.slavovlab.net).

  - **Estimate the reliability of quantification for each protein.**  Estimate the sampling error on a per-protein basis as previously [demonstrated](https://doi.org/10.1101/665307), benchmark protein quantification or estimate the reliability of quantification based on the consistency of the quantification of peptides originating from the same protein as demonstrated by Franks and colleagues. If the reliability of quantification is limited by counting noise, improved sampling can increase reliability. However, if the reliability is limited by other factors, such as sample preparation, improved copy-number sampling may have limited benefits.  

  - **Incorporate estimates of reliability in all subsequent analysis.** Data points should be weighted based on their reliability, with weights proportional to the reliability. Use error propagation methods to reflect the noise in the final results. For example, correlations between noisy variables can be divided by the corresponding reliability to estimate the fraction of explained variance independent from the noise.



  &nbsp;


## Media
- [Twitter thread](https://threadreaderapp.com/thread/1298059519201869826.html)

&nbsp;

&nbsp;  

&nbsp;

## About the project

This project on characterizing the isobaric carrier was conducted in the [Slavov Laboratory](http://slavovlab.net) and [SCP Center](http://center.single-cell.net) at [Northeastern University](https://www.northeastern.edu/), and was authored by [Harrison Specht](http://harrisonspecht.com) and [Nikolai Slavov](https://coe.northeastern.edu/people/slavov-nikolai/). Learn more about [single-cell mass-spectrometry analysis](https://scope2.slavovlab.net/mass-spec/single-cell-proteomics).  


This project was supported by funding from the [NIH Director's Award](https://projectreporter.nih.gov/project_info_description.cfm?aid=9167004&icde=31336575) and by an [Allen Distinguished Investigator Award](https://alleninstitute.org/what-we-do/frontiers-group/distinguished-investigators/projects/tracking-proteome-dynamics-single-cells) from the Paul G. Allen Frontiers Group.

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
