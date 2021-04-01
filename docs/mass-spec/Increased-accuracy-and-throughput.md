---
layout: default
title: Improvments
nav_order: 2
permalink: mass-spec/Increased-accuracy-and-throughput
description: "Improvements in the accuracy and throughout of SCoPE2 over SCoPE-MS"
parent: Optimization
---

# Technical improvements

&nbsp;

[Decreasing coisolation]({{site.baseurl}}#decreasing-coisolation){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
[Apex targetting]({{site.baseurl}}#Apex-targetting){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
<!--
[Sample preparation]({{site.baseurl}}#single-cell-sample-preparation){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
-->

&nbsp;

To increase the throughput and quantitative accuracy of single-cell protein analysis by [SCoPE-MS](https://doi.org/10.1101/102681), we introduced many technical improvements in both the sample preparation and in the mass-spectrometry analysis. The [synergistic effect](https://www.biorxiv.org/content/biorxiv/early/2019/12/05/665307/T1.medium.gif) is to increase quantitative accuracy by 4-fold and the throughput of data acquisition about 8-fold. Below, we outline controlled experiments that illustrate the benefits of **individual** improvements. To comprehensively compare the mass-spec data at all levels (including chromatographic separation, precursor abundance, ion isolation, spectral purity, and peptide sequence identification), we include the full [Data-driven Optimization of MS (DO-MS)](https://do-ms.slavovlab.net) reports for each set of experiments. 

&nbsp;


<h2 style="letter-spacing: 2px; font-size: 26px;" id="decreasing-coisolation" >Decreasing co-isolation with narrow isolation windows</h2>

![]({{site.baseurl}}Figures/SCoPE2_Purity_of_MS2_Spectra.png){: width="50%" .center-image}



The degree of coisolation (or purity of MS2 spectra) was evaluated as a function of the width of the isolation window. [MaxQuant](https://www.maxquant.org/) uses the precursor intensity fractions (PIF) to estimate the fraction of the ion intensity originating from the precursor assigned to a peptide spectral match for each MS2 spectrum. The distributions of PIF values for two controlled experiments indicate that the isolation window width used by SCoPE2 (0.7 Th) results in purer spectra compared to the isolation with used with SCoPE-MS. The associated [DO-MS report]({{site.baseurl}}DO-MS_Reports/DO-MS_Report_MS2_isolation_window.html) includes dozens of additional plots demonstrating that these  experiments were well controlled. It also shows that the wider isolation window of 1 Th results in slightly more identified peptides.



&nbsp;

<h2 style="letter-spacing: 2px; font-size: 26px;" id="Apex-targetting" >Improving apex targetting</h2>

![]({{site.baseurl}}Figures/SCoPE-MS__SCoPE2__Apex-offsets.png){: width="42%" }&nbsp;&nbsp;![]({{site.baseurl}}Figures/SCoPE-MS__SCoPE2__PIF.png){: width="42%" }

Distributions of apex offsets for data from [SCoPE-MS](http://scope.slavovlab.net) and from SCoPE2. These offsets are computed by MaxQuant as the time between the apex of the elution peak and the time when the ion is sampled for MS2 analysis. The discreetness of the SCoPE-MS distributions is due to the long duty cycles, and thus low time resolution of sampling the elution peaks. The application of [DO-MS](http://do-ms.slavovlab.net) allowed sampling elution peaks closer to their apecies. This improved apex sampling with SCoPE2 results in much lower coisolation, qantified by the precursor ion fractions (PIF), which MaxQuant computes as estimates for spectral purity, i.e., the fraction of the ion intensity originating from the precursor assigned to a peptide spectral match for each MS2 spectrum.

&nbsp;



<!--
<h2 style="letter-spacing: 2px; font-size: 26px;" id="single-cell-sample-preparation" >Apex </h2>

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

&nbsp;

&nbsp;

&nbsp;
