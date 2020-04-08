# **SCoPE2**

Single-Cell ProtEomics by Mass Spectrometry

<!--![GitHub release](https://img.shields.io/github/release/SlavovLab/DO-MS.svg)-->
![GitHub](https://img.shields.io/github/license/SlavovLab/DO-MS.svg)

* [Website](https://scope2.slavovlab.net)
* [Running](#running)
* [Download data](https://scope2.slavovlab.net/docs/data)
* [Manuscript](https://www.biorxiv.org/content/10.1101/665307v3)

<img src="http://scope2.slavovlab.net/assets/images/SCOPE2-ac.png" width="70%">



### Requirements

This application has been tested on R >= 3.5.0, OSX 10.14 / Windows 7/8/10. R can be downloaded from the main [R Project page](https://www.r-project.org/) or downloaded with the [RStudio Application](https://www.rstudio.com/products/rstudio/download/). All modules are maintained for MaxQuant >= 1.6.0.16.



### Running

Each thematic part of the SCoPE2_analysis.R script is bounded by an R "region". These can be minimized when not in use. You will need to create "figs","dat", and "code" sub-directories.

To execute the ladder_experiment_analysis.R script you will have to change the path to the ladder experiment evidence file provided on MassIVE, as well as install packages ggpubr and reshape2 from CRAN. The ladder_experiment_design.csv is not necessary for script execution.


------------

## About the project

<!--
DO-MS is a project...


The manuscript for this tool is published at the Journal of Proteome Research: [https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00039](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00039)
-->
The manuscript is freely available on bioRxiv: [Specht et al., 2019](https://www.biorxiv.org/content/10.1101/665307v3).

Contact the authors by email: [nslavov\{at\}northeastern.edu](mailto:nslavov@northeastern.edu).

### License

The SCoPE2 code is distributed by an [MIT license](https://github.com/SlavovLab/DO-MS/blob/master/LICENSE).

### Contributing

Please feel free to contribute to this project by opening an issue or pull request.

<!--
### Data
All data used for the manuscript is available on [UCSD's MassIVE Repository](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=ed5a1ab37dc34985bbedbf3d9a945535)
-->

<!--
### Figures/Analysis
Scripts for the figures in the DART-ID manuscript are available in a separate GitHub repository, [https://github.com/SlavovLab/DART-ID_2018](https://github.com/SlavovLab/DART-ID_2018) 
-->

-------------

## Help!

For any bugs, questions, or feature requests, 
please use the [GitHub issue system](https://github.com/SlavovLab/SCoPE2/issues) to contact the developers.
