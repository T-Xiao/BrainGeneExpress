# BrainGeneExpress

An RShiny app for quick analysis of Brain Tumor public datasets 

## Description

* The goal of this app is to help Brain Tumor researchers do quick research of the gene(s) in public Brain Tumor datasets.
* Datasets include: (Datasets can be found using the link in the Data Source section) 
> TCGA GBM Agilent data(585 patients, 17814 genes, clinical factors, 17 patients with Event & Time to event info missing; <br/>
> TCGA GBM U133a Affymetrix data(539 patients, 12042 genes, 11 clinical factors, 16 patients with Event & Time to event info missing); <br/>
> TCGA GBM RNA seq data (172 patients, 20028 genes, 11 clinical factors, 7 patients with Event & Time to event info missing); <br/>
> TCGA LGG  data（525 patients, 12042 genes, 11 clinical factors, 16 patients with Event & Time to event info missing); <br/>
> IvyGAP data (270 patients, 25873 genes, no clinical factors); <br/>
> TCGA - The Cancer Genome Atlas;  GBM - Glioblastoma; LGG - Low Grade Glioma<br/>
* Modules include: 
> (1) Correlation: <br/> 
> > (1.1) One vs Many: Pearson Correlations between one gene and several other genes of choice. Correlation plots will be shown as several one-on-one correlations. <br/>
> > (1.2) One vs All: Correlations of the chosen gene and all genes are calculated and sorted. The n (user can determine n) most correlated genes will be clustered and displayed in heatmap. <br/>
> > (1.3) Two sets: Canonical correlations will be calculated for <br/> 

> (2) Cox Regression:<br/>
>(3) Survival Tree: 



## Getting Started

### Dependencies

* R and RStudio

### Installing

* Download R: Windows -- https://cran.r-project.org/bin/windows/base/. 
* Download R: Mac -- https://cran.r-project.org/bin/macosx/
* Download RStudio Desktop: https://www.rstudio.com/products/rstudio/download/

## Authors

Ting Xiao

## Version History

* Initial Release

## License

TBD

## Data Sources
* [TCGA data from Xenabrowser](https://xenabrowser.net/datapages/?cohort=TCGA%20Glioblastoma%20(GBM)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
* [IvyGAP](https://glioblastoma.alleninstitute.org/static/download.html)
