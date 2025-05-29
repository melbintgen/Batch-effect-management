# Batch-effect-management

**Authors: Yiwen (Eva) Wang**

**Tutors: Xiaochen Zhang**

| Audience      | Prerequisites | Duration    |
| ------------- | ------------- | ----------- |
| Biologists    | [Intro to R](https://melbintgen.github.io/intro-to-r/intro_r_biologists.html)          |~ 2.5 hours    |
| Biologists    | [Intro to Experimental Designs](https://github.com/melbintgen/intro-to-experimental-design)          |~ 3.0 hours    |
| Biologists    | [Intro to Linear Models](https://melbintgen.github.io/intro-to-linear-models/linear_models.html)          |~ 3.0 hours    |

### Description

This repository includes materials for our workshop 'Batch Effect Management'. This workshop introduces commonly encountered sources of batch effects, batch x treatment designs and the scale of batch influence. We will discuss the suitable applications and limitations of existing methods through illustrative case studies. Practical guidelines will also be provided for preprocessed input data, including batch-effect detection and management, and evaluation of method effectiveness through visual and numerical approaches. While our examples are based on microbiome data, the concepts presented in the workshop are applicable to all types of omics data.

### Installation Requirements

Install R first, then RStudio. Download the most recent version of R and RStudio using the links below:
- [R](https://cran.r-project.org/) (Preferably R version > 4.0)
- [RStudio](https://posit.co/download/rstudio-desktop/#download)

Install the R packages.
Type the R command lines:
``` 
# CRAN packages
cran.pkgs <- c('pheatmap', 'vegan', 'ruv', 'ggplot2', 
               'performance', 'gridExtra')

install.packages(cran.pkgs)

# Bioconductor packages
bioc.pkgs <- c('mixOmics', 'sva', 'limma', 'Biobase', 'metagenomeSeq', 
               'PLSDAbatch', 'TreeSummarizedExperiment')

if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(bioc.pkgs)  

# Test if packages have been installed
sapply(c(cran.pkgs, bioc.pkgs), require, character.only = TRUE)

```


### Material

[Click here](https://melbintgen.github.io/intro-to-linear-models/linear-models-master/Linear_model_slides.pdf) to access the slides.

[Click here](https://melbintgen.github.io/intro-to-linear-models/linear-models-master/linear_models.R) to access the Rscript.

[Click here](https://melbintgen.github.io/intro-to-linear-models/linear_models.html) to access the HTML workshop document.


### References
[1] Wang, Y., & Lê Cao, K. A. (2020). [Managing batch effects in microbiome data.](https://academic.oup.com/bib/article/21/6/1954/5643537) Briefings in bioinformatics, 21(6), 1954-1970.
[2] Wang, Y., & Lê Cao, K. A. (2023). [PLSDA-batch: a multivariate framework to correct for batch effects in microbiome data.](https://academic.oup.com/bib/article/24/2/bbac622/6991121) Briefings in Bioinformatics, 24(2), bbac622.
