# Rscripts for genomics analysis
## 1) gemma_output_results.R
## Overview
This R script extract results from gemma outputs, generates QQ plots and Manhattan plots, extracts significant (SNPs), and computes relevant statistics from log files.

## Usage
Rscript  gemma_output_results.R <P_wald_threshold>
### Requirements
Make sure you have R installed on your system. Additionally, the script requires the following R packages, which will be installed automatically if not found:
- qqman
- extrafont
## 2) MDS.plots.from_plink.R
### Usage 
Rscript MDS.plots.from_plink.R plink_output.mds
#### Pefrom MDS plots from plink.mds file output after preforming IBD analysis please read from https://zzz.bwh.harvard.edu/plink/strat.shtml
