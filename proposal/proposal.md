Machine learning using RNAseq data
================
Alsu

``` r
library(tidyverse)
library(broom)
```

## 1. Introduction

I would like to use the RNAseq dataset for machine learning, to see if
it can predict Pakinson’s disease (PD) based on the transcriptome
profile (in this data set, from the substantia nigra).

The human samples used to create the dataset were provided by the
Parkinson’s UK Brain Bank (Imperial College London, London, UK). A total
of 8 PD (6 males and 2 females; mean age 81.125 and SD of 5.693) and 6
CTRL (1 male and 5 females; mean age 78.166 and SD of 10.381) samples.
All PD cases included in the study are classified at Braak LB Stage 6.
The study adhered to the principles of the Declaration of Helsinki of
1964 and its subsequent amendments. The Ethics Committee of the Oasi
Research Institute—IRCCS of Troina (Italy) granted approval for the
protocol on 5 April 2022 (approval code: 2022/04/05/CE-IRCCS-OASI/52).

RNA was extracted from samples formalin-fixed paraffin-embedded
slide-mounted sections, stored at -80 °C, followed by the RNA sequencing
and data normalization.

The data then was accessed at ArrayExpress
(<https://www.ebi.ac.uk/biostudies/arrayexpress>) for this project.

## 2. Data

Submitter: Genomix4Life Giovanna Ventola

Email: <info@genomix4life.com>

Role: submitter

Affiliation: Genomix4Life s.r.l.

“Gene expression profiling of substantia nigra: a post-mortem study in
Parkinson’s disease and control subjects.” BioStudies, E-MTAB-13295,
2024,
<https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13295>.
Accessed 13 February 2024.

Accosiated article: Salemi M, Ravo M, Lanza G, Schillaci FA, Ventola GM,
Marchese G, Salluzzo MG, Cappelletti G, Ferri R. Gene Expression
Profiling of Post Mortem Midbrain of Parkinson’s Disease Patients and
Healthy Controls. International Journal of Molecular Sciences. 2024;
25(2):707. <https://doi.org/10.3390/ijms25020707>

``` r
setwd("/cloud/project/data/")  

normalized.counts.ibd <- read.table(file="Normalized_Count.csv",
                             sep="",
                             header=T,
                             fill=T,
                             check.names=F)

writeLines(sprintf("%i %s", c(dim(normalized.counts.ibd)[1], dim(normalized.counts.ibd)[2]), c("rows corresponding to transcript IDs", "columns corresponding to samples")))
```

    ## 16148 rows corresponding to transcript IDs
    ## 15 columns corresponding to samples

``` r
# read in the sdrf file (information file that came with dataset, from the original study)
setwd("/cloud/project/data/") 

samp.info.ibd <- read.table(file="E-MTAB-13295.sdrf.txt", sep="\t", header=T, fill=T, check.names=F)

sprintf("There are %i rows, corresponding to the samples", dim(samp.info.ibd)[1])
```

    ## [1] "There are 28 rows, corresponding to the samples"

## 3. Ethics review

The data set is relatively small, including only 6 controls and 8 PD
samples. That could be a potential limitation for training the
algorythm. However, since RNAseq datasets are large due to amount of
genes, that should balance up for the limited computation power of my
laptop and the computing time.

This data set is publicly available, therefore I think that there are no
ethical problems of using it. Additionally, if I achieve intresting
results, it could be potentially of interest to the scientific
community. The data set does not provide any personal information on the
patients from whome it was collected, so there are no negative
consequences for the people. If anything, there could be potential
positive consequences as machine learning for the prediction of certain
diseases based on the RNAseq is very promising area of research.

## 4. Data analysis plan

Supervised machine learning would be performed, potentially on full data
set, or by dividing this data set into two parts (for the learning and
for the testing). Differential expression would be potentially performed
prior to supervised machine learning, to narrow down the number of genes
which can be relevant for the learning process. Alternatively, whole
data set would be used and dimensionality reduction would be performed
prior to the machine learning process.
