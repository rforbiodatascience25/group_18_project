# Group 18 Final Project

This is the final project of the group #18 for the DTU course 22160 (R for Bio Data Science).

## Project Contributors

|   username  | student ID |    student name    |
|:-----------:|:----------:|:------------------:|
|   arwynnnn  |   s243548  | Kacper Maciejewski |
| SofiRusso09 |   s252054  |     Sofia Russo    |
|    alewoz   |   s253713  | Aleksandra Wozniak |
|  loayzapre  |   s252608  |   Gabriel Loayza   |
|   DagmarG   |   s252256  |   Dagmar Geevers   |

## Presentation

**View the final presentation [here](https://raw.githack.com/rforbiodatascience25/group_18_project/refs/heads/main/doc/presentation.html).**

## Data source

> This project directly recreates the following study: *[The Epigenetic Modifiers HDAC2 and HDAC7 Inversely Associate with Cancer Stemness and Immunity in Solid Tumors](https://www.mdpi.com/1422-0067/25/14/7841)*, Int. J. Mol. Sci. 2024, 25(14), 7841

### 1. Stemness Index (mRNAsi)  

The cancer stemness index values were retrieved directly from the supplementary materials of the original publication introducing the mRNAsi metric:

Table S1: Stemness Indices Derived for All PanCancer 33 TCGA Cohort, Related to Figure 1; Malta et al., *[Machine Learning Identifies Stemness Features Associated with Oncogenic Dedifferentiation](https://www.cell.com/cell/fulltext/S0092-8674(18)30358-1)*, Cell, 2018. (PanCanAtlas project)

This file was manually downloaded from the article’s supplementary data repository because automatic retrieval is prohibited by the publisher license granted to the DTU access key.

Direct download link cannot be provided by us because one will not propagate the authorization.

### 2. Expression & Clinical Data

Gene expression data (RNA-seq v2 mRNA) and clinical annotations were downloaded automatically from *[cBioPortal](https://docs.cbioportal.org/web-api-and-clients/#r-clients)* using the R API client

## Project organisation

![project organisation chart](https://r4bds.github.io/images/viz_bio_data_science_project_organisation_qmd.png)
