# Repository name: WNV-seroprevalence
## Description:
This repository contains datasets and R code related to the analysis presented in the publication titled "Quantifying the West Nile virus circulation in the avian host population in Northern Italy". The publication describes a modelling framework used to estimate the West Nile virus (WNV) incidence in the corvid population in the Emilia-Romagna region, which enabled us to eventually assess the fraction of corvids that present anti-WNV antibodies at the end of each epidemiological season.

## Datasets:
1. *corvid_data.xlsx*: number of corvids that were sampled and tested for WNV RNA in organs in each specific subregion and time point. LATIN_NAME: species; WEEK: sampling occurred in WEEK or WEEK+1; YEAR: year of sampling; SUBREGION: location of the trap; TOTAL: number of sampled corvids; WNV: number of WNV-positive corvids.  
2. *corvids_towns.xlsx*: number of corvids sampled in each municipality from 2013 to 2022. TOWN: municipality; TOTAL: number of sampled corvids; WNV: number of WNV-positive corvids; INCIDENCE: WNV/TOTAL; id: unique code assigned to each municipality; long - lat: longitude and latitude of the municipality's centroid; PROVINCE: province of the municipality; SUBREGION: subregion containing the municipality.
3. *mosquito_data.xlsx*: number of mosquito pools tested for WNV RNA presence. WEEK: sampling occurred in WEEK or WEEK+1; YEAR: year of sampling; WNV: number of WNV-positive pools; TOTAL: number of sampled pools.
4. *serology_data.xlsx*: number of corvids, sampled in subregion 2, that were tested for WNV antibody presence. WEEK: sampling occurred in WEEK or WEEK+1; YEAR: year of sampling; TOTAL: number of sampled corvids; IMMUNE: number of corvids classified positive for anti-WNV in serum.  

## R Code:
*main.R*: this R script contains the code used for data preprocessing, bayesian inference, output analysis, and visualization. It is extensively commented to aid understanding and replication of the analyses conducted in the publication.

## Usage: 
### Clone the repository:
```git clone https://github.com/dnaxel/WNV-seroprevalence.git```

### Install required packages: 
Ensure that you have the necessary R packages installed. You can install them using the following command:
```install.packages(c("openxlsx", "RColorBrewer", "ggmap", "osmdata", "raster", "viridis", "deSolve", "mvtnorm", "ggplot2", "cowplot", "latex2exp"))```
### Run the analysis script:
Open analysis_code.R in your preferred R environment (e.g., RStudio) and run the script to replicate the analyses presented in the publication.

## Citation:
If you use the datasets or code provided in this repository for your research, please cite the publication:

Alex De Nardi, Giovanni Marini, Ilaria Dorigatti, Roberto Rosà, Marco Tamba, Luca Gelmini, Alice Prosperi, Francesco Menegale, Piero Poletti, Mattia Calzolari, Andrea Pugliese. 2024. Quantifying the West Nile virus circulation in the avian host population in Northern Italy.

## Contact:
For any questions or inquiries regarding the datasets or code, please contact Alex De Nardi at alex.denardi@unitn.it.

## Acknowledgements:
This study was partially funded by EU grant 874850 MOOD (catalogued as MOOD 087) and by EU funding within the MUR PNRR Extended Partnership initiative on Emerging Infectious Diseases (Project no. PE00000007, INF-ACT) and the VERDI project (EU grant 101045989). The contents of this publication are the sole responsibility of the authors and don't necessarily reflect the views of the European Commission. This work was supported by the Italian Ministry of University and Research (MUR) through the PRIN 2020 project (No. 2020JLWP23) “Integrated Mathematical Approaches to Socio–Epidemiological Dynamics” (CUP: E15F21005420006).

