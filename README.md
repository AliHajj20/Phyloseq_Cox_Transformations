# Phyloseq_Cox_Transformations
This project aims to evaluate the impact of different microbiome abundance transformation formats on the performance and outputs of an elastic-net penalized Cox proportional hazards model, implemented on a Phyloseq objects.
The workflow is based on, and extends, the methodology proposed by Pujolassos et al. (2024) for microbiome survival analysis.

Install Required Packages To facilitate reproducibility, this project relies on a Docker image that includes many of the required dependencies. Instructions for installing and running the Docker image can be found here: https://microbiome.github.io/OMA/docs/devel/pages/session_info.html#sec-docker-image
Additionally, users should install any remaining R packages listed at the beginning of the script if they are not already installed.

--

## Project Context  
Master’s Thesis Research   
Master’s Program: Molecular Biology and Genetics  
University: University of Pavia  
Location: Turku Data Science Group, Turku, Finland  

--

## Data Description

- **Data Type:** Survival data from a microbiome study  
- **Source:** Public dataset  
- **Data File:** biom.Rdata (included in this repository)  
- **Data Structure:** TreeSummarizedExperiment object  
- **Preprocessing:** No additional data preprocessing is required beyond what is implemented in the analysis scripts.

--

## Methods:
Statistical Model: Cox Proportional Hazards Model  
Evaluation Metric: Harrell's Concordance Index (C-index)  

--

## Analysis Goal

This project builds on and extends the supplementary analysis from the article  
“Microbiome compositional data analysis for survival studies”  
([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11044448/) | [code](https://zenodo.org/records/10554488))  

After independently reproducing the supplementary analysis, the goal of this project is to evaluate the impact of different microbiome abundance transformation formats on the
performance of the Cox proportional hazards model, a workflow(similar to the one in the
original study) was implemented using Phyloseq objects, consistent with the data structure
employed in the original study by [Pujolassos et al. (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11044448/).

### Specifically, this project investigates:
1- The impact of the choice of data structure used to store microbiome data
2- The impact of the choice of statistical transformation of the abundance table

--

## Key differences from the original study
### Statistical transformations:
1- Original study: All pairwise log ratio transformation(APLR)
2- This Project: All pairwise log ratio transformation(APLR) / Isometric log-ratio (ILR) transformation / Additive log-ratio (ALR) transformation / Centered log-ratio (CLR) transformation / Relative abundance / No transformation (raw compositional data)
### functions used to implement the model:
1- Original study: coda_coxnet
2- This Project: the custom function which my colleague and I built: abund_coxnet2
But for APLR transformation I used coda_coxnet

--

This repository includes a custom function, "abund_coxnet2", which performs multiple microbiome abundance transformations (e.g., ILR, ALR, CLR, and relative abundance) based on user input. Then, does Cox’s proportional hazard regression model implementation on the Phyloseq data structure.
This function is defined in a separate script and needs to be sourced before use:

source("abund_coxnet2.R")

For the purpose of reproducing this analysis, first run the source code "abund_coxnet2 so it's saved in the environment.

--

performance assessed using the C-index.
The results obtained represent:
i) signature plot that represents the selected taxa with their respective coefficient
ii) risk score plot that represents the risk score of each taxa
iii) Kaplan–Meier survival can be generated separately using the plot_survival() function, Kaplan–Meier survival curve illustrating the probability of survival among participants, stratified by their microbial risk score. Participants are divided into two groups: those with a high microbial risk score (red) and those with a low microbial risk score (blue). The x-axis represents time, while the y-axis shows the survival probability
