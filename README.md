# Phyloseq_Cox_Transformations
This project aims to evaluate the impact of different microbiome abundance transformation formats on the performance and outputs of an elastic-net penalized Cox proportional hazards model, implemented on a Phyloseq objects.
The workflow is based on, and extends, the methodology proposed by Pujolassos et al. (2024) for microbiome survival analysis.

Install Required Packages To facilitate reproducibility, this project relies on a Docker image that includes many of the required dependencies. Instructions for installing and running the Docker image can be found here: https://microbiome.github.io/OMA/docs/devel/pages/session_info.html#sec-docker-image
Additionally, users should install any remaining R packages listed at the beginning of the script if they are not already installed.

--



--
The data on which the analysis is performed can be found in this repository "biom.Rdata".
2-Source Custom Function:
The script includes a custom function abund_coxnet2, which implements a penalized Cox regression model on different data transformations (e.g., ILR, ALR, CLR, and relative abundance). This function is defined in a separate script and needs to be sourced before use:

source("abund_coxnet2.R")
