# MoTrPAC-Pass1B

Template for processing MoTrPAC data internally. The directory structure is loosely based on the [MoTrPAC Proteomics PNNL Prototype](https://github.com/MoTrPAC/motrpac-proteomics-pnnl-prototype).

# Installation

```{r}
library(remotes)
remotes::install_github("vladpetyuk/PlexedPiper", build_vignettes=FALSE)
```

# Usage

The pipeline uses [PlexedPiper](https://github.com/vladpetyuk/PlexedPiper) and requires a connection to the DMS.

```
library(PlexedPiper)

is_PNNL_DMS_connection_successful()
```

The pipeline consists of two scripts: one for global proteomics, and another for PTM-omics (phosphoproteomics, acetylomics, ubiquitinomics). 

## Command line usage

The pipeline can be called from the command line.

```
Rscript process_global_data.R -d 3606 -o data/test_data_global &
```

```
Rscript process_PTM_data.R -d 3625 -g data/test_data_global/results_ratio.txt -o data/test_data_phospho &
```
