
# fishN
### _Random encounter modelling as a viable method to estimate absolute abundance of reef fish_ 
_Methods in Ecology and Evolution_

**Authors**: Cher F Y Chow, Viviana Brambilla, Garrett J Fundakowski, Joshua S Madin, Tiago A Marques, Nina M D Schiettekatte, Andrew S Hoey, Maria Dornelas

**Corresponding author**: Please direct any inquiries or questions to **Cher** via [Github](https://github.com/cherfychow) or [email](mailto:cher.fyc@gmail.com)    
  
This repository contains the data and code used for the analysis and figure generation for _Random encounter modelling: a viable alternative to MaxN for estimating abundance of reef fish from remote underwater videos_. The project implements random encounter modelling with fish video timestamp data and compares it against UVC and MaxN surveys conducted at the same time. 
  
**Disclaimer**: The code may undergo revisions and changes following releases.

## Requirements
All files were written in R 4.3.3. We recommend executing the repository using this version along with the relevant versions of the packages below.  

**Package dependencies:**  
The package dependencies are declared in each file as `require(package name here)`  

- `ape` (5.7.1)
- `beepr` (1.3)
- `ggdist` (3.2.0)
- `iNEXT` (3.0.1)
- `nlstools` (2.0.0)
- `patchwork` (1.1.2)
- `rfishbase` (4.0.0)
- `SpadeR` (0.1.1)
- `tidyverse` (2.0.0)
- `vegan` (2.6.4)
- `worrms` (0.4.3)

## Repository structure
**`data`** : all the data files needed to replicate the analysis workflow  
**`R`** : Contains analysis scripts separated by analysis component/primary data used. Not all scripts are directly used in the pipeline, such as reference scripts used for cleaning. Also includes figure generation code with each analysis component.

## Data
- **`data_LW_LL.csv`** :  length-weight and length-length conversions (if needed)
- **`data_LWab.csv`** :  parameters used in length-weight biomass calculations
- **`data_MaxN.csv`** :  counts from MaxN calculations (downstream output from timestamp data)
- **`ruv_himb_pilot.csv`** :  fish timestamp data
- **`traits_group.csv`** :  aggregation trait classification data
- **`uvc_himb.csv`** :  UVC survey data

## Analysis

- **`01_dataclean.R`** code used to clean up raw entered UVC and RUV data, primarily checking for non-sensicals and taxonomy spelling. Does not need to be run but kept as a record.
- **`02_rest_spsize.R`** Timestamp handling and iterative fitting of random encounter staying time models on fish timestamp data. One model for each species each size class each site.
- **`03_rest_nosize.R`**: The same as above, but not size-class aggregated models. Just one model for each species each site.
- **`04_maxN.R`**: subsampling calculations of MaxN values from timestamp data (so I don't have to process twice)
- **`05_compare_species.R`**: Comparative analyses on species detection/richness criteria, e.g. composition PCoA, species accumulations, species rarefaction.
- **`06_compare_metrics.R`**: Calculating standing stock biomass and biodiversity metrics
- **`07_compare_abund.R`**: Pairwise multiple linear regression on abundance/count estimates.


## Figure generation
Figure code is largely self explanatory. Because there are several small analyses in this paper, figure code is presented at the end of each chunk in the analysis files.
