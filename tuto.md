# eDNA Metabarcoding 

Disclaimer: The code used in this github page is provided by Oceanogen in the Online Course Series. This is just my rendition of the script for easier review in the future. 

### 1. Installing package 
Installing package ` usually is only done once at the beginning of every project. But you would need to activate the package using the `library()` function. 

```sh
install.packages("dplyr")
library(dplyr)
```

```sh
#Instalasi DADA2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.16")

#Instalasi phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

#Aktivasi library
library(dada2)
#Aktivasi dada2
packageVersion("dada2")

#Instalasi library lainnya
install.packages("stringr")

library("phyloseq")
library("Biostrings")
library("ggplot2")
library("dplyr")
library("tidyr")
library("tibble")
library("readxl")
library("readr")
library("stringr")
library("writexl")
```

DADA2 is yadayadayada 
