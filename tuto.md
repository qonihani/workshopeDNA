# eDNA Metabarcoding 

Disclaimer: The code used in this github page is provided by Oceanogen in the Online Course Series. This is just my rendition of the script for easier review in the future. 

### 1. Installing package 
Installing package `install.packages()` (note the s in the word 'packages') is usually only done once at the beginning of every project. But you would need to reactivate the package using the `library()` function everytime you open the project. 

```sh
install.packages("dplyr")
library(dplyr)
```
Package names need to be enclosed in quotation marks (" ") for the `install.packages()` function while it doesn't matter if you use it or not in the `library()` function. Below are the script for installing packages required for metabacoding analysis:

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
