library("phyloseq")     
library("ggplot2")      
library("readxl")       
library("dplyr")        
library("tibble")       
library("ggpubr")       

#Setting directory file dan input data
setwd("E:/Workshop data/Phyloseq")

Abundance <- read_excel("Abundance.xlsx")
Sample <- read_excel("Sample.xlsx")
Taxa <- read_excel("Taxa.xlsx")

#Buat phyloseq dataset
otu_mat <- Abundance
tax_mat <- Taxa
samples_df <- Sample

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("ASVNumber") 

tax_mat <- tax_mat %>%
  tibble::column_to_rownames("ASVNumber")

samples_df <- samples_df %>%
  tibble::column_to_rownames("Sample Id")

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

Workshop <- phyloseq(OTU, TAX, samples)

Workshop

#Pilih taxa dan visualisasi relative abundance di tingkat phylum
Workshop <- Workshop %>% subset_taxa(Phylum   == "Arthropoda" | Phylum   == "Bacillariophyta" | Phylum   == "Chlorophyta" | Phylum   == "Rotifera" | Phylum   == "Proteobacteria")

Workshop_Phylum <- tax_glom(Workshop, taxrank = "Phylum")
subset_taxa(Workshop_Phylum, !is.na(Phylum) & !Phylum %in% c("NA"))
table(tax_table(Workshop_Phylum)[, "Phylum"], exclude = NULL)

Relabund_Phylum <- transform_sample_counts(Workshop_Phylum, function(x) x / sum(x)*100)
otu_table(Relabund_Phylum)[1:5, 1:5]

library("randomcoloR")
colors <- distinctColorPalette(k = 5, altCol = FALSE, runTsne = FALSE)

plot_bar(Relabund_Phylum, fill="Phylum") + 
  geom_bar(aes(fill = Phylum, ), stat="identity", position="stack") +
  labs(x = "Sample Id", y = "Reads Proportion") +
  scale_fill_manual(values = colors)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))


