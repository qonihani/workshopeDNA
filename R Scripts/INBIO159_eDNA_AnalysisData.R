# #####################################################
# INBIO159 - Metabarcoding dan DNA lingkungan (eDNA) untuk organisme eukariota: Pengenalan dan analisis bioinformatika
# Andhika Prima Prasetyo
# andh009@brin.go.id
# Mohon cantumkan persantunan jika menggunakan bagian dari paparan ini: “Prasetyo (2023) INBIO159 - Metabarcoding dan DNA lingkungan (eDNA) untuk organisme eukariota: Pengenalan dan analisis bioinformatika”
# #####################################################


# Data analysis
# In this section we will briefly discuss a few basic types of multivariate data analysis and data visualisation which are often used in metabarcoding studies. It is important to remember that there are several different ways to plot and analyse your data and we have presented only a few examples here. It might be that for your own data set and questions it would be better to approach your analysis in a different way.

# Analisis data
# Pada bagian ini kita akan membahas secara singkat beberapa jenis dasar analisis data multivariat dan visualisasi data yang sering digunakan dalam studi metabarcoding. Penting untuk diingat bahwa ada beberapa cara yang berbeda untuk memplot dan menganalisis data Anda dan kami hanya menyajikan beberapa contoh di sini. Mungkin saja untuk kumpulan data dan pertanyaan Anda sendiri, akan lebih baik untuk melakukan pendekatan analisis dengan cara yang berbeda.

# If you are concerned about contaminants in your own analysis then it could be useful to explore the R package decontam (https://benjjneb.github.io/decontam/vignettes/decontam_intro.html) which is designed to identify and filter contaminants
# Jika Anda khawatir tentang kontaminan dalam analisis Anda sendiri, maka mungkin berguna untuk mengeksplorasi paket R decontam (https://benjjneb.github.io/decontam/vignettes/decontam_intro.html) yang dirancang untuk mengidentifikasi dan menyaring kontaminan

# Install packages if you haven't installed
# Instal paket jika Anda belum menginstal
install.packages("dplyr")     # To manipulate dataframes
install.packages("readxl")    # To read Excel files into R
install.packages("ggplot2")   # for high quality graphics

# Load packages
# Load paket
library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names

# Read the data and create phyloseq objects
# Three tables are needed create in Excel, namel ASVs, Taxonomy and Samples

# Membaca data dan membuat objek phyloseq
# Tiga tabel yang diperlukan untuk dibuat di Excel, yaitu ASV, Taksonomi, dan Sampel

asv_mat<- read_excel("Tab_ASV.xlsx")
tax_mat<- read_excel("Tab_TAXA.xlsx")
samples_df <- read_excel("Tab_SAMPLE.xlsx")

# Phyloseq objects need to have row.names
# define the row names from the asv column

# Objek Phyloseq harus memiliki row.names
# mendefinisikan nama baris dari kolom asv

asv_mat <- asv_mat %>%
  tibble::column_to_rownames("ASV") 

tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV")

samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") 

# Transform into matrixes otu and tax tables (sample table can be left as data frame)
# Transformasi ke dalam matriks otu dan tabel pajak (tabel contoh dapat dibiarkan sebagai bingkai data)
asv_mat <- as.matrix(asv_mat)
tax_mat <- as.matrix(tax_mat)

# Transform to phyloseq objects
# Mengubah ke objek phyloseq
ASV = otu_table(asv_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

ps <- phyloseq(ASV, TAX, samples)
ps

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 53 taxa and 2 samples ]
# sample_data() Sample Data:       [ 2 samples by 3 sample variables ]
# tax_table()   Taxonomy Table:    [ 53 taxa by 7 taxonomic ranks ]

# Visualize data
# Visualisasi data
sample_names(ps)
rank_names(ps)
sample_variables(ps)

# Normalize number of reads in each sample using median sequencing depth.
# Menormalkan jumlah pembacaan dalam setiap sampel menggunakan kedalaman pengurutan rata-rata.
total = median(sample_sums(ps))
standf = function(x, t=total) round(t * (x / sum(x)))
ps = transform_sample_counts(ps, standf)

# Bar graphs
# Basic bar graph based on Division

# Grafik batang
# Grafik batang dasar berdasarkan Ordo
plot_bar(ps, fill = "Order")

# Make the bargraph nicer by removing ASVs boundaries. This is done by adding ggplot2 modifier.
# Membuat grafik batang menjadi lebih bagus dengan menghapus batas-batas ASV. Hal ini dilakukan dengan menambahkan pengubah ggplot2.
plot_bar(ps, fill = "Order") + 
  geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")

# Heatmaps
# A basic heatmap using the default parameters.
# Heatmaps
# Heatmaps dasar yang menggunakan parameter default.
plot_heatmap(ps, method = "NMDS", distance = "bray")

# It is very very cluttered. It is better to only consider the most abundant ASVs for heatmaps. For example one can only take ASVs that represent at least 5% of reads in at least one sample. Remember we normalized all the samples to median number of reads (total). We are left with only 13 ASVs which makes the reading much more easy.

# Ini sangat berantakan. Lebih baik hanya mempertimbangkan ASV yang paling banyak untuk Heatmaps. Sebagai contoh, kita hanya dapat mengambil ASV yang mewakili setidaknya 5% pembacaan dalam setidaknya satu sampel. Ingatlah bahwa kami menormalkan semua sampel ke jumlah rata-rata pembacaan (total). Kita hanya memiliki 13 ASV yang membuat pembacaan menjadi lebih mudah.
ps_abund <- filter_taxa(ps, function(x) sum(x > total*0.05) > 0, TRUE)
ps_abund

plot_heatmap(ps_abund, method = "NMDS", distance = "bray")

# Alpha diversity
# Plot Chao1 richness estimator and Shannon diversity estimator.
# Keanekaragaman alfa
# Plot penaksir kekayaan Chao1 dan penaksir keragaman Shannon.
plot_richness(ps, measures=c("Chao1", "Shannon"))

# Regroup together samples from the same fraction.
# Kumpulkan kembali sampel dari fraksi yang sama.
plot_richness(ps, measures=c("Chao1", "Shannon"), x="location", color="site")

# Ordination
# Do multivariate analysis based on Bray-Curtis distance and NMDS ordination.
# Ordination
# Lakukan analisis multivariat berdasarkan jarak Bray-Curtis dan ordinasi NMDS.
ps.ord <- ordinate(ps, "NMDS", "bray")

# Plot ASVs
# Visualisasi ASVs
plot_ordination(ps, ps.ord, type="taxa", color="Order", shape= "Family", 
                title="ASVs")

# A bit confusing, so make it more easy to visualize by breaking according to taxonomic division.
# Agak membingungkan, jadi buatlah lebih mudah untuk memvisualisasikannya dengan memecahnya menurut pembagian taksonomi.
plot_ordination(ps, ps.ord, type="taxa", color="Order", 
                title="ASVs", label="Order") + 
                facet_wrap(~Family, 3)

# Now display samples and enlarge the points to make it more easy to read.
# Sekarang tampilkan contoh dan perbesar titik-titik agar lebih mudah dibaca.
plot_ordination(ps, ps.ord, type="samples", color="site", 
                shape="location", title="Samples") + geom_point(size=3)

# Display both samples and ASVs but in 2 different panels.
# Menampilkan sampel dan ASV, tetapi dalam 2 panel yang berbeda.
plot_ordination(ps, ps.ord, type="split", color="Order", 
                shape="location", title="biplot", label = "site") +  
                geom_point(size=3)

# Save RData to be load in the future
# Menyimpan variable dan data untuk bisa dimuat kembali
save.image("INBIO159_eDNA_AnalysisData.RData")
