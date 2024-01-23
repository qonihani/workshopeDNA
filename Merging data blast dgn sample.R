#Import hasil taxonomy secara manual dari "Import Dataset" di Environment

Hasil_blast = taxonomic_assignment

#Merging data 
Taxonomy = merge(seqtab.nochim_trans, taxonomic_assignment, by.x = "ASVNumber", by.y = "ASV", all.x=TRUE)
head(Taxonomy)

#Import CO1_taxonomy_Blast.txt scr manual di environment 

Referensi_taksa = CO1_taxonomy_Blast

#Merging data lagi 
Taxonomy = merge(Taxonomy, Referensi_taksa, by.x = "Assignment", by.y = "Accession", all.x=TRUE)

Taxonomy = separate(Taxonomy , col = "Id", into =
                      c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep = ";")

write_xlsx(Taxonomy,"C:/Users/Asus/Documents/workshopeDNA/Data/Hasil.xlsx")

#Taxa of interest
Arthropoda
Bacillariophyta
Chlorophyta
Rotifera
Proteobacteria

