#BLAST on Windows di command prompt/powershell
#BLASTN on Windows

#--BLASTN--

makeblastdb -in CO1_Blast.fasta -dbtype nucl -out CO1_Blast.fasta -title "COI_Database"

blastdbcmd -db CO1_Blast.fasta -info

blastn -query ASV_Sequence.fasta -db CO1_Blast.fasta -out taxonomic_assignment.txt -max_target_seqs 1 -perc_identity 97 -outfmt "6 qseqid sseqid"

Merging Taxonomic result with dada2 result in Rstudio

Taxonomy = merge(seqtab.nochim_trans, taxonomic_assignment, by.x = "ASVNumber", by.y = "V1", all.x=TRUE)
head(Taxonomy)
Taxonomy = merge(Taxonomy, CO1_taxonomy, by.x = "V2", by.y = "Accession", all.x=TRUE)

Taxonomy = separate(Taxonomy , col = "Id", into =
                      c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";")

write_xlsx(Taxonomy,"E:/2k23/Lab - Oceanogen Onlince Course Series 2 Batch #1/Workshop data/Hasil.xlsx")

#Taxa of interest
Arthropoda
Bacillariophyta
Chlorophyta
Rotifera
Proteobacteria

