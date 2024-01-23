#Running DADA2
#setting directory file
#! jangan lupa di akhir path dikasih string "/" biar workshopeDNA ga kebaca sbg nama 
path <- setwd("C:/Users/Asus/Documents/workshopeDNA/")
list.files(path)

#Nama forward and reverse fastq memiliki format: SAMPLENAME_1.fq and SAMPLENAME_2.fq
#Mendefinisikan bahwa sebelum pattern _1.fq adalah forward dan sebelum _2.fq adalah reverse
fnFs = sort(list.files(path, pattern="_1.fq", full.names = TRUE))
fnRs = sort(list.files(path, pattern="_2.fq", full.names = TRUE))

#Mengekstrak nama sampel dari format: SAMPLENAME_XXX.fastq
sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Melihat kualitas profil FnFs (Forward)
plotQualityProfile(fnFs[1:6])

#Melihat kualitas profil FnRs (Reverse)
plotQualityProfile(fnRs[1:6])

#Menempatkan file yang akan difilter ke dalam sub-directory 
filtFs = file.path(path, "filtered", paste0(sample.names, "_F_filt.fq.gz"))
filtRs = file.path(path, "filtered", paste0(sample.names, "_R_filt.fq.gz"))
names(filtFs) = sample.names
names(filtRs) = sample.names

#Memfilter FnFs dan FnRs
#!! ini yg truncLen=c diganti jadi truncLen=c(200,180) di recording 
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 16, trimRight = 20,truncLen=c(160,180),
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=FALSE) 
head(out)

#Melihat sebaran tingkat error Forward
errF <- learnErrors(filtFs, multithread=TRUE)

#Melihat sebaran tingkat error Reverse
errR <- learnErrors(filtRs, multithread=TRUE)

#Menghasil Amplicon Single Variant dari Forward
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

#Menghasil Amplicon Single Variant dari Reverse
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Melihat ASV dari Forward
dadaFs[[1]]

#Melihat ASV dari Reverse
dadaRs[[1]]

#Menyatukan data Forward dan Reverse
mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap = 4)
head(mergers[[1]])

#Membuat tabel sequence
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Melihat distribusi dari panjang sequence
table(nchar(getSequences(seqtab)))

#Menghapus chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#Melihat perbandingan setelah penghapusan chimeras
sum(seqtab.nochim)/sum(seqtab)


#Meresume hasil filtrasi data sequence
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


write_xlsx(data.frame(track),"C:/Users/Asus/Documents/workshopeDNA/Data/Denoised.xlsx")

write_xlsx(data.frame(seqtab.nochim),"C:/Users/Asus/Documents/workshopeDNA/Data/seqtabnochim.xlsx")

seqtab.nochim_trans <- as.data.frame(t(seqtab.nochim)) %>% 
  rownames_to_column(var = "sequence") %>% 
  rowid_to_column(var = "ASVNumber") %>% 
  mutate(ASVNumber = sprintf("asv%04d", ASVNumber)) %>% 
  mutate(sequence = str_replace_all(sequence, "(-|\\.)", ""))

write_xlsx(seqtab.nochim_trans,"C:/Users/Asus/Documents/workshopeDNA/Data/seqtab.nochim_trans.xlsx")

#Menghasilkan fasta file dari ASV sequence
df <- seqtab.nochim_trans
seq_out <- Biostrings::DNAStringSet(df$sequence)

names(seq_out) <- df$ASVNumber

Biostrings::writeXStringSet(seq_out, str_c(path , "Sequences.fasta"), 
                            compress = FALSE, width = 20000)


#Pilihan untuk running RDP classifier di DADA2
#just an alternative cause the process time could take days...
#taxa <- assignTaxonomy(seqtab.nochim,"C:/Users/Asus/Documents/workshopeDNA/Data/COI_database_RDP Classifier.fasta", multithread=TRUE)




