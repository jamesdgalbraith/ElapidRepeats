library(tidyverse)
library(GenomicRanges)
library(BSgenome)
library(plyranges)

genome_dir <- "~/Genomes/Reptiles/"

genome_table <- read_tsv("snake_genomes.tsv", col_names = c("species_name", "genome_name"))

repeat_name <- "Laticauda_colubrina_ltr-1_family-332#DNA_TcMar-Tc2"

i=5

blast_search <- read_tsv(system(paste0("blastn -query ", repeat_name, ".fasta -db ", genome_dir, "/", genome_table$species_name[i], "/", genome_table$genome_name[i], " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen\" -num_threads 12"), intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "qlen", "slen"))

flank5 = 1000
flank3 = 1000

filtered_blast_search <- blast_search %>%
  dplyr::arrange(-bitscore) %>%
  # dplyr::filter(pident >= 90) %>%
  dplyr::filter(length >= 500) %>%
  dplyr::slice(1:30) %>%
  dplyr::mutate(start = ifelse(sstart < send, sstart, send),
                end = ifelse(sstart > send, sstart, send),
                strand = ifelse(sstart < send, "+", "-")) %>%
  dplyr::rename(seqnames = sseqid) %>%
  dplyr::select(seqnames, start, end, strand, slen) %>%
  base::unique() %>%
  dplyr::mutate(start = ifelse(strand == "+", start - flank5, start - flank3),
                end = ifelse(strand == "+", end + flank3, end + flank5),
                start = ifelse(start < 1, 1, start),
                end = ifelse(end > slen, slen, end))

filtered_blast_ranges <- plyranges::as_granges(filtered_blast_search) %>%
  reduce_ranges_directed()


if(length(filtered_blast_ranges) > 30){
  filtered_blast_ranges <- filtered_blast_ranges[1:30]
  }
  
genome_seq <- readDNAStringSet(paste0(genome_dir, "/", genome_table$species_name[i], "/", genome_table$genome_name[i]))
names(genome_seq) <- sub( " .*", "", names(genome_seq))

consensus_seq <- readDNAStringSet(paste0(repeat_name, ".fasta"))

filtered_blast_seq <- getSeq(genome_seq, filtered_blast_ranges)

names(filtered_blast_seq) <- paste0(seqnames(filtered_blast_ranges), ":", ranges(filtered_blast_ranges), "(", strand(filtered_blast_ranges), ")")

filtered_blast_seq <- c(consensus_seq, filtered_blast_seq)

writeXStringSet(filtered_blast_seq, "temp.fa")
  
system(paste0("mafft --thread 12 temp.fa > ", repeat_name, "_aln.fasta"))

