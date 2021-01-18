library(tidyverse)
library(GenomicRanges)
library(BSgenome)
library(plyranges)

species_list <- read_tsv("~/Genomes/Birds/bird_genomes.tsv", col_names = c("species_name", "genome_name"))

repeat_name <- "Phylloscopus_trochilus_acredula_family004041"

i=100

blast_search <- read_tsv(system(paste0("blastn -query second_rnd/", repeat_name, ".fasta -db ~/Genomes/Birds/", species_list$species_name[i], "/", species_list$genome_name[i], " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen\" -num_threads 12"), intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "qlen", "slen"))

flank5 = 2000
flank3 = 20

filtered_blast_search <- blast_search %>%
  dplyr::arrange(-bitscore) %>%
  dplyr::filter(pident >= 90) %>%
  dplyr::filter(length >= 1000) %>%
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
  
genome_seq <- readDNAStringSet(paste0("~/Genomes/Birds/", species_list$species_name[i], "/", species_list$genome_name[i]))
names(genome_seq) <- sub( " .*", "", names(genome_seq))

consensus_seq <- readDNAStringSet(paste0("second_rnd/", repeat_name, ".fasta"))

filtered_blast_seq <- getSeq(genome_seq, filtered_blast_ranges)

names(filtered_blast_seq) <- paste0(seqnames(filtered_blast_ranges), ":", ranges(filtered_blast_ranges), "(", strand(filtered_blast_ranges), ")")

filtered_blast_seq <- c(consensus_seq, filtered_blast_seq)

writeXStringSet(filtered_blast_seq, "second_rnd/temp.fa")
  
system(paste0("mafft --localpair --adjustdirection --thread 12 second_rnd/temp.fa > second_rnd/", repeat_name, "_aln.fasta"))

