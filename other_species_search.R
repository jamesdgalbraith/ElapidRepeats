# read in libraries, hiding messages
suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

# set species and genome names
species_name <- "Laticauda_colubrina"
genome_name <- "latCor_1.0.fasta"

genome_dir <- "~/Genomes/Reptiles/"

genome_path <- paste0(genome_dir, "/", species_name, "/", genome_name)

genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
gc()
names(genome_seq) <- sub(" .*", "", names(genome_seq))

genome_fai <- tibble(seqnames = names(genome_seq), scaffold_length = as.double(width(genome_seq)))

# set and read in rm fasta and index
rm_path <- paste0("RepeatModeler/", species_name, "-families.fa")

rm_seq <- Biostrings::readDNAStringSet(filepath = rm_path)
gc()
names(rm_seq) <- sub(" .*", "", names(rm_seq))

rm_fai <- tibble(seqnames = names(rm_seq), scaffold_length = as.double(width(rm_seq)))

self_blast_out <- readr::read_tsv(system(paste0("blastn -evalue 1e-50 -query RepeatModeler/", species_name, "-families.fa  -subject RepeatModeler/", species_name, "-families.fa -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
  filter(qseqid != seqnames)

redundant <- self_blast_out %>%
  filter(qcovs > 75, qlen < slen, pident > 90) %>%
  dplyr::select(qseqid) %>%
  base::unique()

non_redundant_ranges <- rm_fai %>%
  filter(!seqnames %in% redundant$qseqid, scaffold_length > 200) %>%
  mutate(start = 1, end = scaffold_length) %>%
  as_granges()

non_redundant_seq <- Biostrings::getSeq(rm_seq, non_redundant_ranges)
names(non_redundant_seq) <- seqnames(non_redundant_ranges)
Biostrings::writeXStringSet(x = non_redundant_seq, filepath = paste0("RepeatModeler/", species_name, "-nonredundant-families.fa"))

blast_out <- readr::read_tsv(system(paste0("blastn -num_threads 12 -query RepeatModeler/", species_name, "-nonredundant-families.fa -db ", genome_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))

filtered_blast_out <- blast_out %>%
  filter(pident >= 80, qcovs >= 20)

filtered_blast_out_counts <- tibble(qseqid = names(table(filtered_blast_out$qseqid)), count = table(filtered_blast_out$qseqid)) %>%
  filter(count > 20)

further_blast_ranges <- filtered_blast_out_counts %>%
  dplyr::rename(seqnames = qseqid) %>%
  inner_join(rm_fai) %>%
  mutate(start = 1, end = scaffold_length) %>%
  dplyr::select(seqnames, start, end) %>%
  plyranges::as_granges()

s_species_name <- "Laticauda_laticaudata"
s_genome_name <- "latLat_1.0.fasta"

subject_blast_out <- readr::read_tsv(system(paste0("blastn -num_threads 12 -dust yes -query RepeatModeler/", species_name, "-nonredundant-families.fa  -db ", genome_dir, "/", s_species_name, "/", s_genome_name, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))

filtered_subject_blast_out <- subject_blast_out %>%
  filter(pident >= 80, qcovs >= 20)

filtered_subject_blast_out_counts <- tibble(qseqid = names(table(filtered_subject_blast_out$qseqid)), count = table(filtered_subject_blast_out$qseqid))  %>%
  filter(count > 20)

filtered_blast_out_counts %>%
  filter(!qseqid %in% filtered_subject_blast_out_counts$qseqid) %>%
  arrange(-count) %>%
  dplyr::rename(seqnames = qseqid) %>%
  inner_join(rm_fai)

