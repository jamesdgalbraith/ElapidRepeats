# read in libraries, hiding messages
suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

# set genome path
genome_path <- "~/Genomes/Reptiles/Aipysurus_laevis/kmer_49.pilon_x2.sorted.fasta"

# read in genome
genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
gc()
names(genome_seq) <- sub(" .*", "", names(genome_seq))

# create tibble of genome sequence names and lengths
genome_fai <- tibble(seqnames = names(genome_seq), scaffold_length = as.double(width(genome_seq)))

# set path to query repeat
query_repeat <- "Aipysurus_laevis_family148002#LINE_CR1.fasta"

# read in query repeat
query_seq <- Biostrings::readDNAStringSet(filepath = query_repeat)

# set flanks either side
flank5 <- 500
flank3 <- 500

# blast search for repear
blast_out <- read_tsv(system(paste0("blastn -num_threads 4 -query ", query_repeat, " -db ", genome_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
  mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"), # determine strand of hit
         start = case_when(sstart < send ~ sstart, send < sstart ~ send), # determine stranded start
         end = case_when(sstart > send ~ sstart, send > sstart ~ send)) %>% # determine stranded end
  arrange(-bitscore) # sort by bit score

# filter blast search based on covergae of an identidy to query
blast_out <- blast_out %>%
  filter(length > qlen * 0.5, pident > 97)

# reduce number of hits if > 20
if(nrow(blast_out) > 20){
  blast_out <- blast_out %>%
    dplyr::slice(1:20)
}

# extend flanks and create granges object
blast_out_ranges <- blast_out %>%
  mutate(start = case_when(strand == "+" ~ start - flank5, strand == "-" ~ start - flank3), # adjust start
         end = case_when(strand == "-" ~ end + flank5, strand == "+" ~ end + flank3), # adjust end
         start = case_when(start <= 1 ~ 1, start > 1 ~ start), # correct start for less than 1
         end = case_when(end > slen ~ slen, end <= slen ~ end)) %>% # correct end for over scaffold length
  dplyr::select(seqnames, start, end, strand) %>% # select necessary data
  as_granges() # convert to granges

# get sequences from genome
blast_out_seq <- Biostrings::getSeq(genome_seq, blast_out_ranges)

# name sequences using scaffold_name:cordinates(strand) format
names(blast_out_seq) <- paste0(seqnames(blast_out_ranges), ":", ranges(blast_out_ranges), "(", strand(blast_out_ranges), ")")

# combine sequences with original query sequence
blast_out_seq <- c(query_seq, blast_out_seq)

# write compiled sequences to temporary fasta file
Biostrings::writeXStringSet(x = blast_out_seq, filepath = paste0("temp.fa"), append = F)

# create alignment of sequences, file named after original query with "/" character replaced with "_"
system(paste0("mafft --localpair --thread 4 temp.fa > ", sub("/", "_", names(query_seq)), "_aligned.fa"))