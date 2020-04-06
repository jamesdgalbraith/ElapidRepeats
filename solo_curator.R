# read in libraries, hiding messages
suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

# set species and genome names
species_name <- "Laticauda_colubrina"
genome_name <- "latCor_1.0.fasta"

# set/create output directory
if(!dir.exists(paths = paste0("CARP/", species_name, "/curation/"))){
  dir.create(path = paste0("CARP/", species_name, "/curation/"), recursive = T)
}
out_dir <- paste0("CARP/", species_name, "/curation/")

# set genome and carp paths
genome_path <- paste0("~/Genomes/Reptiles/", species_name, "/", genome_name)
carp_path <- paste0("CARP/", species_name, "/", species_name, "_Denovo_TE_Library.fasta")

# read in genome
genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
gc()
names(genome_seq) <- sub(pattern = " .*", replacement = "", x = names(genome_seq))

# read in carp file
carp_seq <- Biostrings::readDNAStringSet(filepath = carp_path)
gc()
names(carp_seq) <- sub(pattern = " .*", replacement = "", x = names(carp_seq))

# create tibble of genome sequence names and lengths
genome_fai <- tibble(seqnames = names(genome_seq), scaffold_length = as.double(width(genome_seq)))

# create tibble of carp sequence names and lengths
carp_fai <- tibble(seqnames = names(carp_seq), start = 1, end = as.double(width(carp_seq)))

# find repeats of interets
harbingers_tbl <- carp_fai %>%
  filter(grepl(pattern = "harbinger", x = seqnames, ignore.case = T), end > 1000) %>%
  arrange(-end)

# set flanks either side
flank5 <- 3500
flank3 <- 3500

# convert repeats of interest to ranges
harbingers_ranges <- harbingers_tbl %>%
  plyranges::as_granges()

# get and name harbingers
harbingers_seq <- Biostrings::getSeq(carp_seq, harbingers_ranges)
names(harbingers_seq) <- seqnames(harbingers_ranges)
Biostrings::writeXStringSet(x = harbingers_seq, filepath = paste0(out_dir, "temp.fa"), append = F)

# compare to self
self_blast_tbl <- read_tsv(system(command = paste0("blastn -query ", out_dir, "temp.fa -subject ", out_dir, "temp.fa -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))

# identify redundant repeats
redundant_tbl <- self_blast_tbl %>%
  dplyr::filter(qseqid != seqnames, pident > 94, qcovs >= 80, qlen < slen) %>%
  dplyr::select(qseqid) %>%
  base::unique()

# remove redundant repeats
nonredundant_harbingers_tbl <- harbingers_tbl %>%
  dplyr::filter(!seqnames %in% redundant_tbl$qseqid)

# convert to ranges
harbingers_ranges <- nonredundant_harbingers_tbl %>%
  plyranges::as_granges()

# get and name nonredundant
harbingers_seq <- Biostrings::getSeq(carp_seq, harbingers_ranges)
names(harbingers_seq) <- seqnames(harbingers_ranges)
Biostrings::writeXStringSet(x = harbingers_seq, filepath = paste0(out_dir, "temp.fa"), append = F)

# blast search for repeat
blast_out <- read_tsv(system(paste0("blastn -num_threads 12 -query ", out_dir, "temp.fa -db ", genome_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))

# manipulate blast output
filtered_blast_out <- blast_out %>%
  dplyr::filter(pident >= 94, qcovs >= 80) %>%
  dplyr::mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"), # determine strand of hit
                start = case_when(strand == "+" ~ sstart, strand == "-" ~ send), # determine stranded start
                end = case_when(strand == "+" ~ send, strand == "-" ~ sstart))

# select and count blast hits
blast_hits <- filtered_blast_out %>%
  dplyr::select(qseqid)

blast_hits <- tibble(seqnames = base::names(base::table(blast_hits$qseqid)), count = as.double(base::table(blast_hits$qseqid))) %>%
  filter(count >= 5) %>%
  inner_join(carp_fai)

i=1
for(i in 1:nrow(blast_hits)){
  
  print(blast_hits$seqnames[i])
  # set query repeat
  query_repeat <- blast_hits$seqnames[i]
  
  # create ranges of query sequence
  query_repeat_ranges <- blast_hits %>%
    filter(seqnames == blast_hits$seqnames[i]) %>%
    plyranges::as_granges()
  
  # get and name query sequence
  query_repeat_seq <- Biostrings::getSeq(carp_seq, query_repeat_ranges)
  names(query_repeat_seq) <- query_repeat
  
  # subest base on query repeat
  query_blast_subset <- filtered_blast_out %>%
    dplyr::filter(qseqid == query_repeat) %>%
    dplyr::arrange(-bitscore)
  
  
  # reduce number of hits if > 20
  if(base::nrow(query_blast_subset) > 20){
    query_blast_subset <- query_blast_subset %>%
      dplyr::slice(1:20)
  }
  
  # adjust and create granges object
  query_blast_subset_ranges <- query_blast_subset %>%
    dplyr::select(seqnames, start, end, strand, slen) %>%
    dplyr::mutate(start = case_when(strand == "+" ~ start - flank5, strand == "-" ~ start - flank3),
                  end = case_when(strand == "+" ~ end + flank3, strand == "-" ~ end + flank5),
                  start = case_when(start > 0 ~ start, start <= 0 ~ 1),
                  end = case_when(end <= slen ~ end, end > slen ~ slen)) %>%
    plyranges::as_granges()
  
  # get sequences from genome
  query_blast_subset_seq <- Biostrings::getSeq(genome_seq, query_blast_subset_ranges)
  
  # name sequences using scaffold_name:cordinates(strand) format
  names(query_blast_subset_seq) <- paste0(seqnames(query_blast_subset_ranges), ":", ranges(query_blast_subset_ranges), "(", strand(query_blast_subset_ranges), ")")
  
  # combine sequences with original query sequence
  query_blast_subset_seq <- c(query_repeat_seq, query_blast_subset_seq)
  
  # write compiled sequences to temporary fasta file
  Biostrings::writeXStringSet(x = query_blast_subset_seq, filepath = paste0(out_dir, "temp.fa"), append = F)
  
  # create alignment of sequences, file named after original query with "/" character replaced with "_"
  system(paste0("mafft --localpair --thread 12 ", out_dir, "temp.fa > ", out_dir, species_name, "_", query_repeat, "_aligned.fa"))
}