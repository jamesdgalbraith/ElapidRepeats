suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

# set output folder
genome_dir <- "~/Genomes/Reptiles/"

genome_table <- read_tsv("snake_genomes.tsv", col_names = c("species_name", "genome_name"))

# set variables
species_name <- genome_table$species_name[j]
genome_name <- genome_table$genome_name[j]
print(species_name)
num_threads <- 4

# set and read in genome and index
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

DNA_types <- rm_fai %>%
  filter(grepl("DNA", seqnames)) %>%
  filter(scaffold_length > 500) %>%
  dplyr::select(seqnames) %>%
  mutate(seqnames = sub(".*#", "", seqnames))

DNA_types_table <- tibble(names = names(table(DNA_types$seqnames)), count = table(DNA_types$seqnames))

LINE_types <- rm_fai %>%
  filter(grepl("#LINE", seqnames)) %>%
  filter(scaffold_length > 1000) %>%
  dplyr::select(seqnames) %>%
  mutate(seqnames = sub(".*#", "", seqnames))

LINE_types_table <- tibble(names = names(table(LINE_types$seqnames)), count = table(LINE_types$seqnames))

LTR_types <- rm_fai %>%
  filter(grepl("#LTR", seqnames)) %>%
  filter(scaffold_length > 1000) %>%
  dplyr::select(seqnames) %>%
  mutate(seqnames = sub(".*#", "", seqnames))

LTR_types_table <- tibble(names = names(table(LTR_types$seqnames)), count = table(LTR_types$seqnames))

### Method to extract desired repeats
rm_to_blast <- rm_fai %>%
  filter(grepl("#DNA", seqnames)) %>%
  filter(scaffold_length > 500)
  
rm_to_blast_ranges <- rm_to_blast %>%
  mutate(start = 1, end = scaffold_length, repeat_full_name = seqnames) %>%
  plyranges::as_granges()

rm_to_blast_seq <- Biostrings::getSeq(rm_seq, rm_to_blast_ranges)

names(rm_to_blast_seq) <- seqnames(rm_to_blast_ranges)

if(!dir.exists(paste0("RepeatModeler/", species_name, "/curation"))){dir.create(paste0("RepeatModeler/", species_name, "/curation"), recursive = T)}

Biostrings::writeXStringSet(x = rm_to_blast_seq, filepath = paste0("RepeatModeler/", species_name, "/curation/temp.fa"), append = F)

self_blast <- readr::read_tsv(system(paste0("blastn -evalue 1e-50 -query RepeatModeler/", species_name, "/curation/temp.fa  -subject RepeatModeler/", species_name, "/curation/temp.fa -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))

redundant <- self_blast %>%
  filter(qseqid != seqnames) %>%
  filter(pident >= 94) %>%
  filter(qcovs > 80) %>%
  filter(qlen < slen) %>%
  select(qseqid) %>%
  base::unique()

rm_to_blast <- rm_to_blast %>%
  filter(!seqnames %in% redundant$qseqid) %>%
  mutate(start = 1, end = scaffold_length, repeat_full_name = seqnames) %>%
  plyranges::as_granges()

rm_to_blast_seq <- Biostrings::getSeq(rm_seq, rm_to_blast_ranges)

names(rm_to_blast_seq) <- seqnames(rm_to_blast_ranges)

Biostrings::writeXStringSet(x = rm_to_blast_seq, filepath = paste0("RepeatModeler/", species_name, "/curation/temp.fa"), append = F)

rm_blast <- read_tsv(system(paste0("blastn -evalue 1e-50 -num_threads ", num_threads, " -query RepeatModeler/", species_name, "/curation/temp.fa  -db ", genome_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
  mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
         start = case_when(sstart < send ~ sstart, send < sstart ~ send),
         end = case_when(sstart > send ~ sstart, send > sstart ~ send))
i=1
for(i in 1:length(rm_to_blast_seq)){
  
  print(rm_to_blast_ranges$repeat_full_name[i])
  
  rm_blast_bed <- rm_blast %>%
    filter(qseqid == rm_to_blast_ranges$repeat_full_name[i]) %>%
    filter(length > qlen * 0.5) %>%
    filter(pident > 90) %>%
    arrange(-bitscore)
  
  # if(grepl("#LINE/L1", rm_to_blast_ranges$repeat_full_name[i])){
    # flank5 <- 4500
    # flank3 <- 1500
  # } else {
    flank5 <- 2000
    flank3 <- 2000
  # }
  
  if(nrow(rm_blast_bed) < 3){
    next
  } else if(nrow(rm_blast_bed) > 20){
    rm_blast_bed <- rm_blast_bed %>%
      dplyr::slice(1:20)
  }
  
  rm_blast_bed <- rm_blast_bed %>%
    mutate(start = case_when(strand == "+" ~ start - flank5, strand == "-" ~ start - flank3),
           end = case_when(strand == "-" ~ end + flank5, strand == "+" ~ end + flank3),
           start = case_when(start <= 1 ~ 1, start > 1 ~ start),
           end = case_when(end > slen ~ slen, end <= slen ~ end)) %>%
    dplyr::select(seqnames, start, end, strand) %>%
    as_granges()
  
  rm_blast_seq <- Biostrings::getSeq(genome_seq, rm_blast_bed)
  
  names(rm_blast_seq) <- paste0(seqnames(rm_blast_bed), ":", ranges(rm_blast_bed), "(", strand(rm_blast_bed), ")")
  
  rm_blast_seq <- c(rm_to_blast_seq[i], rm_blast_seq)
  
  Biostrings::writeXStringSet(x = rm_blast_seq, filepath = paste0("RepeatModeler/", species_name, "/curation/temp.fa"), append = F)
  
  system(paste0("mafft --thread ", num_threads, " RepeatModeler/", species_name, "/curation/temp.fa > RepeatModeler/", species_name, "/curation/", species_name, "_", sub("\\/", "_", rm_to_blast_ranges$repeat_full_name[i]), ".fa"))
  
}
