suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

num_threads <- 12

# set output folder
genome_dir <- "~/Genomes/Reptiles/"

genome_table <- read_tsv("snake_genomes.tsv", col_names = c("species_name", "genome_name"))
j=5

# for(j in c(1:11)){
  
# set variables
species_name <- "Laticauda_colubrina"
genome_name <- "latCor_1.0.fasta"
print(species_name)

# set and read in genome and index
genome_path <- paste0(genome_dir, "/", species_name, "/", genome_name)

genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
gc()
names(genome_seq) <- sub(" .*", "", names(genome_seq))

genome_fai <- tibble(seqnames = names(genome_seq), scaffold_length = as.double(width(genome_seq)))

# set and read in rm fasta and index
search_path <- paste0("RepeatModeler/Laticauda_colubrina/Unknown_curation/Unknown.fasta")

rm_seq <- Biostrings::readDNAStringSet(filepath = search_path)
gc()
names(rm_seq) <- sub("#.*", "", names(rm_seq))

rm_fai <- tibble(seqnames = names(rm_seq), scaffold_length = as.double(width(rm_seq)))

blast <- read_tsv(system(paste0("blastn -query ", search_path, " -db ", genome_path, " -num_threads 12 -outfmt \"6 qseqid qstart qend qlen sseqid sstart send slen length pident bitscore\""), intern = T),
                  col_names = c("qseqid", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "slen", "length", "pident", "bitscore"))
blast

idx <- tibble(seqnames = blast$sseqid, slen = blast$slen) %>%
  base::unique()

flank5 <- 25
flank3 <- 25
i=3
for(i in 1:nrow(rm_fai)){
  
  print(names(rm_seq[i]))
  
  qlen_ <- blast %>% filter(sub("#.*", "", qseqid) == rm_fai$seqnames[i])
  qlen_ <- qlen_$qlen[1]
  
  align_ranges <- blast %>%
    filter(sub("#.*", "", qseqid) == rm_fai$seqnames[i]) %>%
    mutate(strand = ifelse(sstart<send, "+", "-"),
           start = ifelse(sstart<send, sstart, send),
           end = ifelse(sstart>send, sstart, send),
           seqnames = sseqid) %>%
    arrange(-bitscore) %>%
    mutate(start = ifelse(strand == "+", start - flank5, start - flank3),
           end = ifelse(strand == "-", end + flank5, end + flank3)) %>%
    mutate(start = ifelse(start < 1, 1, start), end  = ifelse(end > slen, slen, end)) %>%
    as_granges() %>%
    GenomicRanges::reduce(min.gapwidth = 100L) %>%
    as_tibble() %>%
    mutate(seqnames = as.character(seqnames)) %>%
    inner_join(idx) %>%
    filter(width/qlen_ >= 0.5) %>%
    dplyr::slice(1:30) %>%
    dplyr::select(-width) %>%
    as_granges()
  
  align_seq <- getSeq(genome_seq, align_ranges)
  names(align_seq) <- paste0(seqnames(align_ranges), ":", ranges(align_ranges), "(", strand(align_ranges), ")")
  align_seq <- c(rm_seq[i], align_seq)
  
  writeXStringSet(align_seq, "RepeatModeler/Laticauda_colubrina/Unknown_curation/align_seq.fa")
  system(paste0("mafft --thread 12 --localpair RepeatModeler/Laticauda_colubrina/Unknown_curation/align_seq.fa > RepeatModeler/Laticauda_colubrina/Unknown_curation/", names(rm_seq[i]), "_aligned.fa"))

}
