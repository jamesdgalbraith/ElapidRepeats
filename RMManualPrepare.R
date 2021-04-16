suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

num_threads <- 12

# set output folder
genome_dir <- "~/Genomes/Reptiles/"

genome_table <- read_tsv("snake_genomes.tsv", col_names = c("species_name", "genome_name"))

flank5 <- 5000
flank3 <- 5000
# j=5
for(j in c(1:11)){

# set variables
species_name <- genome_table$species_name[j]
genome_name <- genome_table$genome_name[j]
print(species_name)

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

### Method to extract desired repeats
rm_to_blast <- rm_fai
  
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
  select(qseqid, qlen) %>%
  base::unique()

nonredundant_ranges <- rm_to_blast %>%
  filter(!seqnames %in% redundant$qseqid) %>%
  mutate(start = 1, end = scaffold_length,
         names = ifelse(grepl("ltr", seqnames) & !grepl("LTR", seqnames), sub("#.*", "#Unknown", seqnames), seqnames)) %>%
  as_granges()
  
nonredundant_seq <- getSeq(rm_seq, nonredundant_ranges)
names(nonredundant_seq) <- nonredundant_ranges$names

rm_to_blast_seq <- Biostrings::getSeq(rm_seq, nonredundant_ranges)

names(rm_to_blast_seq) <- nonredundant_ranges$names

Biostrings::writeXStringSet(x = rm_to_blast_seq, filepath = paste0("RepeatModeler/", species_name, "/curation/temp.fa"), append = F)

rm_blast <- read_tsv(system(paste0("blastn -task dc-megablast -evalue 1e-50 -num_threads ", num_threads, " -query RepeatModeler/", species_name, "/curation/temp.fa  -db ", genome_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
  mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
         start = case_when(sstart < send ~ sstart, send < sstart ~ send),
         end = case_when(sstart > send ~ sstart, send > sstart ~ send)) %>%
  filter(length >= 0.4*qlen, pident > 80)

rm_blast_hits <- tibble(qseqid = names(table(rm_blast$qseqid)), n = as.integer(table(rm_blast$qseqid))) %>%
  filter(n > 2)

rm_blast_few_hits <- tibble(qseqid = names(table(rm_blast$qseqid)), n = as.integer(table(rm_blast$qseqid))) %>%
  filter(n <= 2)

for(i in 1:length(rm_blast_hits$qseqid)){
  
  print(rm_blast_hits$qseqid[i])
  
  rm_blast_bed <- rm_blast %>%
    filter(qseqid == rm_blast_hits$qseqid[i]) %>%
    arrange(-bitscore)
  
  if(nrow(rm_blast_bed) > 30){
    rm_blast_bed <- rm_blast_bed %>%
      dplyr::slice(1:30)
  }
  
  rm_blast_bed <- rm_blast_bed %>%
    mutate(start = ifelse(strand == "+", start - flank5, start - flank3),
           end = ifelse(strand == "+", end + flank5, end + flank3),
           start = case_when(start <= 1 ~ 1, start > 1 ~ start),
           end = case_when(end > slen ~ slen, end <= slen ~ end)) %>%
    dplyr::select(seqnames, start, end, strand) %>%
    as_granges() %>%
    reduce()
  
  rm_blast_seq <- Biostrings::getSeq(genome_seq, rm_blast_bed)
  
  names(rm_blast_seq) <- paste0(seqnames(rm_blast_bed), ":", ranges(rm_blast_bed), "(", strand(rm_blast_bed), ")")
  
  Biostrings::writeXStringSet(x = rm_blast_seq, filepath = paste0("RepeatModeler/", species_name, "/curation/temp.fa"), append = F)
  
  repeat_self_blast <- read_tsv(system(paste0("blastn -evalue 1e-50 -query RepeatModeler/", species_name, "/curation/temp.fa  -subject RepeatModeler/", species_name, "/curation/temp.fa -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\" -task dc-megablast"), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))
  
  repeat_self_ranges <- repeat_self_blast %>%
    dplyr::filter(qseqid != seqnames)
  
  if(nrow(repeat_self_ranges) <1){
   next()
  }
  
  repeat_self_ranges <- repeat_self_ranges %>%
    dplyr::group_by(qseqid, seqnames) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(qseqid, qstart, qend) %>%
    dplyr::group_by(qseqid) %>%
    dplyr::mutate(qstart = min(qstart), qend = max(qend)) %>%
    dplyr::ungroup() %>%
    base::unique() %>%
    dplyr::rename(seqnames = qseqid, start = qstart, end = qend) %>%
    plyranges::as_granges()
  
  rm_blast_seq <- getSeq(rm_blast_seq, repeat_self_ranges)
  names(rm_blast_seq) <- paste0(seqnames(repeat_self_ranges), "#(", ranges(repeat_self_ranges), ")")
  
  rm_blast_seq <- c(rm_to_blast_seq[names(rm_to_blast_seq) == rm_blast_hits$qseqid[i]], rm_blast_seq)
  
  Biostrings::writeXStringSet(x = rm_blast_seq, filepath = paste0("RepeatModeler/", species_name, "/curation/temp.fa"), append = F)
  
  system(paste0("mafft --localpair --thread ", num_threads, " RepeatModeler/", species_name, "/curation/temp.fa > RepeatModeler/", species_name, "/curation/", species_name, "_", sub("\\/", "_", rm_blast_hits$qseqid[i]), ".fa"))
  
}

# names(nonredundant_seq) <- paste0(species_name, "_", names(nonredundant_seq))
# writeXStringSet(nonredundant_seq, paste0("RepeatModeler/", species_name, "_nonredundant.fasta"))


}
