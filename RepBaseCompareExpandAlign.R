suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

# set output folder
out_folder <- "RepeatModeler/curation/aligned/"
genome_dir <- "~/Genomes/Reptiles/"

# set variables
flank5 <- 500
flank3 <- 500
species_name <- "Aipysurus_laevis"
genome_name <- "kmer_49.pilon_x2.sorted.fasta"
print(species_name)

# set and read in genome and index
genome_path <- paste0(genome_dir, "/", species_name, "/", genome_name)

genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
gc()
names(genome_seq) <- sub(" .*", "", names(genome_seq))

genome_fai <- tibble(seqnames = names(genome_seq), scaffold_length = as.double(width(genome_seq)))

# set and read in RepeatModeler fasta and index

rm_path <- paste0("RepeatModeler/", species_name, "-families.fa")

rm_seq <- Biostrings::readDNAStringSet(filepath = rm_path)
gc()
names(rm_seq) <- sub(" .*", "", names(rm_seq))

rm_fai <- tibble(seqnames = names(rm_seq), scaffold_length = as.double(width(rm_seq)))


# set repbase path

repbase_path <- "~/Databases/RepBase/RepBase24.07_classed.fasta"

rm_repbase_out <- read_tsv(system(paste0("blastn -evalue 0.00002 -num_threads 11 -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -word_size 10 -dust yes -gapopen 30 -gapextend 6 -query ", rm_path, " -db ", repbase_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))

absent <- rm_fai %>%
  filter(!seqnames %in% rm_repbase_out$qseqid) %>%
  arrange(-scaffold_length)

present <- rm_fai %>%
  select(seqnames) %>%
  filter(seqnames %in% rm_repbase_out$qseqid) %>%
  base::unique()

accurate <- rm_repbase_out %>%
  filter(length > 0.7 * qlen | qcovs > 70 ) %>%
  group_by(qseqid) %>%
  dplyr::slice(1) %>%
  arrange(-qlen)

inaccurate <- rm_repbase_out %>%
  filter(!qseqid %in% accurate$qseqid) %>%
  group_by(qseqid) %>%
  dplyr::slice(1) %>%
  select(qseqid, qlen, sseqid) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(start = 1, seqnames = qseqid, qseqid = sub("/", "_", qseqid)) %>%
  dplyr::rename(end = qlen)

inaccurate_bed <- inaccurate %>%
  dplyr::select(seqnames, start, end) %>%
  plyranges::as_granges()

inaccurate_seq <- Biostrings::getSeq(rm_seq, inaccurate_bed)

names(inaccurate_seq) <- inaccurate$qseqid

Biostrings::writeXStringSet(x = inaccurate_seq, filepath = paste0("RepeatModeler/curation/temp.fa"))

inaccurate_blast <- read_tsv(system(paste0("blastn -num_threads 12 -query RepeatModeler/curation/temp.fa -db ", genome_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
  mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
         start = case_when(sstart < send ~ sstart, send < sstart ~ send),
         end = case_when(sstart > send ~ sstart, send > sstart ~ send)) %>%
  inner_join(genome_fai)

# Need to add step to create folders and extend LINEs, DNAs, LTRs and Unknowns by different amounts


for(i in 1:nrow(inaccurate_blast)){
  
  if(i == 1){
    inaccurate_missing <- inaccurate[0,]
  }
  
  inaccurate_blast_bed <- inaccurate_blast %>%
    filter(qseqid == inaccurate$qseqid[i], length > qlen * 0.8, pident > 90) %>%
    arrange(-bitscore)
  
  if(nrow(inaccurate_blast_bed) < 2){
    inaccurate_missing <- rbind(inaccurate_missing, inaccurate[i,])
    next
  } else if(nrow(inaccurate_blast_bed) > 25){
    inaccurate_blast_bed <- inaccurate_blast_bed %>%
      dplyr::slice(1:25)
  }
  
  inaccurate_blast_bed <- inaccurate_blast_bed %>%
    mutate(start = case_when(strand == "+" ~ start - flank5, strand == "-" ~ start - flank3),
           end = case_when(strand == "-" ~ end + flank5, strand == "+" ~ end + flank3),
           start = case_when(start <= 1 ~ 1, start > 1 ~ start),
           end = case_when(end > scaffold_length ~ scaffold_length, end <= scaffold_length ~ end)) %>%
    dplyr::select(seqnames, start, end, strand) %>%
    as_granges()
  
  inaccurate_blast_seq <- Biostrings::getSeq(genome_seq, inaccurate_blast_bed)
  
  names(inaccurate_blast_seq) <- paste0(seqnames(inaccurate_blast_bed), ":", ranges(inaccurate_blast_bed), "(", strand(inaccurate_blast_bed), ")")
  
  inaccurate_blast_seq <- c(inaccurate_seq[i], inaccurate_blast_seq)
  
  Biostrings::writeXStringSet(x = inaccurate_blast_seq, filepath = paste0("RepeatModeler/curation/temp.fa"), append = F)
  
  system(paste0("mafft --localpair --maxiterate 10 --thread 12 RepeatModeler/curation/temp.fa > ", out_folder, species_name, "_", inaccurate$qseqid[i], ".fa"))
  
  
}

# CARP_cdd <- read_tsv(system(paste0("rpstblastn -num_threads 6 -query ", "CARP/", species_name, "/", species_name, "_Tc1.fa", " -db ", cdd_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))