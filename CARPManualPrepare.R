suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

# set output folder
genome_dir <- "~/Genomes/Reptiles/"

# set variables
species_name <- "Notechis_scutatus"
genome_name <- "TS10Xv2-PRI.fasta"
print(species_name)

# set and read in genome and index
genome_path <- paste0(genome_dir, "/", species_name, "/", genome_name)

genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
gc()
names(genome_seq) <- sub(" .*", "", names(genome_seq))

genome_fai <- tibble(seqnames = names(genome_seq), scaffold_length = as.double(width(genome_seq)))

# set and read in CARP fasta and index

carp_path <- paste0("CARP/", species_name, "/", species_name, "_Denovo_TE_Library.fasta")

carp_seq <- Biostrings::readDNAStringSet(filepath = carp_path)
gc()
names(carp_seq) <- sub(" .*", "", names(carp_seq))

carp_fai <- tibble(seqnames = names(carp_seq), scaffold_length = as.double(width(carp_seq)))

carp_classified <- carp_fai %>%
  filter(grepl(":", seqnames)) %>%
  mutate(repeat_name = sub(".*:", "", seqnames))

# set repbase LINE path
repbase_path <- "~/Databases/RepBase/Separate15_2_20/"

repbase_LINEs_seq <- Biostrings::readDNAStringSet(filepath = paste0(repbase_path, "LINEs.fa"))

repbase_LINEs_info <- read_tsv(names(repbase_LINEs_seq), col_names = c("repeat_name", "class", "species_name"))

names(repbase_LINEs_seq) <- sub("\t.*", "", names(repbase_LINEs_seq))

repbase_LINEs_info$width = width(repbase_LINEs_seq)
  
carp_classified_LINEs <- carp_classified %>%
  inner_join(repbase_LINEs_info) %>%
  mutate(family_name = sub(":.*", "", seqnames), repeat_full_name = paste0(family_name, "#LINE/", class), start = 1, end = scaffold_length) %>%
  filter(scaffold_length > 1000)

carp_classified_LINEs_bed <- carp_classified_LINEs %>%
  dplyr::select(seqnames, start, end, repeat_full_name, family_name) %>%
  plyranges::as_granges()

carp_classified_LINEs_seq <- Biostrings::getSeq(carp_seq, carp_classified_LINEs_bed)

names(carp_classified_LINEs_seq) <- seqnames(carp_classified_LINEs_bed)

Biostrings::writeXStringSet(x = carp_classified_LINEs_seq, filepath = paste0("CARP/", species_name, "/curation/temp.fa"))

carp_repbase_out <- read_tsv(system(paste0("blastn -evalue 0.00002 -num_threads 11 -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -word_size 10 -dust yes -gapopen 30 -gapextend 6 -query CARP/", species_name, "/curation/temp.fa -db ", repbase_path, "LINEs.fa -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("seqnames", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))

carp_repbase_strand <- carp_repbase_out %>%
  dplyr::group_by(seqnames) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+")) %>%
  select(seqnames, strand)

carp_classified_LINE_strand_bed <- carp_classified_LINEs %>%
  inner_join(carp_repbase_strand) %>%
  dplyr::select(seqnames, start, end, repeat_full_name, family_name, strand) %>%
  plyranges::as_granges()
  
carp_classified_LINEs_seq <- Biostrings::getSeq(carp_seq, carp_classified_LINE_strand_bed)

names(carp_classified_LINEs_seq) <- carp_classified_LINE_strand_bed$repeat_full_name

Biostrings::writeXStringSet(x = carp_classified_LINEs_seq, filepath = paste0("CARP/", species_name, "/curation/temp.fa"))

carp_classified_LINEs_blast <- read_tsv(system(paste0("blastn -num_threads 12 -query CARP/", species_name, "/curation/temp.fa -db ", genome_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
  mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
                  start = case_when(sstart < send ~ sstart, send < sstart ~ send),
                  end = case_when(sstart > send ~ sstart, send > sstart ~ send))
 
i=1 

names(carp_classified_LINEs_seq)

for(i in 1:length(carp_classified_LINEs_seq)){

  print(carp_classified_LINE_strand_bed$repeat_full_name[i])
  
    carp_classified_LINEs_blast_bed <- carp_classified_LINEs_blast %>%
    filter(qseqid == carp_classified_LINE_strand_bed$repeat_full_name[i], length > qlen * 0.8, pident > 90) %>%
    arrange(-bitscore)
    
    if(grepl("#LINE/L1", carp_classified_LINE_strand_bed$repeat_full_name[i])){
      flank5 <- 4500
      flank3 <- 2000
    } else {
      flank5 <- 3000
      flank3 <- 2000
    }

  if(nrow(carp_classified_LINEs_blast_bed) < 3){
    next
  } else if(nrow(carp_classified_LINEs_blast_bed) > 5){
    carp_classified_LINEs_blast_bed <- carp_classified_LINEs_blast_bed %>%
      dplyr::slice(1:25)
  }

  carp_classified_LINEs_blast_bed <- carp_classified_LINEs_blast_bed %>%
    mutate(start = case_when(strand == "+" ~ start - flank5, strand == "-" ~ start - flank3),
           end = case_when(strand == "-" ~ end + flank5, strand == "+" ~ end + flank3),
           start = case_when(start <= 1 ~ 1, start > 1 ~ start),
           end = case_when(end > slen ~ slen, end <= slen ~ end)) %>%
    dplyr::select(seqnames, start, end, strand) %>%
    as_granges()

  carp_classified_LINEs_blast_seq <- Biostrings::getSeq(genome_seq, carp_classified_LINEs_blast_bed)

  names(carp_classified_LINEs_blast_seq) <- paste0(seqnames(carp_classified_LINEs_blast_bed), ":", ranges(carp_classified_LINEs_blast_bed), "(", strand(carp_classified_LINEs_blast_bed), ")")

  carp_classified_LINEs_blast_seq <- c(carp_classified_LINEs_seq[i], carp_classified_LINEs_blast_seq)

  Biostrings::writeXStringSet(x = carp_classified_LINEs_blast_seq, filepath = paste0("CARP/", species_name, "/curation/temp.fa"), append = F)

  system(paste0("mafft --maxiterate 5 --thread 12 CARP/", species_name, "/curation/temp.fa > CARP/", species_name, "/curation/", species_name, "_", carp_classified_LINE_strand_bed$family_name[i], ".fa"))

}

# 
# # 
# # rm_repbase_out <- read_tsv(system(paste0("blastn -evalue 0.00002 -num_threads 11 -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -word_size 10 -dust yes -gapopen 30 -gapextend 6 -query ", rm_path, " -db ", repbase_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))
# # 
# # absent <- rm_fai %>%
# #   filter(!seqnames %in% rm_repbase_out$qseqid) %>%
# #   arrange(-scaffold_length)
# # 
# # present <- rm_fai %>%
# #   select(seqnames) %>%
# #   filter(seqnames %in% rm_repbase_out$qseqid) %>%
# #   base::unique()
# # 
# # accurate <- rm_repbase_out %>%
# #   filter(length > 0.7 * qlen | qcovs > 70 ) %>%
# #   group_by(qseqid) %>%
# #   dplyr::slice(1) %>%
# #   arrange(-qlen)
# # 
# # inaccurate <- rm_repbase_out %>%
# #   filter(!qseqid %in% accurate$qseqid) %>%
# #   group_by(qseqid) %>%
# #   dplyr::slice(1) %>%
# #   select(qseqid, qlen, sseqid) %>%
# #   dplyr::ungroup() %>%
# #   dplyr::mutate(start = 1, seqnames = qseqid, qseqid = sub("/", "_", qseqid)) %>%
# #   dplyr::rename(end = qlen)
# # 
# # inaccurate_LINEs <- inaccurate %>%
# #   filter(grepl("#LINE", qseqid)) %>%
# #   select(qseqid) %>%
# #   base::unique()
# # 
# # inaccurate_bed <- inaccurate %>%
# #   dplyr::select(seqnames, start, end) %>%
# #   plyranges::as_granges()
# # 
# # inaccurate_seq <- Biostrings::getSeq(rm_seq, inaccurate_bed)
# # 
# # names(inaccurate_seq) <- inaccurate$qseqid
# # 
# # Biostrings::writeXStringSet(x = inaccurate_seq, filepath = paste0("RepeatModeler/curation/temp.fa"))
# # 
# # inaccurate_blast <- read_tsv(system(paste0("blastn -num_threads 12 -query RepeatModeler/curation/temp.fa -db ", genome_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
# #   mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
# #          start = case_when(sstart < send ~ sstart, send < sstart ~ send),
# #          end = case_when(sstart > send ~ sstart, send > sstart ~ send)) %>%
# #   inner_join(genome_fai)
# # 
# # # Need to add step to create folders and extend LINEs, DNAs, LTRs and Unknowns by different amounts
# # 
# # 
# # 
# # 