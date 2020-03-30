suppressMessages(library(tidyverse))
suppressMessages(library(BSgenome))
suppressMessages(library(plyranges))

##### functions #####
getTrueWidth <- function(df){
  reqCols <- c("abbreviation", "database", "slen", "seqnames") # define which columns to keep/are necessary
  stopifnot(reqCols %in% names(df)) # cancel if columns missing
  gr <-  makeGRangesFromDataFrame(df, keep.extra.columns = TRUE) # convert tibble to GRanges
  grl <- split(gr, f = gr$abbreviation) # split GRanges by domain abbreviation
  tibList <- lapply(grl, function(x){ # create tibble of AA coverage
    ## Make sure you retain every single column that you may possibly need
    tibble(
      seqnames = unique(seqlevels(x)),
      width = IRanges::reduce(x) %>% width %>% sum(), # this step is the reduction
      database = unique(x$database),
      abbreviation = unique(x$abbreviation), 
      slen = unique(x$slen)
    )
  })
  bind_rows(tibList) # convert list into tibble
  
}

##### input #####
cdd_codes <- read_tsv("~/Databases/localrpsb/names_codes_db.tsv", col_names = c("code", "database", "abbreviation"))

# set repbase LINE path and read in sequences
repbase_path <- "~/Databases/RepBase/Separate15_2_20/"
repbase_LINEs_seq <- Biostrings::readDNAStringSet(filepath = paste0(repbase_path, "LINEs.fa"))
repbase_LINEs_info <- read_tsv(names(repbase_LINEs_seq), col_names = c("repeat_name", "class", "species_name"))
names(repbase_LINEs_seq) <- sub("\t.*", "", names(repbase_LINEs_seq))
repbase_LINEs_info$width = width(repbase_LINEs_seq)

# set genome folder
genome_dir <- "~/Genomes/Reptiles/"

# set variables
genome_list <- read_tsv("snake_genomes.tsv", col_names = c("species_name", "genome_name"))
j=2
species_name <- genome_list$species_name[j]
genome_name <- genome_list$genome_name[j]
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

# read in rpstbalstn output of all de novo TEs vs custom database
carp_rps <- read_tsv(paste0("CARP/", species_name, "/", species_name, "_rps.out"), col_names = c("seqnames", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
  filter(evalue <= 0.01) %>%
  separate(sseqid, into = c("x", "y", "code")) %>%
  mutate(strand = case_when(qend > qstart ~ "+",
                            qend < qstart ~ "-"),
         start = case_when(qend > qstart ~ qstart,
                           qend < qstart ~ qend),
         end = case_when(qend < qstart ~ qstart,
                         qend > qstart ~ qend),
         code = as.double(code)) %>%
  select(-x, -y) %>%
  inner_join(cdd_codes) %>%
  mutate(ID = sub("#.*", "", seqnames), ID = sub(":.*", "", ID))

##### Step 1 - Get Unknown LINEs based on presence of RT and EN #####
# Get all potential LINEs
unclassified_rps <- carp_rps %>%
  filter(!grepl(":", seqnames), !grepl("Unclassified", seqnames), qlen > 2300)

# create list for each sequence
listed_unclassified_rps_out <- unclassified_rps %>%
  split(f = .$ID)

# run getTrueWidth on all sequences to collapse domains
t <- lapply(listed_unclassified_rps_out, getTrueWidth)
head(listed_unclassified_rps_out)

# collapse list into tibble, filter out coverage < 0.8
t2 <- t %>%
  bind_rows() %>%
  mutate(coverage = width / slen) %>%
  filter(coverage >= 0.8)

# extract all sequences containing reverse transcriptases
carp_rt <- t2 %>%
  filter(abbreviation %in% c("RT_like", "RVT_1", "TERT", "RT_nLTR_like", "RVT_2", "RVT_3"))

carp_rve <- t2 %>%
  filter(abbreviation %in% c("rve")) %>%
  arrange(-width)

# extract all sequences containing endonucleases
carp_en <- t2 %>%
  filter(abbreviation %in% c("Ape1-like_AP-endo", "Ape2-like_AP-endo", "Exo_endo_phos", "Exo_endo_phos_2", "EEP", "EEP-1", "EEP-2", "R1-I-EN", "L1-EN", "GIY-YIG_PLEs"))

# extract all sequences containing both a reverse transcriptase and an endonuclease, transform into ranges ready format
carp_unclassified_lines <- carp_rps %>%
  filter(seqnames %in% carp_en$seqnames, seqnames %in% carp_rt$seqnames) %>%
  mutate(start = 1, end = qlen) %>%
  select(seqnames, start, end, strand) %>%
  base::unique() %>%
  mutate(new_seqnames = gsub("#.*", "#LINE", seqnames), repeat_full_name = gsub("#.*", "", seqnames), family_name = gsub("#.*", "", seqnames)) %>%
  arrange(seqnames)

if(nrow(carp_unclassified_lines) > 1){
# convert tibble into ranges object
carp_unclassified_lines_ranges <- carp_unclassified_lines %>%
  plyranges::as_granges()

# extract and name all unclassified LINE seqs
carp_unclassified_lines_seq <- Biostrings::getSeq(carp_seq, carp_unclassified_lines_ranges)
names(carp_unclassified_lines_seq) <- carp_unclassified_lines_ranges$new_seqnames
} else {
  carp_unclassified_lines_seq <- NA
}

##### Step 2 - Get all known LINEs based on presence of RT and EN #####
# select classifed carp sequences
carp_classified <- carp_fai %>%
  filter(grepl(":", seqnames), scaffold_length >= 2300) %>%
  mutate(repeat_name = sub(".*:", "", seqnames))

# determine strand of all classified sequences
carp_classified <- carp_rps %>%
  inner_join(carp_classified) %>%
  group_by(seqnames) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(start = 1, end = qlen) %>%
  select(seqnames, start, end, strand, repeat_name)

# extract big classified carp LINEs from carp_classified, create name for sequences
carp_classified_lines <- carp_classified %>%
  inner_join(repbase_LINEs_info) %>%
  mutate(family_name = sub(":.*", "", seqnames), new_seqnames = paste0(family_name, "#LINE/", class), start = 1, repeat_full_name = new_seqnames) %>%
  filter(end > 1000) %>%
  dplyr::select(seqnames, start, end, strand, new_seqnames, repeat_full_name, family_name) 

if(nrow(carp_classified_lines) > 1){
# create ranges object of classified LINEs
carp_classified_lines_ranges <- carp_classified_lines %>%
  plyranges::as_granges()

# get sequence of and name classified LINEs
carp_classified_lines_seq <- Biostrings::getSeq(carp_seq, carp_classified_lines_ranges)
names(carp_classified_lines_seq) <- carp_classified_lines_ranges$repeat_full_name
} else {
  carp_classified_lines_seq <- NA
}

##### Step 3 - compare all LINEs
# create tibble of all LINEs
carp_all_lines <- rbind(carp_classified_lines, carp_unclassified_lines)

if(nrow(carp_all_lines) > 1){

# combine classified and unclassified LINEs
if(!is.na(carp_classified_lines_seq[1]) & !is.na(carp_unclassified_lines_seq[1])){
  # both classified and unclassified present
  carp_all_lines_seq <- c(carp_classified_lines_seq, carp_unclassified_lines_seq)
  
} else if(!is.na(carp_classified_lines_seq[1]) & is.na(carp_unclassified_lines_seq[1])){
  # only classified present
  carp_all_lines_seq <- carp_classified_lines_seq
  
} else if(is.na(carp_classified_lines_seq[1]) & !is.na(carp_unclassified_lines_seq[1])){
  # only unclassified present
  carp_all_lines_seq <- carp_unclassified_lines_seq
  
} else{
  
  next()
  
}

# write combined LINEs seq to file
Biostrings::writeXStringSet(x = carp_all_lines_seq, filepath = paste0("CARP/", species_name, "/curation/temp.fa"))

# create blast database of LINEs
system(paste0("makeblastdb -in CARP/", species_name, "/curation/temp.fa -out CARP/", species_name, "/curation/temp.fa -dbtype nucl"))

# blast all vs all LINEs, remove self hits
carp_line_comparison_out <- read_tsv(system(paste0("blastn -num_threads 12 -query CARP/", species_name, "/curation/temp.fa -db CARP/", species_name, "/curation/temp.fa -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("seqnames", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
  filter(seqnames != sseqid)

# remove blast database
system(paste0("rm CARP/", species_name, "/curation/temp.fa.n*"))

# determine divergent hits
carp_divergent_lines <- carp_all_lines %>%
  filter(!new_seqnames %in% carp_line_comparison_out$seqnames)

# determine redundant hits at high ID
carp_redundant <- carp_line_comparison_out %>%
  filter(qlen < slen) %>%
  filter(qcovs > 50) %>%
  filter(pident > 97) %>%
  dplyr::select(seqnames) %>%
  base::unique()

# extract nonredundant LINEs
carp_line_comparison_nonredundant_out <- carp_line_comparison_out %>%
  filter(!seqnames %in% carp_redundant$seqnames) %>%
  group_by(seqnames) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(start = 1, end = qlen) %>%
  dplyr::select(seqnames, sseqid, start, end)

# split into three groups - already classified, newly classified, still unclassified
already_classified_lines <- carp_line_comparison_nonredundant_out %>%
  filter(grepl("/", seqnames))

still_unclassified_lines <- carp_line_comparison_nonredundant_out %>%
  filter(!grepl("/", seqnames), !grepl("/", sseqid))

step_one_lines_ranges <- rbind(already_classified_lines, still_unclassified_lines) %>%
  as_granges()

step_one_lines_seq <- Biostrings::getSeq(carp_all_lines_seq, step_one_lines_ranges)
names(step_one_lines_seq) <- seqnames(step_one_lines_ranges)

newly_classified_lines <- carp_line_comparison_nonredundant_out %>%
  filter(!grepl("/", seqnames), grepl("/", sseqid)) %>%
  mutate(classification = sub(".*\\/", "", sseqid)) %>%
  mutate(new_seqnames = paste0(seqnames, "/", classification))

if(nrow(newly_classified_lines) > 0){
  newly_classified_lines_ranges <- newly_classified_lines %>%
    as_granges()
  newly_classified_seq <- Biostrings::getSeq(carp_all_lines_seq, newly_classified_lines_ranges)
  names(newly_classified_seq) <- newly_classified_lines_ranges$new_seqnames
  carp_for_curation_seq <- c(step_one_lines_seq, newly_classified_seq)
} else {
  carp_for_curation_seq <- step_one_lines_seq
}
Biostrings::writeXStringSet(x = carp_for_curation_seq, filepath = paste0("CARP/", species_name, "/curation/temp.fa"))

carp_for_curation_genome_blast <- read_tsv(system(paste0("blastn -num_threads 12 -query CARP/", species_name, "/curation/temp.fa -db ", genome_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))

carp_for_curation <- carp_for_curation_genome_blast %>%
  select(qseqid) %>%
  mutate(family_name = sub("#.*", "", qseqid)) %>%
  base::unique()

  for(i in 1:nrow(carp_for_curation)){
    
    print(carp_for_curation$qseqid[i])
    
    carp_classified_LINEs_blast_subset <- carp_for_curation_genome_blast %>%
      filter(qseqid == carp_for_curation$qseqid[i], length > qlen * 0.5, pident > 90) %>%
      arrange(-bitscore)
    
    flank5 <- 2000
    flank3 <- 2000
    
    
    if(nrow(carp_classified_LINEs_blast_subset) < 3){
      next
    } else if(nrow(carp_classified_LINEs_blast_subset) > 20){
      carp_classified_LINEs_blast_subset <- carp_classified_LINEs_blast_subset %>%
        dplyr::slice(1:20)
    }
    
    # carp_classified_LINEs_blast_subset <- 
    carp_classified_LINEs_blast_subset_ranges <- carp_classified_LINEs_blast_subset %>%
      mutate(strand = case_when(sstart < send ~ "+", send < sstart ~ "-")) %>%
      mutate(start = case_when(strand == "+" ~ sstart - flank5, strand == "-" ~ send - flank3),
             end = case_when(strand == "-" ~ sstart + flank3, strand == "+" ~ send + flank5),
             start = case_when(start < 1 ~ 1, start >= 1 ~ start),
             end = case_when(end > slen ~ slen, end <= slen ~ end),
             seqnames = sseqid,
             qseqid = gsub(":.*", "", qseqid)) %>%
      dplyr::select(seqnames, start, end, strand, qseqid) %>%
      plyranges::as_granges()
    
    carp_classified_LINEs_blast_subset_ranges_seq <- Biostrings::getSeq(genome_seq, carp_classified_LINEs_blast_subset_ranges)
    
    names(carp_classified_LINEs_blast_subset_ranges_seq) <- paste0(seqnames(carp_classified_LINEs_blast_subset_ranges), ":", ranges(carp_classified_LINEs_blast_subset_ranges), "(", strand(carp_classified_LINEs_blast_subset_ranges), ")")
    
    carp_classified_LINEs_blast_subset_ranges_seq <- c(carp_for_curation_seq[i], carp_classified_LINEs_blast_subset_ranges_seq)
    
    Biostrings::writeXStringSet(x = carp_classified_LINEs_blast_subset_ranges_seq, filepath = paste0("CARP/", species_name, "/curation/temp.fa"), append = F)
    
    system(paste0("mafft --localpair --thread 12 CARP/", species_name, "/curation/temp.fa > CARP/", species_name, "/curation/", species_name, "_", carp_for_curation$family_name[i], ".fa"))
    
  }
}

system("echo $PATH")
