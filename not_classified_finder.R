suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

# read cdd codes
cdd_codes <- read_tsv("~/Databases/localrpsb/names_codes_db.tsv", col_names = c("code", "database", "abbreviation"))

# set folder
genome_dir <- "~/Genomes/Reptiles/"

genome_table <- read_tsv("snake_genomes.tsv", col_names = c("species_name", "genome_name"))

i=6

# for(i in 2:6){
  # set variables
  species_name <- genome_table$species_name[i]
  genome_name <- genome_table$genome_name[i]
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
  
  carp_classified_ranges <- carp_fai %>%
    filter(grepl(":", seqnames)) %>%
    dplyr::arrange(seqnames) %>%
    dplyr::mutate(start = 1, end = scaffold_length) %>%
    plyranges::as_granges()
  
  carp_classified_seq <- Biostrings::getSeq(carp_seq, carp_classified_ranges)
  
  names(carp_classified_seq) <- seqnames(carp_classified_ranges)
  
  Biostrings::writeXStringSet(x = carp_classified_seq, filepath = paste0("CARP/", species_name, "/curation/carp_classified.fasta"))
  
  carp_not_classified <- carp_fai %>%
    filter(!grepl(":", seqnames), !grepl("#Unclassified", seqnames), scaffold_length > 1000, scaffold_length< 12000) %>%
    dplyr::arrange(seqnames)
  
  carp_not_classified_ranges <- carp_not_classified %>%
    mutate(start = 1, end = scaffold_length) %>%
    plyranges::as_granges()
  
  carp_not_classified_seq <- Biostrings::getSeq(carp_seq, carp_not_classified_ranges)
  
  names(carp_not_classified_seq) <- seqnames(carp_not_classified_ranges)
  
  Biostrings::writeXStringSet(x = carp_not_classified_seq, filepath = paste0("CARP/", species_name, "/curation/carp_not_classified.fasta"))
  
  system(paste0("makeblastdb -dbtype nucl -in CARP/", species_name, "/curation/carp_classified.fasta -out CARP/", species_name, "/curation/carp_classified.fasta"))
  
  carp_classified_vs_non_blast <- read_tsv(system(paste0("blastn -num_threads 12 -query CARP/", species_name, "/curation/carp_not_classified.fasta -db CARP/", species_name, "/curation/carp_classified.fasta -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))
  
  carp_actually_classified <- carp_classified_vs_non_blast %>%
    dplyr::filter(pident >= 94, qlen < slen, qcovs > 90) %>%
    dplyr::select(qseqid, seqnames) %>%
    base::unique()
  
  carp_not_classified <- carp_not_classified %>%
    filter(!seqnames %in% carp_actually_classified$qseqid)
  
  carp_not_classified_ranges <- carp_not_classified %>%
    mutate(start = 1, end = scaffold_length) %>%
    plyranges::as_granges()
  
  carp_not_classified_seq <- Biostrings::getSeq(carp_seq, carp_not_classified_ranges)
  
  names(carp_not_classified_seq) <- seqnames(carp_not_classified_ranges)
  
  Biostrings::writeXStringSet(x = carp_not_classified_seq, filepath = paste0("CARP/", species_name, "/curation/carp_not_classified.fasta"))
  
  system(paste0("makeblastdb -dbtype nucl -in CARP/", species_name, "/curation/carp_not_classified.fasta -out CARP/", species_name, "/curation/carp_not_classified.fasta"))
  
  carp_not_classified_self_blast <- read_tsv(system(paste0("blastn -num_threads 12 -query CARP/", species_name, "/curation/carp_not_classified.fasta -db CARP/", species_name, "/curation/carp_not_classified.fasta -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))
  
  carp_not_classified_duplicate <- carp_not_classified_self_blast %>%
    dplyr::filter(qseqid != seqnames, pident >= 94, qlen < slen, qcovs > 90) %>%
    dplyr::select(qseqid, seqnames) %>%
    base::unique()
  
  carp_not_classified <- carp_not_classified %>%
    filter(!seqnames %in% carp_not_classified_duplicate$qseqid)
  
  carp_not_classified_ranges <- carp_not_classified %>%
    mutate(start = 1, end = scaffold_length) %>%
    plyranges::as_granges()
  
  carp_not_classified_seq <- Biostrings::getSeq(carp_seq, carp_not_classified_ranges)
  
  names(carp_not_classified_seq) <- seqnames(carp_not_classified_ranges)
  
  Biostrings::writeXStringSet(x = carp_not_classified_seq, filepath = paste0("CARP/", species_name, "/curation/carp_not_classified.fasta"))
  
  system(paste0("rpsblastn -num_threads 12 -query CARP/", species_name, "/curation/carp_not_classified.fasta -out CARP/", species_name, "/curation/carp_not_classified_rps.out -db ~/Databases/localrpsb/db/Cdd -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""))

# }

carp_not_classified_rps_out <- read_tsv(paste0("CARP/", species_name, "/curation/carp_not_classified_rps.out"), col_names = c("seqnames", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
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
  inner_join(cdd_codes)

# convert rps_out to list
listed_rps_out <- carp_not_classified_rps_out %>%
  mutate(ID = str_remove_all(seqnames, "\\#.+")) %>% split(f = .$ID)
  
# function to reduce
getTrueWidth <- function(df){
  reqCols <- c("abbreviation", "database", "slen", "seqnames") # define which columns to keep/are necessary
  stopifnot(reqCols %in% names(df)) # cancel if columns missing
  gr <-  makeGRangesFromDataFrame(df, keep.extra.columns = TRUE) # convert tibble to GRanges
  grl <- split(gr, f = gr$abbreviation) # split GRanges by domain abbreviation
  tibList <- lapply(grl, function(x){ # create tibble of AA coverage
    ## Make sure you retain every single column that you may possibly need
    tibble(
      seqnames = unique(seqlevels(x)),
      width = reduce(x) %>% width %>% sum(), # this step is the reduction
      database = unique(x$database),
      abbreviation = unique(x$abbreviation), 
      slen = unique(x$slen)
    )
  })
  bind_rows(tibList) # convert list into tibble
  
}

# example of getTrueWidth
family054177_consensus_out <- carp_not_classified_rps_out %>%
  filter(grepl("family054177_consensus", seqnames))

getTrueWidth(family054177_consensus_out) %>%
  filter(width / 3 > 0.8 * slen)


# family054177_consensus_out %>% 
#   makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
#   split(f = mcols(.)$abbreviation) %>%
#   lapply(function(x){
#     tibble(
#       width = reduce(x) %>% width %>% sum(),
#       database = unique(x$database),
#       abbreviation = unique(x$abbreviation), 
#       slen = unique(x$slen)
#     )
#   }
#   ) %>%
#   bind_rows() %>% 
#   mutate(prop = width / (3*slen))
# 
# carp_not_classified %>%
#   filter(!seqnames %in% carp_not_classified_rps_out$seqnames)
# 
# domain_count <- tibble(abbreviation = names(base::table(carp_not_classified_rps_out$abbreviation)), count = as.integer(base::table(carp_not_classified_rps_out$abbreviation))) %>%
#   arrange(-count)
# 
# GIY_YIG <- carp_not_classified_rps_out %>%
#   filter(abbreviation == "GIY-YIG_PLEs")
# 
# penelopes <- carp_not_classified_rps_out %>% filter(seqnames %in% GIY_YIG$seqnames)
# 
# RNase_HI_RT_Ty1 <- carp_not_classified_rps_out %>%
#   filter(abbreviation == "RNase_HI_RT_Ty1")
# 
# copias <- carp_not_classified_rps_out %>% filter(seqnames %in% RNase_HI_RT_Ty1$seqnames)
# 
# RNase_HI_RT_Ty3 <- carp_not_classified_rps_out %>%
#   filter(abbreviation == "RNase_HI_RT_Ty3")
# 
# gypsys <- carp_not_classified_rps_out %>% filter(seqnames %in% RNase_HI_RT_Ty3$seqnames)
# 
# RNase_HI_RT_Bel <- carp_not_classified_rps_out %>%
#   filter(abbreviation == "RNase_HI_RT_Bel")
# 
# bels <- carp_not_classified_rps_out %>% filter(seqnames %in% RNase_HI_RT_Bel$seqnames)
# 
# RNase_HI_RT_DIRS1 <- carp_not_classified_rps_out %>%
#   filter(abbreviation == "RNase_HI_RT_DIRS1")
# 
# dirs <- carp_not_classified_rps_out %>% filter(seqnames %in% RNase_HI_RT_DIRS1$seqnames)
# 
# rve <- carp_not_classified_rps_out %>%
#   filter(abbreviation == "rve")
# 
# DDE_3 <- carp_not_classified_rps_out %>%
#   filter(abbreviation == "DDE_3")
# 
# retropepsins <- carp_not_classified_rps_out %>%
#   filter(grepl("retropepsin", abbreviation))
# 
# lines_en <- carp_not_classified_rps_out %>%
#   filter(abbreviation == "L1-EN" | abbreviation == "Exo_endo_phos" | abbreviation == "Exo_endo_phos_2")	
# 
# lines_rt <- carp_not_classified_rps_out %>%
#   filter(abbreviation == "RT_nLTR_like" | abbreviation == "RVT_1" | abbreviation == "RVT_3")
# 
# lines <- carp_not_classified_rps_out %>% filter(seqnames %in% lines_rt$seqnames, seqnames %in% lines_en$seqnames)
# 
# lines <- lines %>%
#   select(-code, -pident, -database, -qstart, -qend)
# 
# 
# family000205_out <- carp_not_classified_rps_out %>% 
#   filter(grepl("family000205", seqnames))
# 
# 
