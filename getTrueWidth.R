suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

## Add functions
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
      width = IRanges::reduce(x) %>% width %>% sum(), # this step is the reduction
      database = unique(x$database),
      abbreviation = unique(x$abbreviation), 
      slen = unique(x$slen)
    )
  })
  bind_rows(tibList) # convert list into tibble
  
}

# read cdd codes
cdd_codes <- read_tsv("~/Databases/localrpsb/names_codes_db.tsv", col_names = c("code", "database", "abbreviation"))

# set folder
genome_dir <- "~/Genomes/Reptiles/"

genome_table <- read_tsv("snake_genomes.tsv", col_names = c("species_name", "genome_name"))
i=1
species_name <- genome_table$species_name[i]
genome_name <- genome_table$genome_name[i]
print(species_name)

# rm_rps <- read_tsv(paste0("RepeatModeler/", species_name, "-families_rps.out"), col_names = c("seqnames", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
#   filter(evalue <= 0.01) %>%
#   separate(sseqid, into = c("x", "y", "code")) %>%
#   mutate(strand = case_when(qend > qstart ~ "+",
#                             qend < qstart ~ "-"),
#          start = case_when(qend > qstart ~ qstart,
#                            qend < qstart ~ qend),
#          end = case_when(qend < qstart ~ qstart,
#                          qend > qstart ~ qend),
#          code = as.double(code)) %>%
#   select(-x, -y) %>%
#   inner_join(cdd_codes)

carp_path <- paste0("CARP/", species_name, "/", species_name, "_Denovo_TE_Library.fasta")

carp_seq <- Biostrings::readDNAStringSet(filepath = carp_path)
gc()
names(carp_seq) <- sub(" .*", "", names(carp_seq))

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

# convert rps_out to list
# listed_rm_rps_out <- rm_rps %>%
#   mutate(ID = str_remove_all(seqnames, "\\#.+")) %>% split(f = .$ID)

single_rps <- carp_rps %>% filter(seqnames == "family000201_consensus:Gypsy-8_AMi-I")

listed_carp_rps_out <- carp_rps %>%
  split(f = .$ID)

# rm_rps %>% dplyr::select(seqnames, abbreviation) %>% filter(grepl("Unknown", seqnames))



# t <- lapply(listed_carp_rps_out, getTrueWidth)
t <- purrr::map(listed_carp_rps_out, getTrueWidth)

t2 <- t %>%
  bind_rows() %>%
  mutate(coverage = width / slen) %>%
  filter(coverage > 0.8)
  
t2_giy <- t2 %>%
  filter(abbreviation == "GIY-YIG_PLEs")

t2_rt <- t2 %>%
  filter(abbreviation %in% c("RT_like", "RVT_1", "TERT", "RT_nLTR_like", "RVT_2", "RVT_3"))

t2_en <- t2 %>%
  filter(abbreviation %in% c("Ape1-like_AP-endo", "Ape2-like_AP-endo", "Exo_endo_phos", "Exo_endo_phos_2", "EEP", "EEP-1", "EEP-2", "R1-I-EN", "L1-EN"))

t2 %>%
  filter(seqnames %in% t2_rt$seqnames, seqnames %in% t2_en$seqnames)
  
LINEs <- carp_rps %>%
  filter(seqnames %in% t2_en$seqnames, seqnames %in% t2_rt$seqnames) %>%
  dplyr::select(seqnames, qlen) %>%
  base::unique()

carp_giy <- carp_rps %>%
  filter(abbreviation == "GIY-YIG_PLEs", length / slen > 0.8)

carp_rt <- carp_rps %>%
  filter(abbreviation %in% c("RT_like", "RVT_1", "TERT", "RT_nLTR_like", "RVT_2", "RVT_3")) %>%
  filter(length / slen > 0.8)

carp_en <- carp_rps %>%
  filter(abbreviation %in% c("Ape1-like_AP-endo", "Ape2-like_AP-endo", "Exo_endo_phos", "Exo_endo_phos_2", "EEP", "EEP-1", "EEP-2", "R1-I-EN", "L1-EN"))  %>%
  filter(length / slen > 0.8)

carp_lines <- carp_rps %>%
  filter(seqnames %in% carp_en$seqnames, seqnames %in% carp_rt$seqnames) %>%
  mutate(start = 1, end = qlen) %>%
  select(seqnames, start, end, strand) %>%
  base::unique() %>%
  arrange(seqnames)

carp_lines_ranges <- carp_lines %>%
  plyranges::as_granges()

carp_lines_seq <- Biostrings::getSeq(carp_seq, carp_lines_ranges)

names(carp_lines_seq) <- seqnames(carp_lines_ranges)

Biostrings::writeXStringSet(x = carp_lines_seq, filepath = paste0("CARP/", species_name, "/", species_name, "_Denovo_LINEs.fasta"), append = F)
