suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

# read cdd codes
cdd_codes <- read_tsv("~/Databases/localrpsb/names_codes_db.tsv", col_names = c("code", "database", "abbreviation"))

# set folder
genome_dir <- "~/Genomes/Reptiles/"

genome_table <- read_tsv("snake_genomes.tsv", col_names = c("species_name", "genome_name"))

species_name <- genome_table$species_name[i]
genome_name <- genome_table$genome_name[i]
print(species_name)

rm_rps <- read_tsv(paste0("RepeatModeler/", species_name, "-families_rps.out"), col_names = c("seqnames", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
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
listed_rps_out <- rm_rps %>%
  mutate(ID = str_remove_all(seqnames, "\\#.+")) %>% split(f = .$ID)

rm_rps %>% dplyr::select(seqnames, abbreviation) %>% filter(grepl("Unknown", seqnames))

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

getTrueWidth(rnd_1_family_304_consensus_out) %>%
  filter(width / 3 > 0.8 * slen)
