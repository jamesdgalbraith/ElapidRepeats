library(tidyverse)


search_term <- "curated_LTRs"
query <- "Laticauda_colubrina"
sister <- "Notechis_scutatus"
outgroup <- "Naja_naja"

query_in_query <- read_tsv(paste0("~/ElapidRepeats/RepeatModeler/Laticauda_colubrina/", search_term, "_in_", query, ".out"),
                           col_names = c("qseqid", "sseqid", "qstart", "qend", "qlen", "sstart", "send", "slen", "pident", "length"))
"qseqid sseqid qstart qend qlen sstart send slen pident length"
query_in_sister <- read_tsv(paste0("~/ElapidRepeats/RepeatModeler/Laticauda_colubrina/", search_term, "_in_", sister, ".out"),
                            col_names = c("qseqid", "sseqid", "qstart", "qend", "qlen", "sstart", "send", "slen", "pident", "length"))

query_in_outgroup <- read_tsv(paste0("~/ElapidRepeats/RepeatModeler/Laticauda_colubrina/", search_term, "_in_", outgroup, ".out"),
                              col_names = c("qseqid", "sseqid", "qstart", "qend", "qlen", "sstart", "send", "slen", "pident", "length"))

query_len <- tibble(qseqid = query_in_query$qseqid, qlen = query_in_query$qlen) %>%
  base::unique()

query_in_query_filtered <- query_in_query %>%
  filter(pident >= 85, length >= 100)

query_in_query_filtered_tbl <- tibble(qseqid = names(table(query_in_query_filtered$qseqid)),
                                      n = as.integer(table(query_in_query_filtered$qseqid))) %>%
  inner_join(query_len) %>%
  arrange(-n)

query_in_sister_filtered <- query_in_sister %>%
  filter(pident >= 85, length >= 100)

query_in_sister_filtered_tbl <- tibble(qseqid = names(table(query_in_sister_filtered$qseqid)),
                                       n = as.integer(table(query_in_sister_filtered$qseqid))) %>%
  inner_join(query_len) %>%
  arrange(-n) %>%
  filter(n > 1)

query_in_outgroup_filtered <- query_in_outgroup %>%
  filter(pident >= 85, length >= 100)

query_in_outgroup_filtered_tbl <- tibble(qseqid = names(table(query_in_outgroup_filtered$qseqid)),
                                         n = as.integer(table(query_in_outgroup_filtered$qseqid))) %>%
  inner_join(query_len) %>%
  arrange(-n) %>%
  filter(n > 1)

new_in_query <- query_in_query_filtered_tbl %>%
  filter(!qseqid %in% query_in_sister_filtered_tbl$qseqid,
         !qseqid %in% query_in_outgroup_filtered_tbl$qseqid) %>%
  mutate(sorter1 = sub(".*rnd-", "", sub(".*ltr-1", "0", qseqid)),
         sorter1 = as.integer(sub("_family.*", "", sorter1)),
         sorter2 = as.integer(sub("#.*", "", sub(".*family-", "", qseqid)))) %>%
  arrange(-n)

new_in_siblings <- query_in_query_filtered_tbl %>%
  filter(!qseqid %in% query_in_outgroup_filtered_tbl$qseqid) %>%
  filter(qseqid %in% query_in_sister_filtered_tbl$qseqid)

query_in_outgroup_filtered_tbl %>% arrange(-qlen)

def_new_in_query <- new_in_query %>%
  filter(!qseqid %in% query_in_sister$qseqid,
         !qseqid %in% query_in_outgroup$qseqid) %>%
  filter(qlen >= 1200)

PIF_1 <- query_in_query %>%
  filter(qseqid == new_in_query$qseqid[1]) %>%
  filter(pident >= 85, length >= 100) 

ggplot(PIF_1, aes(pident, length)) + 
  geom_bin2d() + 
  scale_fill_distiller(palette= "Spectral", direction=-1) +
  scale_x_continuous(expand = c(0,0), name = "Pairwise identity (%)") +
  scale_y_continuous(expand = c(0,0), name = "Length") +
  ggtitle(label = PIF_1$qseqid[1])

ggplot(PIF_1, aes(length)) + geom_density()

query_len %>%
  filter(grepl("PIF", qseqid))

long <- PIF_1 %>%
  filter(length >= 5000) %>%
  mutate(seqnames = sseqid,
         start = ifelse(sstart<send, sstart, send),
         end = ifelse(send>sstart, send, sstart),
         strand = ifelse(sstart<send, "+", "-")) %>%
  as_granges()

long_seq <- getSeq(genome_seq, long)

names(long_seq) <- paste0(seqnames(long), ":", ranges(long), "(", strand(long), ")")

writeXStringSet(long_seq, "long_PIF.fasta")

system("makeblastdb -dbtype nucl -in long_PIF.fasta")

long_seq_self_blast <- readr::read_tsv(system(paste0("blastn -evalue 1e-50 -query long_PIF.fasta -db long_PIF.fasta -num_threads 11 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))

long_seq_self_blast_filtered <- long_seq_self_blast %>%
  filter(qseqid != seqnames, length >= 4000)

ggplot(long_seq_self_blast_filtered, aes(pident)) + geom_density()

new_in_query_ranges <- new_in_query %>%
  mutate(seqnames = qseqid, start = 1, end = qlen) %>%
  as_granges()

families <- readDNAStringSet("RepeatModeler/Laticauda_colubrina-families.fa")
names(families) <- sub(" .*", "", names(families))
new_in_query_seq <- getSeq(families, new_in_query_ranges)
names(new_in_query_seq) <- new_in_query_ranges$qseqid
writeXStringSet(new_in_query_seq, "temp_new.fa")

new_in_query_unknown <- new_in_query %>%
  filter(grepl("Unknown", qseqid))


query_in_query


PIF_merged <- query_in_query %>%
  filter(grepl("PIF", qseqid)) %>%
  mutate(seqnames = sseqid,
         start = ifelse(sstart<send, sstart, send),
         end = ifelse(send>sstart, send, sstart),
         strand = ifelse(sstart<send, "+", "-")) %>%
  as_granges() %>%
  reduce_ranges() %>%
  as_tibble()

big_PIF_merged <- PIF_merged %>%
  filter(width > 4000)

big_PIF_locales <- tibble(seqnames = names(table(big_PIF_merged$seqnames)), n = as.integer(table(big_PIF_merged$seqnames))) %>% arrange(-n)
i=2
big_PIF_plotting <- big_PIF_merged %>%
  filter(seqnames == big_PIF_locales$seqnames[i])

ggplot(big_PIF_plotting, aes(start)) + geom_histogram(binwidth = 100000, aes(y=..density..)) + geom_density(aes(colour = "red")) +
  scale_x_continuous(expand = c(0,0), name = "Binned start position of repeat on scaffold (bp)") +
  scale_y_continuous(expand = c(0,0), name = "Number of Harbingers in bin") +
  ggtitle("Distribution of PIF/Harbingers across Laticauda scaffold")

ggplot(big_PIF_plotting, aes(start)) + geom_density() + theme_bw()


