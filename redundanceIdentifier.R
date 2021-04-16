library(tidyverse)
library(BSgenome)
library(plyranges)

# self align
self_blast <- read_tsv(system(paste0("blastn -query ~/Elapid_Repeats/Laticauda_HT/Laticauda_colubrina_curated.fasta -subject ~/Elapid_Repeats/Laticauda_HT/Laticauda_colubrina_curated.fasta -outfmt \"6 qseqid qstart qend qlen sseqid sstart send slen length pident qcovs\""), intern = T),
                       col_names = c("qseqid", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "slen", "length", "pident", "qcovs"))

# identify redundant sequences
redundant <- self_blast %>%
  dplyr::filter(qseqid != sseqid) %>%
  dplyr::filter(qcovs >= 90, pident >= 99) %>%
  dplyr::filter(qlen < slen) %>%
  dplyr::filter(length/qlen > 0.95) %>%
  dplyr::select(qseqid, sseqid) %>%
  base::unique()

# remove nesyed redundant sequences
redundant <- redundant %>%
  filter(!qseqid %in% redundant$sseqid)

# import curated seq
curated_seq <- readDNAStringSet("~/Elapid_Repeats/Laticauda_HT/Laticauda_colubrina_curated.fasta")
names(curated_seq) <- sub(" .*", "", names(curated_seq))

# export non-redundant curated seq
curated_seq <- curated_seq[!names(curated_seq) %in% redundant$qseqid]
writeXStringSet(curated_seq, "~/Elapid_Repeats/Laticauda_HT/Laticauda_colubrina_curated.fasta")
