library(tidyverse)

genome_name <- "latCor_1.0.fasta"

system(paste0("awk '{print $1 \" \" $2 \" \" $3 \" \" $4 \" \" $5 \" \" $6 \" \" $7 \" \" $8 \" \" $9 \" \" $10 \" \" $11 \" \" $12 \" \" $13 \" \" $14 \" \" $15 \" \" $16}' self_masker_out/", genome_name, ".out > self_masker_out/", genome_name, ".out.spaced"))

rm_out <- readr::read_delim("self_masker_out/latCor_1.0.fasta.out.spaced", skip = 3, col_names = F, delim = " ")

colnames(rm_out) <- c("sw", "div", "del", "ins", "qseqid", "qstart", "qend", "qremain", "strand", "sseqid", "class", "begin", "send", "end", "ID", "overlap")

self_rm_fixed <- rm_out %>%
  mutate(strand = ifelse(strand == "+", "+", "-"), 
         sstart = ifelse(!grepl("\\(", begin), begin, end),
         sstart = as.integer(sub("\\)", "", sstart)),
         length = send - sstart + 1) %>%
  dplyr::select(qseqid, qstart, qend, sseqid, sstart, send, strand, class, sw, div, del, ins, overlap, length)

system(paste0("awk '{print $1 \" \" $2 \" \" $3 \" \" $4 \" \" $5 \" \" $6 \" \" $7 \" \" $8 \" \" $9 \" \" $10 \" \" $11 \" \" $12 \" \" $13 \" \" $14 \" \" $15 \" \" $16}' all_masker_out/", genome_name, ".out > all_masker_out/", genome_name, ".out.spaced"))

rm_out <- readr::read_delim("all_masker_out/latCor_1.0.fasta.out.spaced", skip = 3, col_names = F, delim = " ")

colnames(rm_out) <- c("sw", "div", "del", "ins", "qseqid", "qstart", "qend", "qremain", "strand", "sseqid", "class", "begin", "send", "end", "ID", "overlap")

all_rm_fixed <- rm_out %>%
  mutate(strand = ifelse(strand == "+", "+", "-"), 
         sstart = ifelse(!grepl("\\(", begin), begin, end),
         sstart = as.integer(sub("\\)", "", sstart)),
         length = send - sstart + 1) %>%
  dplyr::select(qseqid, qstart, qend, sseqid, sstart, send, strand, class, sw, div, del, ins, overlap, length)

self_rm_fixed %>%
  group_by(sseqid) %>%
  mutate(sseqid_sum = sum(length)) %>%
  ungroup() %>%
  dplyr::select(sseqid, class, sseqid_sum) %>%
  base::unique() %>%
  dplyr::arrange(-sseqid_sum)

most_common <- self_rm_fixed %>%
  group_by(sseqid) %>%
  mutate(sseqid_sum = sum(length)) %>%
  ungroup() %>%
  dplyr::select(sseqid, class, sseqid_sum) %>%
  base::unique() %>%
  dplyr::arrange(-sseqid_sum)

temp <- all_rm_fixed %>%
  filter(sseqid == "Laticauda_colubrina_rnd-5_family-745") %>%
  arrange(sstart)


ggplot(temp, aes(y = 1:nrow(temp), yend = 1:nrow(temp),
                 x = sstart, xend = send)) + geom_segment()

