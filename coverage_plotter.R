library(tidyverse)
library(Biostrings)



# set query and genome
genome_path <- "~/Genomes/Reptiles/Aipysurus_laevis/kmer_49.pilon_x2.sorted.fasta"
query_repeat <- "Aipysurus_laevis_family148002#LINE_CR1.fasta"

# perform blastn search
blastn <- read_tsv(system(paste0("blastn -dust yes -query ", query_repeat, " -db ", genome_path, " -outfmt \"6 sseqid qseqid qstart qend pident qcovs bitscore length mismatch evalue qlen\""), intern = TRUE), col_names = c("seqnames", "qseqid",  "qstart", "qend", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen")) %>%
  as_tibble() %>%
  dplyr::filter(length >= 50) %>%
  mutate(div = 100 - pident, d = mismatch/length, jc_dist = (-3 / 4) * log(1 - (4 * d / 3))) %>%
  arrange(-bitscore)

# determine length of query repeat
query_width <- blastn$qlen[1]

# select near full length repeats
blastn_fl <- blastn %>%
  filter(length >= qlen * 0.9)

# calculate coverage of each bp
bp_coverage <- matrix(rep(0, length(blastn$seqnames)*as.numeric(query_width)), byrow = T, ncol = as.numeric(query_width))
for(k in 1:length(blastn$seqnames)){
  bp_coverage[k,]<-c(rep(0,blastn$qstart[k]-1),rep(1,abs(blastn$qend[k]-blastn$qstart[k])+1), rep(0,as.numeric(query_width)-blastn$qend[k]))
}
# convert matrix to tibble
bp_coverage_tbl <- tibble(bp_x = 1:query_width, bp_y = colSums(bp_coverage)) %>%
  mutate(bp_y_scaled = max(blastn$div) * bp_y / max(bp_y))

# plot coverage
ggplot() +
  geom_segment(data = blastn, aes(x = qstart, xend = qend, y = 100 - pident, yend = 100 - pident), colour = "black") + # plot all repeats
  geom_segment(data = blastn_fl, aes(x = qstart, xend = qend, y = 100 - pident, yend = 100 - pident), colour = "#db7800", size = 1.5) + # plot full length repeats
  geom_line(data = bp_coverage_tbl, aes(bp_x, bp_y_scaled), colour = "#006edb") + # plot coverage data
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . * max(bp_coverage_tbl$bp_y)/3, name = "coverage"), # add second axis
                     limits = c(-(max(blastn$div) *0.01),(max(blastn$div) *1.01)),  # scale for second axis
                     expand = c(0,0)) + # don't expand axis
  scale_x_continuous(expand = c(0,0)) + # don't expand axis
  theme_bw() +
  labs(x = "bp", y = "pairwise divergence from consensus", title = blastn$qseqid[1])

