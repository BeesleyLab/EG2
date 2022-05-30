## code to prepare `sysdata` dataset goes here
library(devtools) ; setwd("/working/lab_jonathb/alexandT/EG2/") ; load_all()
sysdata_dir <- "/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/sysdata/"
GENCODE_dir <- paste0(sysdata_dir, "/output/GENCODE/")
GENCODE_filter <- paste0(GENCODE_dir, "filter.gencode.v34lift37.basic.")
coding_mutations_dir <- paste0(sysdata_dir, "/output/coding_mutations/")

# annotation weights
# manually copied from
# https://docs.google.com/spreadsheets/d/1De71B8qUdNge9jrH65GvryX6eYYHgdHy707HtDxyU2k/edit#gid=1075783341
# > data/weights.tsv

# GENCODE
pcENSGs <- read_tibble(paste0(GENCODE_dir, "proteincoding.gencode.v34lift37.basic.ENSGs.txt"))$V1
TSSs <- import_BED(gzfile(paste0(GENCODE_filter, "tss.bed.gz")), metadata_cols = c("ensg", "symbol", "enst"))
promoters <- import_BED(gzfile(paste0(GENCODE_filter, "promoter.bed.gz")), metadata_cols = "enst")
introns <- import_BED(gzfile(paste0(GENCODE_filter, "intron.bed.gz")), metadata_cols = "enst")
exons <- import_BED(gzfile(paste0(GENCODE_filter, "exon.bed.gz")), metadata_cols = "enst")

# ChrSizes
ChrSizes <- read_tibble(paste0(sysdata_dir, "data/hg/hg19.genome"))
names(ChrSizes) <- c("chrom", "size")

# coding_mutations - deleterious coding variants
missense <- paste0(coding_mutations_dir, "missense.tsv") %>% read_tibble(header = T)
nonsense <- paste0(coding_mutations_dir, "nonsense.tsv") %>% read_tibble(header = T)
splicesite <- paste0(coding_mutations_dir, "splicesite.tsv") %>% read_tibble(header = T)

usethis::use_data(pcENSGs,
                  TSSs,
                  promoters,
                  introns,
                  exons,
                  ChrSizes,
                  missense,
                  nonsense,
                  splicesite,
                  internal = TRUE,
                  overwrite = TRUE)
