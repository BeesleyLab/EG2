EG2 Tour
================
Alexandra Tidd
May 31, 2022

## predict\_target\_genes()

EG2 is a package to predict the target genes of fine-mapped variants of a trait from GWAS.

`predict_target_genes()` is the master, user-facing function of this package. All other functions are helper functions called by `predict_target_genes()`.

``` r
?predict_target_genes()
```

### Installation

In order to install and run this package yourself...

``` r
# 1. Install and load the devtools package
install.packages("devtools")
library(devtools)
# 2. Install and load the EG2 package
install_github("BeesleyLab/EG2")
library(EG2)
```

### Package data

This package uses both reference genomic annotation datasets and user-provided trait-specific datasets. Genomic coordinates use the hg19 build. The package's BED handling follows UCSC BedTools conventions, so it expects 0-based start positions and 1-based end positions.

#### Reference data

##### Internal reference data

Smaller generic reference datasets, including default annotation weights, chromosome sizes, GENCODE annotations and coding mutation annotations (`default_weights`, `ChrSizes`, `TSSs`, `exons`, `introns`, `promoters`, `missense`, `nonsense`, `splicesite`) are stored internally as parsed objects in `R/sysdata.R`. They are accessible when the package is loaded, but not visible due to lazy loading. The scripts used to generate the input files are published at <https://github.com/alextidd/tgp_paper/tree/main/wrangle_package_data/sysdata/code>.

##### External reference data

Large cell type-specific reference data that accompany the EG2 package have been published to OSF (<https://osf.io/254nq/>) and must be downloaded separately. In order to use the package, first download and unzip the data into a local directory...

``` bash
mkdir path/to/EG2_data/
wget -O EG2_data.zip https://osf.io/254nq/download --no-check-certificate
unzip -d path/to/EG2_data/ EG2_data.zip
```

...and then pass the directory to the `reference_panels_dir` argument of the `predict_target_genes()` function.

``` r
predict_target_genes(reference_panels_dir = "path/to/EG2_data/", ...)
```

The scripts to generate these data are published at <https://github.com/alextidd/tgp_paper/tree/main/wrangle_package_data/reference_panels/code>.

#### User-provided data

There is one required user-provided file for the `predict_target_genes()` function: the trait variants (the `variants_file` argument). Known genes for the trait (the `known_genes_file` argument) are needed only if `do_performance = T`. Example files for breast cancer are in `inst/example_data/`.

``` r
system.file("example_data", "variants.bed", package = "EG2")
```

    ## [1] "/mnt/backedup/home/alexandT/R/x86_64-pc-linux-gnu-library/4.0/EG2/example_data/variants.bed"

``` r
system.file("example_data", "known_genes.txt", package = "EG2")
```

    ## [1] "/mnt/backedup/home/alexandT/R/x86_64-pc-linux-gnu-library/4.0/EG2/example_data/known_genes.txt"

The scripts to generate these and all other trait files used in the study are published at <https://github.com/alextidd/tgp_paper/tree/main/wrangle_package_data/traits/code>.

##### Trait variants

The variants file should be a BED file with metadata columns for the variant name and the credible set to which it belongs.

``` r
variants_file <- system.file("example_data", "variants.bed", package = "EG2")
EG2::import_BED(variants_file, metadata_cols = c("variant", "cs"))
```

    ## # A tibble: 5,375 × 5
    ##    chrom    start      end variant     cs         
    ##    <chr>    <int>    <int> <chr>       <chr>      
    ##  1 chr1  10551762 10551763 rs657244    BCAC_FM_1.1
    ##  2 chr1  10563363 10563364 rs202087283 BCAC_FM_1.1
    ##  3 chr1  10564674 10564675 rs2847344   BCAC_FM_1.1
    ##  4 chr1  10566521 10566522 rs617728    BCAC_FM_1.1
    ##  5 chr1  10569000 10569000 rs60354536  BCAC_FM_1.1
    ##  6 chr1  10569257 10569258 rs2480785   BCAC_FM_1.1
    ##  7 chr1  10579544 10579545 rs1411402   BCAC_FM_1.1
    ##  8 chr1  10580890 10580891 rs2483677   BCAC_FM_1.1
    ##  9 chr1  10581050 10581051 rs2506885   BCAC_FM_1.1
    ## 10 chr1  10581657 10581658 rs2056417   BCAC_FM_1.1
    ## # … with 5,365 more rows

##### Trait known genes

The known genes file should be a text file with a single column of known gene symbols. These symbols must be GENCODE-compatible.

``` r
known_genes_file <- system.file("example_data", "known_genes.txt", package = "EG2")
read.delim(known_genes_file, header = F)$V1
```

    ##  [1] "AKT1"     "ARID1A"   "ATM"      "BRCA1"    "BRCA2"    "CBFB"    
    ##  [7] "CDH1"     "CDKN1B"   "CHEK2"    "CTCF"     "ERBB2"    "ESR1"    
    ## [13] "FGFR2"    "FOXA1"    "GATA3"    "GPS2"     "HS6ST1"   "KMT2C"   
    ## [19] "KRAS"     "LRRC37A3" "MAP2K4"   "MAP3K1"   "NCOR1"    "NF1"     
    ## [25] "NUP93"    "PALB2"    "PIK3CA"   "PTEN"     "RB1"      "RUNX1"   
    ## [31] "SF3B1"    "STK11"    "TBX3"     "TP53"     "ZFP36L1"

##### Alternative Weights

The full annotation weights and descriptions are stored as a TSV file in `data/weights.tsv`. This can be accessed as a dataframe when `EG2` is loaded.

``` r
default_weights
```

A file of alternative weights for annotations can be passed to `weights_file`. The file must be tab-delimited and contain `annotation` and `weight` columns. If an annotation is missing from the `weights_file`, it will be weighted 0.

``` r
predict_target_genes(weights_file = "path/to/alternative/weights.tsv", ...)
```

### Running predict\_target\_genes()

To run `predict_target_genes()`...

``` r
annotations <- predict_target_genes(
  trait = "BC_Michailidou2017_FM",
  variants_file = "example_data/variants.bed",
  known_genes_file = "example_data/known_genes.txt",
  reference_panels_dir = "path/to/EG2_data/"
  )
```

Unless an `out_dir` argument is passed, the results will be saved to "out/${trait}/${celltypes}/". If a `sub_dir` is passed, then run results will be saved to a subdirectory below this. If `do_timestamp = T`, then the run results will be saved to a time-stamped subdirectory.

If you are calling `predict_target_genes()` repeatedly in the same session, you can load the large reference objects `H3K27ac` and `HiChIP` into the global environment once, and then pass them to the function pre-loaded. This prevents redundant re-loading with each call to `predict_target_genes()`.

``` r
# 1. load large reference panel objects
HiChIP <- readRDS("path/to/EG2_data/HiChIP.rds")
H3K27ac <- readRDS("path/to/EG2_data/H3K27ac.rds")
# 2. pass objects to predict_target_genes for a quicker runtime
annotations <- predict_target_genes(
  HiChIP = HiChIP,
  H3K27ac = H3K27ac,
  ...
  )
```

This package was written using package development conventions from <https://r-pkgs.org/>.
