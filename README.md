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

### Package data

This package uses both reference genomic annotation datasets and user-provided trait-specific datasets. Genomic coordinates use the hg19 build. The package's BED handling follows UCSC's `bedtools` conventions, so it expects 0-based start positions and 1-based end positions.

#### Reference data

##### Internal reference data

Smaller generic reference datasets, including chromosome sizes, GENCODE annotations and coding mutation annotations (`ChrSizes`, `TSSs`, `exons`, `introns`, `promoters`, `missense`, `nonsense`, `splicesite`) are stored internally as parsed objects in `R/sysdata.R`. They are accessible when the package is loaded, but not visible due to lazy loading. The scripts used to generate the input files are published at <https://github.com/alextidd/tgp_paper/tree/main/wrangle_package_data/sysdata/code>.

The weightings of annotations used to generate pair scores are stored as a TSV file in data/weights.tsv.

##### External reference data

Large cell type-specific reference data that accompany the EG2 package have been published to OSF (<https://osf.io/czyvr/>) and must be downloaded separately. In order to use the package, first download and unzip the data into a local directory...

``` bash
mkdir path/to/EG2_data/
wget -O EG2_data.zip https://osf.io/czyvr/download --no-check-certificate
unzip -d path/to/EG2_data/ EG2_data.zip
```

...and then pass the directory to the `reference_panels_dir` argument of the `predict_target_genes()` function.

``` r
predict_target_genes(reference_panels_dir = "path/to/EG2_data/", ...)
```

The scripts used to generate these data are published at <https://github.com/alextidd/tgp_paper/tree/main/wrangle_package_data/reference_panels/code>.

#### User-provided data

There is one required user-provided file for the `predict_target_genes()` function: the trait variants (the `variants_file` argument). Known genes for the trait (the `known_genes_file` argument) are needed only if `do_performance = T`. Example files are published at <https://github.com/alextidd/tgp_paper/tree/main/wrangle_package_data/traits/output> for several tested traits. The scripts used to generate all trait files are published at <https://github.com/alextidd/tgp_paper/tree/main/wrangle_package_data/traits/code>.

##### Trait variants

The variants file should be a BED file with metadata columns for the variant name and the credible set to which it belongs.

``` bash
head path/to/trait/variants.bed
chr1    10551762        10551763        rs657244        BCAC_FM_1.1
chr1    10563363        10563364        rs202087283     BCAC_FM_1.1
chr1    10564674        10564675        rs2847344       BCAC_FM_1.1
chr1    10566521        10566522        rs617728        BCAC_FM_1.1
chr1    10569000        10569000        rs60354536      BCAC_FM_1.1
```

##### Trait known genes

The known genes file should be a text file with a single column of known gene symbols. These symbols must be GENCODE-compatible.

``` bash
head path/to/trait/known_genes.txt
AKT1
ARID1A
ATM
BRCA1
BRCA2
```

### Running predict\_target\_genes()

To run `predict_target_genes()`...

``` r
annotations <- predict_target_genes(
  trait = "BC_Michailidou2017_FM",
  variants_file = "path/to/trait/variants.bed",
  known_genes_file = "path/to/trait/known_genes.txt",
  reference_panels_dir = "path/to/EG2_data/"
  )
```

Unless an `out_dir` argument is passed, the results will be saved to "out/${trait}/${celltypes}/". If `sub_dir` is passed, then run results will be saved to a subdirectory below this. If `do_timestamp = T`, then the run results will be saved to a time-stamped subdirectory.

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
