# binmapr: A R package from binmap analysis

## Installation

```
devtools::install("xuzhougeng/binmapr")
```

## How to use it?

```
library(binmapr)
```

load the data

```
marker <- readLines("marker.txt")

info <- getAdFromVcf("snp_flt.recode.vcf.gz", keep = marker)

AD <- info$AD
CHROM <- unique(info$CHROM)

GT <- callGtFromAd(AD$AD)
```

filter marker with high missing ratio

```
miss_ratio <- rowSums(is.na(GT)) / ncol(GT)
GT_flt <- GT[miss_ratio < 0.05, ]
```

window genotype 

```
CHROM <- CHROM[1:8]

geno <- batchCallGeno(GT_flt, CHROM = CHROM, outdir = ".")
```

## Reference

[binmap](https://mp.weixin.qq.com/s/x6zRylSiPn0LmNtmEGbdAA)

