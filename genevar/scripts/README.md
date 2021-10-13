## The scripts in this directory generate the whole-genome data as an input for the server.R and ui.R scripts.

### Trace

```
snakemake --cores 13 -p
```

R/Bioconductor packages used include:
- dplyr
- rtracklayer
- BiocManager
- GenomicRanges


