# Mapping mRNAs to MAGs/Genomes

Software to use dependent on data quality and RNA read length - mine is only 50bp - this is a "small RNA sequence" but does not actually come from a "small RNA"

Going to attempt using BBMap's "[seal.sh](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/seal-guide/)" for mapping and generation of RPKM values.

Looks like `seal.sh` can map and produce RPKM values and also split the reads into bins? - How to format reference for these cases?

For mapping to a genome assembly, split into bins with genes?

```bash
seal.sh in=reads.fq *.fa stats=sealstats.txt rpkm=sealrpkm.txt ambig=random
```

For splitting reads into bins?

```bash
seal.sh in=reads.fq *.fa pattern=out_%.fq outu=unmapped.fq ambig=all refnames=t
```


