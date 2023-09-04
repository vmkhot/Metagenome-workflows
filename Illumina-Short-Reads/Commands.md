# Supplementary Methods

Commands, bioinformatic programs and their parameters used for metagenomic aand metatranscriptomic data analysis

## Metagenomics

### Quality Control

**FastQC**

```bash
fastqc *_R1_001.fastq.gz *_R2_001.fastq.gz -o ./
```

**bbmap and samtools**

```bash
#Step 1: Trim off last bases
bbduk.sh in=$sample_R1_001.fastq.gz in2=$sample_R2_001.fastq.gz out=$sample_R1_001.lastbase_rm.fastq out2=$sample_R2_001.lastbase_rm.fastq ftm=5

#Step 2: trim off the partial adapter
bbduk.sh in=$sample_R1_001.lastbase_rm.fastq in2=$sample_R2_001.lastbase_rm.fastq out=$sample_R1_001.adapter_rm.fastq out2=$sample_R2_001.adapter_rm.fastq ref=~/data/Programs/Metagenomics/bbmap/resources/adapters.fa tbo tpe k=23 mink=11 hdist=1 ktrim=r

# Step 3: filter out contaminates
bbduk.sh in=$sample_R1_001.adapter_rm.fastq in2=$sample_R2_001.adapter_rm.fastq out=$sample_R1_001.dec.fastq out2=$sample_R2_001.dec.fastq outm=$sample_contaminates.fq ref=~/data/Programs/Metagenomics/bbmap/resources/phix_adapters.fa.gz k=31 hdist=1 stats=$sample_stats.txt

# Step 4: clipping off the low quality ends and low complexity regions (high entropy; AAAAAAA)
bbduk.sh in=$sample_R1_001.dec.fastq in2=$sample_R2_001.dec.fastq out=$sample_R1_001.qc.fastq out2=$sample_R2_001.qc.fastq qtrim=rl trimq=15 minlength=30 entropy=0.5
```

## Assembly

**MegaHit**

```bash
#coassembly
megahit -1 CAT_ALL_R1_001.qc.fastq -2 CAT_ALL_R2_001.qc.fastq -o output/dir

#individual assemblies
megahit -1 sample_R1_001.qc.fastq -2 sample_R2_001.qc.fastq -o output/dir/S52

#assemblies by day

```

**MetaViralSpades**

```bash
spades.py --metaviral -t 30 --pe1-1 CAT_S64_R1_001.qc.fastq --pe1-2 CAT_S64_R2_001.qc.fastq --pe2-1 CAT_S65_R1_001.qc.fastq --pe2-2 CAT_S65_R2_001.qc.fastq --pe3-1 CAT_S66_R1_001.qc.fastq --pe3-2 CAT_S66_R2_001.qc.fastq --pe4-1 CAT_S70_R1_001.qc.fastq --pe4-2 CAT_S70_R2_001.qc.fastq -o ./S64_65_66_70_viralSpades/
```

**MetaQUAST**

```bash
metaquast.py -m 500 -t 20 --fast -o metaquast -l subsample_viral contigs.fasta
```

## Mapping

**bbmap & samtools & metabat2**

```bash
bbmap.sh ref= in=sample_R1_001.qc.fastq in2=sample_R2_001.qc.fastq out=./sample_bbmap.sam covstats=./sample_bbmap_covstats.txt scafstats=./sample_bbmap_scafstats.txt threads=20 minid=0.95 ambiguous=toss

samtools view -b ./sample_bbmap.sam | samtools sort -o ./sample_bbmap_sorted.bam

samtools index ./sample_bbmap_sorted.bam

jgi_summarize_bam_contig_depths --outputDepth ./depth.txt ./*.bam
```

## Binning

Binning from coassembly into Metagenome-Assembled-Genomes (MAGs)

```bash
metabat2 -i final.contigs.fa.gz -o /Bin -a depth.txt --unbinned -t 30 -v
```

Individual binning of Ca. Phormidium alkaliphilum MAGs

**MetaQUAST**

```bash

```

## MAG Quality

**CheckM2**

```bash

```

## Taxonomy

**Phyloflash**

```bash
phyloFlash.pl -dbhome silva/138.1 -lib samples  -read1 sample_R1_001.qc.fastq -read2 sample_R2_001.qc.fastq -readlength 150
```

**GTDB-tk**

Genome Taxonomy Database release 202 was used to classify MAG taxonomy

```bash
gtdbtk classify_wf -x fa --cpus 20 --genome_dir ./allBins/ --out_dir ./ --pplacer_cpus 20
```

## MAG Annotation

**MetaERG2**

```bash
apptainer exec -B /work/ebg_lab/ ~/metaerg_latest.sif metaerg --database_dir /referenceDatabases/metaerg_database/ --path_to_signalp /referenceDatabases/metaerg_database/ --path_to_tmhmm /referenceDatabases/metaerg_database/ --contig_file ./allBins/ --rename_genomes --rename_contigs --cpus 40 --file_extension .fna
```

## MAG and Contig Coverage and Relative Abundance

**CoverM**

```bash
```

## CRISPR Identification

**minced**

```bash
```

## Viral Sequence Identification

**Virsorter2**

```bash
```

## Viral Host Predictions

**BLASTn**

```bash
```

**CRISPR Spacers**

```bash
```

**tetraucleotide frequency profiles**

```bash
```

**iphop**

```bash
```

## Viral Taxonomy

An in-house script was used to predict viral taxonomy against the latest viral reference sequence database (v?). Script can be found [here]()

```bash
```

The resulting network was visualized in Cytoscape

## Viral Sequence Annotation

**DRAM-V**

```bash
```
