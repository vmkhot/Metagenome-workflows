# Supplementary Methods

Commands, bioinformatic programs and their parameters used for metagenomic aand metatranscriptomic data analysis

## Metagenomics Analysis

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

### Assembly

**MegaHit**

```bash
#coassembly
megahit -1 CAT_ALL_R1_001.qc.fastq -2 CAT_ALL_R2_001.qc.fastq -o output/dir

#individual assemblies
megahit -1 sample_R1_001.qc.fastq -2 sample_R2_001.qc.fastq -o output/dir/S52

```

**MetaViralSpades**

```bash
spades.py --metaviral -t 30 --pe1-1 CAT_S64_R1_001.qc.fastq --pe1-2 CAT_S64_R2_001.qc.fastq --pe2-1 CAT_S65_R1_001.qc.fastq --pe2-2 CAT_S65_R2_001.qc.fastq --pe3-1 CAT_S66_R1_001.qc.fastq --pe3-2 CAT_S66_R2_001.qc.fastq --pe4-1 CAT_S70_R1_001.qc.fastq --pe4-2 CAT_S70_R2_001.qc.fastq -o ./S64_65_66_70_viralSpades/
```

**MetaQUAST**

```bash
metaquast.py -m 500 -t 20 --fast -o metaquast -l viral contigs.fasta
```

### Mapping

**bbmap & samtools & metabat2**

```bash
bbmap.sh ref= in=sample_R1_001.qc.fastq in2=sample_R2_001.qc.fastq out=./sample_bbmap.sam covstats=./sample_bbmap_covstats.txt scafstats=./sample_bbmap_scafstats.txt threads=20 minid=0.95 ambiguous=toss

samtools view -b ./sample_bbmap.sam | samtools sort -o ./sample_bbmap_sorted.bam

samtools index ./sample_bbmap_sorted.bam

jgi_summarize_bam_contig_depths --outputDepth ./depth.txt ./*.bam
```

### Binning

Binning from coassembly into Metagenome-Assembled-Genomes (MAGs)

```bash
metabat2 -i final.contigs.fa.gz -o /Bin -a depth.txt --unbinned -t 30 -v
```

Individual binning of Ca. Sodalinema alkaliphilum MAGs

**MetaQUAST**

```bash
metaquast.py -m 500 -t 20 -o metaquast_full -l day0,day2,day4,day6,day8,day9 day0/final.contigs.fa day2/final.contigs.fa day4/final.contigs.fa day6/final.contigs.fa day8/final.contigs.fa day9/final.contigs.fa -r phormidium_ref.fasta
```

### MAG Quality

**CheckM2**

```bash
checkm2 predict -t 30 -x fa --input ./allBins/bins/ --output-directory ./Checkm2
```


### Taxonomy

**Phyloflash**

```bash
phyloFlash.pl -dbhome silva/138.1 -lib samples  -read1 sample_R1_001.qc.fastq -read2 sample_R2_001.qc.fastq -readlength 150
```

**GTDB-tk**

Genome Taxonomy Database release 202 was used to classify MAG taxonomy

```bash
gtdbtk classify_wf -x fa --cpus 20 --genome_dir ./allBins/ --out_dir ./ --pplacer_cpus 20
```

### MAG Annotation

**MetaERG2**

```bash
apptainer exec -B /work/ebg_lab/ ~/metaerg_latest.sif metaerg --database_dir /referenceDatabases/metaerg_database/ --path_to_signalp /referenceDatabases/metaerg_database/ --path_to_tmhmm /referenceDatabases/metaerg_database/ --contig_file ./allBins/ --rename_genomes --rename_contigs --cpus 40 --file_extension .fna
```

### MAG and Contig Coverage and Relative Abundance

**CoverM**

```bash
coverm genome -v -x fa -t 20 --methods trimmed_mean -1 *_R1_001.qc.fastq -2 *_R2_001.qc.fastq --genome-fasta-directory bins --bam-file-cache-directory ./coverm_bam -o out_genome_coverage.tsv

#for relative abundance
coverm genome -v -x fa -t 20 --methods relative_abundance --bam-files *.bam -s '~' -o output_relative_abundances.tsv
```

### CRISPR Identification

**minced**

```bash
minced -spacers final.contigs.fa assembly.crisprs assembly.gff
```

### Viral Sequence Identification

**Virsorter2**

```bash
virsorter config --set HMMSEARCH_THREADS=30

virsorter run --keep-original-seq --min-length 1000 --min-score 0.5 -w ./pass1 -i ../final.contigs.fa.gz --include-groups dsDNAphage,ssDNA -j 30 all
```

### Viral Host Predictions

**BLASTn**

```bash
blastn -db final.contigs.fa  -query viral_genomes.fna -out blastn.out
```

**CRISPR Spacers**

```bash
blastn -db viral_genomes.fna  -query ../assembly_spacers.fa -out blastn.out -task "blastn-short" -outfmt 6 -dust no
```
CRISPR-Blastn results were filtered by length >= 25 bp, percent identity >= 95% and e-value <= 1E-5 

**tetraucleotide frequency profiles**

Using an in-house script

**iphop**

```bash
iphop predict --fa_file ../viral_contigs_renamed.singleline.fna --db_dir /work/ebg_lab/referenceDatabases/iphop/db/ --out_dir iphop_output/
```

### Viral Taxonomy

An in-house script was used to predict viral taxonomy against the latest viral reference sequence database (v?). Script can be found [here]()

```bash
blastp -db cat_all_viruses2blast.faa -query cat_all_viruses2blast.faa -out VCproteins_refseq_selfblastp_hsp1.out -outfmt 6 -num_threads 20 -max_hsps 1
```
Blastp results were filtered by e-value <= 1E-5 and percent identity >= 40%

The resulting network was visualized in Cytoscape

### Viral Sequence Annotation

**VirSorter2 + DRAM-V**

```bash
virsorter run --prep-for-dramv --keep-original-seq --min-length 1000 --min-score 0.5 -w ./pass1 -i viral_genomes.fasta --include-groups dsDNAphage,ssDNA -j 50 all

#Annotate
DRAM-v.py annotate -i final-viral-combined-for-dramv.fa -v viral-affi-contigs-for-dramv.tab -o ./dramv/ --threads 30
```

## Metatranscriptome Analysis

### Quality Control

```bash
fastqc -o ./ -t 30 ../raw_fastq/*.fastq.gz

multiqc -o ./ -n rawReads ./
```

Last base and quality trimming using bbduk

```bash
#trim off last base (round to 50bp)
bbduk.sh in=R1_001.fastq.gz in2=R2_001.fastq.gz out=R1_001.lastbase_rm.fastq.gz out2=R2_001.lastbase_rm.fastq.gz ftm=5

#quality filtering (for small RNA, minimum read length of 10 and min quality of 15)
bbduk.sh in=R1_001.lastbase_rm.fastq.gz in2=R2_001.lastbase_rm.fastq.gz  out=R1_001.qc.fastq.gz out2=R2_001.qc.fastq.gz qtrim=rl trimq=15 minlength=10
```

**Sorting rRNA using SortMeRNA**

```bash
sortmerna --ref referenceDatabases/sortmerna_db/smr_v4.3_default_db.fasta \
       --workdir ./sortmerna/ \
       --reads R1_001.qc.fastq.gz --reads R2_001.qc.fastq.gz \
       --aligned rRNA_reads/rRNA_reads.qc \
       --other non_rRNA_reads/non_rRNA_reads.qc \
       --sam --SQ --log --fastx --threads 40 --paired_in
```

### Mapping

Mapping transcriptome reads to nucleotide gene sequences from metagenome

```bash
seal.sh in=reads.fq nucl_seq.fa stats=sealstats.txt rpkm=sealrpkm.txt ambig=all
```

### Differential Expression Analysis

DE analysis was performed using the *DESeq2* pipeline in R. The generic pipeline used to assess quality of raw counts from samples is [here](vmkhot/Metagenome-workflows/Metatranscriptomics/R-scripts/deseq2_script_sample_QC.R) and the time-course script [here](vmkhot/Metagenome-workflows/Metatranscriptomics/R-scripts/deseq2_time_course_script_cyano.R)

