# ASSEMBLY
You are probably looking to assemble full genomes if you are working with nanopore data. There are a number of ways/programs to achieve this depending on if you have long + short reads and whether you're assembling a single genome or metagenome. 

## For Genomes
For single bacterial genomes, I highly recommend assembling your reads with multiple programs/methods and using [Trycycler](https://github.com/rrwick/Trycycler/wiki) to generate a consensus assembly from all of them. 

**Why?** To correct for errors such as missing ends of contigs, spurious contigs, missing plasmids, misassembly on repeat regions etc. Read the Trycycler wiki for more information on these.

While there are many assemblers for single genomes that exist, I've only experimented with success with the following:

### [Canu](https://github.com/marbl/canu)
Canu is an older assembly PIPELINE but still works very well. For nanopore reads, all reads are assumed to be raw and untrimmed. 
```bash
#To download, use the binary release (https://github.com/marbl/canu/releases)

#-useGrid=remote : Required if you are submitting jobs to a job scheduler like SLURM, otherwise the pipeline will stop at each step

#genomeSize=    (required to estimate the genome size)
#-maxInputCoverage=100  : limit input coverage to 100x depth

canu -useGrid=remote -p opitutales -d opitutales genomeSize=4m maxInputCoverage=100 -nanopore ../all_passed_reads.fastq.gz
```
### [Raven](https://github.com/lbcb-sci/raven)
Raven is a relatively new and fast (IMO) assembler. It has GPU integration if you are interested in processing heavy datasets and also includes [Racon](https://github.com/isovic/racon) polishing in its pipeline

```bash
#To download, I used conda
conda install -c bioconda raven-assembler

#To run
#-p: iterations of polishing by racon
raven -t 30 -p 4 ../all_reads_trimmed_HQ.fastq.gz
```
### [Minimap2 + Miniasm](https://github.com/lh3/minimap2)
The "mini" suite of tools is very versatile and work together, therefore I recommend installing them altogether in a single conda environment.
```bash
#To download
conda create -n minisuite
conda activate minisuite
conda install -c bioconda minimap2
conda install -c bioconda miniasm
```

```bash
conda activate minisuite

#First the long reads are mapped to themselves to produce a .paf mapping alignment file
#-x: type of reads to input
minimap2 -x ava-ont ../all_reads_trimmed_HQ.fastq.gz ../all_reads_trimmed_HQ.fastq.gz | gzip > ./minimap.paf.gz

#Use the alignment and reads to assemble into a graph (.gfa)
miniasm -f ../all_reads_trimmed_HQ.fastq.gz ./minimap.paf.gz > miniasm.gfa

#Convert graph to fasta file
awk '/^S/{print ">"$2"\n"$3}' miniasm.gfa > miniasm.fasta
```
### [Flye](https://github.com/fenderglass/Flye)
 

```bash
#To download, I used conda
conda create -n flye
conda activate flye
conda install flye

#--genome-size: estimated genome size for organism or for metagenome (mg size)
#-i: polishing iterations
#-t: threads

flye --nano-raw all_reads_trimmed_HQ.fastq.gz --genome-size 4m --out-dir assembly_fly -i 0 -t 30
```
### [Trycycler](https://github.com/rrwick/Trycycler/wiki)

Even though I've added my scripts here, I HIGHLY recommend reading the [Trycycler](https://github.com/rrwick/Trycycler/wiki) wiki to learn what you're actually doing

```bash
#To download
conda create --n trycycler
conda activate trycycler
conda install trycycler

#to cluster the assemblies

trycycler cluster --assemblies ../*.fasta --reads ../../all_reads_trimmed_HQ.fastq.gz --out_dir ./output

#to reconcile the clusters
trycycler reconcile --reads ../../all_reads_trimmed_HQ.fastq.gz --cluster_dir output/cluster_001

#to run multiple sequence alignment

trycycler msa --cluster_dir output/cluster_001

trycycler partition --reads ../../all_reads_trimmed_HQ.fastq.gz --cluster_dirs output/cluster_001

trycycler consensus --cluster_dir output/cluster_001
```


## For Metagenomes
Here is a useful [comparative paper](https://doi.org/10.1038/s41598-020-70491-3) on nanopore metagenome assemblers.

 ### [Flye](https://github.com/fenderglass/Flye)
 

```bash
#--meta: for metagenomes with uneven coverage
#--genome-size: estimated size for metagenome assembly (mg size)

flye --nano-raw 200214_BR_pass_catall.trimmed.fastq --meta --genome-size 450m --out-dir assembly_fly -i 0 -t 30&
```
There are more programs that now deal with nanopore metagenome assemblies but I have not used them so I'll update this when I try them!

## Hybrid Assembly
Coming soon...


## How good is my assembly?

There's no answer to this. There's only some tools to help you determine whether your assembly is good enough. Refer to [this](https://github.com/vmkhot/Metagenome-workflows/blob/main/Illumina-Short-Reads/Assembly.md#how-good-is-my-assembly). You can run some scripts like *assembly-stats* and [MetaQuast](https://quast.sourceforge.net/metaquast) to get information
### Assembly Statistics
This simple little program will give you essential assembly statistics such as # of contigs, N50, N90. (see below)

```bash
#To download
conda install -c bioconda assembly-stats

#To run, will be written to stdout
assembly-stats ./miniasm.fasta
```

Example stats for ./miniasm.fasta
|||||
|---|---|---|---|
|sum = 5945987, n = 2| ave = 2972993.50| largest = 5938639
|N50 = 5938639 n = 1
N60 = 5938639, n = 1
N70 = 5938639, n = 1
N80 = 5938639, n = 1
N90 = 5938639, n = 1
N100 = 7348, n = 2
N_count = 0
Gaps = 0

### Assembly Graphs
All the nanopore assembly programs produce a .gfa (Assembly graph) file. View your assembly graph (.GFA file) in [Bandage](https://rrwick.github.io/Bandage/) to see where contigs are broken and which contigs are circular.

### Circularity
How to check for circularity

### CheckM2
