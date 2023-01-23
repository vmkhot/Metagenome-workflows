# Polishing
## Polishing with Long-Reads
Typically, [Racon](https://github.com/isovic/racon) is used for polishing using the same long-reads used to create the assembly. The steps are to map with minimap2 and then to .sam file to polish with Racon. We do this 4x, each time with the newer assembly.

Download [Racon](https://anaconda.org/bioconda/racon)
```bash
# to download, I used conda
conda create -n racon
conda activate racon
conda install -c bioconda racon
```
To run long-read polishing
```bash
# long-read mapping
conda activate minisuite
minimap2 -ax map-ont ../02_Assembly/7_final_consensus.fasta ../all_reads_trimmed_HQ.fastq.gz > longreads_consensus_assembly_minimap2.sam
# polish1
conda activate racon
racon -t 20 ../all_reads_trimmed_HQ.fastq.gz longreads_consensus_assembly_minimap2.sam ../02_Assembly/7_final_consensus.fasta  > racon1.fasta

# Mapping2
conda activate minisuite
minimap2 -ax map-ont ./racon1.fasta ../all_reads_trimmed_HQ.fastq.gz > longreads_racon1_minimap2.sam
# polish2
conda activate racon
racon -t 20 ../all_reads_trimmed_HQ.fastq.gz longreads_racon1_minimap2.sam racon1.fasta  > racon2.fasta

# Mapping3
conda activate minisuite
minimap2 -ax map-ont ./racon2.fasta ../all_reads_trimmed_HQ.fastq.gz > longreads_racon2_minimap2.sam
# polish3
conda activate racon
racon -t 20 ../all_reads_trimmed_HQ.fastq.gz longreads_racon2_minimap2.sam racon2.fasta  > racon3.fasta

# Mapping4
conda activate minisuite
minimap2 -ax map-ont ./racon3.fasta ../all_reads_trimmed_HQ.fastq.gz > longreads_racon3_minimap2.sam
# polish4
conda activate racon
racon -t 20 ../all_reads_trimmed_HQ.fastq.gz longreads_racon3_minimap2.sam racon3.fasta  > racon4.fasta
```
## Polishing with Short Reads (Recommended)
Coming soon...

## Correcting Frameshift Errors

**What is a frameshift error?**

This video is [the best explanation](https://www.youtube.com/watch?v=JaW42lROslE&ab_channel=PatrickHaney) I've found on these errors so far. Nanopore reads are prone to errors in regions of the genome that include the same nucleotide (e.g. AAAAA). The most common error are small "indels" (insertion and/or deletion). If one or more nucleotides are inserted or deleted, the amino-acid coding "frame" (codon) is shifted and a different/incorrect amino acid is predicted.

To correct this, the following program will map proteins from closely related organisms back to the raw sequences and correct these frameshifts.

**METHOD**

As I'm trying to correct an Opitutales genome, for this example, I downloaded all the genomes/MAGs from Opitutae order from NCBI as amino acid fasta files (.faa). 

*Download* the genomes you will use for the correction.


*Download* [Diamond](https://github.com/bbuchfink/diamond) - a program very similar to, but much faster than, Blast. It needs to be in your $PATH
```bash
wget http://github.com/bbuchfink/diamond/releases/download/v2.0.14/diamond-linux64.tar.gz

tar xzf diamond-linux64.tar.gz

# add to your path
export PATH=$PATH:/path/to/directory/containing/diamond
```
*Download* [Proovframe](https://github.com/thackl/proovframe) - program for frameshift error correction. It also needs to be in your $PATH
```bash
git clone https://github.com/thackl/proovframe

# add to your path
export PATH=$PATH:/path/to/directory/containing/proovframe

# to initialize these programs in your shell
Restart your terminal
```
Run Proovframe
```bash
# map proteins to reads
 proovframe map -t 50 -a ncbi_dataset/data/Opi_all_proteins.faa -o proovframe/raw-seqs.tsv ../analysis/Consensus_Genome/racon4.fasta
 
# fix frameshifts in reads
proovframe fix -o proovframe/corrected-seqs.fa ../analysis/Consensus_Genome/racon4.fasta proovframe/raw-seqs.tsv
```
## Checking your polished assemblies
To look for "improvement" your polished assemblies, I would recommend running CheckM or [CheckM2](https://github.com/chklovski/CheckM2) to see if the completeness of your genome has increased or whether the contamination has decreased. 

Below is a example of how my genome has changed based on CheckM2. As you can see, the best genome was achieved after long-read polishing and frameshift correction

|Name|Completeness|Contamination
|---|---|---|
1_canu|85.26|4.15|
2_miniasm|71.86|24.82
3_raven|89.69|2.98|
consensus_genome_trycycler|92.17|3.4|
consensus_racon1|87.59|3.65
consensus_racon2|85.37|4.06|
consensus_racon3|80.91|3.68|
consensus_racon4|88.94|5.22|
consensus_racon4_proovframe_corrected_genome|95.86|3.37|




