# Assembly

Assembly is an essential step. This is where DNA sequences are stitched together into longer DNA sequences. 

**How do assemblers work?** There are multiple methods for assembly but Graph-based assembly is the most common method. Rob Edwards has short easy-to-understand YT videos on assembly, so I won't explain each step here. 
Intro: https://www.youtube.com/watch?v=KASvlXYPCBI
Graph-Based: https://www.youtube.com/watch?v=OY9Q_rUCGDw

The assemblers shown here (Megahit and MetaSpades) both are De Bruijn graph assemblers.

#### Co-assemble or not?

In a coassembly, multiple samples are input into a single assembly, compared to individual assemblies, which only has 1 sample per assembly. The decision to make a coassembly or not depends on your research question. Coassemblies do have some advantages: 1) higher read depth, i.e. useful for assembling genomes of lower abundance, 2) single reference assembly to compare all your samples to, and 3) generation of differential coverage profiles for the binner. Good candidates for coassemblies include community time-series and enrichment cultures. 

Conversely, higher complexity communities, e.g. soil samples, where you are interested in the community differences between samples may benefit from individual assemblies instead.

### Megahit Assembly

[Megahit](https://github.com/voutcn/megahit) is a conservative assembler: this means it prefers to keep contigs shorter (fragmented) with less errors (misassembly) rather than assembling longer fragments with more errors. 

```shell
#For co-assembly
#Have to concatenate all QC'ed forward reads from 1 sample into 1 file and reverse into another. 
for j in {52..70}; do cat *S"$j"*_R1_001.qc.fastq > CAT_S"$j"_R1_001.qc.fastq; done;
```

Assembly script

```sh
#!/bin/bash
#SBATCH --job-name=megahit           # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20           # Number of CPU cores per task
#SBATCH --mem=150G                   # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --partition=bigmem			# Partition to bigmem because lots of memory required
#SBATCH --output=megahit%j.log       # Standard output and error log

pwd; hostname; date

#coassembly
megahit -1 CAT_ALL_R1_001.qc.fastq -2 CAT_ALL_R2_001.qc.fastq -o output/dir

#individual assemblies
megahit -1 S52_R1_001.qc.fastq -2 S52_R2_001.qc.fastq -o output/dir/S52
megahit -1 S54_R1_001.qc.fastq -2 S54_R2_001.qc.fastq -o output/dir/S54
megahit -1 S56_R1_001.qc.fastq -2 S56_R2_001.qc.fastq -o output/dir/S56
```

### [MetaSPAdes](https://github.com/ablab/spades#installation)

Metaspades will attempt assembling longer contigs than Megahit

MetaSPAdes (Viral)

This is MetaSPAdes for viral and plasmid reads specifically but the commands are similar for other spades

```shell
#!/bin/bash
#SBATCH --job-name=spades      # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --mem=250G                    # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --partition=bigmem
#SBATCH --output=spades%j.log     # Standard output and error log
pwd; hostname; date

#metaviralspades
#with subsampled at 5%
spades.py --metaviral -t 30 --pe1-1 ../01_qc/subsample_viral_fraction/CAT_S64_R1_001.0.05subsample.qc.fastq --pe1-2 ../01_qc/subsample_viral_fraction/CAT_S64_R2_001.0.05subsample.qc.fastq --pe2-1 ../01_qc/subsample_viral_fraction/CAT_S65_R1_001.0.05subsample.qc.fastq --pe2-2 ../01_qc/subsample_viral_fraction/CAT_S65_R2_001.0.05subsample.qc.fastq --pe3-1 ../01_qc/subsample_viral_fraction/CAT_S66_R1_001.0.05subsample.qc.fastq --pe3-2 ../01_qc/subsample_viral_fraction/CAT_S66_R2_001.0.05subsample.qc.fastq --pe4-1 ../01_qc/subsample_viral_fraction/CAT_S70_R1_001.0.05subsample.qc.fastq --pe4-2 ../01_qc/subsample_viral_fraction/CAT_S70_R2_001.0.05subsample.qc.fastq -o ./S64_65_66_70_viralSpades_subsample/

#all reads
spades.py --metaviral -t 30 --pe1-1 ../01_qc/subsample_viral_fraction/CAT_S64_R1_001.qc.fastq --pe1-2 ../01_qc/subsample_viral_fraction/CAT_S64_R2_001.qc.fastq --pe2-1 ../01_qc/subsample_viral_fraction/CAT_S65_R1_001.qc.fastq --pe2-2 ../01_qc/subsample_viral_fraction/CAT_S65_R2_001.qc.fastq --pe3-1 ../01_qc/subsample_viral_fraction/CAT_S66_R1_001.qc.fastq --pe3-2 ../01_qc/subsample_viral_fraction/CAT_S66_R2_001.qc.fastq --pe4-1 ../01_qc/subsample_viral_fraction/CAT_S70_R1_001.qc.fastq --pe4-2 ../01_qc/subsample_viral_fraction/CAT_S70_R2_001.qc.fastq -o ./S64_65_66_70_viralSpades/
```

### How good is my assembly?

Checking assembly quality with [**MetaQuast**](https://quast.sourceforge.net/metaquast)

--fast option doesn't output graphs and plots and all the junk

```shell
#!/bin/bash
#SBATCH --job-name=metaquast           # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20           # Number of CPU cores per task
#SBATCH --mem=50G                    # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --output=metaquast%j.log       # Standard output and error log

pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate quast

metaquast.py -m 500 -t 20 --fast -o metaquast -l subsample_viral contigs.fasta
```

Parameters to assess a good assembly: 

* high N50 (>2500): N50 is **not** the median length. It's the halfway of the cumulative total of all your contig lengths sorted by length. i.e. Above and below this length are are half your total data in bp.
  * High N50 will ensure that contigs are binned into MAGs well
  * When it comes to interpreting your data/annotations, there is also higher chances of finding full operons together on a single contig if they are long. e.g. contig of length 1000bp is <= 1 gene. length 100,000bp <= 100 genes
* how many total bases are contained in contigs > 10kpb vs. shorter (<1000bp) = how much data is stored in longer contigs?

Lots of short contigs signal assembly issues and many of these may be thrown out during binning and annotation (contigs <= 500bp) as getting meaningful data out of them is difficult. 

**Your assembly may be particularly fragmented if you:** 

* have complex environmental samples in a coassembly
  * May be solved using individual assemblies
* Have high strain diversity
  * May require manually solving assembly using a program like Bandage
* have genomes with many repeat regions like CRISPR spacers. Repeat regions confuse the assembler. Why? https://www.youtube.com/watch?v=SiBgAGHJcWY&t=2s
* if the sequencing depth is too high (< 1500)
  * May be solved by subsampling your reads at 5 or 10%
* Have a viral sequences or plasmids which have irregular k-mers and low-complexity and repeat regions

You can also try MetaSpades to get longer contigs. Spades also has a version for viral sequences and plasmids. MetaSpades is SLOOWWW and memory-heavy for large datasets and may result in some errors as it is not conservative.