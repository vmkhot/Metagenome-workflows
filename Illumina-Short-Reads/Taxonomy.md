# Taxonomy

## Taxonomy using 16s rRNA genes (reads based)
### [Phyloflash](http://hrgv.github.io/phyloFlash/#:~:text=metagenomic%20dataset.-,phyloFlash%20is%20a%20pipeline%20to%20rapidly%20reconstruct%20the%20SSU%20rRNAs,taxa%20(NTUs)%20are%20handled.)

Phyloflash is a useful tool to get a 'first-look' into your community data using 16s rRNA genes. While it is not strictly necessary, I would recommend it to satiate your excitement about getting metagenomic data back!

Phyloflash is a pipeline. It assembles data, finds 16s genes and classifies taxonomy using the [SILVA](https://www.arb-silva.de/) database for each sample. 

```bash
#!/bin/bash
#SBATCH --job-name=phylo      # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --mem=50G                    # Job memory request
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=phylo_test%j.log     # Standard output and error log
pwd; hostname; date

#phyloFlash.pl -dbhome /home/vmkhot/data/Programs/silva/138.1 -lib testy  -read1 ../01_qc/qc_files/Li32276_S52_L001_R1_001.qc.fastq -read2 ../01_qc/qc_files/Li32276_S52_L001_R2_001.qc.fastq -readlength 150 -CPUs 30

for fn in ../01_qc/qc_files/CAT*_R1_001.qc.fastq
do
        base="${fn:0:-16}"
        newname=$(basename $fn .qc.fastq)
        sample="${newname:4:-7}"
    	phyloFlash.pl -dbhome /home/vmkhot/data/Programs/silva/138.1 -lib ${sample}  -read1 ${base}_R1_001.qc.fastq -read2 ${base}_R2_001.qc.fastq -readlength 150

done
```
## Taxonomy of bins/MAGs
For the taxonomy of bins, we typically use the **G**enome **T**axonomy **D**ata**b**ase **T**ool**k**it ([GTDBTk](https://github.com/Ecogenomics/GTDBTk)). This step requires you to have the taxonomy database from [here](https://gtdb.ecogenomic.org/). GTDBTk works by identifying bacterial and arhcaeal single copy genes in your "genome" (bin) using Hidden Markov Models. It then places your "genome" in a tree with reference genomes and identifies which clade and what taxonomic level your genome clusters with. This is a slow process. 

### GTDBTk

```bash
#!/bin/bash
#SBATCH --job-name=gtdb			# Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20           # Number of CPU cores per task
#SBATCH --mem=300G                   # Job memory request
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --partition=bigmem
#SBATCH --output=gtdb%j.log    # Standard output and error log
pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate gtdbtk

export GTDBTK_DATA_PATH=/work/ebg_lab/referenceDatabases/GTDB/release202/ > /home/vmkhot/miniconda3/envs/gtdbtk/etc/conda/activate.d/gtdbtk.sh

gtdbtk classify_wf -x fa --cpus 20 --genome_dir ../../04_binning/allBins/ --out_dir ./ --pplacer_cpus 20
```

The results of this are fairly self-explanatory