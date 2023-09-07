# Binning

In this step, we produce Metagenome-Assembled-Genomes (MAGs) or bins. These are approximate genomes of organisms or populations in your community. The binning algorithm for MetaBat2 uses 2 main methods to group contigs together: 1) k-mer frequency profiles of the contigs and 2) differential coverage. 

```bash
#!/bin/bash
#SBATCH --job-name=metabat2      # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --mem=100G                    # Job memory request
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=metabat2%j.log     # Standard output and error log
pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate metabat2

metabat2 -i ../02_assembly/megahit/final.contigs.fa.gz -o allBins/Bin -a ../03_mapping/depth.txt --unbinned -t 30 -v

metabat2 -o allBins_2/Bin --unbinned -t 30 -v -i ../02_assembly/megahit/final.contigs.fa.gz ../03_mapping/*.bam

```

There are other binners out there, which I have not experimented with like [MaxBin](https://microbiomejournal.biomedcentral.com/articles/10.1186/2049-2618-2-26), [CONCOCT](https://concoct.readthedocs.io/en/latest/). Results from different binners can also be combined and filtered for the best bins using [DASTool.](https://github.com/cmks/DAS_Tool)

### **How good are my MAGs?**

Once you have MAGs, it's once again important to evaluate the quality of your results. For quality control of MAGs, we typically run [CheckM](https://github.com/Ecogenomics/CheckM/wiki). A new version, [CheckM2](https://github.com/chklovski/CheckM2), was recently released but I have not tried it yet. 

CheckM computes completeness, lineage, contamination, genome size, N50 and much more for MAGs. See below for example.

#### CheckM

```bash
#!/bin/bash
#SBATCH --job-name=metabat2      # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --mem=100G                    # Job memory request
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=metabat2%j.log     # Standard output and error log
pwd; hostname; date
source ~/miniconda3/etc/profile.d/conda.sh

conda activate checkm

#compute quality checks on all the bins
checkm lineage_wf -f allBins/CheckM.txt -t 30 -x fa allBins/ allBins/Checkm

echo part1...finished

#to extract the unbinned fraction of the assembly
checkm unbinned -x fa ./allBins/bins ../02_assembly/megahit/final.contigs.fa ./allBins/bins/unbinned.fa unbinned_stats.tsv
```

#### [CheckM2](https://github.com/chklovski/CheckM2)

```bash
conda activate checkm2

checkm2 predict -t 30 -x fa --input ./allBins/bins/ --output-directory ./Checkm2
```

##### Example Output from CheckM

You'll get a simplified but more accurate version of this from CheckM2.

| Bin Id  | Marker lineage                   | # genomes | # markers | # marker sets | 0    | 1    | 2    | 3    | 4    | 5+   | Completeness | Contamination | Strain heterogeneity |
| ------- | -------------------------------- | --------- | --------- | ------------- | ---- | ---- | ---- | ---- | ---- | ---- | ------------ | ------------- | -------------------- |
| Bin.100 | k__Bacteria (UID3187)            | 2258      | 190       | 119           | 19   | 168  | 3    | 0    | 0    | 0    | 92.71        | 1.85          | 0                    |
| Bin.12  | k__Bacteria (UID203)             | 5449      | 104       | 58            | 22   | 82   | 0    | 0    | 0    | 0    | 71.08        | 0             | 0                    |
| Bin.13  | k__Bacteria (UID2569)            | 434       | 278       | 186           | 8    | 270  | 0    | 0    | 0    | 0    | 96.24        | 0             | 0                    |
| Bin.16  | root (UID1)                      | 5656      | 56        | 24            | 55   | 1    | 0    | 0    | 0    | 0    | 4.17         | 0             | 0                    |
| Bin.17  | k__Bacteria (UID203)             | 5449      | 104       | 58            | 8    | 62   | 26   | 8    | 0    | 0    | 86.21        | 18.34         | 52                   |
| Bin.20  | c__Gammaproteobacteria (UID4202) | 67        | 481       | 276           | 25   | 454  | 2    | 0    | 0    | 0    | 97.48        | 0.48          | 0                    |

In the example above, we have 3 'good' bins for sure, Bin.100, Bin.13 and Bin.20. This is based on the high completeness, low contamination and low strain heterogeneity. The threshold values of these are not set in concrete, so you have to assess for yourself whether you think a **bin is worth keeping or not**. Alternatively, you can copy them from literature or other studies. E.g. Bin.12 above is only expected to be 71% complete, however, we might find that it has some exciting taxonomy or belongs to an organism with a reduced genome, in which case, we would want to keep it for downstream analyses. Similarly Bin.16 appears to be very low completeness and hits the "root" markers, rather than kingdom Bacteria - which may indicate that it is a viral bin or something else. Last thing to note is that the values in the Marker Lineage column is not indicative of the taxonomy - simply the taxonomy of marker sets which were hit the most. 

As you can see, a case can be made for keeping very few bins or keeping all the them and this will depend on your research questions and which organisms you're studying. IMHO, it's best to also do figure out the taxonomy for all the bins before deciding to keep a bin or throw it out, since it might be an organism of interest.



### Sequencing Depth and Relative Abundance by MAG/Genome

For this, I have previously used [CoverM](https://github.com/wwood/CoverM/)

```bash
#!/bin/bash
#SBATCH --job-name=coverm      # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --mem=50G                    # Job memory request
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=coverm%j.log     # Standard output and error log
pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate coverm

# COVERAGE OF MAGS (keeps bams)
# list all your forward and reverse reads
coverm genome -v -x fa -t 20 --methods trimmed_mean \
-1 sample1.r1.fq,sample2.r1.fq,sample3.r1.fq -2 sample1.r2.fq,sample2.r2.fq,sample3.r2.fq \
--genome-fasta-directory ./bins --bam-file-cache-directory ./03_mapping/coverm_bam \
-o out_genome_coverage.tsv

# RELATIVE ABUNDANCE OF MAGS (using bam files generated above)
coverm genome -v -x fa -t 20 --methods relative_abundance \
--bam-files sample1.bam,sample2.bam,sample3.bam \
--genome-fasta-directory ./bins/to_keep -o output_relative_abundances.tsv
```