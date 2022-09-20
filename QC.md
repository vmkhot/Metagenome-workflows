## QC

First, you need to check the quality of our sequencing data and do some quality control. This means trimming adapters, low complexity regions and low quality ends. 

We check our raw data using FastQC, do QC with bbmap and check again. 

**Note**: We also evaluate the quality of our results after the major steps in this workflow, such as assembly and binning. This is so that your data remains credible throughout and that your end results are actually meaningful to your research questions. 

### FastQC

Used a for-loop to interate through all the raw sample .fastq files starting with "VK" to not duplicate it.

#### Run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

This will provide an html report of your raw reads. Yes, you must open and check each one. 

```shell
#!/bin/bash
#SBATCH --job-name=fastqc      # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=10            # Number of CPU cores per task
#SBATCH --mem=50G                    # Job memory request
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=fastqc%j.log     # Standard output and error log
pwd; hostname; date

for fn in /work/ebg_lab/gm/Metagenomes_2022/GAPP2/VK*/*R1_001.fastq.gz;
do
    base="${fn:0:-16}"
    echo $base
    ~/data/Programs/Metagenomics/FastQC/fastqc ${base}_R1_001.fastq.gz ${base}_R2_001.fastq.gz -o ./
done
```

What do FastQC results mean?

Check that the reads are paired (R1 and R2 have same # of reads) and other graphs to be "normal"