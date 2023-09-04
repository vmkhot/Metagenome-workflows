# QC

First, you need to check the quality of our sequencing data and do some quality control. This means trimming adapters, low complexity regions and low quality ends. 

We check our raw data using FastQC, do QC with bbmap and check again. 

**Note**: We also evaluate the quality of our results after the major steps in this workflow, such as assembly and binning. This is so that your data remains credible throughout and that your end results are actually meaningful to your research questions. 

### FastQC

Used a for-loop to interate through all the raw sample .fastq files starting with "VK" to not duplicate it.

#### Run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

This will provide an html report of your raw reads. Yes, you must open and check each one. 

```bash
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

### QC with bbmap

Used for-loop to qc all the raw files. Removed all the fastq intermediates afterwards

Run bbduk.sh to trim last bases, adapters, and remove contamination. Option to run fastqc again after to check. 

```bash
#!/bin/bash
#SBATCH --job-name=bbduk      # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --mem=50G                    # Job memory request
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=bbduk%j.log     # Standard output and error log
pwd; hostname; date


for fn in /work/ebg_lab/gm/Metagenomes_2022/GAPP2/VK*/*R1_001.fastq.gz;
do
        base="${fn:0:-16}"
        newname=$(basename $fn .fastq.gz)
        shortname="${newname:0:-7}"
        echo $base
# Step 1 (optional), when your reads are 151 bp instead of 150bp, this command trims off the last base of 151bp reads because that base is usually very low quality.
        bbduk.sh in=${base}_R1_001.fastq.gz in2=${base}_R2_001.fastq.gz out=${shortname}_R1_001.lastbase_rm.fastq out2=${shortname}_R2_001.lastbase_rm.fastq ftm=5

# Step 2: trim off the partial adapter
# ktrim=r trim (right) 3' end
# k=23: kmer to use for trimming
# mink=11: minimum kmer to use at the end of the read
# hdist=1: hamming distance (1 mismatch)
# tbo: trim based on pair overlap
# tpe: trim both reads to same length
        bbduk.sh in=${shortname}_R1_001.lastbase_rm.fastq in2=${shortname}_R2_001.lastbase_rm.fastq out=${shortname}_R1_001.adapter_rm.fastq out2=${shortname}_R2_001.adapter_rm.fastq ref=~/data/Programs/Metagenomics/bbmap/resources/adapters.fa tbo tpe k=23 mink=11 hdist=1 ktrim=r

# Step 3: filter out contaminates
        bbduk.sh in=${shortname}_R1_001.adapter_rm.fastq in2=${shortname}_R2_001.adapter_rm.fastq out=${shortname}_R1_001.dec.fastq out2=${shortname}_R2_001.dec.fastq outm=${shortname}_contaminates.fq ref=~/data/Programs/Metagenomics/bbmap/resources/phix_adapters.fa.gz k=31 hdist=1 stats=${shortname}_stats.txt

# Step 4: clipping off the low quality ends and low complexity regions (high entropy; AAAAAAA)
        bbduk.sh in=${shortname}_R1_001.dec.fastq in2=${shortname}_R2_001.dec.fastq out=${shortname}_R1_001.qc.fastq out2=${shortname}_R2_001.qc.fastq qtrim=rl trimq=15 minlength=30 entropy=0.5
done

echo finished qc

#~~~~~~~~~~~~~~#
#MOVING INTERMEDIATE FILES
#mkdir raw_files
mkdir qc_files
mkdir intermediates #remember to delete this folder once you've verified that your qc'ed reads are all G

#mv *.fastq.gz raw_files
mv *.qc.fastq qc_files
mv *.fastq intermediates
mv *.fq intermediates
```

Sometimes bbmap stuffs up and the reads aren't paired properly. To check:

```bash
#Verify pair
reformat.sh in=S108_L001_R1_001.lastbase_rm.fastq in2=S108_L001_R2_001.lastbase_rm.fastq vpair

#fix paired reads
repair.sh -h
```

After running bbmap, I would recommend running FastQC again on some samples which were not ideal earlier to see if the warnings have been resolved. 