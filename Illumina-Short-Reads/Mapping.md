# Mapping

Mapping reads back onto the assembly essentially means finding out which reads and how many reads from each sample make up each contig. This exercise is essential in assessing sequencing depth (coverage) for contigs and later bins (organisms). Mapping produces a .sam file (**S**equence **A**lignment **M**ap, massive), containing all the information about which/how many/quality of reads have been recruited for each contig.  

Profiles of sequencing depths (differential coverage profiles) also inform the binning algorithm and make it much more accurate. For example: all contigs with similar seq depths in each sample likely belong to the same organism. 

Bowtie2 and bbmap are two very commonly used tools for mapping.

### [bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/)

bbmap does this task very effectively/fast. We can use some parameters to speed up the mapping job:

```sh
ambigious=toss		# if uncertainty about which contig a read belongs to, throw it out
minid=0.95			# minimum 95% match between read sequence and contig sequence
```

After bbmap, we use [samtools](http://www.htslib.org/) to convert the .sam file to a .bam file (**B**inary **A**lignment **M**ap) and index it to a .bai file (**B**inary **A**lignment **I**ndex). .bam file is a binary indexed file (machine readable, not human). FYI You cannot open it in a text editor. It's much more efficient to query and search this binary file than the .sam one. 

Once coverted, .sam files **can and should be deleted**, as they are a massive waste of disk space and you can extract all relevant info from a bam file using samtools anyway. 

**At the end of the mapping**, you should have a .bam and a .bam.bai (index) file for each sample. 

```shell
for fn in ../01_qc/qc_files/CAT*_R1_001.qc.fastq;
do
	base="${fn:0:-16}"
	newname=$(basename $fn .qc.fastq)
	sample="${newname:4:-7}"

bbmap.sh ref= in=${base}_R1_001.qc.fastq in2=${base}_R2_001.qc.fastq out=./{sample}_bbmap.sam covstats=./{sample}_bbmap_covstats.txt scafstats=./{sample}_bbmap_scafstats.txt threads=20 minid=0.95 ambiguous=toss

samtools view -b ./{sample}_bbmap.sam | samtools sort -o ./{sample}_bbmap_sorted.bam

samtools index ./{sample}_bbmap_sorted.bam

mv ./{sample}_bbmap_scafstats.txt stats/
mv ./{sample}_bbmap_covstats.txt stats/
rm ./{sample}_bbmap.sam

done
```

### Sequence Depth

Once you have all your .bam and .bam.bai files for each sample, we can use the following perl script to calculate sequencing depth profiles for each contig across all the samples. This produces a depth.txt file, which will be used by binners, namely Metabat2. This script is part of [MetaBat2](https://bitbucket.org/berkeleylab/metabat/src/master/)

Produce depth.txt file for coverage and binning

```shell
conda activate metabat2

jgi_summarize_bam_contig_depths --outputDepth ./depth.txt ./*.bam
```