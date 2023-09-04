# General notes on metatranscriptomics

## What is metatranscriptomics?

16srRNA = "who is there"

Metagenomics  = "what can they do?"

Metatranscriptomics = "who is what doing?"

### General RNA-seq worflow

1. Extract RNA
2. Fragment and sequence
3. Align reads to known predicted genes to obtain relative expression
4. Readout is relative gene expression

### Microbiome challenges

- lack of a poly-A tail - makes it difficult to isolate bacterial mRNA and results in massive rRNA contamination

- Environmental microbiomes samples lack ref genomes making it difficult to map reads back to their source transcripts

- production of mRNA and protein biosynthesis is not linear so not all mRNA is translated into proteins

## A typical metatranscriptomics analytical pipeline

1. Remove low quality reads
2. Remove rRNA reads
3. Remove host reads (if necessary)
4. Assemble? - the longer the fragments - the better the chance of annotation
5. Identify bacterial transcripts
6. Map to pathways
7. Sample comparisons

How many reads are enough?  
*Could generate a rarefaction curve based on the number of enzymes covered in a metagenome vs the number of reads associated with them*
~25% of the reads will be mRNA

Need to build a workflow and then work with established tools and swap them out as they are updated (same as genomics)

#### 1. Quality Control

Trimmomatic? BWA - Host?BLAT - rRNA?  
Reads < 36bp are discarded

#### 2. Assembly

Trinity?

Does length of the transcripts affect the annotation quality?

#### 3. Annotation methods

Identification of functions can be done with sequence similarity matches using BWA, BLAT & Blast  
BWA and BLAT - require near perfect matches

1. Map to a metagenome/MAGs generated from same sample
2. Map to NCBI reference genomes or genes for overall functional analyses

Could use BlastX to work in "peptide" space and search protein databases  
Very time consuming but could use Diamond and Vsearch

*How many reads are mapped versus unmapped?*

What to use as cut-off for the matches?  
What do e-values look like? - apparently will get very high evalue but that doesn't mean the match is not good? - I don't trust this - e-values comment on how random a match is versus not - so it matters particularly for shorter sequence fragments   

85% ID over 65% of length?

#### 4. Converting Mappings to Expression

Normalizing the expression of the genes - genes have variable lengths - more likely to sample from a longer gene so have to account for this

**RPKM** Reads per kilobase of transcript mapped

RPKM <sub>geneA</sub> = (10<sup>9</sup>)(C<sub>geneA</sub>)/NL

C<sub>geneA</sub> = number of reads mapped to geneA  
N = total number of reads in sample  
L = length of gene in kbp

Bowtie and Cufflinks will do this for you?

#### 5. Taxonomic annotation

Can split reads into bins prior to annotation

Compositional methods: nt frequency, codon bias - kmer profiles, Nearest neighbours

GIST - pipeline - which method to use for taxonomic differentiation?

#### 6. Data Analysis

mRNA expression =/= rRNA abundance

Just because you have a highly abundant organism doesn't mean it's highly active? - this is not intuitive??

##### Visualizing results

Look into more systems based visualizations rather than functional groups? 

i.e. instead of all the transporters together, group genes by which pathways they belong to - might be good to use metaERG for this.

Therefore: what pathways are upregulated 

Protein-protein networks as a scaffold to overlay expression data

1. Gene set enrichment analysis (GSEA) = grouping functionally related genes
2. Network approach = identify groups of genes based on functional interaction and then overlay expression

Piecharts of taxa relative abundance on top of network nodes to show which taxa contribute more or less to the upregulation of which genes

MEGAN - can colour which pathway genes are up or down regulated but is based on KEGG  
KEGG pathways are curated for very few organisms so are they actually representative of the microbes in your sample?

Fold change? vs DMC?

Differentially expressed genes can be subsequently analyzed through GSEA based on hypergeometric distributions. 

**Hypergeometric distributions** = like binomial - likelihood of k successes in n draws without replacement  
