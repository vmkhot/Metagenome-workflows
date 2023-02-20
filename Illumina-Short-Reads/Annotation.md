# Annotation Genomes/MAGs

Functional annotation of genes will give you an idea of what your organism is capable of doing. Annotation typically starts with gene prediction. The most common program being [Prodigal](https://github.com/hyattpd/Prodigal).

Where does annotation go wrong?

Prodigal does not consider about other genomic features such as CRISPRs, tRNAs, rRNA genes, repeat sequences.

Therefore, we have to predict all the other features first and **MASK** the features for prodigal by replacing the region with Ns (NNNN). Metaerg2 takes this approach to annotation. It will also give you GTDB taxonomy per gene annotation so that might give you some idea of what organism the closest relative of the gene lives in.

## [Metaerg 2.0](https://github.com/kinestetika/MetaErg)

Highly recommend reading the metaerg2 readme by [kinestetika](https://github.com/kinestetika) as it has very good explanations of how the annotation in metaerg2 works and what annotation pitfalls it is trying to solve.

### Set up Metaerg2

To download on a HPC cluster - use [apptainer](https://apptainer.org/docs/admin/main/installation.html) (formlerly singularity) to pull a Docker image and convert to a singularity image (.sif file)

```bash
# pull docker image and convert to singularity
apptainer pull docker://kinestetika/metaerg
```

To run metaerg, you have to enter into the apptainer container (like a python environment?). shell command will open a new shell into the container, where metaerg lives and you can interact like a virtual environment.

```bash
apptainer shell ~/metaerg_latest.sif

# to check the download
metaerg -h

# to exit apptainer container 
cntrl + d
```

By default, apptainer only "binds" your $HOME directory. So to access other directories in the filesystem above that you can run the following command.

```bash
apptainer shell -B /work/ebg_lab/ metaerg_latest.sif
```

However, it's **not permanent** so you have to "-B" everytime you run the apptainer container. Haven't figured out how to make this permanent this yet. Something to do with the "export APPTAINER_BIND" environment variable.

### Download the databases

On the HPC cluster, I had to use a sbatch script to do this. Takes 60min+

```bash
 # to download database
 apptainer exec -B /work/ebg_lab/ ~/data/Programs/metaerg_latest.sif metaerg --download_database --database_dir /work/ebg_lab/referenceDatabases metaerg_database/

# to update the database
 apptainer exec -B /work/ebg_lab/ ~/data/Programs/metaerg_latest.sif metaerg --create_database S --database_dir /work/ebg_lab/referenceDatabases metaerg_database/
```

You will also need download the signalP and tmhmm databases from [here](https://services.healthtech.dtu.dk/software.php). Once you have the link in your email, you can download directly to the server using wget. Download to same the folder as above (metaerg_database)

```bash
# download signalP and tmhmm
cd metaerg_database/

wget url/to/signalP.tar.gz
wget url/to/tmhmm.tar.gz
```

### Run metaerg2 from apptainer

To run metaerg/apptainer from slurm batch script on multiple genomes. It takes about 2h per genome with the below options so you might want to split up your genomes/MAGs into multiple directories and run parallely. *You have to use the "apptainer -B" to bind whatever directory you are working from.*

```bash
#!/bin/bash
#SBATCH --job-name=metaerg2_1      
#SBATCH --nodes=1                    
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=40           
#SBATCH --mem=100G                    
#SBATCH --time=24:00:00              
#SBATCH --output=metaerg2_1run%j.log    

pwd; hostname; date

# metaerg options
# --database_dir
# --path_to_signalp
# --path_to_tmhmm
# --contig_file
# --rename_contigs
# --cpus
# --file_extension

apptainer exec -B /work/ebg_lab/ ~/data/Programs/metaerg_latest.sif metaerg --database_dir /work/ebg_lab/referenceDatabases/metaerg_database/ --path_to_signalp /work/ebg_lab/referenceDatabases/metaerg_database/ --path_to_tmhmm /work/ebg_lab/referenceDatabases/metaerg_database/ --contig_file ../genomes2run/genomes_1 --rename_contigs --cpus 40 --file_extension .fna
```

## [DRAM and DRAM-v](https://github.com/WrightonLabCSU/DRAM/wiki)
This is another annotation pipeline. I've only used this for annotating **viral** genomes

### Download the program
Super easy conda installation.
[Read this](https://github.com/WrightonLabCSU/DRAM/wiki#dram-installation)


### Download databases
Downloading and "preparing" databases requires a lot of memory

```bash
#!/bin/bash
#SBATCH --job-name=dramsetup      # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=10            # Number of CPU cores per task
#SBATCH --mem=300G                    # Job memory request
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=dramsetup%j.log     # Standard output and error log
#SBATCH --partition=bigmem
pwd; hostname; date

# to DOWNLOAD the databases for DRAM
DRAM-setup.py prepare_databases --output_dir /work/ebg_lab/referenceDatabases/DRAM  --threads 50

# **if you're using previosuly downloaded databases**

# this is to initialize your DRAM installation with already configured databases
DRAM-setup.py set_database_locations 
        --kofam_hmm_loc /work/ebg_lab/referenceDatabases/DRAM2/kofam_profiles.hmm 
        --kofam_ko_list_loc /work/ebg_lab/referenceDatabases/DRAM2/kofam_ko_list.tsv 
        --uniref_db_loc /work/ebg_lab/referenceDatabases/DRAM2/uniref90.20220122.mmsdb 
        --pfam_db_loc /work/ebg_lab/referenceDatabases/DRAM2/pfam.mmspro 
        --pfam_hmm_dat /work/ebg_lab/referenceDatabases/DRAM2/Pfam-A.hmm.dat.gz 
        --dbcan_db_loc /work/ebg_lab/referenceDatabases/DRAM2/dbCAN-HMMdb-V9.txt 
        --dbcan_fam_activities /work/ebg_lab/referenceDatabases/DRAM2/CAZyDB.07302020.fam-activities.txt
        --vogdb_db_loc /work/ebg_lab/referenceDatabases/DRAM2/vog_latest_hmms.txt 
        --vog_annotations /work/ebg_lab/referenceDatabases/DRAM2/vog_annotations_latest.tsv.gz
        --viral_db_loc /work/ebg_lab/referenceDatabases/DRAM2/refseq_viral.20220122.mmsdb 
        --peptidase_db_loc /work/ebg_lab/referenceDatabases/DRAM2/peptidases.20220122.mmsdb
        --description_db_loc /work/ebg_lab/referenceDatabases/DRAM2/description_db.sqlite 
        --genome_summary_form_loc /work/ebg_lab/referenceDatabases/DRAM2/genome_summary_form.20220122.tsv 
        --module_step_form_loc /work/ebg_lab/referenceDatabases/DRAM2/module_step_form.20220122.tsv 
        --etc_module_database_loc /work/ebg_lab/referenceDatabases/DRAM2/etc_mdoule_database.20220122.tsv 
        --function_heatmap_form_loc /work/ebg_lab/referenceDatabases/DRAM2/function_heatmap_form.20220122.tsv 
        --amg_database_loc /work/ebg_lab/referenceDatabases/DRAM2/amg_database.20220122.tsv
```

### To run DRAM or DRAM-v
```bash
#!/bin/bash
#SBATCH --job-name=dramv      # Job name
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=30            # Number of CPU cores per task
#SBATCH --mem=300G                    # Job memory request
#SBATCH --time=24:00:00              # Time limit hrs:min:sec
#SBATCH --output=dramv%j.log     # Standard output and error log
#SBATCH --partition=bigmem
pwd; hostname; date

source ~/miniconda3/etc/profile.d/conda.sh

conda activate DRAM

#this is prokaryote DRAM trial with "annotate" and "distill"
#Annotate
DRAM.py annotate -i ../../02_assembly/viral_genomes_reassmbly/combined_assembly/combined_Assembly_v2.fasta -o ./dram_combined_assembly/annotation --threads 30

#Distill
#check that trna and rrna files are actually created for each bin before using those options
DRAM.py distill -i ./dram_prok/annotation/annotations.tsv -o ./dram_prok/annotation/genome_summaries --trna_path ./dram_prok/annotation/trnas.tsv --rrna_path ./dram_prok/annotation/rrnas.tsv

```

To run DRAM-v for annotating viral genomes, you will first need to find some viral genomes!

```bash
#this is viral DRAM trial with "annotate" and "distill"
#Annotate
DRAM-v.py annotate -i ../../virsorter2/for_dramv/for-dramv/final-viral-combined-for-dramv.fa -v ../../virsorter2/for_dramv/for-dramv/viral-affi-contigs-for-dramv.tab -o ./dramv/2ndround --threads 30

#Distill
DRAM-v.py distill -i ./dramv/2ndround/annotations.tsv -o ./dramv/2ndround/genome_summaries

tRNAscan-SE -G -o ./dram_prok/annotation/working_dir/bin12/tmp/raw_trnas.txt --thread 20 ./dram_prok/annotation/working_dir/bin12/scaffolds.annotated.fa
```