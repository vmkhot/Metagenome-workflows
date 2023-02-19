# Annotation Genomes/MAGs

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