# Annotation Genomes/MAGs

The first step in annotation is gene prediction. but now pipelines exist which do gene prediction followed by annotating the predicted genes

## [Metaerg 2.0](https://github.com/kinestetika/MetaErg)

To download on a HPC cluster - use [apptainer](https://apptainer.org/docs/admin/main/installation.html) (formlerly singularity) to pull a Docker image and convert to a singularity image (.sif file)

```bash
# pull docker image and convert to singularity
apptainer pull docker://kinestetika/metaerg
```

To run metaerg, you have to enter into the apptainer container (like an environment?). shell command will open a new shell into the container, where metaerg lives and you can interact like a virtual environment.

```bash
apptainer shell ~/metaerg_latest.sif

# to check the download
metaerg -h

# to exit 
cntrl + d
```

By default, apptainer only "binds" your $HOME directory. So to access other directories in the filesystem above that you can run the following command.

It's **not permanent** so you have to "-B" everytime you run the apptainer container. Haven't figured out how to do this yet. Something to do with the "export APPTAINER_BIND" environment variable.

```bash
apptainer shell -B /work/ebg_lab/ metaerg_latest.sif
```

Then have to downlaod the databases. On the HPC cluster, I had to use a sbatch script to do this. Takes 30min+

```bash
 # to download database
 apptainer exec -B /work/ebg_lab/ ~/data/Programs/metaerg_latest.sif metaerg --download_database --database_dir /work/ebg_lab/referenceDatabases metaerg_database/

# to update the database
 apptainer exec -B /work/ebg_lab/ ~/data/Programs/metaerg_latest.sif metaerg --create_database S --database_dir /work/ebg_lab/referenceDatabases metaerg_database/
```

To run metaerg/apptainer from slurm batch script

```bash
#!/bin/bash
#SBATCH --job-name=metaerg2_db_dl      
#SBATCH --nodes=1                    
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=10            
#SBATCH --mem=100G                    
#SBATCH --time=10:00:00              
#SBATCH --output=metaerg2_db_dl%j.log
#SBATCH --partition=bigmem,cpu2019,cpu2021
pwd; hostname; date


apptainer exec -B /work/ebg_lab/ ~/data/Programs/metaerg_latest.sif metaerg --download_database --database_dir /work/ebg_lab/referenceDatabases/metaerg_database/
```
