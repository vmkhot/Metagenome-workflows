A template:  

```bash

#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=<threads>
#SBATCH --partition=<partitions,listed,here>
#SBATCH --mem=<>
#SBATCH --time=<hh:mm:ss>
#SBATCH --job-name=<job_name_id>
#SBATCH --output=<outdir>/<tool>.slurm.out.%j
#SBATCH --error=<outdir>/<tool>.slurm.err.%j

# activate the conda environment

source /path/to/miniconda3/etc/profile.d/conda.sh && conda activate <tool_env>

my commands here

```

An Example:

```bash

#!/bin/bash
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=short,standard
#SBATCH --mem=1G
#SBATCH --time=2:00:00
#SBATCH --job-name=sequencing_QC
#SBATCH --output=10_nanoplot/nanoplot.slurm.out.%j
#SBATCH --error=10_nanoplot/nanoplot.slurm.err.%j

# First, we check the quality of the reads using nanoplot
# activate the conda environment containing nanoplot

source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh && conda activate nanoplot_v1.41.3

```
