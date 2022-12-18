#!/bin/bash
#SBATCH --job-name=downloadFromSRA               # Job name
#SBATCH --mail-type=END,FAIL                     # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=oadebali@sabanciuniv.edu     # Where to send mail	
#SBATCH --ntasks=1                               # Run on a single CPU
#SBATCH --mem=1gb                                # Job memory request
#SBATCH --time=02:00:00                          # Time limit hrs:min:sec
#SBATCH --partition=genomics                     # Partition name
#SBATCH --output=download_from_sra_%j.log        # Standard output and error log
pwd; hostname; date

module load sratoolkit/3.0.0

echo "Downloading SRA files"

fastq-dump --gzip -O ../../resources/samples SRR12775235
fastq-dump --gzip -O ../../resources/samples SRR12775234

date