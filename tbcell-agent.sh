#!/bin/bash
#SBATCH -N 1 --tasks-per-node=1
#SBATCH -t 23:59:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=b.d.vanschaik@amc.uva.nl

module load stopos
module load python
module load python/3.5.0
module load samtools
module load R

export STOPOS_POOL=d8c24f78f9772cbdff54cf62

# Start analysis
python2 StoposGetTokens.py ${STOPOS_POOL} 

echo "FINISHED"
