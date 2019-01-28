#!/bin/bash
#SBATCH -N 1 --tasks-per-node=1
#SBATCH -t 23:59:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=b.d.vanschaik@amc.uva.nl

module load stopos

STOPOS_POOL=d8c24f78f9772cbdff54cf62

# Start analysis
nohup python StoposGetTokens.py ${STOPOS_POOL} 

echo "FINISHED"
