#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=3
#SBATCH --threads-per-core=1
#SBATCH -A b196                                                                          
#SBATCH -p skylake      
#SBATCH --time 10:00:00
#SBATCH --job-name=WEST_parall
#SBATCH -o log_MHDG_notIso_WEST_parall_%j.log
#SBATCH -e err_MHDG_NotIso_WEST_parall_%j.err
export OMP_NUM_THREADS=1

#mpirun -n 10 ./MHDG-NGammaTiTeNeutral-parall-2D ./Meshes/West_NoHole_Nel29960_P4 Sol2D_West_NoHole_Nel29960_P4_DPe0.200E+02_DPai0.314E+06_DPae0.105E+08
#mpirun -n 10 ./MHDG-NGammaTiTeNeutral-parall-2D ./Meshes/West_NoHole_Nel29960_P4 Sol2D_West_NoHole_Nel29960_P4_DPe0.650E+01_DPai0.314E+06_DPae0.105E+08
mpirun -n 10 ./MHDG-NGammaTiTeNeutral-parall-2D ./Meshes/West_NoHole_Nel29960_P4 Sol2D_West_NoHole_Nel29960_P4_DPe0.100E+01_DPai0.314E+06_DPae0.105E+08