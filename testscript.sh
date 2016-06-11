# module rm python/2.7.11
module load python/3.3.4
# module rm mpi4py/1.3.1_py-2.7
module load R/3.2.3
module load mpi4py/1.3.1
module list

# mpirun -np 11 python3 ~/bioinformatics/clusterfit.py
# mpirun -np 60 python3 clusterfit.py --input hzar_10loci.csv --locinames hzar_100loci_names.txt --output small_output/testiter1 --repeat 1
mpirun -np 120 python3 clusterfit.py --input ref7big.csv --locinames ref7small.txt --output small_output/oldoutput --repeat 1 --modelIII False --chainlength 1e5