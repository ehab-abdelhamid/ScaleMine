# ScaleMine
Parallel Frequent Subgraph Mining

This is the ScaleMine Source Release README
Last updated for ScaleMine on 9 August, 2016
-----------------------------------------------------------------------------

OVERVIEW:

ScaleMine is a novel parallel frequent subgraph mining system for a single large
graph. ScaleMine introduces a novel two-phase approach. The first phase is 
approximate; it quickly identifies subgraphs that are frequent with high 
probability, while collecting various statistics. The second phase computes the
exact solution by employing the results of the approximation to achieve good 
load balance; prune the search space; generate efficient execution plans; and 
guide intra-task parallelism.

For more details, check our paper:
Ehab Abdelhamid, Ibrahim Abdelaziz, Panos Kalnis, Zuhair Khayyat and Fouad
Jamour, 2016, “ScaleMine: Scalable Parallel Frequent Subgraph Mining in a Single
Large Graph”, Proceedings of the International Conference for High Performance
Computing, Networking, Storage and Analysis (Supercomputing'16)

CONTENTS:

    README ...................  This file
    LICENSE.txt ..............  License file (Open Source)
    makefile .................  build ScaleMine binary files
    Datasets/ ................  Example graphs
    SRC_*/ ...................  Directory containing ScaleMine source files


REQUIREMENTS:

Java JRE v1.6.0 or later

INSTALLATION:

- install MPI and Boost libraries on the target machine
- Uncompress ScaleMine using any compression tool
- Build ScaleMine using the "compile.sh" script file

Running:
On a single machine:
Run the tool using the following command:
mpirun -n N pfsm -file GRAPH_FILE -freq F -threads T

N: the number of MPI computation nodes, make sure that there is at lease one 
for th emaster and one for a worker. Best practice is to have one computation
node at each separate machine, then for each machine set a number of parallel
threads.
GRAPH_FILE: the input graph file name, the supported graph format is .lg
F: the user-give support threshold
T: the number of threads per compute node

EXAMPLE:

Use the following command:
mpirun -n 2 pfsm -file ./Datasets/patent_citations.lg -freq 28000 -threads 4

to mine teh patent_citation graph for subgraphs having support larger than or
equal to 28000 using 2 compute nodes; 1 master and 1 worker, the worker has 4
threads.

For running on a supercomputer using SLURM Job Scheduler:
#!/bin/bash
#SBATCH --account=user-account
#SBATCH --job-name=pfsm
#SBATCH --output=/output/mining.out
#SBATCH --error=/output/mining.err
#SBATCH --time=01:00:00
#SBATCH --threads-per-core=1
#SBATCH --nodes=256
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=32
export OMP_NUM_THREADS=32
export MKL_NUM_THREADS=32
export MPICH_MAX_THREAD_SAFETY=multiple
srun --ntasks=256 pfsm -file /Datasets/patent_citations.lg -freq 28000 -threads 32

OUTPUT:

ScaleMine outputs the list of frequents subgraphs on the standard output.
Also, the elapsed time is returned at the end.
