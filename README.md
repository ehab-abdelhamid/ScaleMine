# ScaleMine: Scalable Parallel Frequent Subgraph Mining in a Single Large Graph

## Overview:

ScaleMine is a novel parallel frequent subgraph mining system for a single large
graph. ScaleMine introduces a novel two-phase approach. The first phase is 
approximate; it quickly identifies subgraphs that are frequent with high 
probability, while collecting various statistics. The second phase computes the
exact solution by employing the results of the approximation to achieve good 
load balance; prune the search space; generate efficient execution plans; and 
guide intra-task parallelism.

If you use ScaleMine in your research, please cite our paper:
 ```
@inproceedings{hamid2016scalemine,
  title={ScaleMine: Scalable Parallel Frequent Subgraph Mining in a Single Large Graph},
  author={Ehab Abdelhamid, Ibrahim Abdelaziz, Panos Kalnis, Zuhair Khayyat and Fuad Jamour},
  booktitle={Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis},
  pages={11},
  year={2016},
  organization={ACM}
}
```
## Contents:

    README ...................  This file
    LICENSE.txt ..............  License file (Open Source)
    makefile .................  build ScaleMine binary files
    Datasets/ ................  Example graphs
    SRC_*/ ...................  Directory containing ScaleMine source files


## Dependencies
 There are a few dependencies which must be satisfied in order to compile and run ScaleMine.
 
 * build-essential and g++ (>= 4.4.7) [Required]
   +  Needed for compiling ScaleMine.
 
 * openssh-server [Required]
    + Required to initialize MPI and establish connections among compute nodes.
 
 * MPICH2 [Required]
    + ScaleMine uses MPI for inter-node communication. Open MPI is not tested with ScaleMine.

## Installation:

- install MPI and Boost libraries on the target machine
- Uncompress ScaleMine using any compression tool
- Build ScaleMine using the "compile.sh" script file

## Running:
### Single Machine Mode:
Run ScaleMine using the following command:
```
mpirun -n N pfsm -file GRAPH_FILE -freq F -threads T
```

N: the number of MPI computation nodes, make sure that there is at lease one 
for th emaster and one for a worker. Best practice is to have one computation
node at each separate machine, then for each machine set a number of parallel
threads.

GRAPH_FILE: the input graph file name, the supported graph format is .lg

F: the user-give support threshold

T: the number of threads per compute node

Example:

Use the following command:
```
mpirun -n 2 pfsm -file ./Datasets/patent_citations.lg -freq 28000 -threads 4
```
to mine the patent_citation graph for subgraphs having support larger than or
equal to 28000 using 2 compute nodes; 1 master and 1 worker, the worker has 4
threads.

### Distributed Mode
For running on a supercomputer using SLURM job scheduler:
```
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
```

## Output:

ScaleMine outputs the list of frequents subgraphs on the standard output.
Also, the elapsed time is returned at the end.

## License:
ScaleMine is a free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 2 of the License, or (at your option) any later version.

ScaleMine is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with ScaleMine.  If not, see <http://www.gnu.org/licenses/>.

