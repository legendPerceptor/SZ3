# SZ3 layer-by-layer compression: high performance parallel compression

This repo serves as an add-on to the original SZ3 compression method.
The traditional compression methods require loading the whole dataset into memory,
resulting in out-of-memory (OOM) failures. Also, there was no way to compress a single tensor in parallel.
We develop a layer-by-layer compression technique that allows users to compress very large 3D tensors with
limited amount of memory. It can also scale up to multiple nodes with non-uniform memory access (NUMA)
to compress huge files in short amount of time.

The executable file lies in `tools/sz3-split/sz3_split` after built with CMake. Below is the help information.
Users can get the help information with `sz3_split compress --help`. We will demonstrate some example usage later.

```c++
"Usage: sz3_split (de)compress [options]\n"
"options:  --threads/-t     INT   number of threads, default is 1, for MPI mode, this specifies the number of I/O processes\n"
"          --input/-i       STR   the RAW file/compressed file\n"
"          --output/-o      STR   the compressed file/decompressed file location\n"
"          --help/-h              print this help information\n"
"          --dimension/-d   STR   the data dimension of the file, e.g., 256 256 512\n"
"          --errorbound/-e  FLOAT the error bound to use in compression\n"
"          --float64              the default is float32 for each datapoint, this param changes it to float64\n"
"          --mode           STR   select 'layer' for layer-by-layer compression, 'direct' for direct compression\n"
"          --depth          INT   select the layer depth in layer-by-layer compression\n"
"          --mpi                  use MPI to run multiple processes";
```

For example, we can compress the NYX dataset layer-by-layer although the dimension is only (512, 512, 512). The following command
compresses the file with layers of depth 32 and there is one read process and one write process.

The NYX dataset can be downloaded on [this website](https://sdrbench.github.io/)

```bash
mpirun -n 4 sz3_split compress -i ./temperature.f32 -o ./temperature.f32.szsplit -d 512 512 512 -e 0.01 --mode layer --depth 32 --mpi --threads 2
```

The decompression doesn't have to use the same number of processes but the depth and error bounds need to be set the same.
For example, the command below uses 8 processes and there are 4 I/O processes and 4 compute processes.

```bash
mpirun -n 8 sz3_split decompress -i ./temperature.f32.szsplit -o ./temperature.f32.szsplit.dp -d 512 512 512 -e 0.01 --mode layer --depth 32 --mpi --threads 4
```

If the files are stored in a parallel file system, the performance can usually get better with more I/O processes. We recommend
setting the ratio of I/O processes to compute processes at around 1:8.

Our program can also use multi-threading mode by not setting the `--mpi` parameter. In this mode,
there is always one read thread and one write thread, the `--threads` parameter specifies the number of compute threads.
One special case is `--threads 1`, where read, compute, and write all happen in the main thread and there is no multi-threading.

For example, we can use 6 compute threads to do layer-by-layer compression in a single node, and then decompress the file with 4 threads.

```bash
sz3_split compress -i ./temperature.f32 -o ./temperature.f32.szsplit -d 512 512 512 \
-e 0.01 --mode layer --depth 32 --threads 6
sz3_split decompress -i ./temperature.f32.szsplit -o ./temperature.f32.szsplit.dp -d 512 512 512 \
-e 0.01 --mode layer --depth 32 --threads 4
```

Note that no matter how many threads/processes we use, the compressed file will be the same, and decompression
can always be done by any number of threads/processes, as long as the depth and error bound stay the same.

To utilize multiple nodes to compress a very large file, we usually need to submit a batch job to a supercomputer.
We provide an example to use 8 nodes to compress a 960GB file below. The dimension of the tensor is 10240x7680x1536. The script
is tested on Purdue Anvil machine, but should work on any machines with Slurm systems. Remember to change the file paths, partition, and account information before submitting the job.

Users can download the huge file on [this website](https://klacansky.com/open-scivis-datasets/category-simulation.html)

```bash
#!/bin/bash
#SBATCH --job-name=mpisz3
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:10:00
#SBATCH -A cis220161
#SBATCH -p wholenode
#SBATCH -o ./large-8node-32each-1.o      # Name of stdout output file
#SBATCH -e ./large-8node-32each-1.e      # Name of stderr error file


# Load MPI module
module load openmpi

executable="/home/x-yliu4/sz3-yuanjian/official_sz3/SZ3-Split/build/tools/sz3-split/sz3_split"
datafilename="dns_10240x7680x1536_float64.raw"
datafullpath="/anvil/projects/x-cis220161/datasets/large-single-file-data/RAW-Files/$datafilename"
compressedfile="/anvil/scratch/x-yliu4/multinode/$datafilename.sz3split"
decompressedfile="/anvil/scratch/x-yliu4/multinode/$datafilename.sz3split.dp"
# Run MPI program

echo "!!!compress test started!!!"
date
mpiexec -n 256 $executable compress -i $datafilename -o $compressedfile -d 10240 7680 1536 -e 0.01 --mode layer --depth 4 --threads 32 --mpi
date
echo "!!!decompress test started!!!"
date
mpiexec -n 256 $executable decompress -i $compressedfile -o $decompressedfile -d 10240 7680 1536 -e 0.01 --mode layer --depth 4 --threads 32 --mpi
echo "!!!decompress test finished!!"
date
```

We created Python binding for the program for ease of use in the Python language. Users can expect very similar usage of the above programs in Python. More detailed description can be found in [tools/bindings/README.md](tools/bindings/README.md).
