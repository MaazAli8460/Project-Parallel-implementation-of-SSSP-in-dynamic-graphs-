# Dynamic SSSP MPI

This project implements a parallel Dynamic Single-Source Shortest Path (SSSP) algorithm using MPI and OpenMP. It supports dynamic graph updates (edge insertions and deletions) and computes shortest paths from a source vertex (vertex 0) in a weighted undirected graph. The accuracy of the dynamic SSSP algorithm is evaluated using a path validation method, which verifies the correctness of the shortest path tree by comparing path weights against ground-truth distances computed by Dijkstra's algorithm.

## Table of Contents

* [Features](#features)
* [Prerequisites](#prerequisites)
* [Installation](#installation)
* [Usage](#usage)
* [Input Format](#input-format)
* [Output](#output)
* [Contributing](#contributing)
* [License](#license)

## Features

* **Dynamic Graph Updates:** Supports edge insertions and deletions in a weighted undirected graph.
* **Parallel Processing:** Uses MPI for distributed computing and OpenMP for multi-threaded parallelism within nodes.
* **Graph Partitioning:** Employs METIS to partition the graph for load balancing across MPI ranks.
* **Four-Phase Update Algorithm:**

  * Deletion Phase
  * Insertion Phase
  * Propagation Phase
  * Relaxation Phase
* **Accuracy Evaluation:** Validates shortest paths by tracing paths from each vertex to the source and comparing their weights against Dijkstra's ground-truth distances.
* **Timing Metrics:** Reports detailed timings for each phase, communication, CPU, and wall clock time.
* **Error Handling:** Includes robust checks for memory allocation, file access, and invalid inputs.

## Prerequisites

To compile and run the program, ensure the following are installed:

* **MPI Implementation:** OpenMPI or MPICH (e.g., openmpi).
* **OpenMP:** Supported by the compiler (e.g., gcc with OpenMP support).
* **METIS Library:** For graph partitioning (libmetis).
* **C Compiler:** GCC or compatible (gcc recommended).
* **Operating System:** Linux/Unix-based system.

Dependencies:

* Standard C libraries (libc, libm).
* POSIX-compliant system for timing functions.

Install dependencies on a Debian-based system:

```bash
sudo apt-get update
sudo apt-get install build-essential openmpi-bin libopenmpi-dev libomp-dev libmetis-dev
```

## Installation

1. **Clone or Download the Repository:**

   ```bash
   git clone <repository-url>
   cd dynamic-sssp-mpi
   ```

2. **Compile the Program:**
   Compile the source code using `mpicc` with OpenMP and METIS support:

   ```bash
   mpicc -g -O3 -std=c99 -fopenmp dynamic_sssp_mpi.c -lmetis -lm -o dsssp
   ```

3. **Set Up Hosts File:**
   Create a hosts file specifying the machines and slots for MPI execution. Example:

   ```bash
   cat > hosts << EOL
   pc1 slots=2
   pc2 slots=2
   EOL
   ```

   Replace `pc1` and `pc2` with your machine hostnames and adjust slots based on available cores.

4. **Copy Executable to Remote Nodes (if using multiple machines):**

   ```bash
   scp dsssp mpiuser1@pc2:~/mpi_programs/
   ```

   Ensure the executable is in the same path on all nodes (e.g., `~/mpi_programs/`).

## Usage

Run the program using `mpirun` with the following syntax:

```bash
mpirun --verbose -np <num_processes> --hostfile hosts ./dsssp <edgeFile> <numEdits> <percentAdd>
```

* `<num_processes>`: Number of MPI processes (e.g., 4).
* `<edgeFile>`: Path to the input graph file (see Input Format).
* `<numEdits>`: Total number of edge edits (insertions + deletions).
* `<percentAdd>`: Percentage of edits that are insertions (0–100). The rest are deletions.

Example:

```bash
mpirun --verbose -np 4 --hostfile hosts ./dsssp UD_W.txt 500 50
```

This runs the program with:

* 4 MPI processes.
* Input graph from `UD_W.txt`.
* 500 total edits (250 insertions, 250 deletions, as `percentAdd = 50`).

## Input Format

The input graph file (`<edgeFile>`) is a text file containing edge descriptions in the format:

```
u v w
```

* `u`: Source vertex ID (non-negative integer).
* `v`: Destination vertex ID (non-negative integer).
* `w`: Edge weight (positive integer).

Example (`UD_W.txt`):

```
0 1 5
0 2 3
1 2 4
...
```

Notes:

* The graph is undirected, so an edge `u v w` implies `v u w`.
* Vertex IDs are automatically resized to accommodate the maximum ID in the file.
* Weights must be positive (`w > 0`).

## Output

The program produces two types of output:

### Terminal Output:

* **Edit Summary:** Total edits, insertions, and deletions.
* **Graph Info:** Number of vertices and edges loaded.
* **Phase Timings (on rank 0):** Deletion, insertion, propagation, and relaxation + ghost exchange times.
* **Total Timings (on rank 0):** Communication time, CPU time, wall clock time, and cumulative CPU time across all ranks.
* **Per-Rank Timings:** Communication and CPU time for each rank, along with the hostname.

Example:

```
Total edits: 500  Ins: 250  Del: 250
Loaded 10879 vertices, 39994 edges
Phase timings (rank 0, host pc1):
  Deletion:   0.000 s
  Insertion:  0.001 s
  Propagate:  0.000 s
  Relax+Ghost:0.001 s
Total communication time (rank 0): 0.000 s
Total CPU time (rank 0): 3.641 s
Total wall clock time (rank 0): 3.641 s
Total cumulative CPU time (all ranks): 3.641 s
Rank 0 (host pc1):
  Comm time: 0.000 s
  CPU time:  3.641 s
Rank 1 (host pc1):
  Comm time: 0.000 s
  CPU time:  3.641 s
Rank 2 (host pc2):
  Comm time: 0.000 s
  CPU time:  3.641 s
Rank 3 (host pc2):
  Comm time: 0.000 s
  CPU time:  3.641 s
```

### Log File (`log.txt`):

* **SSSP Tree:** Lists each vertex, its distance from the source, and its parent in the shortest path tree.
* **Valid Paths:** Number of vertices with valid shortest paths.

Example:

```
Final SSSP Tree:
Vertex  Distance  Parent
0       0         -1
1       5         0
2       3         0
...
Vertices with Valid Paths: 10650
```
* **Why This Metric?**:

  * Validates the entire shortest path tree, not just distances, ensuring correctness of both distances and parent pointers.
  * Robust to graph fragmentation, as it evaluates each vertex's path independently.
  * Detects subtle errors (e.g., correct distance but incorrect path), providing a stricter and potentially more accurate assessment.

* If the accuracy is low, inspect `log.txt` to identify vertices with invalid paths, which may indicate issues in the dynamic SSSP algorithm (e.g., relaxation or ghost vertex synchronization).

## Contributing

Contributions are welcome! To contribute:

1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/your-feature`).
3. Commit your changes (`git commit -m 'Add your feature'`).
4. Push to the branch (`git push origin feature/your-feature`).
5.
