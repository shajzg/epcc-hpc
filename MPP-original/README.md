# MSc Message-passing Programming C program to 2-d decomposition parallel percolate

## Requirements

* ICC compiler.
* MPI

### Cirrus

To get the ICC compiler and MPI on Cirrus, run:

```console
$ module load intel-compilers-19
$ module load mpt
```

---

## Program structure

source file: arralloc.c percio.c unirand.c percolate.c

header file: arralloc.h percolate.h 

encapsulation of functionalities:

* arralloc.h and arralloc.c
void *arralloc(size_t size, int ndim, ...); 

* percolate.h and percio.c
void percwrite(char *percfile, int **map, int ncluster,int L);
void percwritedynamic(char *percfile, int **map, int l, int ncluster);

* percolate.h and percolate.c
void initialize_oldmap(int **old, int M,int L);
void initialize_map(int **map, int L, double r, double rho, int *nhole);
void compute_newmap(int **old, int **new, int newval, int oldval, int M, int N, int *nchangelocal);
void test_percolation(int **map, int L, int *perc);

---

## Compilation

Compiling command is included in the file Makefile:

```console
$ make Makefile or make
```

---

## Usage

To run the Percolate program:


Run the percolate executable file:

```console
$ mpirun [-n NUM_PROCESS] ./percolate [-l LENGTH] [-s SEED]
```

For example, to run with default values for the parameters:

```console
$ mpirun -n 16 ./percolate -l 432 -s 1564 
```

The program will output information on its operation. For example:

```
percolate: running on 16 process(es)
percolate: L = 432, rho = 0.411000, seed = 1564
percolate: rho = 0.411000, actual density = 0.409326
percolate: number of changes on step 100 is 33296
percolate: number of changes on step 200 is 24362
percolate: number of changes on step 300 is 17327
percolate: number of changes on step 400 is 16811
percolate: number of changes on step 500 is 15653
percolate: number of changes on step 600 is 12434
percolate: number of changes on step 700 is 7296
percolate: number of changes on step 800 is 4239
percolate: number of changes on step 900 is 2071
percolate: number of changes on step 1000 is 506
 time per step is 0.000148386
percolate: cluster DOES NOT percolate
percwrite: visualising the largest 8 clusters
percwrite: cluster sizes are 35196, 12253, 7583, 6422, 4274, 2736, 2084, 1271
percwrite: opening file <map.pgm>
percwrite: writing data ...
percwrite: ... done
percwrite: file closed
```

Other examples of running the program include:

```console
$ mpirun -n 25 ./percolate -s 1564
$ mpirun -n 4 ./percolate -s 1564
$ mpirun -n 1 ./percolate 
```

### Command-line parameters

| Parameter | Description | Default Value |
| --------- |------------ | ------------- |
| -l | Map width/height (length). length >= 0 | 432 |
| -s | Seed for random number generator (seed). Using the same seed allows the same random numbers to be generated every time the program is run. seed>=0 | 6543 |
| -n | Number of processes to use. size > 0 (Must be typed subjected to MPI) | no default value |
### PGM output files

A Portable Grey Map (PGM) file is output.

This file does not include the halo as the use of a halo is an implementation detail.

This file is plain-text so you can view it as you would any plain-text file e.g.:

```console
$ cat map.pgm
```

PGM files can be viewed graphically using ImageMagick commands as follows.

Cirrus users will need first need to run:

```console
$ module load ImageMagick
```

To view a PPM file, run:

```console
$ display map.pgm
```

