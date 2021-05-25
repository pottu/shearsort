# shearsort
Parallel + distributed shear sort with OpenMPI.

## Build/Run Instructions
To build the program, run:
```
make
```
This should create an executable `shearsort` which can be used with `mpirun` e.g. as:
```
mpirun -n 2 ./shearsort <input-file>
```

## Input
The program reads a given input matrix from a file. The expected format is as follows:
1. A positive integer `n` specifying the matrix size (`n`x`n`), followed by a new line.
2. `n` lines of `n` integers each, separated by a space.

The python script `gen.py` can be used to generate random input files:
```
python3 gen.py <matrix size>
```

Some example input is found under the `input/` directory.

## Options
The program accepts the following option flags:
- `-c` tells the program to check if the solution is correct, outputting the result to stdout.
- `-s` suppresses output of the sorted matrix (default is to print it to stdout).
- `-o <output-file>` tells the program to write the sorted matrix to a given file.

Usage is as follows:
```
shearsort [-cs] [-o <output-file>] <input-file>
```
