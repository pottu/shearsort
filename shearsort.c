#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>

#define ROOT 0

// ---- INPUT / OUTPUT -------------------------------------
void bad_input() {
  printf("Bad input file.\n");
  MPI_Abort(MPI_COMM_WORLD, 1);
  exit(0);
}

int32_t* read_input(int* out_n, char* input_file) {
  FILE* infile = fopen(input_file, "r");
  if (!infile) {
    bad_input();
  }
  if (fscanf(infile, "%d", out_n) != 1) {
    bad_input();
  }
  int n = *out_n;
  int32_t* M = calloc(n*n, sizeof(int32_t));
  for (int i = 0; i < n*n; ++i) {
    if (fscanf(infile, "%d", &M[i]) != 1) {
      bad_input();
    }
  }
  return M;
}

void print_matrix(int n, int32_t* M) {
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      printf("%d ", M[row*n + col]);
    }
    printf("\n");
  }
}

void print_matrix_file(int n, int32_t* M, char* output_file) {
  FILE *outfile = fopen(output_file, "w");
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      fprintf(outfile, "%d ", M[row*n + col]);
    }
    fprintf(outfile, "\n");
  }
  fclose(outfile);
}

// ---- HELPERS --------------------------------------------
bool even(int n) {
  return n % 2 == 0;
}

bool odd(int n) {
  return !even(n);
}

// ---- FUNCTIONALITY --------------------------------------
bool checker(int n, int32_t* M) {
  int32_t prev = M[0];
  for (int row = 0; row < n; row++) {
    if (even(row)) { // Even row: go left to right
      for (int col = 0; col < n; col++) {
        if (prev > M[row*n + col]) {
          return false;
        }
        prev = M[row*n + col]; 
      }
    }
    else { // Odd row: go right to left.
      for (int col = n-1; col >= 0; col--) {
        if (prev > M[row*n + col]) {
          return false;
        }
        prev = M[row*n + col]; 
      }
    }
  }
  return true;
}



int main(int argc, char **argv) {
  if (3 != argc) {
    printf("Usage: shearsort <input_file> <output_file>\n");
    return 0;
  }
  char* input_file = argv[1];
  char* output_file = argv[2];

  MPI_Init(&argc, &argv);
  
  // Find own rank and number of PEs.
  int numPEs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numPEs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int n; // Size of matrix
  int32_t* M = NULL; // Input matrix

  // Root reads input matrix.
  if (rank == ROOT) {
    M = read_input(&n, input_file);
  }

  // Broadcast matrix size to all PEs.
  MPI_Bcast(&n, 1, MPI_INT32_T, ROOT, MPI_COMM_WORLD);
  
  int num_steps = ceil(log2(n)) + 1;
  for (int step = 0; step < num_steps; step++) {
     if (even(step)) {
     }

     else {
     }
  }

  if (rank == ROOT) {
    if (checker(n, M)) {
      printf("Correct!\n");
    } else {
      printf("Incorrect!\n");
    }
  }

  // Clean up.
  if (rank == ROOT) {
    free(M);
  }
  MPI_Finalize();
  return 0;
}
