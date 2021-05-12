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

void print_slice(int n, int num_rows, int32_t* slice, int rank) 
{
  printf("Process: %d\n", rank);
  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0; j < n; ++j) {
      printf("%d ", slice[i*n + j]);
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

int ascending(const void* a, const void* b) {
  int32_t x = *((int32_t*)a);
  int32_t y = *((int32_t*)b);
  return (x > y) - (x < y);
}

int descending(const void* a, const void* b) {
  int32_t x = *((int32_t*)a);
  int32_t y = *((int32_t*)b);
  return (x < y) - (x > y);
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


void sort_column(int col, int w, int h, int32_t* M) {
  bool swapped = true;
  while (swapped) {
    swapped = false;
    for (int i = 0; i < h-1; i++) {
      if (M[i*w + col] > M[(i+1)*w + col]) {
        int32_t tmp = M[i*w + col];
        M[i*w + col] = M[(i+1)*w + col];
        M[(i+1)*w + col] = tmp;
        swapped = true;
      }
    }
  }
}

void sort_columns(int w, int h, int32_t* M)
{
  for (int col = 0; col < w; col++) {
    sort_column(col, w, h, M);
  }
}

void exchange_and_merge(int partner, bool even_rank, int w, int h, int32_t* M)
{
}


// ---- MAIN -----------------------------------------------
int main(int argc, char **argv) 
{
  if (3 != argc) {
    printf("Usage: shearsort <input_file> <output_file>\n");
    return 0;
  }
  char* input_file = argv[1];
  char* output_file = argv[2];

  MPI_Init(&argc, &argv);
  
  // Find own rank and number of PEs.
  int num_PEs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &num_PEs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // Track if own rank is odd or even.
  bool even_rank = rank % 2 == 0; 
  
  // Rank of paired partner process for exchanges.
  int partner = even_rank ? rank + 1 : rank - 1;

  int n; // Size of matrix
  int32_t* M = NULL; // Input matrix

  // Root reads input matrix.
  if (rank == ROOT) {
    M = read_input(&n, input_file);
  }

  // Broadcast matrix size to all PEs.
  MPI_Bcast(&n, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

  const int num_rows = n / num_PEs;

  // Type for rows for handy send/receives.
  MPI_Datatype row_type;
  MPI_Type_contiguous(n, MPI_INT32_T, &row_type);
  MPI_Type_commit(&row_type);

  // Each process' individual slice of rows.
  int32_t *slice = calloc(num_rows * n, sizeof(*slice));


  // ---- Shearsort ----------------------------------------

  // Scatter initial matrix slices.
  MPI_Scatter(M, num_rows, row_type, slice,
              num_rows, row_type, ROOT, MPI_COMM_WORLD);

  print_slice(n, num_rows, slice, rank);
  
  int num_steps = ceil(log2(n)) + 1;
  for (int step = 0; step < num_steps; step++) {
    // Do for d steps..

    for (int row = 0; row < num_rows; row++) {
      // Sort rows.
      if (even(row)) {
        qsort(&slice[row*n], n, sizeof(int32_t), ascending);
      } else {
        qsort(&slice[row*n], n, sizeof(int32_t), descending);
      }
    }

    // IDEA:
    // 1. Sort columns locally,
    // 2. Pairwise exchange columns
    // 3. Merge columns
    // 4. Return new split
    sort_columns(n, num_rows, slice);
    exchange_and_merge(partner, even_rank, n, num_rows, slice);
  }

  print_slice(n, num_rows, slice, rank);


  if (rank == ROOT) {
    if (checker(n, M)) {
      printf("Correct!\n");
    } else {
      printf("Incorrect!\n");
    }
  }

  // Clean up.
  if (rank == ROOT) {
    //print_matrix(n, M);
    free(M);
  }
  MPI_Finalize();
  return 0;
}
