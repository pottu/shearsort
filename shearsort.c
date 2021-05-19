#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#define ROOT 0
#define TAG  0

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

void print_slice(int w, int h, int32_t* slice) 
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  printf("Process: %d\n", rank);
  for (int i = 0; i < h; ++i) {
    for (int j = 0; j < w; ++j) {
      printf("%d ", slice[i*w + j]);
    }
    printf("\n");
  }
  fflush(stdout);
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
bool check_same_elements(int n, int32_t *M, char *input_file)
{
  int32_t *initial = read_input(&n, input_file);
  int32_t *result  = calloc(n*n, sizeof(*result));
  memcpy(result, M, n*n * sizeof(int32_t));

  qsort(initial, n*n, sizeof(int32_t), ascending);
  qsort(result, n*n, sizeof(int32_t), ascending);

  bool ret = true;
  for (int i = 0; i < n*n; i++) {
    if (result[i] != initial[i]) {
      ret = false;
      printf("Different elements detected!\n");
      break;
    }
  }
  free(initial);
  free(result);
  return ret;
}

bool check_sorted(int n, int32_t *M)
{
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

bool checker(int n, int32_t* M, char *input_file) 
{
  return check_sorted(n, M) && check_same_elements(n, M, input_file);
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

long partition(int32_t *data, int n, int col, long left, long right, long pivotIndex) {
  const int32_t pivotValue = data[pivotIndex*n+col];
  int32_t temp = data[pivotIndex*n+col];
  data[pivotIndex*n+col] = data[right*n+col];
  data[right*n+col] = temp;
  long storeIndex = left;

  for (long i = left; i < right; i++) {
    if (data[i*n+col] <= pivotValue) {
      temp = data[i*n+col];
      data[i*n+col] = data[storeIndex*n+col];
      data[storeIndex*n+col] = temp;
      storeIndex++;
    }
  }
  temp = data[storeIndex*n+col];
  data[storeIndex*n+col] = data[right*n+col];
  data[right*n+col] = temp;
  return storeIndex;
}

void quicksort(int32_t *data, int n, int col, long left, long right) {

  if (right > left) {
    long pivotIndex = left + (right - left) / 2;
    long pivotNewIndex = partition(data, n, col, left, right, pivotIndex);
    quicksort(data, n, col, left, pivotNewIndex - 1);
    quicksort(data, n, col, pivotNewIndex + 1, right);
  }
}

void sort_rows(int w, int h, int32_t *M)
{
  for (int row = 0; row < h; row++) {
    if (even(row)) {
      qsort(&M[row*w], w, sizeof(int32_t), ascending);
    } else {
      qsort(&M[row*w], w, sizeof(int32_t), descending);
    }
  }
}

void sort_columns(int w, int h, int32_t* M)
{
  for (int col = 0; col < w; col++) {
    quicksort(M, w, col, 0, h-1);
    //sort_column(col, w, h, M);
  }
}

bool less(int32_t a, int32_t b)
{
  return a < b;
}

bool greater_equal(int32_t a, int32_t b)
{
  return a >= b;
}

void merge(const int32_t *v1, int n1, const int32_t *v2, int n2, 
           int32_t *result, bool (*compare)(int32_t, int32_t)) 
{
  int i = 0;
  int j = 0;
  int k = 0;

  while (i < n1 && j < n2) {
    if (compare(v1[i], v2[j])) {
      result[k++] = v1[i++];
    } else {
      result[k++] = v2[j++];
    }
  }
  if (i == n1) {
    while (j < n2) {
      result[k++] = v2[j++];
    }
  } else {
    while (i < n1) {
      result[k++] = v1[i++];
    }
  }
}


void exchange_and_merge(int partner, bool even_rank, int w, int h, int32_t* M)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // TODO: stack array?
  int32_t* partner_slice = calloc(w * (h/2), sizeof(*partner_slice));

  // Handy types for sends/recvs.
  MPI_Datatype TYPE_ROW_SEND, TYPE_ROW_RECV;
  MPI_Type_vector(h/2, w, w*2, MPI_INT32_T, &TYPE_ROW_SEND);
  MPI_Type_vector(h/2, w, w, MPI_INT32_T, &TYPE_ROW_RECV);
  MPI_Type_commit(&TYPE_ROW_SEND);
  MPI_Type_commit(&TYPE_ROW_RECV);

  // Even ranks send even rows, odd ranks odd rows. ????
  int offset = even_rank ? w : 0;

  // TODO: tags?
  MPI_Sendrecv(
      M + offset,       // const void *sendbuf
      1,                // int sendcount
      TYPE_ROW_SEND,    // MPI_Datatype sendtype
      partner,          // int dest 
      TAG,              // int sendtag
      partner_slice,    // void *recvbuf
      1,                // int recvcount
      TYPE_ROW_RECV,    // MPI_Datatype recvtype
      partner,          // int source
      TAG,              // int recvtag
      MPI_COMM_WORLD,   // MPI_Comm comm
      MPI_STATUS_IGNORE // MPI_Status *status
  );



  int partner_row = 0;
  int32_t merged[w*2];
  for (int row = even_rank ? 0 : 1; row < h; row += 2) {
    size_t s = w * sizeof(int32_t);

    if (even_rank) {
      merge(&M[row * w], w, &partner_slice[partner_row * w], w, merged, less);
    } else {
      merge(&M[row * w], w, &partner_slice[partner_row * w], w, merged, greater_equal);
    }

    if (rank < partner) {
      memcpy(&M[row * w], merged, s);
      memcpy(&partner_slice[partner_row * w], &merged[w], s);
    } else {
      memcpy(&M[row * w], &merged[w], s);
      memcpy(&partner_slice[partner_row * w], merged, s);
    }
    partner_row += 1;
  }

  MPI_Sendrecv_replace(
      partner_slice,    // void *buf, 
      w * (h/2),        // int count, 
      MPI_INT32_T,      // MPI_Datatype datatype,
      partner,          // int dest, 
      TAG,              // int sendtag, 
      partner,          // int source, 
      TAG,              // int recvtag,
      MPI_COMM_WORLD,   // MPI_Comm comm, 
      MPI_STATUS_IGNORE // MPI_Status *status
  );

  partner_row = 0;
  for (int row = even_rank ? 1 : 0; row < h; row += 2) {
    memcpy(&M[row * w], &partner_slice[partner_row * w], w * sizeof(int32_t));
    partner_row += 1;
  }

  // TODO: free types, partner slice
  free(partner_slice);
  MPI_Type_free(&TYPE_ROW_SEND);
  MPI_Type_free(&TYPE_ROW_RECV);
}


void odd_even_sort(int w, int h, int32_t *slice, int rank, int num_PEs, bool even_rank)
{
  for (int i = 0; i < num_PEs; i++) {
    int partner = -1;
    if (i % 2 == 0) {
      partner = even_rank ? rank + 1 : rank - 1;
    } else {
      partner = even_rank ? rank - 1 : rank + 1;
    }
    if (partner >= 0 && partner < num_PEs) {
      exchange_and_merge(partner, even_rank, w, h, slice);
    }
  }
}

// ---- main -----------------------------------------------
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

  int n; // Size of matrix
  int32_t* M = NULL; // Input matrix

  // Root reads input matrix.
  if (rank == ROOT) {
    M = read_input(&n, input_file);
  }

  // Broadcast matrix size to all PEs.
  MPI_Bcast(&n, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

  // Width and height of matrix slices.
  const int w = n / num_PEs;
  const int h = n;

  // Each process' individual slice of columns.
  int32_t *slice = calloc(w * h, sizeof(*slice));
  
  // Type for rows for handy send/receives.
  MPI_Datatype TYPE_TMP, TYPE_COL_SEND, TYPE_COL_RECV;

  MPI_Type_vector(h, 1, h, MPI_INT32_T, &TYPE_TMP);
  MPI_Type_create_resized(TYPE_TMP, 0, sizeof(int32_t), &TYPE_COL_SEND);
  MPI_Type_commit(&TYPE_COL_SEND);
  MPI_Type_free(&TYPE_TMP);

  MPI_Type_vector(h, 1, w, MPI_INT32_T, &TYPE_TMP);
  MPI_Type_create_resized(TYPE_TMP, 0, sizeof(int32_t), &TYPE_COL_RECV);
  MPI_Type_commit(&TYPE_COL_RECV);
  MPI_Type_free(&TYPE_TMP);


  // ---- Shearsort ----------------------------------------

  // Scatter initial matrix column slices.
  MPI_Scatter(
      M,             // const void *sendbuf
      w,             // int sendcount
      TYPE_COL_SEND, // MPI_Datatype sendtype
      slice,         // void *recvbuf
      w,             // int recvcount
      TYPE_COL_RECV, // MPI_Datatype recvtype
      ROOT,          // int root
      MPI_COMM_WORLD // MPI_Comm comm
  );

  // Start timer.
  const double start = MPI_Wtime();

  int num_steps = ceil(log2(n)) + 1;
  for (int step = 0; step < num_steps; step++) {
    // Do for d steps..

    // Sort rows locally.
    sort_rows(w, h, slice);
    // Sort rows globally.
    odd_even_sort(w, h, slice, rank, num_PEs, even_rank);
    
    if (step < num_steps-1) {
      // Sort columns locally.
      sort_columns(w, h, slice);
    }
  }

  // Gather sorted slices.
  MPI_Gather(
      slice,         // const void *sendbuf,
      w,             // int sendcount,
      TYPE_COL_RECV, // MPI_Datatype sendtype,
      M,             // void *recvbuf,
      w,             // int recvcount,
      TYPE_COL_SEND, // MPI_Datatype recvtype,
      ROOT,          // int root,
      MPI_COMM_WORLD // MPI_Comm comm
  );

  // Stop timer.
  const double my_execution_time = MPI_Wtime() - start;

  // Find largest execution time.
  double slowest = 0;
  MPI_Reduce(
      &my_execution_time, 
      &slowest, 
      1,
      MPI_DOUBLE, 
      MPI_MAX, 
      0, 
      MPI_COMM_WORLD
  );

  // Output execution time.
  if (rank == ROOT) {
    printf("%lf\n", slowest);
  }

  // TODO: Check only if user wants to.
  if (rank == ROOT) {
    //printf("Sorted:\n");
    print_matrix_file(n, M, output_file);
    if (checker(n, M, input_file)) {
      printf("Correct!\n");
    } else {
      printf("Incorrect!\n");
    }
  }

  // Clean up.
  free(M);
  free(slice);
  MPI_Type_free(&TYPE_COL_SEND);
  MPI_Type_free(&TYPE_COL_RECV);
  MPI_Finalize();
  return 0;
}
