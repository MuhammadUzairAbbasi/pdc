#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#define INF 1000000

int N;
int *mat;

void abort_with_error_message(char *msg) {
    fprintf(stderr, "%s\n", msg);
    abort();
}

int convert_dimension_2D_1D(int x, int y, int n) {
    return x * n + y;
}

int read_file(char *filename) {
    FILE *inputf = fopen(filename, "r");
    if (inputf == NULL) {
        abort_with_error_message("ERROR OCCURRED WHILE READING INPUT FILE");
    }
    fscanf(inputf, "%d", &N);

    assert(N < (1024 * 1024 * 20));
    mat = (int *) malloc(N * N * sizeof(int));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fscanf(inputf, "%d", &mat[convert_dimension_2D_1D(i, j, N)]);
        }
    }
    fclose(inputf);
    return N;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        abort_with_error_message("Usage: ./program <input_file>");
    }

    char *filename = argv[1];
    int N = read_file(filename);

    MPI_Finalize();
    return 0;
}

