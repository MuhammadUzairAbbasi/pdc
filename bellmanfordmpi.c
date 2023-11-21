#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <mpi.h>

#define INF 1000000

// Define constants and global variables
int N;
int *mat;

// Function to handle errors and abort
void abort_with_error_message(const char *msg) {
    fprintf(stderr, "%s\n", msg);
    exit(EXIT_FAILURE);
}

// Convert 2D coordinates to 1D
int convert_dimension_2D_1D(int x, int y, int n) {
    return x * n + y;
}

// Read input from file
int read_file(const char *filename) {
    FILE *inputf = fopen(filename, "r");
    if (inputf == NULL) {
        abort_with_error_message("ERROR OCCURRED WHILE READING INPUT FILE");
    }
    fscanf(inputf, "%d", &N);
    assert(N < (1024 * 1024 * 20));
    mat = (int *)malloc(N * N * sizeof(int));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fscanf(inputf, "%d", &mat[convert_dimension_2D_1D(i, j, N)]);
        }
    }
    fclose(inputf);
    return 0;
}

// Print results to output file
int print_result(bool has_negative_cycle, int *dist) {
    FILE *outputf = fopen("output.txt", "w");
    if (!has_negative_cycle) {
        for (int i = 0; i < N; i++) {
            if (dist[i] > INF) {
                dist[i] = INF;
            }
            fprintf(outputf, "%d\n", dist[i]);
        }
        fflush(outputf);
    } else {
        fprintf(outputf, "FOUND NEGATIVE CYCLE!\n");
    }
    fclose(outputf);
    return 0;
}

// Bellman-Ford algorithm
void bellman_ford(int my_rank, int p, MPI_Comm comm, int n, int *mat, int *dist, bool *has_negative_cycle) {
    // Implementation of the Bellman-Ford algorithm remains unchanged
    // ...
}

int main(int argc, char **argv) {
    if (argc <= 1) {
        abort_with_error_message("INPUT FILE WAS NOT FOUND!");
    }
    const char *filename = argv[1];

    int *dist;
    bool has_negative_cycle = false;

    MPI_Init(&argc, &argv);
    MPI_Comm comm;

    int p;
    int my_rank;
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &my_rank);

    if (my_rank == 0) {
        assert(read_file(filename) == 0);
        dist = (int *)malloc(sizeof(int) * N);
    }

    double t1, t2;
    MPI_Barrier(comm);
    t1 = MPI_Wtime();

    bellman_ford(my_rank, p, comm, N, mat, dist, &has_negative_cycle);
    MPI_Barrier(comm);

    t2 = MPI_Wtime();

    if (my_rank == 0) {
        printf("Time(s): %.6f\n", t2 - t1);
        print_result(has_negative_cycle, dist);
        free(dist);
        free(mat);
    }
    MPI_Finalize();
    return 0;
}
