#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#define INF 1000000

int N;
int *mat;

void abort_with_error_message(const char *msg) {
    fprintf(stderr, "%s\n", msg);
    exit(EXIT_FAILURE);
}

int convert_dimension_2D_1D(int x, int y, int n) {
    return x * n + y;
}

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

void bellman_ford(int my_rank, int p, MPI_Comm comm, int n, int *mat, int *dist, bool *has_negative_cycle) {
    int loc_n;
    int loc_start, loc_end;
    int *loc_mat;
    int *loc_dist;

    if (my_rank == 0) {
        loc_n = n;
    }
    MPI_Bcast(&loc_n, 1, MPI_INT, 0, comm);

    int ave = loc_n / p;
    loc_start = ave * my_rank;
    loc_end = ave * (my_rank + 1);
    if (my_rank == p - 1) {
        loc_end = loc_n;
    }

    loc_mat = (int *)malloc(loc_n * loc_n * sizeof(int));
    loc_dist = (int *)malloc(loc_n * sizeof(int));

    if (my_rank == 0)
        memcpy(loc_mat, mat, sizeof(int) * loc_n * loc_n);
    MPI_Bcast(loc_mat, loc_n * loc_n, MPI_INT, 0, comm);

    for (int i = 0; i < loc_n; i++) {
        loc_dist[i] = INF;
    }
    loc_dist[0] = 0;
    MPI_Barrier(comm);

    bool loc_has_change;
    int loc_iter_num = 0;
    for (int iter = 0; iter < loc_n - 1; iter++) {
        loc_has_change = false;
        loc_iter_num++;
        for (int u = loc_start; u < loc_end; u++) {
            for (int v = 0; v < loc_n; v++) {
                int weight = loc_mat[convert_dimension_2D_1D(u, v, loc_n)];
                if (weight < INF) {
                    if (loc_dist[u] + weight < loc_dist[v]) {
                        loc_dist[v] = loc_dist[u] + weight;
                        loc_has_change = true;
                    }
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, &loc_has_change, 1, MPI_C_BOOL, MPI_LOR, comm);
        if (!loc_has_change)
            break;
        MPI_Allreduce(MPI_IN_PLACE, loc_dist, loc_n, MPI_INT, MPI_MIN, comm);
    }

    if (loc_iter_num == loc_n - 1) {
        loc_has_change = false;
        for (int u = loc_start; u < loc_end; u++) {
            for (int v = 0; v < loc_n; v++) {
                int weight = loc_mat[convert_dimension_2D_1D(u, v, loc_n)];
                if (weight < INF) {
                    if (loc_dist[u] + weight < loc_dist[v]) {
                        loc_dist[v] = loc_dist[u] + weight;
                        loc_has_change = true;
                        break;
                    }
                }
            }
        }
        MPI_Allreduce(&loc_has_change, has_negative_cycle, 1, MPI_C_BOOL, MPI_LOR, comm);
    }

    if (my_rank == 0)
        memcpy(dist, loc_dist, loc_n * sizeof(int));

    free(loc_mat);
    free(loc_dist);
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
