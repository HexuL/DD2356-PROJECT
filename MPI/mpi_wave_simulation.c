// #include <mpi.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
// #include <stdbool.h>

// // Function prototypes
// void initialize_arrays(double *U, bool *mask, int N, int local_N, int rank, int size);
// void apply_boundary_conditions(double *U, double *xlin, bool *mask, int N, int local_N, double t, int rank, int size);
// void update_arrays(double *U, double *Uprev, double *laplacian, int local_N, int N, double fac);
// void compute_laplacian(double *U, double *laplacian, int local_N, int N, int rank, int size);

// #ifndef M_PI
// #define M_PI 3.14159265358979323846
// #endif

// int main(int argc, char *argv[]) {
//     MPI_Init(&argc, &argv);
//     int rank, size;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     // Simulation parameters
//     int N = 256; // resolution
//     double boxsize = 1.0; // box size
//     double c = 1.0; // wave speed
//     double tEnd = 2.0; // stop time
//     bool plotRealTime = false; // switch for plotting simulation in real time
//     int num_iterations = 6; // number of times to run the simulation
//     double total_computation_time = 0.0; // to accumulate total computation time over iterations

//     // Mesh
//     double dx = boxsize / N;
//     double dt = (sqrt(2) / 2) * dx / c;
//     double fac = dt * dt * c * c / (dx * dx);

//     int local_N = N / size; // Local size of the grid

//     for (int iter = 0; iter < num_iterations; iter++) {
//         double *U = (double *)calloc(local_N * N, sizeof(double));
//         double *Uprev = (double *)calloc(local_N * N, sizeof(double));
//         bool *mask = (bool *)calloc(local_N * N, sizeof(bool));
//         double *laplacian = (double *)calloc(local_N * N, sizeof(double));
//         double *xlin = (double *)malloc(N * sizeof(double));

//         if (U == NULL || Uprev == NULL || mask == NULL || laplacian == NULL || xlin == NULL) {
//             fprintf(stderr, "Error: Unable to allocate memory\n");
//             MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
//         }

//         for (int i = 0; i < N; i++) {
//             xlin[i] = 0.5 * dx + i * dx;
//         }

//         initialize_arrays(U, mask, N, local_N, rank, size);

//         double t = 0.0; // time for this iteration

//         // Start the timer
//         double start_time = MPI_Wtime();

//         // Simulation Main Loop
//         while (t < tEnd) {
//             compute_laplacian(U, laplacian, local_N, N, rank, size);
//             update_arrays(U, Uprev, laplacian, local_N, N, fac);
//             apply_boundary_conditions(U, xlin, mask, N, local_N, t, rank, size);

//             // update time
//             t += dt;

//             MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes after each iteration
//         }

//         // Stop the timer
//         double end_time = MPI_Wtime();
//         double iteration_time = end_time - start_time;

//         // Accumulate total computation time
//         if (rank == 0) {
//             total_computation_time += iteration_time;
//             printf("Iteration %d: %f seconds\n", iter + 1, iteration_time);
//         }

//         // Free allocated memory
//         free(U);
//         free(Uprev);
//         free(mask);
//         free(laplacian);
//         free(xlin);
//     }

//     if (rank == 0) {
//         double average_computation_time = total_computation_time / num_iterations;
//         printf("Average total computation time: %f seconds\n", average_computation_time);
//     }

//     MPI_Finalize();
//     return 0;
// }

// void initialize_arrays(double *U, bool *mask, int N, int local_N, int rank, int size) {
//     for (int i = 0; i < local_N; i++) {
//         for (int j = 0; j < N; j++) {
//             int global_i = rank * local_N + i;
//             int idx = i * N + j;
//             mask[idx] = false;
//             U[idx] = 0.0;

//             if (global_i == 0 || global_i == N-1 || j == 0 || j == N-1) {
//                 mask[idx] = true;
//             }
//             if (global_i >= N/4 && global_i < N*9/32 && j < N-1) {
//                 mask[idx] = true;
//             }
//             if ((j >= N*5/16 && j < N*3/8) || (j >= N*5/8 && j < N*11/16)) {
//                 mask[idx] = false;
//             }
//         }
//     }
// }

// void apply_boundary_conditions(double *U, double *xlin, bool *mask, int N, int local_N, double t, int rank, int size) {
//     for (int i = 0; i < N; i++) {
//         U[i] = sin(20 * M_PI * t) * pow(sin(M_PI * xlin[i]), 2);
//     }
//     for (int i = 0; i < local_N * N; i++) {
//         if (mask[i]) {
//             U[i] = 0.0;
//         }
//     }
// }

// void update_arrays(double *U, double *Uprev, double *laplacian, int local_N, int N, double fac) {
//     for (int i = 0; i < local_N * N; i++) {
//         double Unew = 2 * U[i] - Uprev[i] + fac * laplacian[i];
//         Uprev[i] = U[i];
//         U[i] = Unew;
//     }
// }

// void compute_laplacian(double *U, double *laplacian, int local_N, int N, int rank, int size) {
//     double *U_global = (double *)malloc(N * N * sizeof(double));

//     if (U_global == NULL) {
//         fprintf(stderr, "Error: Unable to allocate memory for U_global\n");
//         MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
//     }

//     MPI_Allgather(U, local_N * N, MPI_DOUBLE, U_global, local_N * N, MPI_DOUBLE, MPI_COMM_WORLD);

//     for (int i = rank * local_N; i < (rank + 1) * local_N; i++) {
//         for (int j = 1; j < N - 1; j++) {
//             if (i > 0 && i < N - 1) {
//                 int idx = (i - rank * local_N) * N + j;
//                 laplacian[idx] = U_global[(i-1) * N + j] + U_global[(i+1) * N + j] + U_global[i * N + (j-1)] + U_global[i * N + (j+1)] - 4 * U_global[i * N + j];
//             }
//         }
//     }

//     free(U_global);
// }
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

// Function prototypes
/**
 * @brief Initialize the arrays U and mask.
 * 
 * @param U The array to initialize.
 * @param mask The mask array to initialize.
 * @param N The total number of grid points.
 * @param local_N The number of grid points per process.
 * @param rank The rank of the current process.
 * @param size The total number of processes.
 */
void initialize_arrays(double *U, bool *mask, int N, int local_N, int rank, int size);

/**
 * @brief Apply boundary conditions to the array U.
 * 
 * @param U The array to apply boundary conditions to.
 * @param xlin The linearized x-coordinate array.
 * @param mask The mask array.
 * @param N The total number of grid points.
 * @param local_N The number of grid points per process.
 * @param t The current simulation time.
 * @param rank The rank of the current process.
 * @param size The total number of processes.
 */
void apply_boundary_conditions(double *U, double *xlin, bool *mask, int N, int local_N, double t, int rank, int size);

/**
 * @brief Update the arrays U and Uprev using the computed laplacian.
 * 
 * @param U The current array of values.
 * @param Uprev The previous array of values.
 * @param laplacian The computed laplacian values.
 * @param local_N The number of grid points per process.
 * @param N The total number of grid points.
 * @param fac The update factor based on time step and grid spacing.
 */
void update_arrays(double *U, double *Uprev, double *laplacian, int local_N, int N, double fac);

/**
 * @brief Compute the laplacian of the array U.
 * 
 * @param U The current array of values.
 * @param laplacian The array to store computed laplacian values.
 * @param local_N The number of grid points per process.
 * @param N The total number of grid points.
 * @param rank The rank of the current process.
 * @param size The total number of processes.
 */
void compute_laplacian(double *U, double *laplacian, int local_N, int N, int rank, int size);

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Main function to run the wave simulation.
 * 
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return int Return status code.
 */
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Simulation parameters
    int N = 256; ///< Resolution
    double boxsize = 1.0; ///< Box size
    double c = 1.0; ///< Wave speed
    double tEnd = 2.0; ///< Stop time
    bool plotRealTime = false; ///< Switch for plotting simulation in real time
    int num_iterations = 6; ///< Number of times to run the simulation
    double total_computation_time = 0.0; ///< To accumulate total computation time over iterations

    // Mesh
    double dx = boxsize / N;
    double dt = (sqrt(2) / 2) * dx / c;
    double fac = dt * dt * c * c / (dx * dx);

    int local_N = N / size; ///< Local size of the grid

    for (int iter = 0; iter < num_iterations; iter++) {
        double *U = (double *)calloc(local_N * N, sizeof(double));
        double *Uprev = (double *)calloc(local_N * N, sizeof(double));
        bool *mask = (bool *)calloc(local_N * N, sizeof(bool));
        double *laplacian = (double *)calloc(local_N * N, sizeof(double));
        double *xlin = (double *)malloc(N * sizeof(double));

        if (U == NULL || Uprev == NULL || mask == NULL || laplacian == NULL || xlin == NULL) {
            fprintf(stderr, "Error: Unable to allocate memory\n");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        for (int i = 0; i < N; i++) {
            xlin[i] = 0.5 * dx + i * dx;
        }

        initialize_arrays(U, mask, N, local_N, rank, size);

        double t = 0.0; ///< Time for this iteration

        // Start the timer
        double start_time = MPI_Wtime();

        // Simulation Main Loop
        while (t < tEnd) {
            compute_laplacian(U, laplacian, local_N, N, rank, size);
            update_arrays(U, Uprev, laplacian, local_N, N, fac);
            apply_boundary_conditions(U, xlin, mask, N, local_N, t, rank, size);

            // Update time
            t += dt;

            MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes after each iteration
        }

        // Stop the timer
        double end_time = MPI_Wtime();
        double iteration_time = end_time - start_time;

        // Accumulate total computation time
        if (rank == 0) {
            total_computation_time += iteration_time;
            printf("Iteration %d: %f seconds\n", iter + 1, iteration_time);
        }

        // Free allocated memory
        free(U);
        free(Uprev);
        free(mask);
        free(laplacian);
        free(xlin);
    }

    if (rank == 0) {
        double average_computation_time = total_computation_time / num_iterations;
        printf("Average total computation time: %f seconds\n", average_computation_time);
    }

    MPI_Finalize();
    return 0;
}

void initialize_arrays(double *U, bool *mask, int N, int local_N, int rank, int size) {
    for (int i = 0; i < local_N; i++) {
        for (int j = 0; j < N; j++) {
            int global_i = rank * local_N + i;
            int idx = i * N + j;
            mask[idx] = false;
            U[idx] = 0.0;

            if (global_i == 0 || global_i == N-1 || j == 0 || j == N-1) {
                mask[idx] = true;
            }
            if (global_i >= N/4 && global_i < N*9/32 && j < N-1) {
                mask[idx] = true;
            }
            if ((j >= N*5/16 && j < N*3/8) || (j >= N*5/8 && j < N*11/16)) {
                mask[idx] = false;
            }
        }
    }
}

void apply_boundary_conditions(double *U, double *xlin, bool *mask, int N, int local_N, double t, int rank, int size) {
    for (int i = 0; i < N; i++) {
        U[i] = sin(20 * M_PI * t) * pow(sin(M_PI * xlin[i]), 2);
    }
    for (int i = 0; i < local_N * N; i++) {
        if (mask[i]) {
            U[i] = 0.0;
        }
    }
}

void update_arrays(double *U, double *Uprev, double *laplacian, int local_N, int N, double fac) {
    for (int i = 0; i < local_N * N; i++) {
        double Unew = 2 * U[i] - Uprev[i] + fac * laplacian[i];
        Uprev[i] = U[i];
        U[i] = Unew;
    }
}

void compute_laplacian(double *U, double *laplacian, int local_N, int N, int rank, int size) {
    double *U_global = (double *)malloc(N * N * sizeof(double));

    if (U_global == NULL) {
        fprintf(stderr, "Error: Unable to allocate memory for U_global\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    MPI_Allgather(U, local_N * N, MPI_DOUBLE, U_global, local_N * N, MPI_DOUBLE, MPI_COMM_WORLD);

    for (int i = rank * local_N; i < (rank + 1) * local_N; i++) {
        for (int j = 1; j < N - 1; j++) {
            if (i > 0 && i < N - 1) {
                int idx = (i - rank * local_N) * N + j;
                laplacian[idx] = U_global[(i-1) * N + j] + U_global[(i+1) * N + j] + U_global[i * N + (j-1)] + U_global[i * N + (j+1)] - 4 * U_global[i * N + j];
            }
        }
    }

    free(U_global);
}
