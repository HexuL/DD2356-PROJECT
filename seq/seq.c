/**
 * @file main.c
 * @brief Main simulation for wave propagation using finite differences.
 * 
 * This program simulates wave propagation in a medium with given boundary conditions.
 * It saves intermediate results to files for visualization with outputdata_show.py.
 * Note: Saving data significantly affects performance, disable it for timing measurements.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

// Function prototypes
void initialize_arrays(double *U, bool *mask, int N);
void apply_boundary_conditions(double *U, double *xlin, bool *mask, int N, double t);
void update_arrays(double *U, double *Uprev, double *laplacian, int N, double fac);
void compute_laplacian(double *U, double *laplacian, int N);
void save_to_file(double *U, bool *mask, int N, const char *filename);

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * Main execution function for the simulation.
 * Initializes arrays, runs the simulation loop, and saves output for visualization.
 */
int main() {
    clock_t start, end;
    double cpu_time_used;
    
    // Simulation parameters
    int N = 256;          ///< Grid resolution
    double boxsize = 1.0; ///< Size of the simulation box
    double c = 1.0;       ///< Wave speed
    double t = 0.0;       ///< Initial time
    double tEnd = 2.0;    ///< End time for the simulation
    bool plotRealTime = true; ///< Flag to control real-time plotting

    // Mesh parameters
    double dx = boxsize / N;
    double dt = (sqrt(2) / 2) * dx / c;
    double fac = dt * dt * c * c / (dx * dx);

    double *U = (double *)calloc(N * N, sizeof(double));
    double *Uprev = (double *)calloc(N * N, sizeof(double));
    bool *mask = (bool *)calloc(N * N, sizeof(bool));
    double *laplacian = (double *)calloc(N * N, sizeof(double));
    double *xlin = (double *)malloc(N * sizeof(double));

    for (int i = 0; i < N; i++) {
        xlin[i] = 0.5 * dx + i * dx;
    }

    initialize_arrays(U, mask, N);

    int outputCount = 1;
    start = clock();
    // Simulation Main Loop
    while (t < tEnd) {
        compute_laplacian(U, laplacian, N);
        update_arrays(U, Uprev, laplacian, N, fac);
        apply_boundary_conditions(U, xlin, mask, N, t);

        // update time
        t += dt;

        printf("Time: %f\n", t);

//         // Create directory if it doesn't exist
//         const char *dirName = "outputdata";
// #ifdef _WIN32
//         mkdir(dirName);
// #else
//         mkdir(dirName, 0777);
// #endif

        // Save state for plotting in real time
        if (plotRealTime || t >= tEnd) {
            char filename[50];
            sprintf(filename, "%s/output_%d.txt", dirName, outputCount);
            save_to_file(U, mask, N, filename);
            outputCount++;
        }

    }
    end = clock();
    // Final save
    save_to_file(U, mask, N, "finitedifference.txt");

    // Free allocated memory
    free(U);
    free(Uprev);
    free(mask);
    free(laplacian);
    free(xlin);

    
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("Time taken: %f seconds\n", cpu_time_used);

    return 0;
}

/**
 * @brief Initializes the simulation arrays and boundary mask.
 * @param U Pointer to the wave field array.
 * @param mask Pointer to the boundary mask array.
 * @param N Size of the grid.
 */ 
void initialize_arrays(double *U, bool *mask, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int idx = i * N + j;
            mask[idx] = false;
            U[idx] = 0.0;

            if (i == 0 || i == N-1 || j == 0 || j == N-1) {
                mask[idx] = true;
            }
            if (i >= N/4 && i < N*9/32 && j < N-1) {
                mask[idx] = true;
            }
            if ((j >= N*5/16 && j < N*3/8) || (j >= N*5/8 && j < N*11/16)) {
                mask[idx] = false;
            }
        }
    }
}

/**
 * @brief Applies boundary conditions to the wave field.
 * @param U Pointer to the wave field array.
 * @param xlin Pointer to the array of x linear space coordinates.
 * @param mask Pointer to the boundary mask array.
 * @param N Size of the grid.
 * @param t Current simulation time.
 */
void apply_boundary_conditions(double *U, double *xlin, bool *mask, int N, double t) {
    for (int i = 0; i < N; i++) {
        U[i] = sin(20 * M_PI * t) * pow(sin(M_PI * xlin[i]), 2);
    }
    for (int i = 0; i < N * N; i++) {
        if (mask[i]) {
            U[i] = 0.0;
        }
    }
}


/**
 * @brief Updates the wave field array for the next time step.
 * @param U Pointer to the current wave field array.
 * @param Uprev Pointer to the previous wave field array.
 * @param laplacian Pointer to the array containing the Laplacian of U.
 * @param N Size of the grid.
 * @param fac Precomputed factor used in the update formula.
 */
void update_arrays(double *U, double *Uprev, double *laplacian, int N, double fac) {
    double *Unew = (double *)malloc(N * N * sizeof(double));
    for (int i = 0; i < N * N; i++) {
        Unew[i] = 2 * U[i] - Uprev[i] + fac * laplacian[i];
        Uprev[i] = U[i];
        U[i] = Unew[i];
    }
    free(Unew);
}

/**
 * @brief Computes the Laplacian of the wave field using finite difference.
 * @param U Pointer to the wave field array.
 * @param laplacian Pointer to the array where the Laplacian will be stored.
 * @param N Size of the grid.
 */
void compute_laplacian(double *U, double *laplacian, int N) {
    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            int idx = i * N + j;
            laplacian[idx] = U[(i-1) * N + j] + U[(i+1) * N + j] + U[i * N + (j-1)] + U[i * N + (j+1)] - 4 * U[idx];
        }
    }
}


/**
 * @brief Saves the current state of the wave field to a file.
 * @param U Pointer to the wave field array.
 * @param mask Pointer to the boundary mask array.
 * @param N Size of the grid.
 * @param filename Name of the file to save the data.
 */
void save_to_file(double *U, bool *mask, int N, const char *filename) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Error opening file!\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int idx = i * N + j;
            if (mask[idx]) {
                fprintf(f, "nan ");
            } else {
                fprintf(f, "%f ", U[idx]);
            }
        }
        fprintf(f, "\n");
    }
    fclose(f);
}