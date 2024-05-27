#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>

void initialize_arrays(double *restrict U, bool *restrict mask, int N);
void apply_boundary_conditions(double *restrict U, double *restrict xlin, bool *restrict mask, int N, double t);
void update_and_compute(double *restrict U, double *restrict Uprev, double *restrict laplacian, double *restrict Unew, int N, double fac);
void save_to_file(double *restrict U, bool *restrict mask, int N, const char *filename);

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main() {
    int N = 256;
    double boxsize = 1.0;
    double c = 1.0;
    double t = 0.0;
    double tEnd = 2.0;
    bool plotRealTime = true;

    double dx = boxsize / N;
    double dt = (sqrt(2) / 2) * dx / c;
    double fac = dt * dt * c * c / (dx * dx);

    double *U = (double *)malloc(N * N * sizeof(double));
    double *Uprev = (double *)malloc(N * N * sizeof(double));
    bool *mask = (bool *)malloc(N * N * sizeof(bool));
    double *laplacian = (double *)malloc(N * N * sizeof(double));
    double *xlin = (double *)malloc(N * sizeof(double));
    double *Unew = (double *)malloc(N * N * sizeof(double));

    for (int i = 0; i < N; i++) {
        xlin[i] = 0.5 * dx + i * dx;
    }

    initialize_arrays(U, mask, N);
    double start_time = omp_get_wtime();
    while (t < tEnd) {
        update_and_compute(U, Uprev, laplacian, Unew, N, fac);
        apply_boundary_conditions(U, xlin, mask, N, t);

        t += dt;

        printf("Time: %f\n", t);
    }

    double end_time = omp_get_wtime(); // 记录结束时间
    double time_spent = end_time - start_time; // 计算运行时间
    printf("Total computation time: %f seconds\n", time_spent);
    // 输出最终结果
    save_to_file(U, mask, N, "output.txt");

    free(U);
    free(Uprev);
    free(mask);
    free(laplacian);
    free(Unew);
    free(xlin);

    return 0;
}

void initialize_arrays(double *restrict U, bool *restrict mask, int N) {
    #pragma omp parallel for collapse(2)
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

void apply_boundary_conditions(double *restrict U, double *restrict xlin, bool *restrict mask, int N, double t) {

    for (int i = 0; i < N; i++) {
        U[i] = sin(20 * M_PI * t) * pow(sin(M_PI * xlin[i]), 2);
    }
    #pragma omp parallel for
    for (int i = 0; i < N * N; i++) {
        if (mask[i]) {
            U[i] = 0.0;
        }
    }
}

void update_and_compute(double *restrict U, double *restrict Uprev, double *restrict laplacian, double *restrict Unew, int N, double fac) {
    #pragma omp parallel
    {
        // 计算拉普拉斯算子
        #pragma omp for collapse(2)
        for (int i = 1; i < N - 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                int idx = i * N + j;
                laplacian[idx] = U[(i-1) * N + j] + U[(i+1) * N + j] + U[i * N + (j-1)] + U[i * N + (j+1)] - 4 * U[idx];
            }
        }

        // 计算 Unew
        #pragma omp for
        for (int i = 0; i < N * N; i++) {
            Unew[i] = 2 * U[i] - Uprev[i] + fac * laplacian[i];
        }

        // 更新 U 和 Uprev
        #pragma omp for
        for (int i = 0; i < N * N; i++) {
            Uprev[i] = U[i];
            U[i] = Unew[i];
        }
    }
}

void save_to_file(double *restrict U, bool *restrict mask, int N, const char *filename) {
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
