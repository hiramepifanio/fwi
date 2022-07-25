#include "alloc.hpp"

float max(float** A, int* shape) {
    float max = A[0][0];
    for (int i = 0; i < shape[0]; i++)
        for (int j = 0; j < shape[1]; j++)
            if (max < A[i][j])
                max = A[i][j];
    
    return max;
}

float min(float** A, int* shape) {
    float min = A[0][0];
    for (int i = 0; i < shape[0]; i++)
        for (int j = 0; j < shape[1]; j++)
            if (min > A[i][j])
                min = A[i][j];
    
    return min;
}

float **laplacian(float **wave, int* wave_shape, float dz, float dx) {
    float **lap = alloc_float2d(0.0, wave_shape);
    int nz = wave_shape[0], nx = wave_shape[1];

    // Derivative in z
    // Center
    for(int i = 2; i < nz - 2; i++)
        for(int j = 0; j < nx; j++)
            lap[i][j] = (- wave[i + 2][j] + 16*wave[i + 1][j] - 30*wave[i][j] + 16*wave[i - 1][j] - wave[i - 2][j])/(12*dz*dz);

    // The 4 outer rows
    for(int j = 0; j < nx; j++){
        lap[0][j] = lap[1][j] = (wave[2][j] - 2*wave[1][j] + wave[0][j])/(dz*dz);
        lap[nz - 2][j] = lap[nz - 1][j] = (wave[nz - 1][j] - 2*wave[nz - 2][j] + wave[nz - 3][j])/(dz*dz);
    }

    // Derivative in x
    // Center
    for(int j = 2; j < nx - 2; j++)
        for(int i = 0; i < nz; i++)
            lap[i][j] += (- wave[i][j + 2] + 16*wave[i][j + 1] - 30*wave[i][j] + 16*wave[i][j - 1] - wave[i][j - 2])/(12*dx*dx);

    // The 4 outer colums
    for(int i = 0; i < nz; i++){
        lap[i][0] = lap[i][1] =  (wave[i][2] - 2*wave[i][1] + wave[i][0])/(dx*dx);
        lap[i][nx - 2] = lap[i][nx - 1] = (wave[i][nx - 1] - 2*wave[i][nx - 2] + wave[i][nx - 3])/(dx*dx);
    }

    return lap;
}

float** vel_dt_2(float** vel, int* vel_shape, float dt) {
    float** con = alloc_float2d(vel_shape);

    for(int i = 0; i < vel_shape[0]; i++)
        for(int j = 0; j < vel_shape[1]; j++)
            con[i][j] = vel[i][j]*vel[i][j]*dt*dt;

    return con;
}