#include "memory.hpp"

float max(float** A, int* shape) {
    float max = A[0][0];
    for (int i = 0; i < shape[0]; i++)
        for (int j = 0; j < shape[1]; j++)
            if (max < A[i][j])
                max = A[i][j];
    
    return max;
}

float** vel_dtp_2(float** vel, int nz, int nx, float dtp) {
    float** con = alloc2Arr(nz, nx);
    for(int i = 0; i < nz; i++)
        for(int j = 0; j < nx; j++)
            con[i][j] = vel[i][j]*vel[i][j]*dtp*dtp;

    return con;
}