#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <fftw3.h>
#include "../lib/fwi.hpp"
#include "../lib/io.hpp"
#include "../lib/memory.hpp"
#include "../lib/grad.hpp"

using namespace std;

int* receptorsz();
int* receptorsx(int nx);
float** velocityMap(int nz, int nx);
float** velocityMap0(int nz, int nx);

int main() {

    // Parameters
    string script = "grads";
    int nz = 201, nx = 251, nt = 300;               // Number of grid points and timesteps
    float dx = 20, dz = 20, dt = 0.01;              // Step size in space and time (m, s)
    int border = 50;                                // Number of points in the border
    int nzb = nz + 2*border, nxb = nx + 2*border;   // Number of grid points with borders
    int souz = 1, soux = nx/2;                      // Source position
    float fp = 10;                                  // Peak frequency
    int sorder = 4/2;                               // Half laplacian's order
    int instant = 100;                              // Instant to visualize field
    float* wavelet;                                     // Ricker wavelet
    int nrecs = 2;
    int* recz;                                      // Receptors grid position
    int* recx;                                      // Receptors grid position
    int** rec;                                      // Receptors grid position
    float** vel;                                    // Velocity map
    float** vel0;                                    // Velocity map

    recz = receptorsz();
    recx = receptorsx(nx);
    rec = new int *[2];
    rec[0] = recz;
    rec[1] = recx;

    wavelet = ricker(fp, dt, nt);                       // Getting Ricker wavelet values
    vel = velocityMap(nz, nx);    
    vel0 = velocityMap0(nz, nx); 

    float** shot = propagator(wavelet, dt, nt, souz, soux, rec, nrecs, vel, dz, dx, nz, nx, border);

    float** grad1 = grad(wavelet, shot, vel0, dt, dz, dx, nt, nz, nx, border, souz, soux, recz, recx, nrecs);
    float** grad2 = grad_recon(wavelet, shot, vel0, dt, dz, dx, nt, nz, nx, border, souz, soux, recz, recx, nrecs);

    int grad_shape[] = {nz, nx};
    save2Arr(grad1, grad_shape, "data/"+script, "grad");
    save2Arr(grad2, grad_shape, "data/"+script, "gradr");

    int shotShape[] = {nt, nrecs};
    saveShot(shot, shotShape,"data/" + script, "shot");

    return 0;
}

int* receptorsz() {
    int s[] = {2};
    int* r = alloc1IntArr(s);
    r[0] = 2;
    r[1] = 2;

    return r;
}

int* receptorsx(int nx) {
    int s[] = {2};
    int* r = alloc1IntArr(s);
    r[0] = nx/15;
    r[1] = 14*nx/15;

    return r;
}

float** velocityMap(int nz, int nx) {
    float** velocityMap = allocateM(nz, nx);

    for(int i = 0; i < nz; i++) 
        for(int j = 0; j < nx; j++) 
            velocityMap[i][j] = 3000;

    return velocityMap;
}

float** velocityMap0(int nz, int nx) {
    float** velocityMap = allocateM(nz, nx);

    for(int i = 0; i < nz; i++) {
        for(int j = 0; j < nx; j++) {
            if (i >= nz/2) velocityMap[i][j] = 4500;
            else velocityMap[i][j] = 3000;
        }
    }

    return velocityMap;
}