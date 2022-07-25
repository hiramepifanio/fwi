#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <fftw3.h>
#include "../lib/fwi.hpp"
#include "../lib/io.hpp"
#include "../lib/memory.hpp"
#include "../lib/grad.hpp"
#include "../lib/ricker.hpp"
#include "../lib/savedata.hpp"
#include "../lib/alloc.hpp"

using namespace std;

int* receptorsz(int nrecs);
int* receptorsx(int nrecs);
int** receptors(int nrecs, int* recz, int* recx);
float** velocityMap(int nz, int nx);
float** velocityMap0(int nz, int nx);

int main() {

    // Parameters
    string script = "grads-overthurst";
    int nz = 187, nx = 801, nt = 1500;               // Number of grid points and timesteps
    float dx = 25, dz = 25, dt = 0.002;              // Step size in space and time (m, s)
    int border = 50;                                // Number of points in the border
    int nzb = nz + 2*border, nxb = nx + 2*border;   // Number of grid points with borders
    int souz = 0, soux = nx/2;                        // Source position
    float fp = 8;                                  // Peak frequency
    int sorder = 2;                               // Half laplacian's order
    float* wavelet;                                     // Ricker wavelet
    int nrecs = nx;
    int* recz;                                      // Receptors grid position
    int* recx;                                      // Receptors grid position
    int** recs;                                      // Receptors grid position
    float** vel;                                    // Velocity map
    float** vel0;                                    // Velocity map

    recz = receptorsz(nrecs);
    recx = receptorsx(nrecs);
    recs = receptors(nrecs, recz, recx);

    int source_shape[] = {2, 1};
    int** source = alloc_int2d(source_shape);
    source[0][0] = souz, source[1][0] = soux;
    save_data_int2d(source, source_shape, "/home/hiram/vscode-workspace/fwi/"+script+"/data", "source");

    int recs_shape[] = {2, nrecs};
    save_data_int2d(recs, recs_shape, "/home/hiram/vscode-workspace/fwi/"+script+"/data", "recs");

    wavelet = ricker(fp, dt, nt);                       // Getting Ricker wavelet values
    vel = velocityMap(nz, nx);    
    vel0 = velocityMap0(nz, nx); 

    int vel_shape[] = {nz, nx};
    save_data_float2d(vel, vel_shape, "/home/hiram/vscode-workspace/fwi/"+script+"/data", "vel");
    save_data_float2d(vel0, vel_shape, "/home/hiram/vscode-workspace/fwi/"+script+"/data", "vel0");

    float** shot = propagator(wavelet, dt, nt, souz, soux, recs, nrecs, vel, dz, dx, nz, nx, border);

    int shot_shape[] = {nt, nrecs};
    save_data_float2d(shot, shot_shape, "/home/hiram/vscode-workspace/fwi/"+script+"/data", "shot");

    //float** grad1 = grad(wavelet, shot, vel0, dt, dz, dx, nt, nz, nx, border, souz, soux, recz, recx, nrecs);
    float** grad2 = grad_recon(wavelet, shot, vel0, dt, dz, dx, nt, nz, nx, border, souz, soux, recz, recx, nrecs);

    int grad_shape[] = {nz, nx};
    //save_data_float2d(grad1, grad_shape, "/home/hiram/vscode-workspace/fwi/"+script+"/data", "grad");
    save_data_float2d(grad2, grad_shape, "/home/hiram/vscode-workspace/fwi/"+script+"/data", "gradr");

    return 0;
}

int* receptorsz(int nrecs) {
    int* r = alloc_int1d(5, nrecs);

    return r;
}

int* receptorsx(int nrecs) {
    int* r = alloc_int1d(nrecs);

    for (int i = 0; i < nrecs; i++)
        r[i] = i;

    return r;
}

int** receptors(int nrecs, int* recz, int* recx) {
    int r_shape[] = {2, nrecs};
    int** r = alloc_int2d(r_shape);

    for (int i = 0; i < nrecs; i++) {
        r[0][i] = recz[i];
        r[1][i] = recx[i];
    }

    return r;
}

float** velocityMap(int nz, int nx) {


    return loadMbin("data/overthrust.bin", nz, nx);;
}

float** velocityMap0(int nz, int nx) {
    int vel_shape[] = {nz, nx};
    float** velocityMap = alloc_float2d(6000, vel_shape);

    return velocityMap;
}