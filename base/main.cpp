#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include "fwi.hpp"
#include "io.hpp"
#include "memory.hpp"
#include <fftw3.h>
//#include "fftw-m/fftw3.h"

using namespace std;

int main() {
    // Parameters
    int nz = 187, nx = 801, nt = 400, ntp; // Number of grid points and timesteps
    float dx = 25, dz = 25, dt = 0.01, dtp; // Step size in space and time (m, s)
    int border = 50; // Number of points in the border
    int nzb = nz + 2*border, nxb = nx + 2*border; // Number of grid points with borders
    int souz = 0, soux = nx/2; // Source position
    int recstep = 4; // Number of points between receptors
    int nrecs = nx/recstep + 1; // Number of receptors
    float fp = 6.0; // Peak frequency
    float **shot; // Keep record of the wave intensity at the receptors
    float *wav, *wavp; // Ricker wavelet
    int **rec; // Receptors grid position
    float **vel; // Velocity map

    rec = receptors(recstep, nz, nx); // Receptors position
    wav = ricker(fp, dt, nt); // Getting Ricker wavelet values
    vel = loadMbin("data/overthrust.bin", nz, nx); // Loading velocity map from binary file

    // Criteria
    dtp = cfl_criteria(vel, nz, nx, dx, dz, dt); // CFL with substitution of dt real to dtp (propagation)
    minimum_wavelenth(vel, nz, nx, fp, dx); // Wavelength

    // New parameters
    ntp *= int(dt/dtp);
    wavp = sincint(wav, nt, dt, dtp);

    
    // Propagation
    // shot = propagator(wav, dt, nt, souz, soux, rec, nrecs, vel, dz, dx, nz, nx, border);
    // saveMbin(shot, nt, nrecs, "data/shot.bin");

    // Free memory
    freeV(wav);
    freeM(vel, nz);
    freeMint(rec, 2);
    //freeM(shot, nt);

    return 0;
}