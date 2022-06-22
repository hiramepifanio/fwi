#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <fftw3.h>
#include "../lib/fwi.hpp"
#include "../lib/io.hpp"
#include "../lib/memory.hpp"

using namespace std;

float** velocityMap(int nz, int nx);

int main() {
    // Parameters
    string script = "singleTrace";
    int nz = 61, nx = 161, nt = 600;                // Number of grid points and timesteps
    float dx = 25, dz = 25, dt = 0.002;             // Step size in space and time (m, s)
    int border = 100;                               // Number of points in the border
    int nzb = nz + 2*border, nxb = nx + 2*border;   // Number of grid points with borders
    int souz = 5, soux = 50;                        // Source position
    int recstep = 1;                                // Number of points between receptors
    int nrecs = nx/recstep + 1;                     // Number of receptors
    int recz = 5;                                   // Receptors' Depth
    float fp = 12;                                  // Peak frequency
    float** shot;                                   // Keep record of the wave intensity at the receptors
    float* wav;                                     // Ricker wavelet
    float* wavp;                                    // Interpolated wavelet for propagation
    int** rec;                                      // Receptors grid position
    float** vel;                                    // Velocity map

    rec = receptors(recstep, nx, recz);             // Receptors position
    wav = ricker(fp, dt, nt);                       // Getting Ricker wavelet values
    vel = velocityMap(nz, nx);      

    int samplingShape[] = {nx};                      
    int* sampling = alloc1IntArr(0, samplingShape);
    sampling[110] = 1;                

    // Propagation
    shot = propagator(wav, dt, nt, souz, soux, rec, nrecs, vel, dz, dx, nz, nx, border);
    //saveMbin(shot, nt, nrecs, "data/gradient/shot.bin");
    int shape[] = {nt, nrecs};
    saveShot(shot, shape, "data/" + script, "shot");

    //Gradient
    float** v0 = valueM(6000, nz, nx);
    int reamos = 5;
    int ntp = nt/reamos;
    int nxa = nx + 2*border, nza = nz + 2*border;

    float** buffer = taper(nz, nx, border);
    float** wave0 = valueM(0, nza, nxa);
    float** wave1 = valueM(0, nza, nxa);
    float** wave2 = valueM(0, nza, nxa);
    float*** d2u = alloc3Arr(0, nz, nx, ntp);
    float** vel0 = extendModel(v0, nz, nx, border);
    float** diff = alloc2Arr(0, nt, nx);

    cout << "Calculating second derivative and difference at surface" << endl;
    // Step 1
    for(int t = 0; t < nt; t++) {
        if(t%(nt/100) == 0){
            cout << "Progress... " << 100*t/nt << "%\n";
        }
        
        // Wave propagation
        float **lap = laplacian(wave1, nzb, nxb, dz, dx);
        for(int i = 0; i < nzb; i++){
            for(int j = 0; j < nxb; j++){
                wave2[i][j] = 2*wave1[i][j] - wave0[i][j] + vel0[i][j]*vel0[i][j]*dt*dt*lap[i][j];
            }
        }
        wave2[souz + border][soux + border] += vel0[souz + border][soux + border]*vel0[souz + border][soux + border]*dt*dt*wav[t];

        if (t%reamos == 0)
            for(int i = 0; i < nz; i++)
                for(int j = 0; j < nx; j++)
                    d2u[i][j][t/reamos] = (wave2[i + border][j + border]-2*wave1[i + border][j + border]+wave0[i + border][j + border]);

        // Buffering
        buffer_wave(buffer, wave2, nzb, nxb);
        buffer_wave(buffer, wave1, nzb, nxb);

        // Roll arrays forward
        copy_M_to(wave1, wave0, nzb, nxb);
        copy_M_to(wave2, wave1, nzb, nxb);

        // Step 2
        for(int i = 0; i < nx; i++)
            diff[t][i] = wave2[border + recz][border + i] - shot[t][i];

    }
    cout << "... Completed." << endl; 

    // int shapeDiff[] = {nt, nx};
    // save2Arr(diff, shapeDiff, "data/" + script, "diff");

    // float** d2ut = alloc2Arr(0, nz, nx);
    // for(int i = 0; i < nz; i++)
    //     for(int j = 0; j < nx; j++)
    //         d2ut[i][j] += d2u[i][j][ntp-1];

    // int shapeD2ut[] = {nz, nx};
    // save2Arr(d2ut, shapeD2ut, "data/" + script, "d2ut"); 

    wave0 = valueM(0, nza, nxa);
    wave1 = valueM(0, nza, nxa);
    wave2 = valueM(0, nza, nxa);
    cout << "Going Backwards." << endl;
    for(int t = 0; t < nt; t++) {
        if(t%(nt/100) == 0){
            cout << "Progress... " << 100*t/nt << "%\n";
        }
        
        // Wave propagation
        float **lap = laplacian(wave1, nzb, nxb, dz, dx);
        for(int i = 0; i < nzb; i++){
            for(int j = 0; j < nxb; j++){
                wave2[i][j] = 2*wave1[i][j] - wave0[i][j] + vel0[i][j]*vel0[i][j]*dt*dt*lap[i][j];
            }
        }
        for(int i = 0; i < nx; i++) {
            wave2[border + recz][border + i] += sampling[i]*vel0[recz + border][i + border]*vel0[recz + border][i + border]*dt*dt*diff[nt - t - 1][i];
        }

        // Buffering
        buffer_wave(buffer, wave2, nzb, nxb);
        buffer_wave(buffer, wave1, nzb, nxb);

        // Roll arrays forward
        copy_M_to(wave1, wave0, nzb, nxb);
        copy_M_to(wave2, wave1, nzb, nxb);

        // Step 3
        if ((nt - 1 - t)%reamos == 0)
            for(int i = 0; i < nz; i++)
                for(int j = 0; j < nx; j++)
                    d2u[i][j][ntp-1-t/reamos] *= wave2[i + border][j + border];

    }
    cout << "... Completed." << endl;
    
    // Steps 4 and 5
    float** grad = alloc2Arr(0, nz, nx);
    for(int t = 0; t < ntp; t ++)
        for(int i = 0; i < nz; i++)
            for(int j = 0; j < nx; j++)
                grad[i][j] += d2u[i][j][t];

    // int shapeD2uSum[] = {nz, nx};
    // save2Arr(grad, shapeD2uSum, "data/" + script, "d2u_sum");  
    
    for(int i = 0; i < nz; i++)
        for(int j = 0; j < nx; j++)
            grad[i][j] *= -reamos*dt;//*2/(v0[i][j]*v0[i][j]*v0[i][j]);

    int shape2[] = {nz, nx};
    save2Arr(grad, shape2, "data/" + script, "grad");

    // Free memory
    freeV(wav);
    freeM(vel, nz);
    freeMint(rec, 2);
    freeM(shot, nt);

    return 0;
}

float** velocityMap(int nz, int nx) {
    float** velocityMap = allocateM(nz, nx);
    for(int i = 0; i < nz; i++) {
        for(int j = 0; j < nx; j++) {
            velocityMap[i][j] = 3000*(1 + 1.0*i/nz);
        }
    }

    return velocityMap;
}