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
    string script = "reconstruction";
    int nz = 201, nx = 251, nt = 300;               // Number of grid points and timesteps
    float dx = 20, dz = 20, dt = 0.01;              // Step size in space and time (m, s)
    int border = 50;                                // Number of points in the border
    int nzb = nz + 2*border, nxb = nx + 2*border;   // Number of grid points with borders
    int souz = 1, soux = nx/2;                      // Source position
    int recstep = 4;                                // Number of points between receptors
    int nrecs = nx/recstep + 1;                     // Number of receptors
    int recz = 2;                                   // Receptors' Depth
    float fp = 10;                                  // Peak frequency
    int sorder = 4/2;                               // Half laplacian's order
    int instant = 100;                              // Instant to visualize field
    float** shotd = alloc2Arr(nt, nrecs);                                  // Direct shot
    float** shotr = alloc2Arr(nt, nrecs);                                  // Reconstructed shot
    float* wav;                                     // Ricker wavelet
    //float* wavp;                                    // Interpolated wavelet for propagation
    int** rec;                                      // Receptors grid position
    float** vel;                                    // Velocity map

    rec = receptors(recstep, nx, recz);             // Receptors position
    wav = ricker(fp, dt, nt);                       // Getting Ricker wavelet values
    vel = velocityMap(nz, nx);                    

    // Propagation
    //int nzb = nz + 2*border, nxb = nx + 2*border; // Full grid with borders
    float **model = extendModel(vel, nz, nx, border); // Velocity model with borders
    float **buffer = taper(nz, nx, border); // Scale factor to buffer the wave at the borders
    float **wave0 = valueM(0.0, nzb, nxb); // Wave in the past
    float **wave1 = valueM(0.0, nzb, nxb); // Wave in the present
    float **wave2 = valueM(0.0, nzb, nxb); // Wave in the future
    float** fieldd = alloc2Arr(0, nz, nx);
    float** fieldr = alloc2Arr(0, nz, nx);
    
    float dtp = cfl_criteria(vel, nz, nx, dx, dz, dt);    // CFL with substitution of dt real to dtp (propagation)

    int ntp = nt*int(dt/dtp);
    float* wavp = sincint(wav, nt, dt, dtp);
    float maxwavp = 0;
    for (int i = 0; i < ntp; i++)
        if (wavp[i] > maxwavp)
            maxwavp = wavp[i];
    cout << "maxwavp " << maxwavp << endl;
    int reamos = int(dt/dtp);
    cout << "Interpolation: " << int(dt/dtp) << endl;

    float** con = alloc2Arr(nz + 2*border, nx + 2*border);
    for(int i = 0; i < nz + 2*border; i++)
        for(int j = 0; j < nx + 2*border; j++)
            con[i][j] = model[i][j]*model[i][j]*dtp*dtp;

    float*** surfacexl = alloc3Arr(ntp, nz, sorder);
    float*** surfacexr = alloc3Arr(ntp, nz, sorder);
    float*** surfacezl = alloc3Arr(ntp, sorder, nx);
    float*** surfacezr = alloc3Arr(ntp, sorder, nx);

    // Main loop
    cout << "Initiating propagator...\n";
    for(int k = 0; k < ntp; k++){
        if(k%(ntp/100) == 0){
            cout << "Progress... " << 100*k/ntp << "%\n";
        }
        
        float **lap = laplacian(wave1, nzb, nxb, dz, dx);

        // Wave propagation
        for(int i = 0; i < nzb; i++){
            for(int j = 0; j < nxb; j++){
                // wave2[i][j] = 2*wave1[i][j] - wave0[i][j] + model[i][j]*model[i][j]*dtp*dtp*lap[i][j];
                wave2[i][j] = 2*wave1[i][j] - wave0[i][j] + con[i][j]*lap[i][j];
            }
        }

        // Input source
        // wave2[souz + border][soux + border] += model[souz + border][soux + border]*model[souz + border][soux + border]*dtp*dtp*wavp[k]/(dx*dz);
        wave2[souz + border][soux + border] += con[souz + border][soux + border]*wavp[k]/(dx*dz);


        for (int i = 0; i < nz; i++)
            for (int j = 0; j < sorder; j++) {
                surfacexl[k][i][j] = wave1[border + i][border + j];
                surfacexr[k][i][j] = wave1[border + i][border + nx - j];
            }

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < sorder; j++) {
                surfacezl[k][j][i] = wave1[border + j][border + i];
                surfacezr[k][j][i] = wave1[border + nz - j][border + i];
            }


        // Collect wave intensity at receptors
        if (k%(ntp/nt) == 0) {
            int t = k/(ntp/nt);
            for(int r = 0; r < nrecs; r++){
                shotd[t][r] = wave1[border + rec[0][r]][border + rec[1][r]];
            }
        }
        
        if (k == instant*reamos)
            for (int i = 0; i < nz; i++)
                for (int j = 0; j < nx; j++)
                    fieldd[i][j] = wave1[border + i][border + j];

        

        freeM(lap, nzb);

        if (k != ntp - 1) {
            // Buffering
            buffer_wave(buffer, wave2, nzb, nxb);
            buffer_wave(buffer, wave1, nzb, nxb);

            // Roll arrays forward
            copy_M_to(wave1, wave0, nzb, nxb);
            copy_M_to(wave2, wave1, nzb, nxb);
        } else {
            copy_M_to(wave2, wave0, nzb, nxb);
        }
    }

    cout << "Propagation completed!" << endl;

    int ta = nt - 1;
    cout << "Initiating reconstruction...\n";
    for(int k = ntp - 1; k > -1; k--){
        if(k%(ntp/100) == 0){
            cout << "Progress... " << 100*(ntp - 1 - k)/ntp << "%\n";
        }
        for (int i = 0; i < nz; i++)
            for (int j = 0; j < sorder; j++) {
                wave1[border + i][border + j] = surfacexl[k][i][j];
                wave1[border + i][border + nx - j] = surfacexr[k][i][j];
            }

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < sorder; j++) {
                wave1[border + j][border + i] = surfacezl[k][j][i];
                wave1[border + nz - j][border + i] = surfacezr[k][j][i];
            }

        float **lap = laplacian(wave1, nzb, nxb, dz, dx);

        // Wave propagation
        for(int i = 0; i < nzb; i++){
            for(int j = 0; j < nxb; j++){
                wave2[i][j] = 2*wave1[i][j] - wave0[i][j] + con[i][j]*lap[i][j];
            }
        }

        // Input source
        wave2[souz + border][soux + border] += con[souz + border][soux + border]*wavp[k]/(dx*dz);

        // Collect wave intensity at receptors
        if (k%(ntp/nt) == 0) {
            //int t = k/(ntp/nt);
            int t = ta;
            for(int r = 0; r < nrecs; r++){
                shotr[t][r] = wave1[border + rec[0][r]][border + rec[1][r]];
            }
            ta -= 1;
        }
        if (k == instant*reamos)
            for (int i = 0; i < nz; i++)
                for (int j = 0; j < nx; j++)
                    fieldr[i][j] = wave1[border + i][border + j];

        

        freeM(lap, nzb);

        // Buffering
        buffer_wave(buffer, wave2, nzb, nxb);
        buffer_wave(buffer, wave1, nzb, nxb);

        // Roll arrays forward
        copy_M_to(wave1, wave0, nzb, nxb);
        copy_M_to(wave2, wave1, nzb, nxb);

    }

    cout << "Propagation completed!" << endl;

    int shape[] = {nz, nx};
    save2Arr(fieldd, shape,"data/" + script, "fieldd");
    save2Arr(fieldr, shape,"data/" + script, "fieldr");
    int shotShape[] = {nt, nrecs};
    saveShot(shotd, shotShape,"data/" + script, "shotd");
    saveShot(shotr, shotShape,"data/" + script, "shotr");

    

    return 0;
}

float** velocityMap(int nz, int nx) {
    float** velocityMap = allocateM(nz, nx);

    for(int i = 0; i < nz; i++) {
        for(int j = 0; j < nx; j++) {
            if (i >= nz/2) velocityMap[i][j] = 5500;
            else velocityMap[i][j] = 3000;
        }
    }

    return velocityMap;
}