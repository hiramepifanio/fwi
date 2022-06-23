#include <iostream>
#include <fstream>
#include "fwi.hpp"
#include "memory.hpp"
#include "math.hpp"
using namespace std;

float** grad(float* wavelet, float** shot, float** vel, float dt, float dz, float dx, int nt, int nz, int nx, int border, int souz, int soux, int* recz, int* recx, int nrecs) {
    int scale = cfl_criteria_scale(vel, nz, nx, dt, dz, dx);
    int ntp = nt*scale;
    float dtp = dt/scale;
    float* wavelet_p = sincint(wavelet, nt, dt, dtp);
    
    float** velb = extend_model(vel, nz, nx, border);
    float** buffer = taper(nz, nx, border);
    
    int nzb = nz + 2*border;
    int nxb = nx + 2*border;
    float** wave0 = alloc2Arr(0, nzb, nxb);
    float** wave1 = alloc2Arr(0, nzb, nxb);
    float** wave2 = alloc2Arr(0, nzb, nxb);
    float*** field = alloc3Arr(0, nz, nx, nt);
    float** diff = alloc2Arr(0, nt, nrecs);
    float** con = vel_dtp_2(velb, nzb, nxb, dtp);

    int t_sample = 0;
    cout << "Going forward." << endl;
    for(int t = 0; t < ntp; t++) {
        if (t%(ntp/100) == 0)
            cout << "Progress... " << 100*t/ntp << "%" << endl;
        
        // Wave propagation
        float **lap = laplacian(wave1, nzb, nxb, dz, dx);
        for (int i = 0; i < nzb; i++)
            for (int j = 0; j < nxb; j++)
                wave2[i][j] = 2*wave1[i][j] - wave0[i][j] + con[i][j]*lap[i][j];
        
        // At source
        wave2[souz + border][soux + border] += con[souz + border][soux + border]*wavelet_p[t]/(dx*dz);

        if (t%scale == 0) {
            for(int i = 0; i < nz; i++)
                for(int j = 0; j < nx; j++)
                    field[i][j][t_sample] = (wave2[i + border][j + border]-2*wave1[i + border][j + border]+wave0[i + border][j + border])/(dtp*dtp);

            for(int i = 0; i < nrecs; i++)
                diff[t][i] = wave2[border + recz[i]][border + recx[i]] - shot[t_sample][i];
            
            t_sample += 1;
        }
        // Buffering
        buffer_wave(buffer, wave2, nzb, nxb);
        buffer_wave(buffer, wave1, nzb, nxb);

        // Roll arrays forward
        copy_M_to(wave1, wave0, nzb, nxb);
        copy_M_to(wave2, wave1, nzb, nxb);
    }
    cout << "Finished." << endl;

    wave0 = alloc2Arr(0, nzb, nxb);
    wave1 = alloc2Arr(0, nzb, nxb);
    wave2 = alloc2Arr(0, nzb, nxb);
    diff = sincint_shot(diff, dt, dtp, nt, ntp, nrecs);

    t_sample = nt - 1;
    cout << "Going backward." << endl;
    for(int t = ntp - 1; t > -1; t--) {
        if ((ntp - 1 - t)%(ntp/100) == 0)
            cout << "Progress... " << 100*(ntp - 1 - t)/ntp << "%" << endl;
        
        // Wave propagation
        float **lap = laplacian(wave1, nzb, nxb, dz, dx);
        for (int i = 0; i < nzb; i++)
            for (int j = 0; j < nxb; j++)
                wave2[i][j] = 2*wave1[i][j] - wave0[i][j] + con[i][j]*lap[i][j];
        
        // At receptors
        for (int i = 0; i < nrecs; i++)
            wave2[border + recz[i]][border + recx[i]] += con[border + recz[i]][border + recx[i]]*diff[t][i];

        if (t%scale == 0) {
            for(int i = 0; i < nz; i++)
                for(int j = 0; j < nx; j++)
                    field[i][j][t_sample] *= wave1[border + i][border + j];
            
            t_sample -= 1;
        }
        // Buffering
        buffer_wave(buffer, wave2, nzb, nxb);
        buffer_wave(buffer, wave1, nzb, nxb);

        // Roll arrays forward
        copy_M_to(wave1, wave0, nzb, nxb);
        copy_M_to(wave2, wave1, nzb, nxb);
    }
    cout << "Finished." << endl;

    float** grad = alloc2Arr(0.0, nz, nx);
    for (int t = 0; t < nt; t++)
        for (int i = 0; i < nz; i++)
            for (int j = 0; j < nx; j++)
                grad[i][j] += field[i][j][t];

    return grad;
}
