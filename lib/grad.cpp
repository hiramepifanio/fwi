#include <iostream>
#include <fstream>
#include "fwi.hpp"
#include "memory.hpp"
#include "io.hpp"
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
        // if (t%(ntp/100) == 0)
        //     cout << "Progress... " << 100*t/ntp << "%" << endl;
        
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
                diff[t_sample][i] = wave2[border + recz[i]][border + recx[i]] - shot[t_sample][i];
            
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

    // float** f = alloc2Arr(nz, nx);
    // for(int i = 0; i < nz; i++)
    //     for(int j = 0; j < nx; j++)
    //         f[i][j] = field[i][j][100];
    // string script = "grads";
    // int field_shape[] = {nz, nx};
    // save2Arr(f, field_shape, "data/"+script, "field");

    // string script = "grads";
    // int diff_shape[] = {nt, nrecs};
    // save2Arr(diff, diff_shape, "data/"+script, "diff");
    
    float cost = 0;
    for (int t = 0; t < nt; t++) 
        for (int r = 0; r < nrecs; r++)
            cost += 0.5*diff[t][r]*diff[t][r];
            
    diff = sincint_shot(diff, dt, dtp, nt, ntp, nrecs);

    // string script = "grads";
    // int diff_shape[] = {ntp, nrecs};
    // save2Arr(diff, diff_shape, "data/"+script, "diffsinc");

    wave0 = alloc2Arr(0, nzb, nxb);
    wave1 = alloc2Arr(0, nzb, nxb);
    wave2 = alloc2Arr(0, nzb, nxb);

    t_sample = nt - 1;
    cout << "Going backward." << endl;
    for(int t = ntp - 1; t > -1; t--) {
        // if ((ntp - 1 - t)%(ntp/100) == 0)
        //     cout << "Progress... " << 100*(ntp - 1 - t)/ntp << "%" << endl;
        
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

    float** f = alloc2Arr(nz, nx);
    for(int i = 0; i < nz; i++)
        for(int j = 0; j < nx; j++)
            f[i][j] = field[i][j][200];
    string script = "grads";
    int field_shape[] = {nz, nx};
    save2Arr(f, field_shape, "data/"+script, "field");

    float** grad = alloc2Arr(0.0, nz, nx);
    for (int t = 0; t < nt; t++)
        for (int i = 0; i < nz; i++)
            for (int j = 0; j < nx; j++)
                grad[i][j] -= field[i][j][t];

    cout << "Cost: " << cost << endl;

    return grad;
}

float** grad_recon(float* wavelet, float** shot, float** vel, float dt, float dz, float dx, int nt, int nz, int nx, int border, int souz, int soux, int* recz, int* recx, int nrecs) {
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

    int sorder = 2;
    float*** surfacexl = alloc3Arr(ntp, nz, sorder);
    float*** surfacexr = alloc3Arr(ntp, nz, sorder);
    float*** surfacezl = alloc3Arr(ntp, sorder, nx);
    float*** surfacezr = alloc3Arr(ntp, sorder, nx);

    int t_sample = 0;
    cout << "Going forward." << endl;
    for(int t = 0; t < ntp; t++) {
        // if (t%(ntp/100) == 0)
        //     cout << "Progress... " << 100*t/ntp << "%" << endl;
        
        // Wave propagation
        float **lap = laplacian(wave1, nzb, nxb, dz, dx);
        for (int i = 0; i < nzb; i++)
            for (int j = 0; j < nxb; j++)
                wave2[i][j] = 2*wave1[i][j] - wave0[i][j] + con[i][j]*lap[i][j];
        
        // At source
        wave2[souz + border][soux + border] += con[souz + border][soux + border]*wavelet_p[t]/(dx*dz);

        // At surfaces
        for (int i = 0; i < nz; i++)
            for (int j = 0; j < sorder; j++) {
                surfacexl[t][i][j] = wave1[border + i][border + j];
                surfacexr[t][i][j] = wave1[border + i][border + nx - j];
            }

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < sorder; j++) {
                surfacezl[t][j][i] = wave1[border + j][border + i];
                surfacezr[t][j][i] = wave1[border + nz - j][border + i];
            }

        if (t%scale == 0) {
            for(int i = 0; i < nz; i++)
                for(int j = 0; j < nx; j++)
                    field[i][j][t_sample] = (wave2[i + border][j + border]-2*wave1[i + border][j + border]+wave0[i + border][j + border])/(dtp*dtp);

            for(int i = 0; i < nrecs; i++)
                diff[t_sample][i] = wave2[border + recz[i]][border + recx[i]] - shot[t_sample][i];
            
            t_sample += 1;
        }
        if (t != ntp - 1) {
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
    cout << "Finished." << endl;
    
    float cost = 0;
    for (int t = 0; t < nt; t++) 
        for (int r = 0; r < nrecs; r++)
            cost += 0.5*diff[t][r]*diff[t][r];

    diff = sincint_shot(diff, dt, dtp, nt, ntp, nrecs);
    wave0 = alloc2Arr(0, nzb, nxb);
    wave1 = alloc2Arr(0, nzb, nxb);
    wave2 = alloc2Arr(0, nzb, nxb);
    float** waver0 = alloc2Arr(0, nzb, nxb);
    float** waver1 = alloc2Arr(0, nzb, nxb);
    float** waver2 = alloc2Arr(0, nzb, nxb);
    float** grad = alloc2Arr(0.0, nz, nx);

    t_sample = nt - 1;
    cout << "Going backward." << endl;
    for(int t = ntp - 1; t > -1; t--) {
        // if ((ntp - 1 - t)%(ntp/100) == 0)
        //     cout << "Progress... " << 100*(ntp - 1 - t)/ntp << "%" << endl;

        for (int i = 0; i < nz; i++)
            for (int j = 0; j < sorder; j++) {
                wave1[border + i][border + j] = surfacexl[t][i][j];
                wave1[border + i][border + nx - j] = surfacexr[t][i][j];
            }

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < sorder; j++) {
                wave1[border + j][border + i] = surfacezl[t][j][i];
                wave1[border + nz - j][border + i] = surfacezr[t][j][i];
            }
        
        // Wave propagation
        float **lap = laplacian(wave1, nzb, nxb, dz, dx);
        for (int i = 0; i < nzb; i++)
            for (int j = 0; j < nxb; j++)
                wave2[i][j] = 2*wave1[i][j] - wave0[i][j] + con[i][j]*lap[i][j];

        wave2[souz + border][soux + border] += con[souz + border][soux + border]*wavelet_p[t]/(dx*dz);

        float **lapr = laplacian(waver1, nzb, nxb, dz, dx);
        for (int i = 0; i < nzb; i++)
            for (int j = 0; j < nxb; j++)
                waver2[i][j] = 2*waver1[i][j] - waver0[i][j] + con[i][j]*lapr[i][j];

        // At receptors
        for (int i = 0; i < nrecs; i++)
            waver2[border + recz[i]][border + recx[i]] += con[border + recz[i]][border + recx[i]]*diff[t][i];

        if (t%scale == 0) {
            
            for(int i = 0; i < nz; i++)
                for(int j = 0; j < nx; j++)
                    grad[i][j] -= waver1[i + border][j + border]*(wave2[i + border][j + border]-2*wave1[i + border][j + border]+wave0[i + border][j + border])/(dtp*dtp);

        }

        // Buffering
        buffer_wave(buffer, wave2, nzb, nxb);
        buffer_wave(buffer, wave1, nzb, nxb);
        buffer_wave(buffer, waver2, nzb, nxb);
        buffer_wave(buffer, waver1, nzb, nxb);

        // Roll arrays forward
        copy_M_to(wave1, wave0, nzb, nxb);
        copy_M_to(wave2, wave1, nzb, nxb);
        copy_M_to(waver1, waver0, nzb, nxb);
        copy_M_to(waver2, waver1, nzb, nxb);
    }
    cout << "Finished." << endl;

    cout << "Cost: " << cost << endl;

    return grad;
}
