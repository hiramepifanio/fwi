#include <iostream>
#include <fftw3.h>
#include "memory.hpp"
#include "io.hpp"
#include <math.h>
#include "math.hpp"
#include "alloc.hpp"

using namespace std;

// // Receptors position
// int **receptors(int step, int nx, int z){
//     int nr = nx/step + 1; // Number of receptors
//     int **rec = allocateMint(2, nr); // Receptors
//     // Filling the positions
//     for(int r = 0; r < nr; r++){
//         rec[0][r] = z;
//         rec[1][r] = r*step;
//     }

//     return rec;
// }

// int **receptors(int step, int nx){
//     return receptors(step, nx, 2);
// }

// CFL criteria
int cfl_criteria_scale(float **vel, int nz, int nx, float dt, float dz, float dx) {
    float dtp = dt;
    int scale = 1;

    // Searching for maximum velocity
    int vel_shape[] = {nz, nx};
    float vel_max = max(vel, vel_shape);

    // Finding dt propagation
    dtp = dt/scale; 
    float CFL = vel_max*dtp*pow(1/(dx*dx) + 1/(dz*dz), 0.5), target = pow(3./4, 0.5);
    while(CFL > target) {
        scale++;
        dtp = dt/scale;
        CFL = vel_max*dtp*pow(1/(dx*dx) + 1/(dz*dz), 0.5);
    }
    
    cout << "The CFL criteria is " << CFL << " and should be equal to or less than " << target << ": ";
    if(CFL <= target)
        cout << "True." << endl;
    else
        cout << "False." << endl;

    cout << "The ratio dt/dtp is " << scale << " with dt = " << dt*1000 << " ms and dtp = " << dtp*1000 << " ms." << endl;

    return scale;
}

float cfl_criteria(float **vel, int n, int m, float dz, float dx, float dt) {
    float vmax = vel[0][0];
    float dtp = dt;
    int reamos = 1;

    // Searching for maximum velocity
    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            if(vel[i][j] > vmax)
                vmax = vel[i][j];

    // Finding dt propagation
    dtp = dt/reamos; 
    float CFL = vmax*dtp*pow(1/(dx*dx) + 1/(dz*dz), 0.5), target = pow(3./4, 0.5);
    while(CFL > target) {
        reamos++;
        dtp = dt/reamos;
        CFL = vmax*dtp*pow(1/(dx*dx) + 1/(dz*dz), 0.5);
    }
    
    cout << "The CFL criteria is " << CFL << " and should be equal to or less than " << target << ": ";
    if(CFL <= target)
        cout << "True." << endl;
    else
        cout << "False." << endl;

    cout << "The ratio dt/dtp is " << reamos << " with dt = " << dt*1000 << " ms and dtp = " << dtp*1000 << " ms." << endl;

    return dtp;
}

// // Minimum wavelength criteria
// void minimum_wavelenth(float **vel, int n, int m, float fp, float dx){
//     float vmin = vel[0][0];
//     for(int i = 0; i < n; i++)
//         for(int j = 0; j < m; j++)
//             if(vel[i][j] < vmin)
//                 vmin = vel[i][j];
//     float minlambda = vmin/(3*fp);
//     cout << "The minimum wavelength is " << minlambda << " and should be equal to or bigger than " << 5*dx << ": ";
//     if(minlambda >= 5*dx)
//         cout << "True." << endl;
//     else
//         cout << "False." << endl;
// }

float **extend_model(float **vel, int nz, int nx, int border) {
    int model_shape[] = {nz + 2*border, nx + 2*border};
    float **model = alloc_float2d(0.0, model_shape);

    // Center
    for(int i = 0; i < nz; i++){
        for(int j = 0; j < nx; j++){
            model[border + i][border + j] = vel[i][j];
        }
    }

    // Sides of the borders
    for(int i = 0; i < border; i++){
        for(int j = 0; j < nz; j++){
            model[border + j][i] = vel[j][0];
            model[border + j][nx + 2*border - 1 - i] = vel[j][nx - 1];
        }
    }

    // Top and bottom
    for(int i = 0; i < nx + 2*border; i++){
        for(int j = 0; j < border; j++){
            model[j][i] = model[border][i];
            model[nz + 2*border - 1 - j][i] = model[nz + border - 1][i];
        }
    }

    return model;
}

float **taper(int nz, int nx, int border){
    int buffer_shape[] = {nz + 2*border, nx + 2*border};
    float **buffer = alloc_float2d(1.0, buffer_shape);
    float damp = 6.5*border; // Buffering parameter

    // Decrease the border
    for(int i = 0; i < border; i++){
        for(int j = 0; j < nx + 2*border; j++){
            buffer[i][j] = buffer[nz + 2*border - 1 - i][j] = 0.5*(1 + cos(M_PI*(border - i)/damp));
        }
    }
    for(int j = 0; j < border; j++){
        for(int i = 0; i < nz + 2*border; i++){
            buffer[i][j] = buffer[i][nx + 2*border - 1 - j] *= 0.5*(1 + cos(M_PI*(border - j)/damp));
        }
    }

    return buffer;
}

void buffer_wave(float **buffer, float **wave, int nz, int nx){
    for(int i = 0; i < nz; i++){
        for(int j = 0; j < nx; j++){
            wave[i][j] *= buffer[i][j];
        }
    }
}

float* sincint(float* y, int Nin, float dt, float dtp) {
    int scale = int(dt/dtp);

    if(scale == 1) 
        return y;

    // FFTW lib works with double values
    double* yDouble = new double [Nin];
    for(int i = 0; i < Nin; i++) {
        yDouble[i] = double(y[i]);
    }

    double yDoubleMean = 0;
    for (int i = 0; i < Nin; i++) {
        double val = yDouble[i];
        if (val >= 0) 
            yDoubleMean += val;
        else
            yDoubleMean += -val;
    }
    yDoubleMean /= Nin;

    // Real data DFT
    fftw_complex* yDFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nin/2 + 1));
    fftw_plan planR2C = fftw_plan_dft_r2c_1d(Nin, yDouble, yDFT, FFTW_ESTIMATE);
    fftw_execute(planR2C);

    // Build the output array DFT
    int Nout = scale*Nin;
    fftw_complex* yOutDFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nout/2 + 1));
    for(int i = 0; i < Nout/2 + 1; i++) {
        if (i < Nin/2 + 1) {
            yOutDFT[i][0] = yDFT[i][0];
            yOutDFT[i][1] = yDFT[i][1];
        }
        else {
            yOutDFT[i][0] = 0.0;
            yOutDFT[i][1] = 0.0;
        }
    }

    // Hermitian data DFT
    double* ypDouble = new double [Nout];
    fftw_plan planC2R = fftw_plan_dft_c2r_1d(Nout, yOutDFT, ypDouble, FFTW_ESTIMATE);
    fftw_execute(planC2R);

    double ypDoubleMean = 0;
    for (int i = 0; i < Nout; i++) {
        double val = ypDouble[i];
        if (val >= 0) 
            ypDoubleMean += val;
        else
            ypDoubleMean += -val;
    }
    ypDoubleMean /= Nout;

    // Transform the output from double to float array
    float* yp = new float [Nout];
    for(int i = 0; i < Nout; i++) {
        yp[i] = (float) (yDoubleMean/ypDoubleMean)*ypDouble[i];
    }

    fftw_destroy_plan(planR2C);
    fftw_destroy_plan(planC2R);

    return yp;
}

float** sincint_shot(float** shot, float dt, float dtp, int nt, int ntp, int nrecs) {
    float** output = alloc2Arr(0.0, ntp, nrecs);
    
    for (int i = 0; i < nrecs; i++) {
        float* shot_rec = alloc1Arr(nt);
        for (int t = 0; t < nt; t++)
            shot_rec[t] = shot[t][i];

        float* sincint_shot_rec = sincint(shot_rec, nt, dt, dtp);
        for (int tp = 0; tp < ntp; tp++) {
            output[tp][i] = sincint_shot_rec[tp];
        }
    }

    return output;
}

float **propagator(float *wav, float dt, int nt, int souz, int soux, int **rec, int nrecs, float **vel, float dz, float dx, int nz, int nx, int border){
    int nzb = nz + 2*border, nxb = nx + 2*border; // Full grid with borders
    float **model = extend_model(vel, nz, nx, border); // Velocity model with borders
    float **buffer = taper(nz, nx, border); // Scale factor to buffer the wave at the borders
    float **wave0 = valueM(0.0, nzb, nxb); // Wave in the past
    float **wave1 = valueM(0.0, nzb, nxb); // Wave in the present
    float **wave2 = valueM(0.0, nzb, nxb); // Wave in the future
    
    float dtp = cfl_criteria(vel, nz, nx, dx, dz, dt);    // CFL with substitution of dt real to dtp (propagation)

    int ntp = nt*int(dt/dtp);
    float* wavp = sincint(wav, nt, dt, dtp);

    cout << "Interpolation: " << int(dt/dtp) << endl;

    // Main loop
    float **shot = valueM(0.0, nt, nrecs); // Information about the intensity of the wave at the receptors
    cout << "Initiating propagator...\n";
    for(int k = 0; k < ntp; k++){
        // if(k%(ntp/100) == 0){
        //     cout << "Progress... " << 100*k/ntp << "%\n";
        // }
        int lap_shape[] = {nzb, nxb};
        float **lap = laplacian(wave1, lap_shape, dz, dx);

        // Wave propagation
        for(int i = 0; i < nzb; i++){
            for(int j = 0; j < nxb; j++){
                wave2[i][j] = 2*wave1[i][j] - wave0[i][j] + model[i][j]*model[i][j]*dtp*dtp*lap[i][j];
            }
        }

        // Input source
        wave2[souz + border][soux + border] += model[souz + border][soux + border]*model[souz + border][soux + border]*dtp*dtp*wavp[k];// /(dx*dz);

        // Collect wave intensity at receptors
        if (k%(ntp/nt) == 0) {
            int t = k/(ntp/nt);
            for(int r = 0; r < nrecs; r++){
                shot[t][r] = wave2[border + rec[0][r]][border + rec[1][r]];
            }
        }

        

        freeM(lap, nzb);

        // Buffering
        buffer_wave(buffer, wave2, nzb, nxb);
        buffer_wave(buffer, wave1, nzb, nxb);

        // Roll arrays forward
        int wave_shape[] = {nzb, nxb};
        transfer_float2d(wave1, wave0, wave_shape);
        transfer_float2d(wave2, wave1, wave_shape);
    }

    cout << "Propagation completed!" << endl;

    // Free memory
    freeM(wave0, nzb);
    freeM(wave1, nzb);
    freeM(wave2, nzb);
    freeM(model, nzb);
    freeM(buffer, nzb);

    return shot;
}