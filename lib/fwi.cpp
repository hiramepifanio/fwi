#include <iostream>
#include <fftw3.h>
#include "memory.cpp"
//#include "io.cpp"
#include <math.h>

using namespace std;

// int **receptors(int step, int nz, int nx);
// float* ricker(float freq, float h, int N);
// float** velocity_map(int Lz, int Lx);
// float cfl_criteria(float **vel, int n, int m, float dz, float dx, float dt);
// void minimum_wavelenth(float **vel, int n, int m, float fp, float dx);
// float **laplacian(float **wave, int nz, int nx, float dz, float dx);
// float **extendModel(float **vel, int nz, int nx, int border);
// float **taper(int nz, int nx, int border);
// void buffer_wave(float **buffer, float **wave, int nz, int nx);
// float **propagator(float *wavelet, float dt, int nt, int souz, int soux, int **rec, int nrecs, float **vel, float dz, float dx, int nz, int nx, int border);

// Receptors position
int **receptors(int step, int nz, int nx){
    int nr = nx/step + 1; // Number of receptors
    int **rec = allocateMint(2, nr); // Receptors
    // Filling the positions
    for(int r = 0; r < nr; r++){
        rec[0][r] = 2; // x
        rec[1][r] = r*step; // z
    }

    return rec;
}

// Returns the Ricker wavelet values
float* ricker(float freq, float h, int N){
    float t0 = 6/(M_PI*freq*sqrt(2));
    float *ric = allocateV(N);
    float t;

    for(int i = 0; i < N; i++){
        t = i*h;
        ric[i] = (1 - 2*M_PI*M_PI*freq*freq*(t - t0)*(t - t0))*exp(-M_PI*M_PI*freq*freq*(t - t0)*(t - t0));//pow(M_E, -M_PI*M_PI*freq*freq*(t - t0)*(t - t0));
    }

    return ric;
}

// Generates the velocity map
float** velocity_map(int Lz, int Lx){
    float **M = allocateM(Lz, Lx);

    for(int i = 0; i < Lz; i++){
        for(int j = 0; j < Lx; j++){
            if(i < 40)
                M[i][j] = 2000.0;
            else
                M[i][j] = 3000.0;
        }
    }

    return M;
}

// CFL criteria
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

    return dtp;
}

// Minimum wavelength criteria
void minimum_wavelenth(float **vel, int n, int m, float fp, float dx){
    float vmin = vel[0][0];
    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            if(vel[i][j] < vmin)
                vmin = vel[i][j];
    float minlambda = vmin/(3*fp);
    cout << "The minimum wavelength is " << minlambda << " and should be equal to or bigger than " << 5*dx << ": ";
    if(minlambda >= 5*dx)
        cout << "True." << endl;
    else
        cout << "False." << endl;
}

// Returns the laplacian of the wavestate
float **laplacian(float **wave, int nz, int nx, float dz, float dx){
    float **lap = valueM(0.0, nz, nx);
    float ddx;

    // Derivative in z
    // Center
    for(int i = 2; i < nz - 2; i++){
        for(int j = 0; j < nx; j++){
            lap[i][j] = (- wave[i + 2][j] + 16*wave[i + 1][j] - 30*wave[i][j] + 16*wave[i - 1][j] - wave[i - 2][j])/(12*dz*dz);
        }
    }
    // The 4 outer rows
    for(int j = 0; j < nx; j++){
        lap[0][j] = lap[1][j] = (wave[2][j] - 2*wave[1][j] + wave[0][j])/(dz*dz);
        lap[nz - 2][j] = lap[nz - 1][j] = (wave[nz - 1][j] - 2*wave[nz - 2][j] + wave[nz - 3][j])/(dz*dz);
    }

    // Derivative in x
    // Center
    for(int j = 2; j < nx - 2; j++){
        for(int i = 0; i < nz; i++){
            lap[i][j] += (- wave[i][j + 2] + 16*wave[i][j + 1] - 30*wave[i][j] + 16*wave[i][j - 1] - wave[i][j - 2])/(12*dx*dx);
        }
    }
    // The 4 outer colums
    for(int i = 0; i < nz; i++){
        ddx = (wave[i][2] - 2*wave[i][1] + wave[i][0])/(dx*dx);
        lap[i][0] += ddx;
        lap[i][1] += ddx;
        ddx = (wave[i][nx - 1] - 2*wave[i][nx - 2] + wave[i][nx - 3])/(dx*dx);
        lap[i][nx - 2] += ddx;
        lap[i][nx - 1] += ddx;
    }

    return lap;
}

float **extendModel(float **vel, int nz, int nx, int border){
    float **model = valueM(0.0, nz + 2*border, nx + 2*border);

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
    float **buffer = valueM(1.0, nz + 2*border, nx + 2*border);
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

float *sincint(float *wav, int Nin, float dt, float dtp) {
    if(int(dt/dtp) == 1) {
        return wav;
    }

    float* wavp = wav;
    // float *wavp;
    // double *wavDouble = new double [Nin];
    // fftw_complex *out;
    // fftw_plan plan;

    // for(int i = 0; i < Nin; i++) {
    //     wavDouble[i] = wav[i];
    // }
    // out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (Nin/2 + 1));

    // plan = fftw_plan_dft_r2c_1d(Nin, wavDouble, out, FFTW_MEASURE);
    
    // fftw_execute(plan);


    // cout << "Foi" << endl;

    return wavp;
}

float **propagator(float *wav, float dt, int nt, int souz, int soux, int **rec, int nrecs, float **vel, float dz, float dx, int nz, int nx, int border){
    int nzb = nz + 2*border, nxb = nx + 2*border; // Full grid with borders
    float **model = extendModel(vel, nz, nx, border); // Velocity model with borders
    float **buffer = taper(nz, nx, border); // Scale factor to buffer the wave at the borders
    float **wave0 = valueM(0.0, nzb, nxb); // Wave in the past
    float **wave1 = valueM(0.0, nzb, nxb); // Wave in the present
    float **wave2 = valueM(0.0, nzb, nxb); // Wave in the future
    float **shot = valueM(0.0, nt, nrecs); // Information about the intensity of the wave at the receptors

    // Main loop
    cout << "Initiating propagator...\n";
    for(int k = 0; k < nt; k++){
        if(k%(nt/100) == 0){
            cout << "Progress... " << 100*k/nt << "%\n";
        }
        float **lap = laplacian(wave1, nzb, nxb, dz, dx);

        // Input source
        wave1[souz + border][soux + border] += vel[souz + border][soux + border]*vel[souz + border][soux + border]*dt*dt*wav[k];

        // Collect wave intensity at receptors
        for(int r = 0; r < nrecs; r++){
            shot[k][r] = wave1[border + rec[0][r]][border + rec[1][r]];
        }

        // Wave propagation
        for(int i = 0; i < nzb; i++){
            for(int j = 0; j < nxb; j++){
                wave2[i][j] = 2*wave1[i][j] - wave0[i][j] + model[i][j]*model[i][j]*dt*dt*lap[i][j];
            }
        }
        freeM(lap, nzb);

        // Buffering
        buffer_wave(buffer, wave2, nzb, nxb);
        buffer_wave(buffer, wave1, nzb, nxb);

        // Roll arrays forward
        copy_M_to(wave1, wave0, nzb, nxb);
        copy_M_to(wave2, wave1, nzb, nxb);
    }

    // Free memory
    freeM(wave0, nzb);
    freeM(wave1, nzb);
    freeM(wave2, nzb);
    freeM(model, nzb);
    freeM(buffer, nzb);

    return shot;
}