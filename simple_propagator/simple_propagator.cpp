#include <iostream>
#include "../lib/alloc.hpp"
#include "../lib/ricker.hpp"
#include "../lib/math.hpp"
#include "../lib/cond.hpp"
#include "../lib/savedata.hpp"
#include "../lib/memory.hpp"

using namespace std;

float** velocity_map(int nz, int nx);

int main() {

    string script = "simple_propagator";
    cout << "\nStarting script: " << script << "." << endl;

    // Parameters
    int nz = 101, nx = 151, nt = 800;       // Number of grid points and timesteps
    float dz = 15, dx = 15, dt = 0.003;     // Step and timestep sizes
    int souz = nz/2, soux = nx/2;           // Source grid position    
    float fp = 5;                           // Ricker wavelet peak frequency
    float* wav = ricker(fp, dt, nt);        // Source input: ricker wavelet
    float** vel = velocity_map(nz, nx);     // Velocity map

    int vel_shape[] = {nz, nx};
    verify_cfl_condition(dz, dx, dt, max(vel, vel_shape));
    verify_spatial_condition(dz, dx, min(vel, vel_shape), fp);

    save_data_float2d(vel, vel_shape, "/home/hiram/vscode-workspace/fwi/"+script+"/data", "velocity_map");
    
    int source_shape[] = {2, 1};
    int** source = alloc_int2d(source_shape);
    source[0][0] = souz, source[1][0] = soux;
    save_data_int2d(source, source_shape, "/home/hiram/vscode-workspace/fwi/"+script+"/data", "source");

    int wave_shape[] = {nz, nx};
    float **wave0 = alloc_float2d(0.0, wave_shape); // Wave in the past
    float **wave1 = alloc_float2d(0.0, wave_shape); // Wave in the present
    float **wave2 = alloc_float2d(0.0, wave_shape); // Wave in the future
    float** con = vel_dt_2(vel, vel_shape, dt);     // vel^2*dt^2
    cout << "Initiating propagator...\n";
    for (int t = 0; t < nt; t++) {
        float **lap = laplacian(wave1, wave_shape, dz, dx);

        // Find new wave values
        for(int i = 0; i < nz; i++)
            for(int j = 0; j < nx; j++)
                wave2[i][j] = 2*wave1[i][j] - wave0[i][j] + con[i][j]*lap[i][j];

        // Set wave value at source
        wave2[souz][soux] += con[souz][soux]*wav[t];

        // Save wavestate
        if (t%60 == 0)
            save_data_float2d(wave1, wave_shape, "/home/hiram/vscode-workspace/fwi/"+script+"/data", "wave_"+to_string(t));
        
        // Prepare to next timestep
        transfer_float2d(wave1, wave0, wave_shape);
        transfer_float2d(wave2, wave1, wave_shape);
    }

    save_data_float2d(wave1, wave_shape, "/home/hiram/vscode-workspace/fwi/"+script+"/data", "wave_800");

    cout << "Script finished: " << script << "." << endl;

    return 0;
}

float** velocity_map(int nz, int nx) {
    int shape[] = {nz, nx};
    float** vel = alloc_float2d(2000, shape);

    for (int i = 70; i < nz; i++)
        for (int j = 0; j < nx; j++)
            vel[i][j] = 3000;

    return vel;
}