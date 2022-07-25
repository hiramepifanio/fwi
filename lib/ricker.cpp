#include <iostream>
#include "memory.hpp"
#include <math.h>

using namespace std;

// Returns the Ricker wavelet values
float* ricker(float freq, float h, int N){
    cout << "Generating ricker wavelet... ";
    float t0 = 6/(M_PI*freq*sqrt(2));
    float *ric = allocateV(N);
    float t;

    for(int i = 0; i < N; i++){
        t = i*h;
        ric[i] = (1 - 2*M_PI*M_PI*freq*freq*(t - t0)*(t - t0))*exp(-M_PI*M_PI*freq*freq*(t - t0)*(t - t0));
    }

    cout << "Ok." << endl;

    return ric;
}