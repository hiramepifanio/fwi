#include <iostream>
#include "math.hpp"
#include <math.h>

using namespace std;

float hefec(float dz, float dx) {
    return 1/pow(1/(dx*dx) + 1/(dz*dz), 0.5);
}

float cfl_condition(float dz, float dx, float vel_max) {
    float h = hefec(dz, dx);
    float cfl_dt = pow(3.0/4, 0.5)*h/vel_max;

    return cfl_dt;
}

bool verify_cfl_condition(float dz, float dx, float dt, float vel_max) {
    float cfl_dt = cfl_condition(dz, dx, vel_max);
    bool isSatisfied = (dt <= cfl_dt);

    cout << "CFL condition: ";
    if (isSatisfied)
        cout << "OK";
    else
        cout << "Failed";
    cout << " - The dt is " << 1000*dt << " ms and should be <= " << 1000*cfl_dt << " ms." << endl;

    return isSatisfied;
}

float spatial_condition(float dz, float dx, float vel_min) {
    float D = 5;
    float h = hefec(dz, dx);
    float f_max = vel_min/(3*D*h);

    return f_max;
}

bool verify_spatial_condition(float dz, float dx, float vel_min, float fp) {
    float f_max = spatial_condition(dz, dx, vel_min);

    bool isSatisfied = (fp <= f_max);

    cout << "Spatial condition: ";
    if (isSatisfied)
        cout << "OK";
    else
        cout << "Failed";
    cout << " - The fp is " << fp << " Hz and should be <= " << f_max << " Hz." << endl;

    return isSatisfied;
}