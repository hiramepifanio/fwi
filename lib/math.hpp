float max(float** A, int* shape);
float min(float** A, int* shape);
float **laplacian(float **wave, int* wave_shape, float dz, float dx);
float** vel_dt_2(float** vel, int* vel_shape, float dt);