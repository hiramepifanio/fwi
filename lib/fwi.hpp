int **receptors(int step, int nx, int z);
int **receptors(int step, int nx);
float* ricker(float freq, float h, int N);
float** velocity_map(int Lz, int Lx);

float cfl_criteria(float **vel, int n, int m, float dz, float dx, float dt);
int cfl_criteria_scale(float **vel, int nz, int nx, float dt, float dz, float dx);

void minimum_wavelenth(float **vel, int n, int m, float fp, float dx);
float **laplacian(float **wave, int nz, int nx, float dz, float dx);
float **extend_model(float **vel, int nz, int nx, int border);
float **taper(int nz, int nx, int border);
void buffer_wave(float **buffer, float **wave, int nz, int nx);
float **propagator(float *wavelet, float dt, int nt, int souz, int soux, int **rec, int nrecs, float **vel, float dz, float dx, int nz, int nx, int border);
float *sincint(float *wav, int Nin, float dt, float dtp);
float** sincint_shot(float** shot, float dt, float dtp, int nt, int ntp, int nrecs);