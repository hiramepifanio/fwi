int* alloc_int1d(int n);
int* alloc_int1d(int value, int n);
float* alloc_float1d(int n);
float* alloc_float1d(float value, int n);

int** alloc_int2d(int* shape);
float** alloc_float2d(int* shape);
float** alloc_float2d(float value, int* shape);