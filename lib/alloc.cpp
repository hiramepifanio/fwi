// Alloc 1D arrays
int* alloc_int1d(int n) {
    int *A = new int [n];

    return A;
}

int* alloc_int1d(int value, int n) {
    int *A = new int [n];

    for (int i = 0; i < n; i++)
        A[i] = value;

    return A;
}

float* alloc_float1d(int n) {
    float *A = new float [n];

    return A;
}

float* alloc_float1d(float value, int n) {
    float *A = new float [n];

    for (int i = 0; i < n; i++)
        A[i] = value;

    return A;
}

// Alloc 2D arrays
int** alloc_int2d(int* shape) {
    int **A = new int *[shape[0]];

    for(int i = 0; i < shape[0]; i++)
        A[i] = alloc_int1d(shape[1]);

    return A;
}

float** alloc_float2d(int* shape) {
    float **A = new float *[shape[0]];

    for(int i = 0; i < shape[0]; i++)
        A[i] = alloc_float1d(shape[1]);

    return A;
}

float** alloc_float2d(float value, int* shape) {
    float **A = new float *[shape[0]];

    for(int i = 0; i < shape[0]; i++)
        A[i] = alloc_float1d(value, shape[1]);
                
    return A;
}

// Alloc 3D arrays
float** alloc_float3d(int* shape) {
    float **A = new float *[shape[0]];


    for(int i = 0; i < shape[0]; i++)
        A[i] = alloc_float1d(shape[1]);

    return A;
}

float** alloc_float3d(float value, int* shape) {
    float **A = new float *[shape[0]];

    for(int i = 0; i < shape[0]; i++)
        A[i] = alloc_float1d(value, shape[1]);
                
    return A;
}