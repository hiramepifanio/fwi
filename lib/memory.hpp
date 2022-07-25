float *allocateV(int numItems);
float **allocateM(int numRows, int numCols);
int **allocateMint(int numRows, int numCols);
float ***allocateT(int numPlanes, int numRows, int numCols);
float **valueM(float value, int numRows, int numCols);

int* alloc1IntArr(int* shape);
int* alloc1IntArr(int value, int* shape);
float* alloc1Arr(int numRows);

float** alloc2Arr(int numRows, int numCols);
float** alloc2Arr(float value, int numRows, int numCols);

float*** alloc3Arr(int numPlanes, int numRows, int numCols);
float*** alloc3Arr(float value, int numPlanes, int numRows, int numCols);

void freeV(float *V);
void freeM(float **M, int numRows);
void freeMint(int **M, int numRows);
void freeT(float ***T, int numPlanes, int numRows);

void transfer_float2d(float **A, float **B, int* shape);