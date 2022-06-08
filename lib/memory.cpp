// Allocate memory for a vector[]
float *allocateV(int numItems){
    float *V = new float [numItems];

    return V;
}

// Allocate memory for a matrix[][]
float **allocateM(int numRows, int numCols){
    float **M = new float *[numRows];

    for(int i = 0; i < numRows; i++){
        M[i] = new float [numCols];
    }

    return M;
}

// Allocate memory for a matrix[][]
int **allocateMint(int numRows, int numCols){
    int **M = new int *[numRows];

    for(int i = 0; i < numRows; i++){
        M[i] = new int [numCols];
    }

    return M;
}

// Allocate memory for a tensor[][][]
float ***allocateT(int numPlanes, int numRows, int numCols){
    float ***T = new float **[numPlanes];

    for(int k = 0; k < numPlanes; k++){
        T[k] = new float *[numRows];
        for(int i = 0; i < numRows; i++){
            T[k][i] = new float [numCols];
        }
    }

    return T;
}

// Allocate memory for a matrix filled with float(value)
float **valueM(float value, int numRows, int numCols){
    float **M = allocateM(numRows, numCols);

    for(int i = 0; i < numRows; i++){
        for(int j = 0; j < numCols; j++){
            M[i][j] = value;
        }
    }

    return M;
}

// Allocate memory for 3-Array filled with float(value)
float** alloc2Arr(int numRows, int numCols){
    float **A = new float *[numRows];

    for(int k = 0; k < numRows; k++){
        A[k] = new float [numCols];
    }

    return A;
}

float** alloc2Arr(float value, int numRows, int numCols){
    float **A = alloc2Arr(numRows, numCols);

    for(int i = 0; i < numRows; i++)
        for(int j = 0; j < numCols; j++)
            A[i][j] = value;
                
    return A;
}

// Allocate memory for 3-Array filled with float(value)
float*** alloc3Arr(int numPlanes, int numRows, int numCols){
    float ***A = new float **[numPlanes];

    for(int k = 0; k < numPlanes; k++){
        A[k] = new float *[numRows];
        for(int i = 0; i < numRows; i++){
            A[k][i] = new float [numCols];
        }
    }

    return A;
}

float*** alloc3Arr(float value, int numPlanes, int numRows, int numCols){
    float ***A = alloc3Arr(numPlanes, numRows, numCols);

    for(int k = 0; k < numPlanes; k++)
        for(int i = 0; i < numRows; i++)
            for(int j = 0; j < numCols; j++)
                A[k][i][j] = value;

    return A;
}

// Deallocate memory for a vector
void freeV(float *V){
    delete[] V;
}

// Deallocate memory for a matrix
void freeM(float **M, int numRows){
    for(int i = 0; i < numRows; i ++){
        delete[] M[i];
    }
    delete[] M;
}

// Deallocate memory for a matrix
void freeMint(int **M, int numRows){
    for(int i = 0; i < numRows; i ++){
        delete[] M[i];
    }
    delete[] M;
}

// Deallocate memory for a tensor
void freeT(float ***T, int numPlanes, int numRows){
    for(int k = 0; k < numPlanes; k++){
        for(int i = 0; i < numRows; i++){
            delete[] T[k][i];
        }
        delete[] T[k];
    }
    delete[] T;
}

void copy_M_to(float **M1, float **M2, int numRows, int numCols){
    for(int i = 0; i < numRows; i++){
        for(int j = 0; j < numCols; j++){
            M2[i][j] = M1[i][j];
        }
    }
}