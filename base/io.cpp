#include <iostream>
#include <fstream>
#include <cassert>
#include "memory.hpp"
using namespace std;

// Save a vector to file
void save_V_to_file(float *V, int N, string file_name){
    fstream file;
    file.open(file_name, ios::out);
    assert(file.is_open());
    for(int k = 0; k < N; k++){
        file << V[k] << endl;
    }
    file.close();
}

// Save a vector to file
void save_M_to_file(float **M, int n, int m, string file_name){
    fstream file;
    file.open(file_name, ios::out);
    assert(file.is_open());
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            file << M[i][j] << " ";
        }
        file << endl;
    }
    file.close();
}

// Save a vector to file
void save_Mint_to_file(int **M, int n, int m, string file_name){
    fstream file;
    file.open(file_name, ios::out);
    assert(file.is_open());
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            file << M[i][j] << " ";
        }
        file << endl;
    }
    file.close();
}

// Save a vector to file
void save_T_to_file(float ***T, int numPlanes, int numRows, int numCols, string file_name){
    fstream file;
    file.open(file_name, ios::out);
    assert(file.is_open());
    for(int k = 0; k < numPlanes; k++){
        for(int i = 0; i < numRows; i++){
            for(int j = 0; j < numCols; j++){
                file << T[k][i][j] << " ";
            }
            file << endl;
        }
    }
    file.close();
}

// Loads a float matrix from binary file
float **loadMbin(string file_name, int n, int m){
    streampos length;
    float *memblock;
    float **M = allocateM(n, m);

    ifstream file (file_name, ios::in | ios::binary | ios::ate);
    if(file.is_open()){
        length = file.tellg();
        memblock = new float [length/sizeof(float)];
        file.seekg(0, ios::beg);
        file.read((char *)memblock, length);
        file.close();
    }

    for(int j = 0; j < m; j++){
        for(int i = 0; i < n; i++){
            M[i][j] = memblock[i + j*n];
        }
    }

    freeV(memblock);

    return M;
}

// Saves a float matrix to binary file
void saveMbin(float **M, int n, int m, string filename){
    ofstream out;

    out.open(filename, ios::out | ios::binary);
    for(int j = 0; j < m; j++){
        for(int i = 0; i < n; i++){
            out.write(reinterpret_cast<const char*>(&M[i][j]), sizeof(float));
        }
    }
    out.close();
}