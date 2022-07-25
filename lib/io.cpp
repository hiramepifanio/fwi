#include <iostream>
#include <fstream>
#include <cassert>
#include "memory.hpp"
using namespace std;

// Write an 1-dimensional array to a currently open file, intended to be an auxiliary function
void write1D(fstream* file, float* A, int n) {
    for(int k = 0; k < n; k++){
        *file << A[k] << endl;
    }
}

void write1D(fstream* file, int* A, int n) {
    for(int k = 0; k < n; k++){
        *file << A[k] << endl;
    }
}


// Save an array to file
void saveArray(float *A, int n, string folderPath, string filename){
    cout << "Saving array to "+folderPath+"/"+filename+"... ";
    fstream file;
    string path = folderPath + "/" + filename;

    file.open(path, ios::out);
    assert(file.is_open());
    write1D(&file, A, n);
    file.close();

    cout << "Donne." << endl;
}

void saveArray(int *A, int n, string folderPath, string filename){
    cout << "Saving array to "+folderPath+"/"+filename+"... ";
    fstream file;
    string path = folderPath + "/" + filename;

    file.open(path, ios::out);
    assert(file.is_open());
    write1D(&file, A, n);
    file.close();

    cout << "Donne." << endl;
}

void saveArray(float **A, int* shape, string folderPath, string filename){
    fstream file;
    string path = folderPath + "/" + filename;

    file.open(path, ios::out);
    assert(file.is_open());
    for(int i = 0; i < shape[0]; i++){
        write1D(&file, A[i], shape[1]);
    }
    file.close();
}

void saveArray(int **A, int* shape, string folderPath, string filename){
    fstream file;
    string path = folderPath + "/" + filename;

    file.open(path, ios::out);
    assert(file.is_open());
    for(int i = 0; i < shape[0]; i++){
        write1D(&file, A[i], shape[1]);
    }
    file.close();
}

void saveArray(float ***A, int* shape, string folderPath, string filename){
    fstream file;
    string path = folderPath + "/" + filename;

    file.open(path, ios::out);
    assert(file.is_open());
    for(int i = 0; i < shape[0]; i++){
        for(int j = 0; j < shape[1]; j++){
            write1D(&file, A[i][j], shape[2]);
        }
    }
    file.close();
}

// Saves an array to binary file
void saveArrayBin(float **A, int* shape, string folderPath, string filename){
    ofstream out;
    string path = folderPath + "/" + filename;

    out.open(path, ios::out | ios::binary);
    for(int i = 0; i < shape[1]; i++){
        for(int j = 0; j < shape[0]; j++){
            out.write(reinterpret_cast<const char*>(&A[j][i]), sizeof(float));
        }
    }
    out.close();
}

void save2ArrBin(float **A, int* shape, string folderPath, string filename){
    ofstream out;
    string path = folderPath + "/" + filename;

    out.open(path, ios::out | ios::binary);
    for(int i = 0; i < shape[0]; i++){
        for(int j = 0; j < shape[1]; j++){
            out.write(reinterpret_cast<const char*>(&A[i][j]), sizeof(float));
        }
    }
    out.close();
}

void saveArray2FileConfig(int* shape, int dimensions, string folderPath, string filename) {
    saveArray(shape, dimensions, folderPath, filename + ".config");
}

void save1Array(float* A, int* shape, string folderPath, string filename) {
    saveArray(A, shape[0], folderPath, filename + ".txt");
    saveArray2FileConfig(shape, 1, folderPath, filename);
}

void saveShot(float** shot, int* shape, string folderPath, string filename) {
    saveArrayBin(shot, shape, folderPath, filename + ".bin");
    saveArray2FileConfig(shape, 2, folderPath, filename);
}

void save2Arr(float** A, int* shape, string folderPath, string filename) {
    save2ArrBin(A, shape, folderPath, filename + ".bin");
    saveArray2FileConfig(shape, 2, folderPath, filename);
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