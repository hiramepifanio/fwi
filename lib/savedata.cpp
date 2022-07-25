#include <iostream>
#include <fstream>
#include <cassert>
#include "memory.hpp"
using namespace std;

// Write 1d
void write_float1d(fstream* file, float* A, int n) {
    for(int k = 0; k < n; k++){
        *file << A[k] << endl;
    }
}

void write_int1d(fstream* file, int* A, int n) {
    for(int k = 0; k < n; k++){
        *file << A[k] << endl;
    }
}

// Save array 1d
void save_array_int1d(int *A, int n, string folderPath, string filename){
    fstream file;
    string path = folderPath + "/" + filename;

    file.open(path, ios::out);
    assert(file.is_open());
    write_int1d(&file, A, n);
    file.close();
}

// Save array 2d
void save_array_int2d(int **A, int* shape, string folderPath, string filename){
    fstream file;
    string path = folderPath + "/" + filename;

    file.open(path, ios::out);
    for(int i = 0; i < shape[0]; i++){
        write_int1d(&file, A[i], shape[1]);
    }
    file.close();
}

void save_array_float2d(float **A, int* shape, string folderPath, string filename){
    fstream file;
    string path = folderPath + "/" + filename;

    file.open(path, ios::out);
    for(int i = 0; i < shape[0]; i++){
        write_float1d(&file, A[i], shape[1]);
    }
    file.close();
}

// Save data
void save_data_int1d(int* A, int len, string folderPath, string filename) {
    save_array_int1d(A, len, folderPath, filename + ".data");
}

void save_data_int2d(int** A, int* shape, string folderPath, string filename) {
    cout << "Saving data to "+folderPath+"/"+filename+"... ";

    save_array_int2d(A, shape, folderPath, filename + ".data");
    save_array_int1d(shape, 2, folderPath, filename + ".shp");

    cout << "OK." << endl;
}

void save_data_float2d(float** A, int* shape, string folderPath, string filename) {
    cout << "Saving data to "+folderPath+"/"+filename+"... ";
    
    save_array_float2d(A, shape, folderPath, filename + ".data");
    save_array_int1d(shape, 2, folderPath, filename + ".shp");

    cout << "OK." << endl;
}