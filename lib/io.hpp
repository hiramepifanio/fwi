void save_V_to_file(float *V, int N, std::string file_name);
void save_M_to_file(float **M, int n, int m, std::string file_name);
void save_Mint_to_file(int **M, int n, int m, std::string file_name);
void save_T_to_file(float ***T, int numPlanes, int numRows, int numCols, std::string file_name);
float **loadMbin(std::string file_name, int n, int m);

void saveArray(float *A, int n, std::string folderPath, std::string filename);
void save1Array(float* A, int* shape, std::string folderPath, std::string filename);
void saveArrayBin(float **A, int* shape, std::string folderPath, std::string filename);
void save2Arr(float** A, int* shape, std::string folderPath, std::string filename);
void saveShot(float** shot, int* shape, std::string folderPath, std::string filename);