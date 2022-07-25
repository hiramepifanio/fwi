#include <iostream>
#include "../lib/fwi.hpp"
#include "../lib/io.hpp"
#include "../lib/ricker.hpp"

using namespace std;

int main() {

    float fp = 13;
    float dt = 0.002;
    int nt = 201;

    float* wav = ricker(fp, dt, nt);
    int shape[] = {nt};
    save1Array(wav, shape, "/home/hiram/vscode-workspace/fwi/ricker_wavelet/data", "ricker_wavelet");

    return 0;
}