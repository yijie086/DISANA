#include <string>
#include <iostream>
#include <chrono> // 用于计时
#include "RHipoDS.hxx"
#include "TCanvas.h"
#include <optional>

using namespace ROOT;
using namespace ROOT::RDF;

void Run103ECALfiducialcut(const float beam_energy, const std::string& FilePath);

int main(int argc, char **argv) {

    auto start_time = std::chrono::high_resolution_clock::now();

    if(argc < 3){
        std::cout << "ERROR: WRONG INPUT! usage: ./main101 <beam_energy> <hipo_file_folder_path>\n";
        return 1;
    }else{
        std::cout << "Opening files in " << argv[2] << std::endl;
        std::cout << "Beam energy: " << argv[1] << std::endl;
    }

    float beam_energy = std::stof(argv[1]);
    Run103ECALfiducialcut(beam_energy, argv[2]);

    auto end_time = std::chrono::high_resolution_clock::now();


    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    std::cout << "Total runtime: " << duration << " seconds" << std::endl;

    return 0;
}