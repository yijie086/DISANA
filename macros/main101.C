#include <string>
#include <iostream>
#include "RHipoDS.hxx"
#include "TCanvas.h"
#include <optional>

using namespace ROOT;
using namespace ROOT::RDF;

void Run101(const std::string& FilePath);

int main(int argc, char **argv) {
    //ROOT::EnableImplicitMT(2);
    if(argc < 2){
        std::cout << "Please specify a HIPO data file on the command line. (Only one file.) \n";
        return 1;
    }else{
        std::cout << "Opening file " << argv[1] << std::endl;
    }

    Run101(argv[1]);
}