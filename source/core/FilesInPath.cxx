#include "FilesInPath.h"
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

std::vector<std::string> GetHipoFilesInPath(const std::string& directory) {
    std::vector<std::string> files;
    for (const auto& entry : fs::recursive_directory_iterator(directory)) {
        if (entry.path().extension() == ".hipo") {
            files.push_back(entry.path().string());
            std::cout << "Found file: " << entry.path() << std::endl;
        }
    }
    std::cout << "************ " << files.size() << " Files Found ************" << std::endl;
    return files;
}