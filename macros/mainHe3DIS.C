#include <chrono>
#include <iostream>
#include <string>

void RunHe3DISAnalysis(const std::string& inputFile, int nfile, int nthreads);

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "Usage: ./AnalysisHe3DIS <path_to_hipo_files> "
                 "<number_of_files> <number_of_threads>\n";
    return 1;
  }

  const auto start = std::chrono::high_resolution_clock::now();
  RunHe3DISAnalysis(argv[1], std::stoi(argv[2]), std::stoi(argv[3]));
  const auto end = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> elapsed = end - start;
  std::cout << "Total running time: " << elapsed.count() << " seconds\n";
  return 0;
}
