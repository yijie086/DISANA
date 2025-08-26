#include <iostream>
#include <string>
#include <thread>
#include "TError.h"
void show_runtime_bar(double seconds) {
  const int max_width = 30;      // total width of the bar
  const double max_time = 10.0;  // assume 10 seconds = full bar

  int bar_width = static_cast<int>((seconds / max_time) * max_width);
  if (bar_width > max_width) bar_width = max_width;

  std::cout << "[";
  for (int i = 0; i < bar_width; ++i) std::cout << "█";
  for (int i = bar_width; i < max_width; ++i) std::cout << " ";
  std::cout << "]" << std::endl;
}

void RunPhiAnalysis(const std::string& inputFile, int nfile, int nthreads);
int main(int argc, char* argv[]) {

  // Check if the correct number of arguments were passed
  if (argc !=4) {
    std::cerr << "Usage: ./AnalysisDVCS <path_to_hipo_file> <number_of_files>" << std::endl;
    std::cerr << "Example: ./AnalysisDVCS /..pathtohipofiles/ 1 OR CHOOSE LARGE NUMBER TO SELECT ALL THE FILES IN THE DIR" << std::endl;
    return 1;
  }
  auto start = std::chrono::high_resolution_clock::now();  
  // Pass the command-line argument to FirstDVCSAnalysis
  std::string inputFilePath = argv[1];  // argv[1] is the first argument passed to the program
  RunPhiAnalysis(inputFilePath, std::stoi(argv[2]), std::stoi(argv[3]));
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;

  std::cout << "Total running time: " << elapsed.count() << " seconds\n";
  show_runtime_bar(elapsed.count());
  return 0;
}
