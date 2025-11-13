// AnalysisDVCS.cpp
#include <iostream>
#include <string>
#include <thread>
#include <chrono>
#include <cstdlib>
#include <vector>

// ROOT (you already included TError.h in your macro file; not strictly needed here)
#include "TError.h"

// ----------------------------------------------------------------------
// Progress bar for total runtime
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

// ----------------------------------------------------------------------
// Your macro function (defined elsewhere). Keep the order exactly as below.
void RunPhiAnalysis(const std::string& inputDir, int nfile, int nthreads,
                    const std::string& outputDir,
                    const std::string dataconfig,
                    bool IsMC, bool IsreprocRootFile,
                    bool IsInbending /* you can ignore and derive from dataconfig if you like */,
                    bool IsMinimalBook);

// ----------------------------------------------------------------------
static void print_usage(const char* prog) {
  std::cout <<
R"(Usage:
  # Flag form (recommended)
  )" << prog << R"( -i <inputDir> [-n <nfile>] [-t <nthreads>] [-o <outputDir>]
                [-c <dataconfig>] [--mc] [--reproc] [--inbending] [--minimal]

  # Positional form (quick & dirty; anything omitted uses defaults)
  )" << prog << R"( [inputDir] [nfile] [nthreads] [outputDir] [dataconfig] [mc] [reproc] [inbend] [minimal]

Defaults:
  inputDir="."   nfile=-1 (all)   nthreads=0 (auto)
  outputDir="./out"   dataconfig="rgasp18_outb"
  mc/reproc/inbending/minimal = false

Examples:
  )" << prog << R"( -i /data/phi/rgasp18_outb -n 200 -t 8 -o results -c rgasp18_outb
  )" << prog << R"( /data 100 0 ./out rgasp18_outb 0 0 1 0
  )" << std::endl;
}

// Parse an integer from a C-string, returning true on success.
static bool parse_int(const char* s, int& out) {
  if (!s) return false;
  char* end = nullptr;
  long v = std::strtol(s, &end, 10);
  if (end == s || *end != '\0') return false;
  out = static_cast<int>(v);
  return true;
}

// Parse a bool from "0/1", "true/false", "on/off", "yes/no"
static bool parse_bool(const char* s, bool& out) {
  if (!s) return false;
  std::string v(s);
  for (auto& ch : v) ch = static_cast<char>(std::tolower(ch));
  if (v == "1" || v == "true" || v == "on" || v == "yes") { out = true; return true; }
  if (v == "0" || v == "false" || v == "off" || v == "no") { out = false; return true; }
  return false;
}

int main(int argc, char* argv[]) {
  // -------- Defaults that match your macro's intended behavior
  std::string inputDir   = ".";
  int         nfile      = -1;          // -1 => "all files" convention
  int         nthreads   = 0;           // 0  => auto (macro will pick)
  std::string outputDir  = "./out";
  std::string dataconfig = "rgasp18_outb";
  bool IsMC              = false;
  bool IsreprocRootFile  = false;
  bool IsInbending       = false;
  bool IsMinimalBook     = false;

  // -------- If no args or help requested, show usage
  if (argc == 1) {
    print_usage(argv[0]);
    // not an error; user just asked for help / no args
    return 0;
  }
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "-h" || a == "--help") {
      print_usage(argv[0]);
      return 0;
    }
  }

  // -------- First try flag-style parsing
  // Recognized flags:
  // -i/--input, -n/--nfile, -t/--threads, -o/--out, -c/--config
  // --mc, --reproc, --inbending, --minimal (boolean toggles)
  bool saw_flag = false;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a.rfind("-", 0) != 0) continue;  // not a flag
    saw_flag = true;

    if (a == "-i" || a == "--input") {
      if (i + 1 >= argc) { std::cerr << "ERROR: " << a << " requires a value\n"; return 2; }
      inputDir = argv[++i];
    } else if (a == "-n" || a == "--nfile") {
      if (i + 1 >= argc) { std::cerr << "ERROR: " << a << " requires an integer\n"; return 2; }
      if (!parse_int(argv[++i], nfile)) { std::cerr << "ERROR: invalid nfile\n"; return 2; }
    } else if (a == "-t" || a == "--threads") {
      if (i + 1 >= argc) { std::cerr << "ERROR: " << a << " requires an integer\n"; return 2; }
      if (!parse_int(argv[++i], nthreads)) { std::cerr << "ERROR: invalid nthreads\n"; return 2; }
    } else if (a == "-o" || a == "--out") {
      if (i + 1 >= argc) { std::cerr << "ERROR: " << a << " requires a value\n"; return 2; }
      outputDir = argv[++i];
    } else if (a == "-c" || a == "--config") {
      if (i + 1 >= argc) { std::cerr << "ERROR: " << a << " requires a value\n"; return 2; }
      dataconfig = argv[++i];
    } else if (a == "--mc") {
      IsMC = true;
    } else if (a == "--reproc") {
      IsreprocRootFile = true;
    } else if (a == "--inbending") {
      IsInbending = true;
    } else if (a == "--minimal") {
      IsMinimalBook = true;
    } else {
      std::cerr << "ERROR: unknown option: " << a << "\n";
      print_usage(argv[0]);
      return 2;
    }
  }

  // -------- If the user didn’t use flags, accept positional args:
  // [inputDir] [nfile] [nthreads] [outputDir] [dataconfig] [mc] [reproc] [inbend] [minimal]
  if (!saw_flag) {
    int pos = 1;
    if (pos < argc) inputDir = argv[pos++];
    if (pos < argc) { if (!parse_int(argv[pos++], nfile)) { std::cerr << "ERROR: nfile must be int\n"; return 2; } }
    if (pos < argc) { if (!parse_int(argv[pos++], nthreads)) { std::cerr << "ERROR: nthreads must be int\n"; return 2; } }
    if (pos < argc) outputDir = argv[pos++];
    if (pos < argc) dataconfig = argv[pos++];
    if (pos < argc) { if (!parse_bool(argv[pos++], IsMC)) { std::cerr << "ERROR: mc must be bool (0/1,true/false)\n"; return 2; } }
    if (pos < argc) { if (!parse_bool(argv[pos++], IsreprocRootFile)) { std::cerr << "ERROR: reproc must be bool\n"; return 2; } }
    if (pos < argc) { if (!parse_bool(argv[pos++], IsInbending)) { std::cerr << "ERROR: inbend must be bool\n"; return 2; } }
    if (pos < argc) { if (!parse_bool(argv[pos++], IsMinimalBook)) { std::cerr << "ERROR: minimal must be bool\n"; return 2; } }
    if (pos < argc) {
      std::cerr << "ERROR: too many positional arguments.\n";
      print_usage(argv[0]);
      return 2;
    }
  }

  // -------- Sanity checks
  if (nthreads < 0) {
    std::cerr << "WARN: nthreads < 0; using 0 (auto)\n";
    nthreads = 0;
  }

  // -------- Show the resolved configuration
  std::cout << "------------------------------------------------------------\n";
  std::cout << " inputDir    : " << inputDir << "\n";
  std::cout << " nfile       : " << nfile << "\n";
  std::cout << " nthreads    : " << nthreads << " (0 => auto)\n";
  std::cout << " outputDir   : " << outputDir << "\n";
  std::cout << " dataconfig  : " << dataconfig << "\n";
  std::cout << " flags       : mc=" << (IsMC ? "1" : "0")
            << " reproc=" << (IsreprocRootFile ? "1" : "0")
            << " inbending=" << (IsInbending ? "1" : "0")
            << " minimal=" << (IsMinimalBook ? "1" : "0") << "\n";
  std::cout << "------------------------------------------------------------\n";

  // -------- Time the run
  auto start = std::chrono::high_resolution_clock::now();

  RunPhiAnalysis(inputDir, nfile, nthreads,
                 outputDir, dataconfig,
                 IsMC, IsreprocRootFile, IsInbending, IsMinimalBook);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;

  std::cout << "Total running time: " << elapsed.count() << " seconds\n";
  show_runtime_bar(elapsed.count());
  return 0;
}
