// Run with: root -l 'lund2rdf.cpp("dvcsgen*.dat","out.root")'
#include <TFile.h>
#include <TTree.h>
#include <glob.h>  // POSIX glob for pattern matching

#include <ROOT/RDataFrame.hxx>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Split a line by whitespace
static std::vector<std::string> split_ws(const std::string& line) {
  std::istringstream iss(line);
  std::vector<std::string> toks;
  std::string t;
  while (iss >> t) toks.push_back(t);
  return toks;
}

// Trim leading/trailing spaces
static std::string trim(const std::string& s) {
  size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  size_t b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

// Return basename without path.


static std::string basename_of(const std::string& path) {
    size_t pos = path.find_last_of("/\\");
    return (pos == std::string::npos) ? path : path.substr(pos + 1);
}

// Expand a glob pattern into a list of file names
static std::vector<std::string> glob_files(const char* pattern) {
  glob_t results;
  std::vector<std::string> files;
  if (glob(pattern, 0, nullptr, &results) == 0) {
    for (size_t i = 0; i < results.gl_pathc; ++i) {
      files.emplace_back(results.gl_pathv[i]);
    }
  }
  globfree(&results);
  return files;
}

// Main function: read all files matching pattern, write a single ROOT file
void lund2rdf(const char* inPattern = "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/dvcsgen/script_no_rad/*.dat",
              const char* outPath = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/aaogen/no_rad_gen/no_rad_gen_10312025.root") {

                //const char* inPattern = "/w/hallb-scshelf2102/clas12/singh/Softwares/Generators/dvcsgen/script_rad/*.dat",
             // const char* outPath = "/w/hallb-scshelf2102/clas12/singh/Softwares/DISANA_main/data_processed/sims/aaogen/rad_gen/rad100M_new.root"
  bool rga_sp2018 = false;
  // Prepare output ROOT file
  TFile fout(outPath, "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "ERROR: cannot create output file: " << outPath << std::endl;
    return;
  }
  std::unordered_set<std::string> exclude_names = {};
  // Exact basenames to exclude:
  if (rga_sp2018) {
    exclude_names = {"nor095M1.dat", "nor088M1.dat", "nor015M1.dat", "nor022M1.dat", "nor028M1.dat", "nor033M1.dat", "nor036M1.dat",
                     "nor037M1.dat", "nor040M1.dat", "nor042M1.dat", "nor047M1.dat", "nor050M1.dat", "nor059M1.dat", "nor065M1.dat",
                     "nor068M1.dat", "nor069M1.dat", "nor075M1.dat", "nor078M1.dat", "nor080M1.dat", "nor092M1.dat"};
  } else {
    exclude_names = {""};
  }

  // Branch containers
  std::vector<int> v_pid;
  std::vector<float> v_px, v_py, v_pz;
  bool hasphoto = false;

  TTree tree("MC", "MC particles from LUND (dvcsgen)");
  tree.Branch("MC_Particle_pid", &v_pid);
  tree.Branch("MC_Particle_px", &v_px);
  tree.Branch("MC_Particle_py", &v_py);
  tree.Branch("MC_Particle_pz", &v_pz);

  long long evCount = 0;

  // Expand glob pattern
  auto files = glob_files(inPattern);
  if (files.empty()) {
    std::cerr << "No files match pattern: " << inPattern << std::endl;
    return;
  }

  // Loop over files
  for (auto& fname : files) {
    std::ifstream fin(fname);

    if (!fin.is_open()) {
      std::cerr << "WARNING: cannot open file: " << fname << std::endl;
      continue;
    }
    const std::string base = basename_of(fname);
    if (exclude_names.count(base)) {
      std::cout << "Skipping " << fname << " (in exclude list)\n";
      continue;
    }
    std::cout << "Processing " << fname << " ..." << std::endl;

    std::string line;
    while (std::getline(fin, line)) {
      line = trim(line);
      if (line.empty()) continue;

      auto toks = split_ws(line);
      if (toks.empty()) continue;

      int nPart = -1;
      try {
        nPart = std::stoi(toks[0]);  // Npart
      } catch (...) {
        continue;  // not an event header
      }
      if (nPart <= 0) continue;

      v_pid.clear();
      v_px.clear();
      v_py.clear();
      v_pz.clear();
      v_pid.reserve(nPart);
      v_px.reserve(nPart);
      v_py.reserve(nPart);
      v_pz.reserve(nPart);

      for (int i = 0; i < nPart; ++i) {
        std::string pline;
        if (!std::getline(fin, pline)) break;
        pline = trim(pline);
        if (pline.empty()) {
          --i;
          continue;
        }

        auto ptoks = split_ws(pline);
        if (ptoks.size() < 9) continue;

        try {
          int pid = std::stoi(ptoks[3]);
          float px = std::stof(ptoks[6]);
          float py = std::stof(ptoks[7]);
          float pz = std::stof(ptoks[8]);
          if (pid == 22 && hasphoto) continue;  // only keep first photon
          if (pid == 22) hasphoto = true;

          v_pid.push_back(pid);
          v_px.push_back(px);
          v_py.push_back(py);
          v_pz.push_back(pz);
        } catch (...) {
          // skip bad line
        }
      }
      hasphoto = false;
      tree.Fill();
      ++evCount;
    }
  }

  fout.Write();
  fout.Close();

  std::cout << "Done. Wrote " << evCount << " events to " << outPath << std::endl;
}
