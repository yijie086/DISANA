# **DISANA**

---

## **Directory Structure**

```
DISANA/
├── macros/
│   ├── Run101.C
│   ├── main101.C
├── source/
│   ├── Cuts/
│   │   ├── EventCut.h
│   │   ├── EventCut.cxx
│   ├── ParticleInformation/
│   │   ├── RECParticle.h
│   │   ├── RECParticle.cxx
│   ├── Math/
│   │   ├── RECParticleKinematic.h
│   │   ├── RECParticleKinematic.cxx
│   ├── core/
│   │   ├── FilesInPath.h
│   │   ├── FilesInPath.cxx
│   ├── DrawHist/
│       ├── DrawAndSave.h
│       ├── DrawAndSave.cxx
```

---

## **About `macros`**

### **1. macros/Run101.C**

* **Functionality**:

  * This is the main analysis script. It loads HIPO data files, applies selection cuts, computes particle kinematic variables, generates histograms, and saves results.
* **Main Workflow**:

  1. Calls `GetHipoFilesInPath` to retrieve all `.hipo` files in a specified directory (including recursively within subdirectories).
  2. Loads data using `ROOT::RDataFrame`.
  3. Defines selection criteria for particles (e.g., electrons, protons, and photons; customizable, see `EventCut.cxx` and `EventCut.h`).
  4. Computes particle kinematic variables (e.g., `theta`, `phi`, `p`) and attaches them to the `ROOT::RDF::RNode`.
  5. Calls functions in `DrawAndSave.cxx` to generate and save histograms in `.png` and `.root` formats.
  6. Saves the filtered data as a snapshot into a ROOT file.

---

### **2. macros/main101.C**

* **Functionality**:

  * This is the program’s entry point. It parses command-line arguments and calls the `Run101` function.
* **Main Workflow**:

  1. Checks if a path to the HIPO data directory is provided via command-line arguments.
  2. Calls the `Run101` function to perform data analysis.

---

## **About `source`**

### **1. Cuts/EventCut.h and Cuts/EventCut.cxx**

* **Functionality**:

  * Defines a general event selector `EventCut`.
  * Used to filter events based on particle charge, momentum, vertex position, etc.
* **Main Workflow**:

  1. Provides functions to configure cut criteria (e.g., `SetChargeCut`, `SetMomentumCut`, etc.).
  2. Fully compatible with `ROOT::RDataFrame`.
* **Dependencies**:

  * No direct external dependencies, but works with the internal data structure.

---

### **2. ParticleInformation/RECParticle.h and RECParticle.cxx**

* **Functionality**:

  * Provides access to particle information, including momentum components, vertex position, charge, etc.
  * Includes utility functions to filter specific particles and extract their properties.
* **Main Workflow**:

  1. Defines data structures for particles (e.g., `RECParticle::All()`, `RECParticle::Extend()`).
  2. Provides `get_RECParticle_float_var` and `get_RECParticle_int_var` for extracting specific attributes.
  3. `RECParticle::Extend()` defines extended variables like `theta`, `phi`, and `p`.

---

### **3. Math/RECParticleKinematic.h and Math/RECParticleKinematic.cxx**

* **Functionality**:

  * Provides functions to compute particle kinematic variables such as `theta`, `phi`, and `p`.
* **Main Workflow**:

  1. Defines `RECParticletheta`: computes the polar angle.
  2. Defines `RECParticlephi`: computes the azimuthal angle.
  3. Defines `RECParticleP`: computes the momentum magnitude.

---

### **4. core/FilesInPath.h and core/FilesInPath.cxx**

* **Functionality**:

  * Provides the utility function `GetHipoFilesInPath` to recursively search for `.hipo` files in a given directory.
* **Main Workflow**:

  1. Uses C++17 `std::filesystem` to traverse directories.
  2. Filters files with the `.hipo` extension.
  3. Returns a list of valid file paths.
* **Dependencies**:

  * Requires C++17 standard library.

---

### **5. DrawHist/DrawAndSave.h and DrawHist/DrawAndSave.cxx**

* **Functionality**:

  * Provides utility functions to draw and save histograms.
* **Main Workflow**:

  1. Defines `DrawAndSaveParticleHistograms` to generate `theta`, `phi`, and `p` histograms for particles.
  2. Dynamically generates histogram names and titles.
  3. Uses `Draw1DHist` and `Save1DHist` to draw and save histograms as `.png` and `.root`.

---

## **Usage**

### **1. Build the Project**

On Jefferson Lab ifarm, compile the project by executing:

```bash
module load clas12
```

```bash
git clone https://github.com/yijie086/DISANA
```

```bash
cd DISANA
mkdir build
cd build
cmake ..
make
```

Or with multithreading:

```bash
make -j8
```

### **2. Run the Program**

Run the compiled executable `main101` and specify the directory containing `.hipo` files:

```bash
./main101 /path/to/hipo/files
```

### **3. Output**

* Filtered data is saved to `snapshot.root`.
* Histograms are saved to `test.root`.
* Histogram images will be exported as `.png` files.

---

### **Notes**

1. To modify selection criteria, you don’t need to edit `EventCut` directly. Refer to how it is used in `Run101.C`.
2. To add new particle properties or kinematic variables, extend `RECParticle` and `RECParticleKinematic` accordingly.
