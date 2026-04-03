cmake_minimum_required(VERSION 3.10)
project(AnalysisFramework)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# --- Environment checks ---
foreach(VAR CLAS12ROOT HIPO ROOTSYS QADB)
    if(NOT DEFINED ENV{${VAR}})
        message(FATAL_ERROR "${VAR} environment variable is not set!")
    endif()
    message(STATUS "${VAR} = $ENV{${VAR}}")
endforeach()

# --- HIPO path (priority: -DHIPO_DIR=... > $HIPO env var > hardcoded fallback) ---
if(DEFINED HIPO_DIR AND NOT HIPO_DIR STREQUAL "")
    message(STATUS "Using HIPO_DIR from command line: ${HIPO_DIR}")
elseif(DEFINED ENV{HIPO} AND NOT "$ENV{HIPO}" STREQUAL "")
    set(HIPO_DIR "$ENV{HIPO}")
    message(STATUS "Using HIPO from environment: ${HIPO_DIR}")
else()
    set(HIPO_DIR "/u/home/singh/Softwares/hipo/install")
    message(STATUS "Using HIPO fallback: ${HIPO_DIR}")
endif()

# --- ROOT config ---
execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND root-config --libs   OUTPUT_VARIABLE ROOT_LIBS   OUTPUT_STRIP_TRAILING_WHITESPACE)

# --- Include directories ---
include_directories($ENV{CLAS12ROOT}/Clas12Root)
include_directories($ENV{CLAS12ROOT}/Clas12Banks)
include_directories($ENV{ROOTSYS}/include)
include_directories(${HIPO_DIR}/include)
include_directories(${HIPO_DIR}/include/hipo4)
include_directories($ENV{QADB}/srcC/include)
include_directories($ENV{QADB}/srcC/rapidjson/include)

# --- Find CLAS12ROOT libraries (lib64 for cvmfs installs) ---
find_library(CLAS12ROOT_LIB
    NAMES Clas12Root
    HINTS $ENV{CLAS12ROOT}/lib $ENV{CLAS12ROOT}/lib64
    REQUIRED
)
find_library(CLAS12BANKS_LIB
    NAMES Clas12Banks
    HINTS $ENV{CLAS12ROOT}/lib $ENV{CLAS12ROOT}/lib64
    REQUIRED
)
find_library(HIPO4_LIB
    NAMES hipo4
    HINTS ${HIPO_DIR}/lib ${HIPO_DIR}/lib64
    REQUIRED
)
find_library(HIPODATAFRAME_LIB
    NAMES HipoDataFrame
    HINTS ${HIPO_DIR}/lib ${HIPO_DIR}/lib64
    REQUIRED
)

message(STATUS "Clas12Root    lib: ${CLAS12ROOT_LIB}")
message(STATUS "Clas12Banks   lib: ${CLAS12BANKS_LIB}")
message(STATUS "hipo4         lib: ${HIPO4_LIB}")
message(STATUS "HipoDataFrame lib: ${HIPODATAFRAME_LIB}")
message(STATUS "ROOT Libraries:    ${ROOT_LIBS}")

# --- Shared source files ---
set(COMMON_SOURCES
    DreamAN/core/AnalysisTask.cxx
    DreamAN/core/EventProcessor.cxx
    DreamAN/core/AnalysisTaskManager.cxx
    DreamAN/core/Events.cxx
    DreamAN/core/Columns.cxx
    DreamAN/ParticleInformation/RECParticle.cxx
    DreamAN/ParticleInformation/RECTraj.cxx
    DreamAN/ParticleInformation/RECTrack.cxx
    DreamAN/ParticleInformation/RECCalorimeter.cxx
    DreamAN/ParticleInformation/RECForwardTagger.cxx
    DreamAN/Cuts/EventCut.cxx
    DreamAN/Cuts/TrackCut.cxx
    DreamAN/Cuts/QADBCuts.cxx
    DreamAN/Correction/MomentumCorrection.cxx
    DreamAN/Math/RECParticleKinematic.cxx
    DreamAN/Math/MathKinematicVariable.cxx
    DreamAN/Math/ParticleMassTable.cxx
)

set(COMMON_LIBS
    ${ROOT_LIBS}
    ${CLAS12ROOT_LIB}
    ${CLAS12BANKS_LIB}
    ${HIPO4_LIB}
    ${HIPODATAFRAME_LIB}
    pthread
)

# --- AnalysisDVCS ---
add_executable(AnalysisDVCS
    macros/mainDVCS.C
    macros/RunDVCSAnalysis.C
    ${COMMON_SOURCES}
    DreamAN/core/DVCSAnalysis.cxx
)
target_link_libraries(AnalysisDVCS PRIVATE ${COMMON_LIBS})

# --- AnalysisPhi ---
add_executable(AnalysisPhi
    macros/mainPhi.C
    macros/RunPhiAnalysis.C
    ${COMMON_SOURCES}
    DreamAN/core/PhiAnalysis.cxx
)
target_link_libraries(AnalysisPhi PRIVATE ${COMMON_LIBS})