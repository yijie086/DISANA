cmake_minimum_required(VERSION 3.10)
project(AnalysisFramework LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -------- Environment checks --------
if(NOT DEFINED ENV{CLAS12ROOT})
  message(FATAL_ERROR "CLAS12ROOT environment variable not set!")
endif()

if(NOT DEFINED ENV{HIPO})
  message(FATAL_ERROR "HIPO environment variable not set!")
endif()

# -------- ROOT via root-config --------
execute_process(
  COMMAND root-config --cflags
  OUTPUT_VARIABLE ROOT_CFLAGS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
  COMMAND root-config --libs
  OUTPUT_VARIABLE ROOT_LIBS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

separate_arguments(ROOT_CFLAGS_LIST NATIVE_COMMAND "${ROOT_CFLAGS}")
separate_arguments(ROOT_LIBS_LIST   NATIVE_COMMAND "${ROOT_LIBS}")

message(STATUS "ROOT_CFLAGS: ${ROOT_CFLAGS}")
message(STATUS "ROOT_LIBS:   ${ROOT_LIBS}")
message(STATUS "CLAS12ROOT:  $ENV{CLAS12ROOT}")
message(STATUS "HIPO:       $ENV{HIPO}")

# -------- Include dirs (ONLY directories!) --------
include_directories(
  $ENV{CLAS12ROOT}/Clas12Root
  $ENV{CLAS12ROOT}/Clas12Banks
  $ENV{ROOTSYS}/include
  $ENV{HIPO}/include
  $ENV{HIPO}/include/hipo4
  $ENV{QADB}/srcC/include
  $ENV{QADB}/srcC/rapidjson/include
)

# -------- Find CLAS12ROOT libs robustly --------
set(CLAS12_LIBDIR "$ENV{CLAS12ROOT}/lib64")

find_library(CLAS12ROOT_LIB  NAMES Clas12Root  PATHS "${CLAS12_LIBDIR}" NO_DEFAULT_PATH)
find_library(CLAS12BANKS_LIB NAMES Clas12Banks PATHS "${CLAS12_LIBDIR}" NO_DEFAULT_PATH)

if(NOT CLAS12ROOT_LIB)
  message(FATAL_ERROR "Could not find libClas12Root in ${CLAS12_LIBDIR}")
endif()
if(NOT CLAS12BANKS_LIB)
  message(FATAL_ERROR "Could not find libClas12Banks in ${CLAS12_LIBDIR}")
endif()

message(STATUS "Found Clas12Root:  ${CLAS12ROOT_LIB}")
message(STATUS "Found Clas12Banks: ${CLAS12BANKS_LIB}")

# -------- Find HIPO libs robustly (CVMFS or local) --------
set(HIPO_PREFIX "$ENV{HIPO}")

# try lib64 first, then lib
set(HIPO_LIBDIR "${HIPO_PREFIX}/lib64")
if(NOT EXISTS "${HIPO_LIBDIR}")
  set(HIPO_LIBDIR "${HIPO_PREFIX}/lib")
endif()

find_library(HIPO4_LIB   NAMES hipo4         PATHS "${HIPO_LIBDIR}" NO_DEFAULT_PATH)
find_library(HIPODF_LIB  NAMES HipoDataFrame PATHS "${HIPO_LIBDIR}" NO_DEFAULT_PATH)

if(NOT HIPO4_LIB)
  message(FATAL_ERROR "Could not find libhipo4 in ${HIPO_LIBDIR} (HIPO_PREFIX=${HIPO_PREFIX})")
endif()
if(NOT HIPODF_LIB)
  message(FATAL_ERROR "Could not find libHipoDataFrame in ${HIPO_LIBDIR} (HIPO_PREFIX=${HIPO_PREFIX})")
endif()

message(STATUS "HIPO_LIBDIR: ${HIPO_LIBDIR}")
message(STATUS "Found hipo4:         ${HIPO4_LIB}")
message(STATUS "Found HipoDataFrame: ${HIPODF_LIB}")

# -------- RPATH so executables run without extra env vars --------
set(CMAKE_BUILD_RPATH "${CLAS12_LIBDIR};${HIPO_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# -------- AnalysisDVCS --------
add_executable(AnalysisDVCS
  macros/mainDVCS.C
  macros/RunDVCSAnalysis.C

  DreamAN/core/AnalysisTask.cxx
  DreamAN/core/EventProcessor.cxx
  DreamAN/core/AnalysisTaskManager.cxx
  DreamAN/core/HipoToRootConverter.cxx
  DreamAN/core/Events.cxx
  DreamAN/ParticleInformation/RECParticle.cxx
  DreamAN/ParticleInformation/RECTraj.cxx
  DreamAN/ParticleInformation/RECTrack.cxx
  DreamAN/ParticleInformation/RECCalorimeter.cxx
  DreamAN/ParticleInformation/RECForwardTagger.cxx
  DreamAN/core/Columns.cxx
  DreamAN/Cuts/EventCut.cxx
  DreamAN/Cuts/TrackCut.cxx
  DreamAN/Cuts/QADBCuts.cxx
  DreamAN/Correction/MomentumCorrection.cxx
  DreamAN/Math/RECParticleKinematic.cxx
  DreamAN/Math/MathKinematicVariable.cxx
  DreamAN/Math/ParticleMassTable.cxx
  DreamAN/core/DVCSAnalysis.cxx
)

target_compile_options(AnalysisDVCS PRIVATE ${ROOT_CFLAGS_LIST})

target_link_libraries(AnalysisDVCS PRIVATE
  ${ROOT_LIBS_LIST}
  pthread
  ${CLAS12ROOT_LIB}
  ${CLAS12BANKS_LIB}
  ${HIPO4_LIB}
  ${HIPODF_LIB}
)

# -------- AnalysisPhi --------
add_executable(AnalysisPhi
  macros/mainPhi.C
  macros/RunPhiAnalysis.C

  DreamAN/core/AnalysisTask.cxx
  DreamAN/core/EventProcessor.cxx
  DreamAN/core/AnalysisTaskManager.cxx
  DreamAN/core/HipoToRootConverter.cxx
  DreamAN/core/Events.cxx
  DreamAN/ParticleInformation/RECParticle.cxx
  DreamAN/ParticleInformation/RECTraj.cxx
  DreamAN/ParticleInformation/RECTrack.cxx
  DreamAN/ParticleInformation/RECCalorimeter.cxx
  DreamAN/ParticleInformation/RECForwardTagger.cxx
  DreamAN/core/Columns.cxx
  DreamAN/Cuts/EventCut.cxx
  DreamAN/Cuts/TrackCut.cxx
  DreamAN/Cuts/QADBCuts.cxx
  DreamAN/Correction/MomentumCorrection.cxx
  DreamAN/Math/RECParticleKinematic.cxx
  DreamAN/Math/MathKinematicVariable.cxx
  DreamAN/Math/ParticleMassTable.cxx
  DreamAN/core/PhiAnalysis.cxx
)

target_compile_options(AnalysisPhi PRIVATE ${ROOT_CFLAGS_LIST})

target_link_libraries(AnalysisPhi PRIVATE
  ${ROOT_LIBS_LIST}
  pthread
  ${CLAS12ROOT_LIB}
  ${CLAS12BANKS_LIB}
  ${HIPO4_LIB}
  ${HIPODF_LIB}
)