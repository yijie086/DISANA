cmake_minimum_required(VERSION 3.10)
project(AnalysisFramework)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Check if CLAS12ROOT environment is set
if(NOT DEFINED ENV{CLAS12ROOT})
    message(FATAL_ERROR "CLAS12ROOT environment variable not set!")
endif()

# Set the HIPO4 path
# set(HIPO4 "/u/scigroup/cvmfs/hallb/clas12/sw/almalinux9-gcc11/local/hipo/4.2.0")

# Include ROOT headers using the flags from root-config (if necessary)
execute_process(
    COMMAND root-config --cflags
    OUTPUT_VARIABLE ROOT_CFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# You can also use the ROOT library via 'root-config'
execute_process(
    COMMAND root-config --libs
    OUTPUT_VARIABLE ROOT_LIBS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Include ROOT headers
include_directories(${ROOT_CFLAGS})

# Include CLAS12ROOT headers
include_directories($ENV{CLAS12ROOT}/Clas12Root)
include_directories($ENV{CLAS12ROOT}/Clas12Banks)
include_directories($ENV{ROOTSYS}/include)
include_directories($ENV{HIPO}/include/hipo4)


# Include HIPO4 headers
include_directories($ENV{HIPO}/include)

# Link the libraries
link_directories($ENV{CLAS12ROOT}/lib)  # Make sure this is correct
link_directories($ENV{HIPO}/lib)


# Add the path to Clas12Root and Clas12Banks libraries
link_directories($ENV{CLAS12ROOT}/lib)

# Add your source files bhawani's analysis
add_executable(AnalysisDVCS
    macros/mainDVCS.C
    macros/RunDVCSAnalysis.C

    # main core methods and classes
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
    DreamAN/Correction/MomentumCorrection.cxx
    DreamAN/Math/RECParticleKinematic.cxx
    DreamAN/Math/MathKinematicVariable.cxx
    DreamAN/Math/ParticleMassTable.cxx

    #analysis related classes
    DreamAN/core/DVCSAnalysis.cxx
)

# Link against ROOT libraries, CLAS12ROOT libraries, and HIPO4
target_link_libraries(AnalysisDVCS
    ${ROOT_LIBS}
    pthread
    Clas12Root  # Link against the precompiled Clas12Root library
    Clas12Banks
    hipo4
    HipoDataFrame
)


#PhiAanylsis
add_executable(AnalysisPhi
    macros/mainPhi.C
    macros/RunPhiAnalysis.C

    # main core methods and classes
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
    DreamAN/Correction/MomentumCorrection.cxx
    DreamAN/Math/RECParticleKinematic.cxx
    DreamAN/Math/MathKinematicVariable.cxx
    DreamAN/Math/ParticleMassTable.cxx

    #analysis related classes
    DreamAN/core/PhiAnalysis.cxx
)

# Link against ROOT libraries, CLAS12ROOT libraries, and HIPO4
target_link_libraries(AnalysisPhi
    PRIVATE
        ${ROOT_LIBS}
        Clas12Root
        Clas12Banks
        hipo4
        HipoDataFrame
        pthread
)



# Debugging info (optional)
message(STATUS "ROOT Libraries: ${ROOT_LIBS}")
