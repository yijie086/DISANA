#ifndef EVENTPROCESSOR_H
#define EVENTPROCESSOR_H

#include <string>
#include "AnalysisTaskManager.h"
#include "Events.h"

class EventProcessor {
public:
    EventProcessor(AnalysisTaskManager& taskMgr,const std::string& inputDirectory, const std::string& OuptpuDirectory, bool fIsReprocessRootFile, const std::string& fInputROOTtreeName, const std::string& fInputROOTfileName, int nfiles, const int nthreads);
    void ProcessEvents();

private:
    Events evt;
    AnalysisTaskManager& tasks;
};

#endif