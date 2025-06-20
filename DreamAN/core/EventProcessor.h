#ifndef EVENTPROCESSOR_H
#define EVENTPROCESSOR_H

#include <string>
#include "AnalysisTaskManager.h"
#include "Events.h"

class EventProcessor {
public:
    EventProcessor(AnalysisTaskManager& taskMgr,const std::string& inputDirectory, bool fIsReprocessRootFile, const std::string& fInputROOTtreeName, const std::string& fInputROOTfileName, int nfiles);
    void ProcessEvents();

private:
    Events evt;
    AnalysisTaskManager& tasks;
};

#endif