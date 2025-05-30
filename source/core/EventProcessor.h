#ifndef EVENTPROCESSOR_H
#define EVENTPROCESSOR_H

#include <string>
#include "AnalysisTaskManager.h"
#include "Events.h"

class EventProcessor {
public:
    EventProcessor(const std::string& inputDirectory, AnalysisTaskManager& taskMgr);
    void ProcessEvents();

private:
    Events evt;
    AnalysisTaskManager& tasks;
};

#endif