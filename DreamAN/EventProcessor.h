#ifndef EVENTPROCESSOR_H
#define EVENTPROCESSOR_H

#include "Event.h"
#include "AnalysisTaskManager.h"

class EventProcessor {
public:
    EventProcessor(const std::string& inputFile, AnalysisTaskManager& taskManager);
    void ProcessEvents();

private:
    Event evt;
    AnalysisTaskManager& tasks; // Reference to TaskManager
};

#endif // EVENTPROCESSOR_H
