#include "EventProcessor.h"
#include <iostream>

EventProcessor::EventProcessor(const std::string& inputFile, AnalysisTaskManager& taskManager)
    : evt(inputFile), tasks(taskManager) {}

void EventProcessor::ProcessEvents() {
        if (!evt.isOpen()) {
            std::cerr << "Error: Could not open file." << std::endl;
            return;
        }
    
    // Initialize output objects by calling UserCreateOutputObjects()
     tasks.UserCreateOutputObjects();
        int count = 0;
        while (evt.getNextEvent() && count < 2) {
            if (evt.isValid()) {
                tasks.ExecuteTasks(evt); // Calls all registered tasks
            } else {
                std::cerr << "No valid event retrieved." << std::endl;
            }
            ++count;
        }
        // Save all the output at the end
        tasks.SaveOutput();
        std::cout << "Finished processing all events." << std::endl;
}
    

