#include "FirstAnalysis.h"
#include "DVCSAnalysis.h"
#include "EventProcessor.h"
#include "AnalysisTaskManager.h"
#include <iostream>

void RunDVCSAnalysis(const std::string& inputFile) {
    // Check if the correct number of arguments were passed
    if (inputFile== "") {
        std::cerr << "Usage: " << inputFile << " <input_file>" << std::endl;
        return ;
    }

    DVCSAnalysis userTask; // Create user-defined task

    AnalysisTaskManager taskManager;
    taskManager.AddTask(userTask);
    //taskManager.AddTask(new UserTask2());

    // Run event processing
    EventProcessor processor(inputFile, taskManager);
    processor.ProcessEvents(); // Start event loop
    return;
}
