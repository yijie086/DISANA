#include "EventProcessor.h"

#include <iostream>

EventProcessor::EventProcessor(const std::string& inputDirectory, AnalysisTaskManager& taskMgr) : evt(inputDirectory), tasks(taskMgr) {}

void EventProcessor::ProcessEvents() {
  auto dfOpt = evt.getNode();
  if (!dfOpt.has_value()) {
    std::cerr << "EventProcessor: No input data found." << std::endl;
    return;
  }

  std::cout << "[EventProcessor] DataFrame initialized with " << evt.getFileCount() << " file(s)" << std::endl;

  ROOT::RDF::RNode df = dfOpt.value();

  tasks.UserCreateOutputObjects();
  tasks.Execute(df);
  tasks.SaveOutput();

  std::cout << "[EventProcessor] Finished processing all events." << std::endl;
}