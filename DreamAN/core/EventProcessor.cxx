#include "EventProcessor.h"

#include <iostream>

EventProcessor::EventProcessor(AnalysisTaskManager& taskMgr, const std::string& inputDirectory, bool fIsReprocessRootFile, const std::string& fInputROOTtreeName, const std::string& fInputROOTfileName, int nfiles ) : evt(inputDirectory,fIsReprocessRootFile, fInputROOTtreeName, fInputROOTfileName, nfiles), tasks(taskMgr) {}

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