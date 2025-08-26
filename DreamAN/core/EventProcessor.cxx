#include "EventProcessor.h"

#include <iostream>

EventProcessor::EventProcessor(AnalysisTaskManager& taskMgr, const std::string& inputDirectory,const std::string& OuptpuDirectory, bool fIsReprocessRootFile, const std::string& fInputROOTtreeName, const std::string& fInputROOTfileName, int nfiles, const int nthreads ) : evt(inputDirectory, OuptpuDirectory,fIsReprocessRootFile, fInputROOTtreeName, fInputROOTfileName, nfiles, nthreads), tasks(taskMgr) {}

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