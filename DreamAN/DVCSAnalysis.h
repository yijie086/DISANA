#ifndef DVCSANALYSIS_H
#define DVCSANALYSIS_H
// this task related includes

#include "AnalysisTask.h"
#include "Event.h"
#include "Tracks.h"
#include "EventCuts.h"
#include "TrackCuts.h"
#include <TObject.h>
#include <iostream>
#include "AnalysisTaskManager.h"  // Add the include for the complete definition


class DVCSAnalysis : public AnalysisTask {
public:
  DVCSAnalysis(); // Declare the constructor
  DVCSAnalysis(bool isMC); // Declare the constructor
  virtual ~DVCSAnalysis();
  virtual void UserCreateOutputObjects() ;//override;
  void UserExec(Event &event) ;//override;
  virtual void Terminate(Option_t *) {};

  void SetEventCuts(EventCuts *evtCuts) { fEventCuts = evtCuts; };
  void SetPhotonCuts(TrackCuts *trkCuts) {fTrackCutsPhoton = trkCuts;};
  void SetElectronCuts(TrackCuts *trkCuts) {fTrackCutsElectron = trkCuts;};
  void SetProtonCuts(TrackCuts *trkCuts) {fTrackCutsProton = trkCuts;};


private:
  DVCSAnalysis(const DVCSAnalysis &task);
  DVCSAnalysis &operator=(const DVCSAnalysis &task);
  bool fisLightWeight;  
  bool fIsMC;           
  Event *fEvent;        
  Tracks *fTrack;               
  EventCuts *fEventCuts; 
  TrackCuts *fTrackCutsPhoton; 
  TrackCuts *fTrackCutsElectron; 
  TrackCuts *fTrackCutsProton; 

  TH1F* fHistMomentum;
  TTree* fTree;
  float fMomentum;

  void FillTree();
  void FillHistogram();
};
#endif // DVCSANALYSIS_H
