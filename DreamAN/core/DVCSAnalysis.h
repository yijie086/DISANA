#ifndef DVCSANALYSIS_H
#define DVCSANALYSIS_H

#include <TH1F.h>
#include <ROOT/RDF/RInterface.hxx>
#include <optional> 
#include "./../Cuts/EventCut.h"
#include "./../Cuts/TrackCut.h"
#include "../ParticleInformation/RECParticle.h"
#include "../ParticleInformation/RECTraj.h"
#include "../ParticleInformation/RECCalorimeter.h"
#include "../Cuts/ElectronCut.h"
#include "../DrawHist/DrawAndSave.h"
#include "../core/FilesInPath.h"
#include "../core/Columns.h"
#include "../Cuts/EventCut.h"
#include "../Cuts/TrackCut.h"
#include "../Math/RECParticleKinematic.h"
#include "../Math/ParticleMassTable.h"
#include "AnalysisTask.h"

class DVCSAnalysis : public AnalysisTask {
 public:
  DVCSAnalysis();
  virtual ~DVCSAnalysis();

  void UserCreateOutputObjects() override;
  void UserExec(ROOT::RDF::RNode &df) override;
  void SaveOutput() override;
  void SetTrackCuts(TrackCut *trkCuts) { fTrackCuts = trkCuts; };
  void SetPhotonCuts(EventCut *trkCuts) { fTrackCutsPhoton = trkCuts; };
  void SetElectronCuts(EventCut *trkCuts) { fTrackCutsElectron = trkCuts; };
  void SetProtonCuts(EventCut *trkCuts) { fTrackCutsProton = trkCuts; };
  void SetBeamEnergy(float beam_energy) { fbeam_energy = beam_energy; };


  void SetOutputFile(TFile* file) override;
  void SetOutputDir(const std::string& dir) override;


 private:
  std::optional<ROOT::RDF::RNode> dfSelected;
  std::string fOutputDir;
  float fbeam_energy = 10.6;
  TH1F *fHistPhotonP = nullptr;
  TrackCut *fTrackCuts = nullptr;
  EventCut *fTrackCutsPhoton = nullptr;
  EventCut *fTrackCutsElectron = nullptr;
  EventCut *fTrackCutsProton = nullptr;

  TFile* fOutFile = nullptr; // Output file pointer set by manager
};

#endif
