#ifndef DVCSANALYSIS_H
#define DVCSANALYSIS_H

#include <TH1F.h>

#include <ROOT/RDF/RInterface.hxx>
#include <optional>

#include "../Cuts/EventCut.h"
#include "../Cuts/TrackCut.h"
#include "../Math/ParticleMassTable.h"
#include "../Math/RECParticleKinematic.h"
#include "../ParticleInformation/RECCalorimeter.h"
#include "../ParticleInformation/RECParticle.h"
#include "../ParticleInformation/RECTrack.h"
#include "../ParticleInformation/RECTraj.h"
#include "../ParticleInformation/RECForwardTagger.h"
#include "../core/Columns.h"
#include "./../Cuts/EventCut.h"
#include "./../Cuts/TrackCut.h"
#include "AnalysisTask.h"

class DVCSAnalysis : public AnalysisTask {
 public:
  DVCSAnalysis(bool IsMC = false, bool IsReproc = false);
  virtual ~DVCSAnalysis();

  void UserCreateOutputObjects() override;
  void UserExec(ROOT::RDF::RNode &df) override;
  void SaveOutput() override;
  // Bad for raw pointer setup
  void SetTrackCuts(std::shared_ptr<TrackCut> cuts) { fTrackCuts = std::move(cuts); };

  void SetPhotonCuts(EventCut *trkCuts) { fTrackCutsPhoton = trkCuts; };
  void SetElectronCuts(EventCut *trkCuts) { fTrackCutsElectron = trkCuts; };
  void SetProtonCuts(EventCut *trkCuts) { fTrackCutsProton = trkCuts; };
  void SetDoFiducialCut(bool cut) { fFiducialCut = cut; };
  void SetBeamEnergy(float beam_energy) { fbeam_energy = beam_energy; };

  void SetOutputFile(TFile *file) override;
  void SetOutputDir(const std::string &dir) override;

 private:
  bool IsMC = false;
  bool IsReproc = false;  // Flag to indicate if fiducial cut is applied
  bool fFiducialCut = false;  // Flag to indicate if fiducial cut is applied
  std::optional<ROOT::RDF::RNode> dfSelected;
  std::optional<ROOT::RDF::RNode> dfSelected_after;  // DataFrame after fiducial cuts
  std::string fOutputDir;
  float fbeam_energy = 10.6;
  TH1F *fHistPhotonP = nullptr;
  std::shared_ptr<TrackCut> fTrackCuts;
  std::shared_ptr<TrackCut> fTrackCutsNoFid;
  std::shared_ptr<TrackCut> fTrackCutsWithFid;
  EventCut *fTrackCutsPhoton = nullptr;
  EventCut *fTrackCutsElectron = nullptr;
  EventCut *fTrackCutsProton = nullptr;

  TFile *fOutFile = nullptr;  // Output file pointer set by manager
};

#endif
