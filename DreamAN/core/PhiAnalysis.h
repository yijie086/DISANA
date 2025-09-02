#ifndef PHIANALYSIS_H
#define PHIANALYSIS_H

#include <TH1F.h>

#include <ROOT/RDF/RInterface.hxx>
#include <optional>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include "../Cuts/EventCut.h"
#include "../Cuts/TrackCut.h"
#include "../Correction/MomentumCorrection.h"
#include "../Math/ParticleMassTable.h"
#include "../Math/RECParticleKinematic.h"
#include "../ParticleInformation/RECCalorimeter.h"
#include "../ParticleInformation/RECParticle.h"
#include "../ParticleInformation/RECTrack.h"
#include "../ParticleInformation/RECTraj.h"
#include "../ParticleInformation/RECForwardTagger.h"
#include "../core/Columns.h"
#include "AnalysisTask.h"

class PhiAnalysis : public AnalysisTask {
 public:
  PhiAnalysis(bool IsMC = false, bool IsReproc = false, bool IsMinBook = false);
  virtual ~PhiAnalysis();

  void UserCreateOutputObjects() override;
  void UserExec(ROOT::RDF::RNode &df) override;
  void SaveOutput() override;
  void SetMaxEvents(size_t n) { fMaxEvents = n; }
  // Bad for raw pointer setup
  void SetTrackCuts(std::shared_ptr<TrackCut> cuts) { fTrackCuts = std::move(cuts); };

  void SetEventCuts(EventCut *evtCuts) { fEventCuts = evtCuts; };
  void SetDoInvMassCut(bool cut) { fDoInvMassCut = cut; };
 
  void SetDoFiducialCut(bool cut) { fFiducialCut = cut; };

  void SetBeamEnergy(float beam_energy) { fbeam_energy = beam_energy; };

  void SetOutputDir(const std::string &dir) override;

  void SetFTonConfig(bool config) { fFTonConfig = config; }

  void SetDoMomentumCorrection(bool do_correction) { fDoMomentumCorrection = do_correction; }
  void SetMomentumCorrection(std::shared_ptr<MomentumCorrection> corr) { fMomCorr = std::move(corr); }



 private:
  bool IsMC = false;
  bool fDoInvMassCut = false;  // Flag to indicate if invMass cut is applied
  bool IsMinBooking = false;  // reduces the output to minimum only after fiducial 
  bool IsReproc = false;  // Flag to indicate if fiducial cut is applied
  bool fFiducialCut = false;  // Flag to indicate if fiducial cut is applied
  bool fFTonConfig = true;
  bool fDoMomentumCorrection = false;  // Flag to indicate if momentum correction is applied
  size_t fMaxEvents{0}; // Maximum number of events to process, 0 means no limit

  std::optional<ROOT::RDF::RNode> dforginal;

  std::optional<ROOT::RDF::RNode> dfSelected;
  std::optional<ROOT::RDF::RNode> dfSelected_afterFid;  // DataFrame after fiducial cuts
  std::optional<ROOT::RDF::RNode> dfSelected_afterFid_afterCorr;  // DataFrame after fiducial cuts and momentum correction
  std::string fOutputDir;
  
  float fbeam_energy = 10.6;
  
  TH1F *fHistPhotonP = nullptr;

  std::shared_ptr<TrackCut> fTrackCuts;
  EventCut *fEventCuts = nullptr;
  std::shared_ptr<TrackCut> fTrackCutsNoFid;
  std::shared_ptr<TrackCut> fTrackCutsWithFid;

  std::shared_ptr<MomentumCorrection> fMomCorr = nullptr;  // Pointer to momentum correction object
};

#endif
