#include "TrackCuts.h"
// Default constructor
TrackCuts::TrackCuts()
    : // fMCHists(nullptr),
      fMinimalBooking(false), fHists(nullptr), fPmin(0.0), fPmax(100.0),
      fcutP(0), fVtXmin(-100.0), fVtXmax(100.0), fVtYmin(-100.0),
      fVtYmax(100.0), fVtZmin(-100.0), fVtZmax(100.0), fCharge(0),
      fChiSquared(0), fcutPID(0) {}

TrackCuts::TrackCuts(const TrackCuts &cuts)
    : // fMCHists(cuts.fMCHists),
      fMinimalBooking(cuts.fMinimalBooking), fHists(cuts.fHists),
      fPmin(cuts.fPmin), fPmax(cuts.fPmax), fcutP(cuts.fcutP),
      fcutPID(cuts.fcutPID), fCharge(cuts.fCharge),
      fChiSquared(cuts.fChiSquared) {}

TrackCuts &TrackCuts::operator=(const TrackCuts &cuts) {
  if (this == &cuts) {
    return *this;
  }
  {
    // this->fMCHists = cuts.fMCHists;
    this->fMinimalBooking = cuts.fMinimalBooking;
    this->fHists = cuts.fHists;
    this->fPmin = cuts.fPmin;
    this->fPmax = cuts.fPmax;
    this->fcutP = cuts.fcutP;
    this->fcutPID = cuts.fcutPID;
    this->fPdg = cuts.fPdg;
    this->fChiSquared = cuts.fChiSquared;
  }
  return *this; 
}

TrackCuts::~TrackCuts() {
  // if (fMCHists) {
  // delete fMCHists;
  //}
  if (fHists) {
    delete fHists;
  } else {
    // if (!fMinimalBooking)
    //  fHists->FillTrackCounter(0);
  }
}

bool TrackCuts::isSelected(Tracks *track) {
  if (!track) {
    std::cerr << "No Input Track received" << std::endl;
  }
  bool pass = true;
  if (pass && fcutP) {
    if (!PassesMomentumCut(track)) {
      pass = false;
    }
  }
  if (pass && fcutPID) {
    if (!PIDCuts(track)) {
      pass = false;
    } else {
      // we may wanna have some histograms related to cuts here
    }
    track->SetUse(pass); //This is actaully important since we need to set the

    // check for PDG code
    if (!(pass && PassesPIDCut(track))) {
      pass = false;
    }
  }
  return pass;
}

bool TrackCuts::TrackingCuts(Tracks *track) {
  bool pass = true;
  // detector track related cuts
  bool DetStatus = true; // should be accessed from the particle banks
  if (!DetStatus) {
    pass = false;
  } else {
  }
  return pass;
}

bool TrackCuts::PIDCuts(Tracks *track) {
  // We would need some detector response
  bool pass = true;

  // detector response check
  bool ECalisthere = false;
  bool HTCCisthere = false;

  if (pass) {
    float ChiSquared = 3.0; // needed from the track write a function here
    if (TMath::Abs(ChiSquared) > fChiSquared) {
      pass = false;
    } else {
      // if (!fMinimalBooking)
      //   fHists;
    }
  }
  return pass;
}

/// photon for study of DVCS or other events with dummy values
TrackCuts *TrackCuts::PhotonCuts() {
  TrackCuts *trackCuts = new TrackCuts();
  trackCuts->SetChargeCut(0);
  trackCuts->SetMomentumCut(0.0, 10.0); // yijie please retune these
  trackCuts->SetVtxX(-10, 10);
  trackCuts->SetVtxY(-10, 10);
  trackCuts->SetVtxZ(-10, 10);
  trackCuts->SetPIDCut(0, -999);
  return trackCuts;
}

/// electron for study of DVCS or other events with dummy values
TrackCuts *TrackCuts::ElectronCuts() {
  TrackCuts *trackCuts = new TrackCuts();
  trackCuts->SetChargeCut(1);
  trackCuts->SetMomentumCut(0.0, 10.0); // yijie please retune these
  trackCuts->SetVtxX(-10, 10);
  trackCuts->SetVtxY(-10, 10);
  trackCuts->SetVtxZ(-10, 10);
  trackCuts->SetPIDCut(11, 3.0);
  return trackCuts;
}

/// proton for study of DVCS or other events with dummy values
TrackCuts *TrackCuts::ProtonCuts() {
  TrackCuts *trackCuts = new TrackCuts();
  trackCuts->SetChargeCut(1);
  trackCuts->SetMomentumCut(0.0, 10.0); // yijie please retune these
  trackCuts->SetVtxX(-10, 10);
  trackCuts->SetVtxY(-10, 10);
  trackCuts->SetVtxZ(-10, 10);
  trackCuts->SetPIDCut(2212, 3.0);
  return trackCuts;
}

//// booking histograms here
void TrackCuts::BookQA(Tracks *Track) {
  if (!fMinimalBooking) {
    // it might be usefull for us to create a class where we can store these
    // variable already, thinking about the resuing these class after skimming
    // and not sure if this would be possible we should be able to figure out

    float Px = Track->GetPx();
    float Py = Track->GetPy();
    float Pz = Track->GetPz();

    /// we might need to store more variable separately but it seems there are
    /// already some defined in the clas12reader

    float p = Track->GetP();
    float Vz = Track->GetVz(); // not sure if its the correct one, check required
    float beta =Track->GetBeta(); // not sure if its the correct one, check required

    float phi = 0; // not sure if phi is required but lets keep might need for
                   // the fiducial cuts, need to find out how to compute

    for (int i = 0; i < 2; ++i) {
       if (i == 0 || (i == 1 && Track->UseParticle())) {

      fHists->FillphiCut(i, phi);
      fHists->FillpCut(i, p);
      fHists->FillVzCut(i, p, Vz);
      // fHists->FillBeta(i, p, beta);
      //  we may wanna include more histograms here

       }
    }
  } else {
     if (Track->UseParticle()){
    fHists->FillpCut(1, Track->GetP());
     }
    return;
  }
}

// we need more such functions for the track cuts
//  Momentum cut
bool TrackCuts::PassesMomentumCut(Tracks *track) const {
  float p = track->GetP();
  return (p >= fPmin && p <= fPmax);
}

// Vertex Z position cut
bool TrackCuts::PassesVertexCut(Tracks *track) const {
  return (track->GetVz() >= fVtZmin && track->GetVz() <= fVtZmax);
}

// PID cut
bool TrackCuts::PassesPIDCut(Tracks *track) const {
  return (track->GetPDGCode() == fPdg);
}

// Charge cut
bool TrackCuts::PassesChargeCut(Tracks *track) const {
  return (track->GetCharge() == fCharge);
}
