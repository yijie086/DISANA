#ifndef TRACKCUTS_H_
#define TRACKCUTS_H_
#include "Event.h"
#include "Tracks.h"
#include "TrackHists.h"

class TrackCuts {
 public:
  TrackCuts();
  TrackCuts(const TrackCuts& cuts);
  TrackCuts &operator=(const TrackCuts& cuts);
  virtual ~TrackCuts();

  static TrackCuts *PhotonCuts();
  static TrackCuts *ElectronCuts();
  static TrackCuts *ProtonCuts();

 
    void SetMomentumCut(float minP, float maxP){
        fPmin = minP;
        fPmax = maxP;
        fcutP = true;
    };
    void SetVtxX(float minVx, float maxVx){
        fVtXmin = minVx;
        fVtXmax = maxVx;
    };
    void SetVtxY(float minVy, float maxVy){
        fVtYmin = minVy;
        fVtYmax = maxVy;
    };
    void SetVtxZ(float minVz, float maxVz){
        fVtZmin = minVz;
        fVtZmax = maxVz;
    };
    void SetPIDCut(int pid = -999, float chiSquaredVal =-999){
        fPdg = pid;
        fChiSquared = chiSquaredVal; // this is the chiSquared for the particle/// looking at DVCS Analysis note
        fcutPID = kTRUE;
    };
    void SetChargeCut(int charge){
        fCharge = charge;
        };

    void SetMinimalBook(bool book ){
        fMinimalBooking = book;
    }
    
    // Check if a track passes all cuts
    bool isSelected(Tracks *track);
    bool TrackingCuts(Tracks *track);
    bool PIDCuts(Tracks *track);
    void BookQA(Tracks *track);

private:
    bool fMinimalBooking;               
    int  fPdg;
    bool  fcutPID;
    int   fCharge;
    float fPmin;
    float fPmax;
    float fcutP;
    float fVtXmin;
    float fVtXmax;
    float fVtYmin;
    float fVtYmax;
    float fVtZmin;
    float fVtZmax;
    float fChiSquared;

    TrackHists *fHists;    

    

    // Individual cut functions
    bool PassesMomentumCut(Tracks* track) const;
    bool PassesVertexCut(Tracks*track) const;
    bool PassesPIDCut(Tracks* track) const;
    bool PassesChargeCut(Tracks*track) const;
};

#endif // TRACKCUTS_H
