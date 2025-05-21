#ifndef TRACKHISS_H
#define TRACKHISTS_H

 #include "Rtypes.h"
 #include "TH1F.h"
 #include "TH2F.h"
 #include "TH3F.h"
 #include "TList.h"
 #include "TProfile.h"
 
 class TrackHists {
  public:

   TrackHists();
   TrackHists(bool VzDist, float pmin = 0.0, float pmax = 12);
   TrackHists(TString MinimalBooking);
   virtual ~TrackHists();
   void FillConfig(int iBin, float val) {
     if (!fMinimalBooking)
       fConfig->Fill(iBin, val);
   }
   ;
   void FillTrackCounter(int iBin) {
     if (!fMinimalBooking)
       fCutCounter->Fill(iBin);
   }
   ;
   void FillpCut(int i, float p) {
     if (!fMinimalBooking)
       fpDist[i]->Fill(p);
   }
   ;   
   void FillphiCut(int i, float phi) {
     if (!fMinimalBooking)
       fphiDist[i]->Fill(phi);
   }
   ;
   void FillVzCut(int i, float p, float vz) {
     if (!fMinimalBooking)
       fVz[i]->Fill(p, vz);
   }
   ;

   void FillTrackChiSquare(int i, float p, float chi2) {
     if (!fMinimalBooking)
       fTrackChi2[i]->Fill(p, chi2);
   }
 
   void SetName(TString name) {
     fHistList->SetName(name.Data());
   }
   ;
   TList *GetHistList() {
     return fHistList;
   }
   ;
  private:
   TrackHists &operator=(const TrackHists &obj);
   TrackHists(const TrackHists&);
   bool fMinimalBooking;     
   float fpmin; 
   float fpmax;                           
   TList *fHistList;  
   TList *fTrackCutQA[2];           
   TProfile *fConfig;       
   TH1F *fCutCounter;        
   TH1F *fpDist[2];          
   TH1F *fetaDist[2];        
   TH1F *fphiDist[2];        
   TH2F *fVz[2];             
   TH2F *fTrackChi2[2];      
 };
 
 #endif /* TRACKHISTS_H_ */