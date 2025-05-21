 #ifndef TRACKS_H_
 #define TRACKS_H_
 
 #include "TVector3.h"
 #include "clas12reader.h"  // Include CLAS12 reader

 using namespace clas12;
 class Tracks {

  public:
   Tracks(const int part = 1);
   Tracks(const Tracks& part);
   Tracks &operator=(const Tracks &obj);
   void SetTrack( region_particle *track);
   virtual ~Tracks();


   void SetMomentum(unsigned int iEntr, TVector3 mom) {
     SetMomentum(iEntr, mom.X(), mom.Y(), mom.Z());
   }
   void SetMomentum(unsigned int iEntr, float px, float py, float pz) {
     if (iEntr < fP.size()) {
       fP.at(iEntr).SetXYZ(px, py, pz);
     } else {
      std::cerr<<"Trying to set momentum for a too small vector!"<<std::endl;
     }
   }
   ;
   
   TVector3 GetMomentum(unsigned int iEntr = 0) const {
     if (iEntr > fP.size()) {
       std::cout << "Trying to get a momentum which is out of bounds. iEntr = "
                 << iEntr << std::endl;
     } else {
       return fP.at(iEntr);
     }
     return TVector3(-999,-999,-999);
   }
   ;
   std::vector<TVector3> GetMomenta() const {
     return fP;
   }
   float GetP() const {
     return GetMomentum().Mag();
   }
   float GetPx() const {
     return GetMomentum().Px();
   }
   float GetPy() const {
     return GetMomentum().Py();
   }
   float GetPz() const {
     return GetMomentum().Pz();
   }

   void SetIDTracks(int idTracks) {
     fIDTracks.push_back(idTracks);
   }
   ;
   std::vector<int> GetIDTracks() const {
     return fIDTracks;
   }
   ;
   void SetCharge(int charge) {
     fCharge=charge;
   }
   ;
   int GetCharge() const {
     return fCharge;
   }
   ;

   void SetInvMass(float invMass) {
     fInvMass = invMass;
   }
   float GetInvMass() const {
     return fInvMass;
   }

   void SetPDGCode(int pdgCode) {
     fPDGCode = pdgCode;
   }
   ;
   int GetPDGCode() const {
     return fPDGCode;
   }
   ;
   void SetID(int ID) {
     fID = ID;
   }
   int GetID() const {
     return fID;
   }
   void SetEvtNumber(int evtNumb) {
     fEvtNumber = evtNumb;
   }
   ;
   int GetEvtNumber() const {
     return fEvtNumber;
   }
   ;
   bool IsSet() const {
     return fIsSet;
   }
   ;
   void SetUse(bool use) {
     fUse = use;
   }
   ;
   void SetVz(float Vz) { fZVtx = Vz; };
   float GetVz()  const { return fZVtx; };

   
   bool UseParticle() const {
     return fUse;
   }
   ;

   float GetBeta() const;

   TString ClassName() {
     return "Tracks";
   }
   ;
 


  private:
  void Reset();
  void SetTrackingInformation() ;
   bool fIsReset;
   std::vector<TVector3> fP;
   std::vector<int> fIDTracks;

   // Ask Richard for the particle object if this makes any good or there are other type of object one can use. to me it looks like there aren't better options from the clas12reader method.
   region_particle * fParticle;
   float fInvMass;
   float fZVtx;
   float fChi2;
 // pdg code as set by the track cuts, used for invariant mass calculation/mc matching in v0s
   int fPDGCode;

 // pdg code as obtained from the MC for this particle
   int fID;
   int fCharge;
   int fEvtNumber;
   bool fUse;    //passes cuts
   bool fIsSet;  //has all the attributes set properly
  private:
 };
 
 #endif /* TRACKS_H_*/