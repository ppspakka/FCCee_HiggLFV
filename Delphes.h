//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar  7 00:20:23 2025 by ROOT version 6.30/04
// from TTree Delphes/Analysis tree
// found on file: /work/project/physics/psriling/public/delphes_output_10.root
//////////////////////////////////////////////////////////

#ifndef Delphes_h
#define Delphes_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"

class Delphes {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxEvent = 1;
   static constexpr Int_t kMaxWeight = 2;
   static constexpr Int_t kMaxParticle = 10946;
   static constexpr Int_t kMaxGenMissingET = 1;
   static constexpr Int_t kMaxEFlowTrack = 1011;
   static constexpr Int_t kMaxEFlowPhoton = 1102;
   static constexpr Int_t kMaxEFlowNeutralHadron = 742;
   static constexpr Int_t kMaxPhoton = 11;
   static constexpr Int_t kMaxElectron = 3;
   static constexpr Int_t kMaxMuon = 4;
   static constexpr Int_t kMaxJet = 26;
   static constexpr Int_t kMaxPhotonNoIso = 122;
   static constexpr Int_t kMaxElectronNoIso = 6;
   static constexpr Int_t kMaxMuonNoIso = 5;
   static constexpr Int_t kMaxJetNoIso = 26;
   static constexpr Int_t kMaxGenJet02 = 13;
   static constexpr Int_t kMaxGenJet04 = 26;
   static constexpr Int_t kMaxGenJet08 = 33;
   static constexpr Int_t kMaxGenJet15 = 17;
   static constexpr Int_t kMaxParticleFlowJet02 = 12;
   static constexpr Int_t kMaxParticleFlowJet04 = 26;
   static constexpr Int_t kMaxParticleFlowJet08 = 34;
   static constexpr Int_t kMaxParticleFlowJet15 = 16;
   static constexpr Int_t kMaxCaloJet02 = 10;
   static constexpr Int_t kMaxCaloJet04 = 17;
   static constexpr Int_t kMaxCaloJet08 = 33;
   static constexpr Int_t kMaxCaloJet15 = 15;
   static constexpr Int_t kMaxTrackJet02 = 10;
   static constexpr Int_t kMaxTrackJet04 = 10;
   static constexpr Int_t kMaxTrackJet08 = 20;
   static constexpr Int_t kMaxTrackJet15 = 13;
   static constexpr Int_t kMaxMissingET = 1;
   static constexpr Int_t kMaxScalarHT = 1;

   // Declaration of leaf types
   Int_t           Event_;
   UInt_t          Event_fUniqueID[kMaxEvent];   //[Event_]
   UInt_t          Event_fBits[kMaxEvent];   //[Event_]
   Long64_t        Event_Number[kMaxEvent];   //[Event_]
   Float_t         Event_ReadTime[kMaxEvent];   //[Event_]
   Float_t         Event_ProcTime[kMaxEvent];   //[Event_]
   Int_t           Event_ProcessID[kMaxEvent];   //[Event_]
   Int_t           Event_MPI[kMaxEvent];   //[Event_]
   Float_t         Event_Weight[kMaxEvent];   //[Event_]
   Float_t         Event_CrossSection[kMaxEvent];   //[Event_]
   Float_t         Event_CrossSectionError[kMaxEvent];   //[Event_]
   Float_t         Event_Scale[kMaxEvent];   //[Event_]
   Float_t         Event_AlphaQED[kMaxEvent];   //[Event_]
   Float_t         Event_AlphaQCD[kMaxEvent];   //[Event_]
   Int_t           Event_ID1[kMaxEvent];   //[Event_]
   Int_t           Event_ID2[kMaxEvent];   //[Event_]
   Float_t         Event_X1[kMaxEvent];   //[Event_]
   Float_t         Event_X2[kMaxEvent];   //[Event_]
   Float_t         Event_ScalePDF[kMaxEvent];   //[Event_]
   Float_t         Event_PDF1[kMaxEvent];   //[Event_]
   Float_t         Event_PDF2[kMaxEvent];   //[Event_]
   Int_t           Event_size;
   Int_t           Weight_;
   UInt_t          Weight_fUniqueID[kMaxWeight];   //[Weight_]
   UInt_t          Weight_fBits[kMaxWeight];   //[Weight_]
   Float_t         Weight_Weight[kMaxWeight];   //[Weight_]
   Int_t           Weight_size;
   Int_t           Particle_;
   UInt_t          Particle_fUniqueID[kMaxParticle];   //[Particle_]
   UInt_t          Particle_fBits[kMaxParticle];   //[Particle_]
   Int_t           Particle_PID[kMaxParticle];   //[Particle_]
   Int_t           Particle_Status[kMaxParticle];   //[Particle_]
   Int_t           Particle_IsPU[kMaxParticle];   //[Particle_]
   Int_t           Particle_M1[kMaxParticle];   //[Particle_]
   Int_t           Particle_M2[kMaxParticle];   //[Particle_]
   Int_t           Particle_D1[kMaxParticle];   //[Particle_]
   Int_t           Particle_D2[kMaxParticle];   //[Particle_]
   Int_t           Particle_Charge[kMaxParticle];   //[Particle_]
   Float_t         Particle_Mass[kMaxParticle];   //[Particle_]
   Float_t         Particle_E[kMaxParticle];   //[Particle_]
   Float_t         Particle_Px[kMaxParticle];   //[Particle_]
   Float_t         Particle_Py[kMaxParticle];   //[Particle_]
   Float_t         Particle_Pz[kMaxParticle];   //[Particle_]
   Float_t         Particle_P[kMaxParticle];   //[Particle_]
   Float_t         Particle_PT[kMaxParticle];   //[Particle_]
   Float_t         Particle_Eta[kMaxParticle];   //[Particle_]
   Float_t         Particle_Phi[kMaxParticle];   //[Particle_]
   Float_t         Particle_Rapidity[kMaxParticle];   //[Particle_]
   Float_t         Particle_T[kMaxParticle];   //[Particle_]
   Float_t         Particle_X[kMaxParticle];   //[Particle_]
   Float_t         Particle_Y[kMaxParticle];   //[Particle_]
   Float_t         Particle_Z[kMaxParticle];   //[Particle_]
   Int_t           Particle_size;
   Int_t           GenMissingET_;
   UInt_t          GenMissingET_fUniqueID[kMaxGenMissingET];   //[GenMissingET_]
   UInt_t          GenMissingET_fBits[kMaxGenMissingET];   //[GenMissingET_]
   Float_t         GenMissingET_MET[kMaxGenMissingET];   //[GenMissingET_]
   Float_t         GenMissingET_Eta[kMaxGenMissingET];   //[GenMissingET_]
   Float_t         GenMissingET_Phi[kMaxGenMissingET];   //[GenMissingET_]
   Int_t           GenMissingET_size;
   Int_t           EFlowTrack_;
   UInt_t          EFlowTrack_fUniqueID[kMaxEFlowTrack];   //[EFlowTrack_]
   UInt_t          EFlowTrack_fBits[kMaxEFlowTrack];   //[EFlowTrack_]
   Int_t           EFlowTrack_PID[kMaxEFlowTrack];   //[EFlowTrack_]
   Int_t           EFlowTrack_Charge[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_P[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_PT[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Eta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Phi[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_CtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_C[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Mass[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_EtaOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_PhiOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_T[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_X[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Y[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Z[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_TOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_XOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_YOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ZOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Xd[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Yd[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Zd[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_L[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_D0[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_DZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Nclusters[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_dNdx[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorP[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPT[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPhi[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorCtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorT[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorDZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorC[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0Phi[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0C[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0DZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0CtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPhiC[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPhiDZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPhiCtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorCDZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorCCtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorDZCtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   TRef            EFlowTrack_Particle[kMaxEFlowTrack];
   Int_t           EFlowTrack_VertexIndex[kMaxEFlowTrack];   //[EFlowTrack_]
   Int_t           EFlowTrack_size;
   Int_t           EFlowPhoton_;
   UInt_t          EFlowPhoton_fUniqueID[kMaxEFlowPhoton];   //[EFlowPhoton_]
   UInt_t          EFlowPhoton_fBits[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_ET[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_Eta[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_Phi[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_E[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_T[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Int_t           EFlowPhoton_NTimeHits[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_Eem[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_Ehad[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_Edges[kMaxEFlowPhoton][4];   //[EFlowPhoton_]
   TRefArray       EFlowPhoton_Particles[kMaxEFlowPhoton];
   Int_t           EFlowPhoton_size;
   Int_t           EFlowNeutralHadron_;
   UInt_t          EFlowNeutralHadron_fUniqueID[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   UInt_t          EFlowNeutralHadron_fBits[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_ET[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_Eta[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_Phi[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_E[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_T[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Int_t           EFlowNeutralHadron_NTimeHits[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_Eem[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_Ehad[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_Edges[kMaxEFlowNeutralHadron][4];   //[EFlowNeutralHadron_]
   TRefArray       EFlowNeutralHadron_Particles[kMaxEFlowNeutralHadron];
   Int_t           EFlowNeutralHadron_size;
   Int_t           Photon_;
   UInt_t          Photon_fUniqueID[kMaxPhoton];   //[Photon_]
   UInt_t          Photon_fBits[kMaxPhoton];   //[Photon_]
   Float_t         Photon_PT[kMaxPhoton];   //[Photon_]
   Float_t         Photon_Eta[kMaxPhoton];   //[Photon_]
   Float_t         Photon_Phi[kMaxPhoton];   //[Photon_]
   Float_t         Photon_E[kMaxPhoton];   //[Photon_]
   Float_t         Photon_T[kMaxPhoton];   //[Photon_]
   Float_t         Photon_EhadOverEem[kMaxPhoton];   //[Photon_]
   TRefArray       Photon_Particles[kMaxPhoton];
   Float_t         Photon_IsolationVar[kMaxPhoton];   //[Photon_]
   Float_t         Photon_IsolationVarRhoCorr[kMaxPhoton];   //[Photon_]
   Float_t         Photon_SumPtCharged[kMaxPhoton];   //[Photon_]
   Float_t         Photon_SumPtNeutral[kMaxPhoton];   //[Photon_]
   Float_t         Photon_SumPtChargedPU[kMaxPhoton];   //[Photon_]
   Float_t         Photon_SumPt[kMaxPhoton];   //[Photon_]
   Int_t           Photon_Status[kMaxPhoton];   //[Photon_]
   Int_t           Photon_size;
   Int_t           Electron_;
   UInt_t          Electron_fUniqueID[kMaxElectron];   //[Electron_]
   UInt_t          Electron_fBits[kMaxElectron];   //[Electron_]
   Float_t         Electron_PT[kMaxElectron];   //[Electron_]
   Float_t         Electron_Eta[kMaxElectron];   //[Electron_]
   Float_t         Electron_Phi[kMaxElectron];   //[Electron_]
   Float_t         Electron_T[kMaxElectron];   //[Electron_]
   Int_t           Electron_Charge[kMaxElectron];   //[Electron_]
   Float_t         Electron_EhadOverEem[kMaxElectron];   //[Electron_]
   TRef            Electron_Particle[kMaxElectron];
   Float_t         Electron_IsolationVar[kMaxElectron];   //[Electron_]
   Float_t         Electron_IsolationVarRhoCorr[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPtCharged[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPtNeutral[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPtChargedPU[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPt[kMaxElectron];   //[Electron_]
   Float_t         Electron_D0[kMaxElectron];   //[Electron_]
   Float_t         Electron_DZ[kMaxElectron];   //[Electron_]
   Float_t         Electron_ErrorD0[kMaxElectron];   //[Electron_]
   Float_t         Electron_ErrorDZ[kMaxElectron];   //[Electron_]
   Int_t           Electron_size;
   Int_t           Muon_;
   UInt_t          Muon_fUniqueID[kMaxMuon];   //[Muon_]
   UInt_t          Muon_fBits[kMaxMuon];   //[Muon_]
   Float_t         Muon_PT[kMaxMuon];   //[Muon_]
   Float_t         Muon_Eta[kMaxMuon];   //[Muon_]
   Float_t         Muon_Phi[kMaxMuon];   //[Muon_]
   Float_t         Muon_T[kMaxMuon];   //[Muon_]
   Int_t           Muon_Charge[kMaxMuon];   //[Muon_]
   TRef            Muon_Particle[kMaxMuon];
   Float_t         Muon_IsolationVar[kMaxMuon];   //[Muon_]
   Float_t         Muon_IsolationVarRhoCorr[kMaxMuon];   //[Muon_]
   Float_t         Muon_SumPtCharged[kMaxMuon];   //[Muon_]
   Float_t         Muon_SumPtNeutral[kMaxMuon];   //[Muon_]
   Float_t         Muon_SumPtChargedPU[kMaxMuon];   //[Muon_]
   Float_t         Muon_SumPt[kMaxMuon];   //[Muon_]
   Float_t         Muon_D0[kMaxMuon];   //[Muon_]
   Float_t         Muon_DZ[kMaxMuon];   //[Muon_]
   Float_t         Muon_ErrorD0[kMaxMuon];   //[Muon_]
   Float_t         Muon_ErrorDZ[kMaxMuon];   //[Muon_]
   Int_t           Muon_size;
   Int_t           Jet_;
   UInt_t          Jet_fUniqueID[kMaxJet];   //[Jet_]
   UInt_t          Jet_fBits[kMaxJet];   //[Jet_]
   Float_t         Jet_PT[kMaxJet];   //[Jet_]
   Float_t         Jet_Eta[kMaxJet];   //[Jet_]
   Float_t         Jet_Phi[kMaxJet];   //[Jet_]
   Float_t         Jet_T[kMaxJet];   //[Jet_]
   Float_t         Jet_Mass[kMaxJet];   //[Jet_]
   Float_t         Jet_DeltaEta[kMaxJet];   //[Jet_]
   Float_t         Jet_DeltaPhi[kMaxJet];   //[Jet_]
   UInt_t          Jet_Flavor[kMaxJet];   //[Jet_]
   UInt_t          Jet_FlavorAlgo[kMaxJet];   //[Jet_]
   UInt_t          Jet_FlavorPhys[kMaxJet];   //[Jet_]
   UInt_t          Jet_BTag[kMaxJet];   //[Jet_]
   UInt_t          Jet_BTagAlgo[kMaxJet];   //[Jet_]
   UInt_t          Jet_BTagPhys[kMaxJet];   //[Jet_]
   UInt_t          Jet_TauTag[kMaxJet];   //[Jet_]
   Float_t         Jet_TauWeight[kMaxJet];   //[Jet_]
   Int_t           Jet_Charge[kMaxJet];   //[Jet_]
   Float_t         Jet_EhadOverEem[kMaxJet];   //[Jet_]
   Int_t           Jet_NCharged[kMaxJet];   //[Jet_]
   Int_t           Jet_NNeutrals[kMaxJet];   //[Jet_]
   Float_t         Jet_NeutralEnergyFraction[kMaxJet];   //[Jet_]
   Float_t         Jet_ChargedEnergyFraction[kMaxJet];   //[Jet_]
   Float_t         Jet_Beta[kMaxJet];   //[Jet_]
   Float_t         Jet_BetaStar[kMaxJet];   //[Jet_]
   Float_t         Jet_MeanSqDeltaR[kMaxJet];   //[Jet_]
   Float_t         Jet_PTD[kMaxJet];   //[Jet_]
   Float_t         Jet_FracPt[kMaxJet][5];   //[Jet_]
   Float_t         Jet_Tau[kMaxJet][5];   //[Jet_]
   TLorentzVector  Jet_SoftDroppedJet[kMaxJet];
   TLorentzVector  Jet_SoftDroppedSubJet1[kMaxJet];
   TLorentzVector  Jet_SoftDroppedSubJet2[kMaxJet];
   TLorentzVector  Jet_TrimmedP4[5][kMaxJet];
   TLorentzVector  Jet_PrunedP4[5][kMaxJet];
   TLorentzVector  Jet_SoftDroppedP4[5][kMaxJet];
   Int_t           Jet_NSubJetsTrimmed[kMaxJet];   //[Jet_]
   Int_t           Jet_NSubJetsPruned[kMaxJet];   //[Jet_]
   Int_t           Jet_NSubJetsSoftDropped[kMaxJet];   //[Jet_]
   Double_t        Jet_ExclYmerge23[kMaxJet];   //[Jet_]
   Double_t        Jet_ExclYmerge34[kMaxJet];   //[Jet_]
   Double_t        Jet_ExclYmerge45[kMaxJet];   //[Jet_]
   Double_t        Jet_ExclYmerge56[kMaxJet];   //[Jet_]
   TRefArray       Jet_Constituents[kMaxJet];
   TRefArray       Jet_Particles[kMaxJet];
   TLorentzVector  Jet_Area[kMaxJet];
   Int_t           Jet_size;
   Int_t           PhotonNoIso_;
   UInt_t          PhotonNoIso_fUniqueID[kMaxPhotonNoIso];   //[PhotonNoIso_]
   UInt_t          PhotonNoIso_fBits[kMaxPhotonNoIso];   //[PhotonNoIso_]
   Float_t         PhotonNoIso_PT[kMaxPhotonNoIso];   //[PhotonNoIso_]
   Float_t         PhotonNoIso_Eta[kMaxPhotonNoIso];   //[PhotonNoIso_]
   Float_t         PhotonNoIso_Phi[kMaxPhotonNoIso];   //[PhotonNoIso_]
   Float_t         PhotonNoIso_E[kMaxPhotonNoIso];   //[PhotonNoIso_]
   Float_t         PhotonNoIso_T[kMaxPhotonNoIso];   //[PhotonNoIso_]
   Float_t         PhotonNoIso_EhadOverEem[kMaxPhotonNoIso];   //[PhotonNoIso_]
   TRefArray       PhotonNoIso_Particles[kMaxPhotonNoIso];
   Float_t         PhotonNoIso_IsolationVar[kMaxPhotonNoIso];   //[PhotonNoIso_]
   Float_t         PhotonNoIso_IsolationVarRhoCorr[kMaxPhotonNoIso];   //[PhotonNoIso_]
   Float_t         PhotonNoIso_SumPtCharged[kMaxPhotonNoIso];   //[PhotonNoIso_]
   Float_t         PhotonNoIso_SumPtNeutral[kMaxPhotonNoIso];   //[PhotonNoIso_]
   Float_t         PhotonNoIso_SumPtChargedPU[kMaxPhotonNoIso];   //[PhotonNoIso_]
   Float_t         PhotonNoIso_SumPt[kMaxPhotonNoIso];   //[PhotonNoIso_]
   Int_t           PhotonNoIso_Status[kMaxPhotonNoIso];   //[PhotonNoIso_]
   Int_t           PhotonNoIso_size;
   Int_t           ElectronNoIso_;
   UInt_t          ElectronNoIso_fUniqueID[kMaxElectronNoIso];   //[ElectronNoIso_]
   UInt_t          ElectronNoIso_fBits[kMaxElectronNoIso];   //[ElectronNoIso_]
   Float_t         ElectronNoIso_PT[kMaxElectronNoIso];   //[ElectronNoIso_]
   Float_t         ElectronNoIso_Eta[kMaxElectronNoIso];   //[ElectronNoIso_]
   Float_t         ElectronNoIso_Phi[kMaxElectronNoIso];   //[ElectronNoIso_]
   Float_t         ElectronNoIso_T[kMaxElectronNoIso];   //[ElectronNoIso_]
   Int_t           ElectronNoIso_Charge[kMaxElectronNoIso];   //[ElectronNoIso_]
   Float_t         ElectronNoIso_EhadOverEem[kMaxElectronNoIso];   //[ElectronNoIso_]
   TRef            ElectronNoIso_Particle[kMaxElectronNoIso];
   Float_t         ElectronNoIso_IsolationVar[kMaxElectronNoIso];   //[ElectronNoIso_]
   Float_t         ElectronNoIso_IsolationVarRhoCorr[kMaxElectronNoIso];   //[ElectronNoIso_]
   Float_t         ElectronNoIso_SumPtCharged[kMaxElectronNoIso];   //[ElectronNoIso_]
   Float_t         ElectronNoIso_SumPtNeutral[kMaxElectronNoIso];   //[ElectronNoIso_]
   Float_t         ElectronNoIso_SumPtChargedPU[kMaxElectronNoIso];   //[ElectronNoIso_]
   Float_t         ElectronNoIso_SumPt[kMaxElectronNoIso];   //[ElectronNoIso_]
   Float_t         ElectronNoIso_D0[kMaxElectronNoIso];   //[ElectronNoIso_]
   Float_t         ElectronNoIso_DZ[kMaxElectronNoIso];   //[ElectronNoIso_]
   Float_t         ElectronNoIso_ErrorD0[kMaxElectronNoIso];   //[ElectronNoIso_]
   Float_t         ElectronNoIso_ErrorDZ[kMaxElectronNoIso];   //[ElectronNoIso_]
   Int_t           ElectronNoIso_size;
   Int_t           MuonNoIso_;
   UInt_t          MuonNoIso_fUniqueID[kMaxMuonNoIso];   //[MuonNoIso_]
   UInt_t          MuonNoIso_fBits[kMaxMuonNoIso];   //[MuonNoIso_]
   Float_t         MuonNoIso_PT[kMaxMuonNoIso];   //[MuonNoIso_]
   Float_t         MuonNoIso_Eta[kMaxMuonNoIso];   //[MuonNoIso_]
   Float_t         MuonNoIso_Phi[kMaxMuonNoIso];   //[MuonNoIso_]
   Float_t         MuonNoIso_T[kMaxMuonNoIso];   //[MuonNoIso_]
   Int_t           MuonNoIso_Charge[kMaxMuonNoIso];   //[MuonNoIso_]
   TRef            MuonNoIso_Particle[kMaxMuonNoIso];
   Float_t         MuonNoIso_IsolationVar[kMaxMuonNoIso];   //[MuonNoIso_]
   Float_t         MuonNoIso_IsolationVarRhoCorr[kMaxMuonNoIso];   //[MuonNoIso_]
   Float_t         MuonNoIso_SumPtCharged[kMaxMuonNoIso];   //[MuonNoIso_]
   Float_t         MuonNoIso_SumPtNeutral[kMaxMuonNoIso];   //[MuonNoIso_]
   Float_t         MuonNoIso_SumPtChargedPU[kMaxMuonNoIso];   //[MuonNoIso_]
   Float_t         MuonNoIso_SumPt[kMaxMuonNoIso];   //[MuonNoIso_]
   Float_t         MuonNoIso_D0[kMaxMuonNoIso];   //[MuonNoIso_]
   Float_t         MuonNoIso_DZ[kMaxMuonNoIso];   //[MuonNoIso_]
   Float_t         MuonNoIso_ErrorD0[kMaxMuonNoIso];   //[MuonNoIso_]
   Float_t         MuonNoIso_ErrorDZ[kMaxMuonNoIso];   //[MuonNoIso_]
   Int_t           MuonNoIso_size;
   Int_t           JetNoIso_;
   UInt_t          JetNoIso_fUniqueID[kMaxJetNoIso];   //[JetNoIso_]
   UInt_t          JetNoIso_fBits[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_PT[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_Eta[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_Phi[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_T[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_Mass[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_DeltaEta[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_DeltaPhi[kMaxJetNoIso];   //[JetNoIso_]
   UInt_t          JetNoIso_Flavor[kMaxJetNoIso];   //[JetNoIso_]
   UInt_t          JetNoIso_FlavorAlgo[kMaxJetNoIso];   //[JetNoIso_]
   UInt_t          JetNoIso_FlavorPhys[kMaxJetNoIso];   //[JetNoIso_]
   UInt_t          JetNoIso_BTag[kMaxJetNoIso];   //[JetNoIso_]
   UInt_t          JetNoIso_BTagAlgo[kMaxJetNoIso];   //[JetNoIso_]
   UInt_t          JetNoIso_BTagPhys[kMaxJetNoIso];   //[JetNoIso_]
   UInt_t          JetNoIso_TauTag[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_TauWeight[kMaxJetNoIso];   //[JetNoIso_]
   Int_t           JetNoIso_Charge[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_EhadOverEem[kMaxJetNoIso];   //[JetNoIso_]
   Int_t           JetNoIso_NCharged[kMaxJetNoIso];   //[JetNoIso_]
   Int_t           JetNoIso_NNeutrals[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_NeutralEnergyFraction[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_ChargedEnergyFraction[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_Beta[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_BetaStar[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_MeanSqDeltaR[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_PTD[kMaxJetNoIso];   //[JetNoIso_]
   Float_t         JetNoIso_FracPt[kMaxJetNoIso][5];   //[JetNoIso_]
   Float_t         JetNoIso_Tau[kMaxJetNoIso][5];   //[JetNoIso_]
   TLorentzVector  JetNoIso_SoftDroppedJet[kMaxJetNoIso];
   TLorentzVector  JetNoIso_SoftDroppedSubJet1[kMaxJetNoIso];
   TLorentzVector  JetNoIso_SoftDroppedSubJet2[kMaxJetNoIso];
   TLorentzVector  JetNoIso_TrimmedP4[5][kMaxJetNoIso];
   TLorentzVector  JetNoIso_PrunedP4[5][kMaxJetNoIso];
   TLorentzVector  JetNoIso_SoftDroppedP4[5][kMaxJetNoIso];
   Int_t           JetNoIso_NSubJetsTrimmed[kMaxJetNoIso];   //[JetNoIso_]
   Int_t           JetNoIso_NSubJetsPruned[kMaxJetNoIso];   //[JetNoIso_]
   Int_t           JetNoIso_NSubJetsSoftDropped[kMaxJetNoIso];   //[JetNoIso_]
   Double_t        JetNoIso_ExclYmerge23[kMaxJetNoIso];   //[JetNoIso_]
   Double_t        JetNoIso_ExclYmerge34[kMaxJetNoIso];   //[JetNoIso_]
   Double_t        JetNoIso_ExclYmerge45[kMaxJetNoIso];   //[JetNoIso_]
   Double_t        JetNoIso_ExclYmerge56[kMaxJetNoIso];   //[JetNoIso_]
   TRefArray       JetNoIso_Constituents[kMaxJetNoIso];
   TRefArray       JetNoIso_Particles[kMaxJetNoIso];
   TLorentzVector  JetNoIso_Area[kMaxJetNoIso];
   Int_t           JetNoIso_size;
   Int_t           GenJet02_;
   UInt_t          GenJet02_fUniqueID[kMaxGenJet02];   //[GenJet02_]
   UInt_t          GenJet02_fBits[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_PT[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_Eta[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_Phi[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_T[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_Mass[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_DeltaEta[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_DeltaPhi[kMaxGenJet02];   //[GenJet02_]
   UInt_t          GenJet02_Flavor[kMaxGenJet02];   //[GenJet02_]
   UInt_t          GenJet02_FlavorAlgo[kMaxGenJet02];   //[GenJet02_]
   UInt_t          GenJet02_FlavorPhys[kMaxGenJet02];   //[GenJet02_]
   UInt_t          GenJet02_BTag[kMaxGenJet02];   //[GenJet02_]
   UInt_t          GenJet02_BTagAlgo[kMaxGenJet02];   //[GenJet02_]
   UInt_t          GenJet02_BTagPhys[kMaxGenJet02];   //[GenJet02_]
   UInt_t          GenJet02_TauTag[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_TauWeight[kMaxGenJet02];   //[GenJet02_]
   Int_t           GenJet02_Charge[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_EhadOverEem[kMaxGenJet02];   //[GenJet02_]
   Int_t           GenJet02_NCharged[kMaxGenJet02];   //[GenJet02_]
   Int_t           GenJet02_NNeutrals[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_NeutralEnergyFraction[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_ChargedEnergyFraction[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_Beta[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_BetaStar[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_MeanSqDeltaR[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_PTD[kMaxGenJet02];   //[GenJet02_]
   Float_t         GenJet02_FracPt[kMaxGenJet02][5];   //[GenJet02_]
   Float_t         GenJet02_Tau[kMaxGenJet02][5];   //[GenJet02_]
   TLorentzVector  GenJet02_SoftDroppedJet[kMaxGenJet02];
   TLorentzVector  GenJet02_SoftDroppedSubJet1[kMaxGenJet02];
   TLorentzVector  GenJet02_SoftDroppedSubJet2[kMaxGenJet02];
   TLorentzVector  GenJet02_TrimmedP4[5][kMaxGenJet02];
   TLorentzVector  GenJet02_PrunedP4[5][kMaxGenJet02];
   TLorentzVector  GenJet02_SoftDroppedP4[5][kMaxGenJet02];
   Int_t           GenJet02_NSubJetsTrimmed[kMaxGenJet02];   //[GenJet02_]
   Int_t           GenJet02_NSubJetsPruned[kMaxGenJet02];   //[GenJet02_]
   Int_t           GenJet02_NSubJetsSoftDropped[kMaxGenJet02];   //[GenJet02_]
   Double_t        GenJet02_ExclYmerge23[kMaxGenJet02];   //[GenJet02_]
   Double_t        GenJet02_ExclYmerge34[kMaxGenJet02];   //[GenJet02_]
   Double_t        GenJet02_ExclYmerge45[kMaxGenJet02];   //[GenJet02_]
   Double_t        GenJet02_ExclYmerge56[kMaxGenJet02];   //[GenJet02_]
   TRefArray       GenJet02_Constituents[kMaxGenJet02];
   TRefArray       GenJet02_Particles[kMaxGenJet02];
   TLorentzVector  GenJet02_Area[kMaxGenJet02];
   Int_t           GenJet02_size;
   Int_t           GenJet04_;
   UInt_t          GenJet04_fUniqueID[kMaxGenJet04];   //[GenJet04_]
   UInt_t          GenJet04_fBits[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_PT[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_Eta[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_Phi[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_T[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_Mass[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_DeltaEta[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_DeltaPhi[kMaxGenJet04];   //[GenJet04_]
   UInt_t          GenJet04_Flavor[kMaxGenJet04];   //[GenJet04_]
   UInt_t          GenJet04_FlavorAlgo[kMaxGenJet04];   //[GenJet04_]
   UInt_t          GenJet04_FlavorPhys[kMaxGenJet04];   //[GenJet04_]
   UInt_t          GenJet04_BTag[kMaxGenJet04];   //[GenJet04_]
   UInt_t          GenJet04_BTagAlgo[kMaxGenJet04];   //[GenJet04_]
   UInt_t          GenJet04_BTagPhys[kMaxGenJet04];   //[GenJet04_]
   UInt_t          GenJet04_TauTag[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_TauWeight[kMaxGenJet04];   //[GenJet04_]
   Int_t           GenJet04_Charge[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_EhadOverEem[kMaxGenJet04];   //[GenJet04_]
   Int_t           GenJet04_NCharged[kMaxGenJet04];   //[GenJet04_]
   Int_t           GenJet04_NNeutrals[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_NeutralEnergyFraction[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_ChargedEnergyFraction[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_Beta[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_BetaStar[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_MeanSqDeltaR[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_PTD[kMaxGenJet04];   //[GenJet04_]
   Float_t         GenJet04_FracPt[kMaxGenJet04][5];   //[GenJet04_]
   Float_t         GenJet04_Tau[kMaxGenJet04][5];   //[GenJet04_]
   TLorentzVector  GenJet04_SoftDroppedJet[kMaxGenJet04];
   TLorentzVector  GenJet04_SoftDroppedSubJet1[kMaxGenJet04];
   TLorentzVector  GenJet04_SoftDroppedSubJet2[kMaxGenJet04];
   TLorentzVector  GenJet04_TrimmedP4[5][kMaxGenJet04];
   TLorentzVector  GenJet04_PrunedP4[5][kMaxGenJet04];
   TLorentzVector  GenJet04_SoftDroppedP4[5][kMaxGenJet04];
   Int_t           GenJet04_NSubJetsTrimmed[kMaxGenJet04];   //[GenJet04_]
   Int_t           GenJet04_NSubJetsPruned[kMaxGenJet04];   //[GenJet04_]
   Int_t           GenJet04_NSubJetsSoftDropped[kMaxGenJet04];   //[GenJet04_]
   Double_t        GenJet04_ExclYmerge23[kMaxGenJet04];   //[GenJet04_]
   Double_t        GenJet04_ExclYmerge34[kMaxGenJet04];   //[GenJet04_]
   Double_t        GenJet04_ExclYmerge45[kMaxGenJet04];   //[GenJet04_]
   Double_t        GenJet04_ExclYmerge56[kMaxGenJet04];   //[GenJet04_]
   TRefArray       GenJet04_Constituents[kMaxGenJet04];
   TRefArray       GenJet04_Particles[kMaxGenJet04];
   TLorentzVector  GenJet04_Area[kMaxGenJet04];
   Int_t           GenJet04_size;
   Int_t           GenJet08_;
   UInt_t          GenJet08_fUniqueID[kMaxGenJet08];   //[GenJet08_]
   UInt_t          GenJet08_fBits[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_PT[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_Eta[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_Phi[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_T[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_Mass[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_DeltaEta[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_DeltaPhi[kMaxGenJet08];   //[GenJet08_]
   UInt_t          GenJet08_Flavor[kMaxGenJet08];   //[GenJet08_]
   UInt_t          GenJet08_FlavorAlgo[kMaxGenJet08];   //[GenJet08_]
   UInt_t          GenJet08_FlavorPhys[kMaxGenJet08];   //[GenJet08_]
   UInt_t          GenJet08_BTag[kMaxGenJet08];   //[GenJet08_]
   UInt_t          GenJet08_BTagAlgo[kMaxGenJet08];   //[GenJet08_]
   UInt_t          GenJet08_BTagPhys[kMaxGenJet08];   //[GenJet08_]
   UInt_t          GenJet08_TauTag[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_TauWeight[kMaxGenJet08];   //[GenJet08_]
   Int_t           GenJet08_Charge[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_EhadOverEem[kMaxGenJet08];   //[GenJet08_]
   Int_t           GenJet08_NCharged[kMaxGenJet08];   //[GenJet08_]
   Int_t           GenJet08_NNeutrals[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_NeutralEnergyFraction[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_ChargedEnergyFraction[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_Beta[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_BetaStar[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_MeanSqDeltaR[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_PTD[kMaxGenJet08];   //[GenJet08_]
   Float_t         GenJet08_FracPt[kMaxGenJet08][5];   //[GenJet08_]
   Float_t         GenJet08_Tau[kMaxGenJet08][5];   //[GenJet08_]
   TLorentzVector  GenJet08_SoftDroppedJet[kMaxGenJet08];
   TLorentzVector  GenJet08_SoftDroppedSubJet1[kMaxGenJet08];
   TLorentzVector  GenJet08_SoftDroppedSubJet2[kMaxGenJet08];
   TLorentzVector  GenJet08_TrimmedP4[5][kMaxGenJet08];
   TLorentzVector  GenJet08_PrunedP4[5][kMaxGenJet08];
   TLorentzVector  GenJet08_SoftDroppedP4[5][kMaxGenJet08];
   Int_t           GenJet08_NSubJetsTrimmed[kMaxGenJet08];   //[GenJet08_]
   Int_t           GenJet08_NSubJetsPruned[kMaxGenJet08];   //[GenJet08_]
   Int_t           GenJet08_NSubJetsSoftDropped[kMaxGenJet08];   //[GenJet08_]
   Double_t        GenJet08_ExclYmerge23[kMaxGenJet08];   //[GenJet08_]
   Double_t        GenJet08_ExclYmerge34[kMaxGenJet08];   //[GenJet08_]
   Double_t        GenJet08_ExclYmerge45[kMaxGenJet08];   //[GenJet08_]
   Double_t        GenJet08_ExclYmerge56[kMaxGenJet08];   //[GenJet08_]
   TRefArray       GenJet08_Constituents[kMaxGenJet08];
   TRefArray       GenJet08_Particles[kMaxGenJet08];
   TLorentzVector  GenJet08_Area[kMaxGenJet08];
   Int_t           GenJet08_size;
   Int_t           GenJet15_;
   UInt_t          GenJet15_fUniqueID[kMaxGenJet15];   //[GenJet15_]
   UInt_t          GenJet15_fBits[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_PT[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_Eta[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_Phi[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_T[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_Mass[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_DeltaEta[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_DeltaPhi[kMaxGenJet15];   //[GenJet15_]
   UInt_t          GenJet15_Flavor[kMaxGenJet15];   //[GenJet15_]
   UInt_t          GenJet15_FlavorAlgo[kMaxGenJet15];   //[GenJet15_]
   UInt_t          GenJet15_FlavorPhys[kMaxGenJet15];   //[GenJet15_]
   UInt_t          GenJet15_BTag[kMaxGenJet15];   //[GenJet15_]
   UInt_t          GenJet15_BTagAlgo[kMaxGenJet15];   //[GenJet15_]
   UInt_t          GenJet15_BTagPhys[kMaxGenJet15];   //[GenJet15_]
   UInt_t          GenJet15_TauTag[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_TauWeight[kMaxGenJet15];   //[GenJet15_]
   Int_t           GenJet15_Charge[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_EhadOverEem[kMaxGenJet15];   //[GenJet15_]
   Int_t           GenJet15_NCharged[kMaxGenJet15];   //[GenJet15_]
   Int_t           GenJet15_NNeutrals[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_NeutralEnergyFraction[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_ChargedEnergyFraction[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_Beta[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_BetaStar[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_MeanSqDeltaR[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_PTD[kMaxGenJet15];   //[GenJet15_]
   Float_t         GenJet15_FracPt[kMaxGenJet15][5];   //[GenJet15_]
   Float_t         GenJet15_Tau[kMaxGenJet15][5];   //[GenJet15_]
   TLorentzVector  GenJet15_SoftDroppedJet[kMaxGenJet15];
   TLorentzVector  GenJet15_SoftDroppedSubJet1[kMaxGenJet15];
   TLorentzVector  GenJet15_SoftDroppedSubJet2[kMaxGenJet15];
   TLorentzVector  GenJet15_TrimmedP4[5][kMaxGenJet15];
   TLorentzVector  GenJet15_PrunedP4[5][kMaxGenJet15];
   TLorentzVector  GenJet15_SoftDroppedP4[5][kMaxGenJet15];
   Int_t           GenJet15_NSubJetsTrimmed[kMaxGenJet15];   //[GenJet15_]
   Int_t           GenJet15_NSubJetsPruned[kMaxGenJet15];   //[GenJet15_]
   Int_t           GenJet15_NSubJetsSoftDropped[kMaxGenJet15];   //[GenJet15_]
   Double_t        GenJet15_ExclYmerge23[kMaxGenJet15];   //[GenJet15_]
   Double_t        GenJet15_ExclYmerge34[kMaxGenJet15];   //[GenJet15_]
   Double_t        GenJet15_ExclYmerge45[kMaxGenJet15];   //[GenJet15_]
   Double_t        GenJet15_ExclYmerge56[kMaxGenJet15];   //[GenJet15_]
   TRefArray       GenJet15_Constituents[kMaxGenJet15];
   TRefArray       GenJet15_Particles[kMaxGenJet15];
   TLorentzVector  GenJet15_Area[kMaxGenJet15];
   Int_t           GenJet15_size;
   Int_t           ParticleFlowJet02_;
   UInt_t          ParticleFlowJet02_fUniqueID[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   UInt_t          ParticleFlowJet02_fBits[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_PT[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_Eta[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_Phi[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_T[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_Mass[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_DeltaEta[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_DeltaPhi[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   UInt_t          ParticleFlowJet02_Flavor[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   UInt_t          ParticleFlowJet02_FlavorAlgo[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   UInt_t          ParticleFlowJet02_FlavorPhys[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   UInt_t          ParticleFlowJet02_BTag[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   UInt_t          ParticleFlowJet02_BTagAlgo[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   UInt_t          ParticleFlowJet02_BTagPhys[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   UInt_t          ParticleFlowJet02_TauTag[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_TauWeight[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Int_t           ParticleFlowJet02_Charge[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_EhadOverEem[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Int_t           ParticleFlowJet02_NCharged[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Int_t           ParticleFlowJet02_NNeutrals[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_NeutralEnergyFraction[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_ChargedEnergyFraction[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_Beta[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_BetaStar[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_MeanSqDeltaR[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_PTD[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_FracPt[kMaxParticleFlowJet02][5];   //[ParticleFlowJet02_]
   Float_t         ParticleFlowJet02_Tau[kMaxParticleFlowJet02][5];   //[ParticleFlowJet02_]
   TLorentzVector  ParticleFlowJet02_SoftDroppedJet[kMaxParticleFlowJet02];
   TLorentzVector  ParticleFlowJet02_SoftDroppedSubJet1[kMaxParticleFlowJet02];
   TLorentzVector  ParticleFlowJet02_SoftDroppedSubJet2[kMaxParticleFlowJet02];
   TLorentzVector  ParticleFlowJet02_TrimmedP4[5][kMaxParticleFlowJet02];
   TLorentzVector  ParticleFlowJet02_PrunedP4[5][kMaxParticleFlowJet02];
   TLorentzVector  ParticleFlowJet02_SoftDroppedP4[5][kMaxParticleFlowJet02];
   Int_t           ParticleFlowJet02_NSubJetsTrimmed[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Int_t           ParticleFlowJet02_NSubJetsPruned[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Int_t           ParticleFlowJet02_NSubJetsSoftDropped[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Double_t        ParticleFlowJet02_ExclYmerge23[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Double_t        ParticleFlowJet02_ExclYmerge34[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Double_t        ParticleFlowJet02_ExclYmerge45[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   Double_t        ParticleFlowJet02_ExclYmerge56[kMaxParticleFlowJet02];   //[ParticleFlowJet02_]
   TRefArray       ParticleFlowJet02_Constituents[kMaxParticleFlowJet02];
   TRefArray       ParticleFlowJet02_Particles[kMaxParticleFlowJet02];
   TLorentzVector  ParticleFlowJet02_Area[kMaxParticleFlowJet02];
   Int_t           ParticleFlowJet02_size;
   Int_t           ParticleFlowJet04_;
   UInt_t          ParticleFlowJet04_fUniqueID[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   UInt_t          ParticleFlowJet04_fBits[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_PT[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_Eta[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_Phi[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_T[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_Mass[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_DeltaEta[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_DeltaPhi[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   UInt_t          ParticleFlowJet04_Flavor[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   UInt_t          ParticleFlowJet04_FlavorAlgo[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   UInt_t          ParticleFlowJet04_FlavorPhys[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   UInt_t          ParticleFlowJet04_BTag[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   UInt_t          ParticleFlowJet04_BTagAlgo[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   UInt_t          ParticleFlowJet04_BTagPhys[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   UInt_t          ParticleFlowJet04_TauTag[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_TauWeight[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Int_t           ParticleFlowJet04_Charge[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_EhadOverEem[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Int_t           ParticleFlowJet04_NCharged[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Int_t           ParticleFlowJet04_NNeutrals[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_NeutralEnergyFraction[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_ChargedEnergyFraction[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_Beta[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_BetaStar[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_MeanSqDeltaR[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_PTD[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_FracPt[kMaxParticleFlowJet04][5];   //[ParticleFlowJet04_]
   Float_t         ParticleFlowJet04_Tau[kMaxParticleFlowJet04][5];   //[ParticleFlowJet04_]
   TLorentzVector  ParticleFlowJet04_SoftDroppedJet[kMaxParticleFlowJet04];
   TLorentzVector  ParticleFlowJet04_SoftDroppedSubJet1[kMaxParticleFlowJet04];
   TLorentzVector  ParticleFlowJet04_SoftDroppedSubJet2[kMaxParticleFlowJet04];
   TLorentzVector  ParticleFlowJet04_TrimmedP4[5][kMaxParticleFlowJet04];
   TLorentzVector  ParticleFlowJet04_PrunedP4[5][kMaxParticleFlowJet04];
   TLorentzVector  ParticleFlowJet04_SoftDroppedP4[5][kMaxParticleFlowJet04];
   Int_t           ParticleFlowJet04_NSubJetsTrimmed[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Int_t           ParticleFlowJet04_NSubJetsPruned[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Int_t           ParticleFlowJet04_NSubJetsSoftDropped[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Double_t        ParticleFlowJet04_ExclYmerge23[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Double_t        ParticleFlowJet04_ExclYmerge34[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Double_t        ParticleFlowJet04_ExclYmerge45[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   Double_t        ParticleFlowJet04_ExclYmerge56[kMaxParticleFlowJet04];   //[ParticleFlowJet04_]
   TRefArray       ParticleFlowJet04_Constituents[kMaxParticleFlowJet04];
   TRefArray       ParticleFlowJet04_Particles[kMaxParticleFlowJet04];
   TLorentzVector  ParticleFlowJet04_Area[kMaxParticleFlowJet04];
   Int_t           ParticleFlowJet04_size;
   Int_t           ParticleFlowJet08_;
   UInt_t          ParticleFlowJet08_fUniqueID[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   UInt_t          ParticleFlowJet08_fBits[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_PT[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_Eta[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_Phi[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_T[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_Mass[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_DeltaEta[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_DeltaPhi[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   UInt_t          ParticleFlowJet08_Flavor[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   UInt_t          ParticleFlowJet08_FlavorAlgo[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   UInt_t          ParticleFlowJet08_FlavorPhys[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   UInt_t          ParticleFlowJet08_BTag[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   UInt_t          ParticleFlowJet08_BTagAlgo[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   UInt_t          ParticleFlowJet08_BTagPhys[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   UInt_t          ParticleFlowJet08_TauTag[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_TauWeight[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Int_t           ParticleFlowJet08_Charge[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_EhadOverEem[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Int_t           ParticleFlowJet08_NCharged[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Int_t           ParticleFlowJet08_NNeutrals[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_NeutralEnergyFraction[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_ChargedEnergyFraction[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_Beta[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_BetaStar[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_MeanSqDeltaR[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_PTD[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_FracPt[kMaxParticleFlowJet08][5];   //[ParticleFlowJet08_]
   Float_t         ParticleFlowJet08_Tau[kMaxParticleFlowJet08][5];   //[ParticleFlowJet08_]
   TLorentzVector  ParticleFlowJet08_SoftDroppedJet[kMaxParticleFlowJet08];
   TLorentzVector  ParticleFlowJet08_SoftDroppedSubJet1[kMaxParticleFlowJet08];
   TLorentzVector  ParticleFlowJet08_SoftDroppedSubJet2[kMaxParticleFlowJet08];
   TLorentzVector  ParticleFlowJet08_TrimmedP4[5][kMaxParticleFlowJet08];
   TLorentzVector  ParticleFlowJet08_PrunedP4[5][kMaxParticleFlowJet08];
   TLorentzVector  ParticleFlowJet08_SoftDroppedP4[5][kMaxParticleFlowJet08];
   Int_t           ParticleFlowJet08_NSubJetsTrimmed[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Int_t           ParticleFlowJet08_NSubJetsPruned[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Int_t           ParticleFlowJet08_NSubJetsSoftDropped[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Double_t        ParticleFlowJet08_ExclYmerge23[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Double_t        ParticleFlowJet08_ExclYmerge34[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Double_t        ParticleFlowJet08_ExclYmerge45[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   Double_t        ParticleFlowJet08_ExclYmerge56[kMaxParticleFlowJet08];   //[ParticleFlowJet08_]
   TRefArray       ParticleFlowJet08_Constituents[kMaxParticleFlowJet08];
   TRefArray       ParticleFlowJet08_Particles[kMaxParticleFlowJet08];
   TLorentzVector  ParticleFlowJet08_Area[kMaxParticleFlowJet08];
   Int_t           ParticleFlowJet08_size;
   Int_t           ParticleFlowJet15_;
   UInt_t          ParticleFlowJet15_fUniqueID[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   UInt_t          ParticleFlowJet15_fBits[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_PT[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_Eta[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_Phi[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_T[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_Mass[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_DeltaEta[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_DeltaPhi[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   UInt_t          ParticleFlowJet15_Flavor[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   UInt_t          ParticleFlowJet15_FlavorAlgo[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   UInt_t          ParticleFlowJet15_FlavorPhys[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   UInt_t          ParticleFlowJet15_BTag[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   UInt_t          ParticleFlowJet15_BTagAlgo[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   UInt_t          ParticleFlowJet15_BTagPhys[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   UInt_t          ParticleFlowJet15_TauTag[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_TauWeight[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Int_t           ParticleFlowJet15_Charge[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_EhadOverEem[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Int_t           ParticleFlowJet15_NCharged[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Int_t           ParticleFlowJet15_NNeutrals[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_NeutralEnergyFraction[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_ChargedEnergyFraction[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_Beta[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_BetaStar[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_MeanSqDeltaR[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_PTD[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_FracPt[kMaxParticleFlowJet15][5];   //[ParticleFlowJet15_]
   Float_t         ParticleFlowJet15_Tau[kMaxParticleFlowJet15][5];   //[ParticleFlowJet15_]
   TLorentzVector  ParticleFlowJet15_SoftDroppedJet[kMaxParticleFlowJet15];
   TLorentzVector  ParticleFlowJet15_SoftDroppedSubJet1[kMaxParticleFlowJet15];
   TLorentzVector  ParticleFlowJet15_SoftDroppedSubJet2[kMaxParticleFlowJet15];
   TLorentzVector  ParticleFlowJet15_TrimmedP4[5][kMaxParticleFlowJet15];
   TLorentzVector  ParticleFlowJet15_PrunedP4[5][kMaxParticleFlowJet15];
   TLorentzVector  ParticleFlowJet15_SoftDroppedP4[5][kMaxParticleFlowJet15];
   Int_t           ParticleFlowJet15_NSubJetsTrimmed[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Int_t           ParticleFlowJet15_NSubJetsPruned[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Int_t           ParticleFlowJet15_NSubJetsSoftDropped[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Double_t        ParticleFlowJet15_ExclYmerge23[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Double_t        ParticleFlowJet15_ExclYmerge34[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Double_t        ParticleFlowJet15_ExclYmerge45[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   Double_t        ParticleFlowJet15_ExclYmerge56[kMaxParticleFlowJet15];   //[ParticleFlowJet15_]
   TRefArray       ParticleFlowJet15_Constituents[kMaxParticleFlowJet15];
   TRefArray       ParticleFlowJet15_Particles[kMaxParticleFlowJet15];
   TLorentzVector  ParticleFlowJet15_Area[kMaxParticleFlowJet15];
   Int_t           ParticleFlowJet15_size;
   Int_t           CaloJet02_;
   UInt_t          CaloJet02_fUniqueID[kMaxCaloJet02];   //[CaloJet02_]
   UInt_t          CaloJet02_fBits[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_PT[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_Eta[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_Phi[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_T[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_Mass[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_DeltaEta[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_DeltaPhi[kMaxCaloJet02];   //[CaloJet02_]
   UInt_t          CaloJet02_Flavor[kMaxCaloJet02];   //[CaloJet02_]
   UInt_t          CaloJet02_FlavorAlgo[kMaxCaloJet02];   //[CaloJet02_]
   UInt_t          CaloJet02_FlavorPhys[kMaxCaloJet02];   //[CaloJet02_]
   UInt_t          CaloJet02_BTag[kMaxCaloJet02];   //[CaloJet02_]
   UInt_t          CaloJet02_BTagAlgo[kMaxCaloJet02];   //[CaloJet02_]
   UInt_t          CaloJet02_BTagPhys[kMaxCaloJet02];   //[CaloJet02_]
   UInt_t          CaloJet02_TauTag[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_TauWeight[kMaxCaloJet02];   //[CaloJet02_]
   Int_t           CaloJet02_Charge[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_EhadOverEem[kMaxCaloJet02];   //[CaloJet02_]
   Int_t           CaloJet02_NCharged[kMaxCaloJet02];   //[CaloJet02_]
   Int_t           CaloJet02_NNeutrals[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_NeutralEnergyFraction[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_ChargedEnergyFraction[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_Beta[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_BetaStar[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_MeanSqDeltaR[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_PTD[kMaxCaloJet02];   //[CaloJet02_]
   Float_t         CaloJet02_FracPt[kMaxCaloJet02][5];   //[CaloJet02_]
   Float_t         CaloJet02_Tau[kMaxCaloJet02][5];   //[CaloJet02_]
   TLorentzVector  CaloJet02_SoftDroppedJet[kMaxCaloJet02];
   TLorentzVector  CaloJet02_SoftDroppedSubJet1[kMaxCaloJet02];
   TLorentzVector  CaloJet02_SoftDroppedSubJet2[kMaxCaloJet02];
   TLorentzVector  CaloJet02_TrimmedP4[5][kMaxCaloJet02];
   TLorentzVector  CaloJet02_PrunedP4[5][kMaxCaloJet02];
   TLorentzVector  CaloJet02_SoftDroppedP4[5][kMaxCaloJet02];
   Int_t           CaloJet02_NSubJetsTrimmed[kMaxCaloJet02];   //[CaloJet02_]
   Int_t           CaloJet02_NSubJetsPruned[kMaxCaloJet02];   //[CaloJet02_]
   Int_t           CaloJet02_NSubJetsSoftDropped[kMaxCaloJet02];   //[CaloJet02_]
   Double_t        CaloJet02_ExclYmerge23[kMaxCaloJet02];   //[CaloJet02_]
   Double_t        CaloJet02_ExclYmerge34[kMaxCaloJet02];   //[CaloJet02_]
   Double_t        CaloJet02_ExclYmerge45[kMaxCaloJet02];   //[CaloJet02_]
   Double_t        CaloJet02_ExclYmerge56[kMaxCaloJet02];   //[CaloJet02_]
   TRefArray       CaloJet02_Constituents[kMaxCaloJet02];
   TRefArray       CaloJet02_Particles[kMaxCaloJet02];
   TLorentzVector  CaloJet02_Area[kMaxCaloJet02];
   Int_t           CaloJet02_size;
   Int_t           CaloJet04_;
   UInt_t          CaloJet04_fUniqueID[kMaxCaloJet04];   //[CaloJet04_]
   UInt_t          CaloJet04_fBits[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_PT[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_Eta[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_Phi[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_T[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_Mass[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_DeltaEta[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_DeltaPhi[kMaxCaloJet04];   //[CaloJet04_]
   UInt_t          CaloJet04_Flavor[kMaxCaloJet04];   //[CaloJet04_]
   UInt_t          CaloJet04_FlavorAlgo[kMaxCaloJet04];   //[CaloJet04_]
   UInt_t          CaloJet04_FlavorPhys[kMaxCaloJet04];   //[CaloJet04_]
   UInt_t          CaloJet04_BTag[kMaxCaloJet04];   //[CaloJet04_]
   UInt_t          CaloJet04_BTagAlgo[kMaxCaloJet04];   //[CaloJet04_]
   UInt_t          CaloJet04_BTagPhys[kMaxCaloJet04];   //[CaloJet04_]
   UInt_t          CaloJet04_TauTag[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_TauWeight[kMaxCaloJet04];   //[CaloJet04_]
   Int_t           CaloJet04_Charge[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_EhadOverEem[kMaxCaloJet04];   //[CaloJet04_]
   Int_t           CaloJet04_NCharged[kMaxCaloJet04];   //[CaloJet04_]
   Int_t           CaloJet04_NNeutrals[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_NeutralEnergyFraction[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_ChargedEnergyFraction[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_Beta[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_BetaStar[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_MeanSqDeltaR[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_PTD[kMaxCaloJet04];   //[CaloJet04_]
   Float_t         CaloJet04_FracPt[kMaxCaloJet04][5];   //[CaloJet04_]
   Float_t         CaloJet04_Tau[kMaxCaloJet04][5];   //[CaloJet04_]
   TLorentzVector  CaloJet04_SoftDroppedJet[kMaxCaloJet04];
   TLorentzVector  CaloJet04_SoftDroppedSubJet1[kMaxCaloJet04];
   TLorentzVector  CaloJet04_SoftDroppedSubJet2[kMaxCaloJet04];
   TLorentzVector  CaloJet04_TrimmedP4[5][kMaxCaloJet04];
   TLorentzVector  CaloJet04_PrunedP4[5][kMaxCaloJet04];
   TLorentzVector  CaloJet04_SoftDroppedP4[5][kMaxCaloJet04];
   Int_t           CaloJet04_NSubJetsTrimmed[kMaxCaloJet04];   //[CaloJet04_]
   Int_t           CaloJet04_NSubJetsPruned[kMaxCaloJet04];   //[CaloJet04_]
   Int_t           CaloJet04_NSubJetsSoftDropped[kMaxCaloJet04];   //[CaloJet04_]
   Double_t        CaloJet04_ExclYmerge23[kMaxCaloJet04];   //[CaloJet04_]
   Double_t        CaloJet04_ExclYmerge34[kMaxCaloJet04];   //[CaloJet04_]
   Double_t        CaloJet04_ExclYmerge45[kMaxCaloJet04];   //[CaloJet04_]
   Double_t        CaloJet04_ExclYmerge56[kMaxCaloJet04];   //[CaloJet04_]
   TRefArray       CaloJet04_Constituents[kMaxCaloJet04];
   TRefArray       CaloJet04_Particles[kMaxCaloJet04];
   TLorentzVector  CaloJet04_Area[kMaxCaloJet04];
   Int_t           CaloJet04_size;
   Int_t           CaloJet08_;
   UInt_t          CaloJet08_fUniqueID[kMaxCaloJet08];   //[CaloJet08_]
   UInt_t          CaloJet08_fBits[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_PT[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_Eta[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_Phi[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_T[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_Mass[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_DeltaEta[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_DeltaPhi[kMaxCaloJet08];   //[CaloJet08_]
   UInt_t          CaloJet08_Flavor[kMaxCaloJet08];   //[CaloJet08_]
   UInt_t          CaloJet08_FlavorAlgo[kMaxCaloJet08];   //[CaloJet08_]
   UInt_t          CaloJet08_FlavorPhys[kMaxCaloJet08];   //[CaloJet08_]
   UInt_t          CaloJet08_BTag[kMaxCaloJet08];   //[CaloJet08_]
   UInt_t          CaloJet08_BTagAlgo[kMaxCaloJet08];   //[CaloJet08_]
   UInt_t          CaloJet08_BTagPhys[kMaxCaloJet08];   //[CaloJet08_]
   UInt_t          CaloJet08_TauTag[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_TauWeight[kMaxCaloJet08];   //[CaloJet08_]
   Int_t           CaloJet08_Charge[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_EhadOverEem[kMaxCaloJet08];   //[CaloJet08_]
   Int_t           CaloJet08_NCharged[kMaxCaloJet08];   //[CaloJet08_]
   Int_t           CaloJet08_NNeutrals[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_NeutralEnergyFraction[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_ChargedEnergyFraction[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_Beta[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_BetaStar[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_MeanSqDeltaR[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_PTD[kMaxCaloJet08];   //[CaloJet08_]
   Float_t         CaloJet08_FracPt[kMaxCaloJet08][5];   //[CaloJet08_]
   Float_t         CaloJet08_Tau[kMaxCaloJet08][5];   //[CaloJet08_]
   TLorentzVector  CaloJet08_SoftDroppedJet[kMaxCaloJet08];
   TLorentzVector  CaloJet08_SoftDroppedSubJet1[kMaxCaloJet08];
   TLorentzVector  CaloJet08_SoftDroppedSubJet2[kMaxCaloJet08];
   TLorentzVector  CaloJet08_TrimmedP4[5][kMaxCaloJet08];
   TLorentzVector  CaloJet08_PrunedP4[5][kMaxCaloJet08];
   TLorentzVector  CaloJet08_SoftDroppedP4[5][kMaxCaloJet08];
   Int_t           CaloJet08_NSubJetsTrimmed[kMaxCaloJet08];   //[CaloJet08_]
   Int_t           CaloJet08_NSubJetsPruned[kMaxCaloJet08];   //[CaloJet08_]
   Int_t           CaloJet08_NSubJetsSoftDropped[kMaxCaloJet08];   //[CaloJet08_]
   Double_t        CaloJet08_ExclYmerge23[kMaxCaloJet08];   //[CaloJet08_]
   Double_t        CaloJet08_ExclYmerge34[kMaxCaloJet08];   //[CaloJet08_]
   Double_t        CaloJet08_ExclYmerge45[kMaxCaloJet08];   //[CaloJet08_]
   Double_t        CaloJet08_ExclYmerge56[kMaxCaloJet08];   //[CaloJet08_]
   TRefArray       CaloJet08_Constituents[kMaxCaloJet08];
   TRefArray       CaloJet08_Particles[kMaxCaloJet08];
   TLorentzVector  CaloJet08_Area[kMaxCaloJet08];
   Int_t           CaloJet08_size;
   Int_t           CaloJet15_;
   UInt_t          CaloJet15_fUniqueID[kMaxCaloJet15];   //[CaloJet15_]
   UInt_t          CaloJet15_fBits[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_PT[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_Eta[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_Phi[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_T[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_Mass[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_DeltaEta[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_DeltaPhi[kMaxCaloJet15];   //[CaloJet15_]
   UInt_t          CaloJet15_Flavor[kMaxCaloJet15];   //[CaloJet15_]
   UInt_t          CaloJet15_FlavorAlgo[kMaxCaloJet15];   //[CaloJet15_]
   UInt_t          CaloJet15_FlavorPhys[kMaxCaloJet15];   //[CaloJet15_]
   UInt_t          CaloJet15_BTag[kMaxCaloJet15];   //[CaloJet15_]
   UInt_t          CaloJet15_BTagAlgo[kMaxCaloJet15];   //[CaloJet15_]
   UInt_t          CaloJet15_BTagPhys[kMaxCaloJet15];   //[CaloJet15_]
   UInt_t          CaloJet15_TauTag[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_TauWeight[kMaxCaloJet15];   //[CaloJet15_]
   Int_t           CaloJet15_Charge[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_EhadOverEem[kMaxCaloJet15];   //[CaloJet15_]
   Int_t           CaloJet15_NCharged[kMaxCaloJet15];   //[CaloJet15_]
   Int_t           CaloJet15_NNeutrals[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_NeutralEnergyFraction[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_ChargedEnergyFraction[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_Beta[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_BetaStar[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_MeanSqDeltaR[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_PTD[kMaxCaloJet15];   //[CaloJet15_]
   Float_t         CaloJet15_FracPt[kMaxCaloJet15][5];   //[CaloJet15_]
   Float_t         CaloJet15_Tau[kMaxCaloJet15][5];   //[CaloJet15_]
   TLorentzVector  CaloJet15_SoftDroppedJet[kMaxCaloJet15];
   TLorentzVector  CaloJet15_SoftDroppedSubJet1[kMaxCaloJet15];
   TLorentzVector  CaloJet15_SoftDroppedSubJet2[kMaxCaloJet15];
   TLorentzVector  CaloJet15_TrimmedP4[5][kMaxCaloJet15];
   TLorentzVector  CaloJet15_PrunedP4[5][kMaxCaloJet15];
   TLorentzVector  CaloJet15_SoftDroppedP4[5][kMaxCaloJet15];
   Int_t           CaloJet15_NSubJetsTrimmed[kMaxCaloJet15];   //[CaloJet15_]
   Int_t           CaloJet15_NSubJetsPruned[kMaxCaloJet15];   //[CaloJet15_]
   Int_t           CaloJet15_NSubJetsSoftDropped[kMaxCaloJet15];   //[CaloJet15_]
   Double_t        CaloJet15_ExclYmerge23[kMaxCaloJet15];   //[CaloJet15_]
   Double_t        CaloJet15_ExclYmerge34[kMaxCaloJet15];   //[CaloJet15_]
   Double_t        CaloJet15_ExclYmerge45[kMaxCaloJet15];   //[CaloJet15_]
   Double_t        CaloJet15_ExclYmerge56[kMaxCaloJet15];   //[CaloJet15_]
   TRefArray       CaloJet15_Constituents[kMaxCaloJet15];
   TRefArray       CaloJet15_Particles[kMaxCaloJet15];
   TLorentzVector  CaloJet15_Area[kMaxCaloJet15];
   Int_t           CaloJet15_size;
   Int_t           TrackJet02_;
   UInt_t          TrackJet02_fUniqueID[kMaxTrackJet02];   //[TrackJet02_]
   UInt_t          TrackJet02_fBits[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_PT[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_Eta[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_Phi[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_T[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_Mass[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_DeltaEta[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_DeltaPhi[kMaxTrackJet02];   //[TrackJet02_]
   UInt_t          TrackJet02_Flavor[kMaxTrackJet02];   //[TrackJet02_]
   UInt_t          TrackJet02_FlavorAlgo[kMaxTrackJet02];   //[TrackJet02_]
   UInt_t          TrackJet02_FlavorPhys[kMaxTrackJet02];   //[TrackJet02_]
   UInt_t          TrackJet02_BTag[kMaxTrackJet02];   //[TrackJet02_]
   UInt_t          TrackJet02_BTagAlgo[kMaxTrackJet02];   //[TrackJet02_]
   UInt_t          TrackJet02_BTagPhys[kMaxTrackJet02];   //[TrackJet02_]
   UInt_t          TrackJet02_TauTag[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_TauWeight[kMaxTrackJet02];   //[TrackJet02_]
   Int_t           TrackJet02_Charge[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_EhadOverEem[kMaxTrackJet02];   //[TrackJet02_]
   Int_t           TrackJet02_NCharged[kMaxTrackJet02];   //[TrackJet02_]
   Int_t           TrackJet02_NNeutrals[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_NeutralEnergyFraction[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_ChargedEnergyFraction[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_Beta[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_BetaStar[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_MeanSqDeltaR[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_PTD[kMaxTrackJet02];   //[TrackJet02_]
   Float_t         TrackJet02_FracPt[kMaxTrackJet02][5];   //[TrackJet02_]
   Float_t         TrackJet02_Tau[kMaxTrackJet02][5];   //[TrackJet02_]
   TLorentzVector  TrackJet02_SoftDroppedJet[kMaxTrackJet02];
   TLorentzVector  TrackJet02_SoftDroppedSubJet1[kMaxTrackJet02];
   TLorentzVector  TrackJet02_SoftDroppedSubJet2[kMaxTrackJet02];
   TLorentzVector  TrackJet02_TrimmedP4[5][kMaxTrackJet02];
   TLorentzVector  TrackJet02_PrunedP4[5][kMaxTrackJet02];
   TLorentzVector  TrackJet02_SoftDroppedP4[5][kMaxTrackJet02];
   Int_t           TrackJet02_NSubJetsTrimmed[kMaxTrackJet02];   //[TrackJet02_]
   Int_t           TrackJet02_NSubJetsPruned[kMaxTrackJet02];   //[TrackJet02_]
   Int_t           TrackJet02_NSubJetsSoftDropped[kMaxTrackJet02];   //[TrackJet02_]
   Double_t        TrackJet02_ExclYmerge23[kMaxTrackJet02];   //[TrackJet02_]
   Double_t        TrackJet02_ExclYmerge34[kMaxTrackJet02];   //[TrackJet02_]
   Double_t        TrackJet02_ExclYmerge45[kMaxTrackJet02];   //[TrackJet02_]
   Double_t        TrackJet02_ExclYmerge56[kMaxTrackJet02];   //[TrackJet02_]
   TRefArray       TrackJet02_Constituents[kMaxTrackJet02];
   TRefArray       TrackJet02_Particles[kMaxTrackJet02];
   TLorentzVector  TrackJet02_Area[kMaxTrackJet02];
   Int_t           TrackJet02_size;
   Int_t           TrackJet04_;
   UInt_t          TrackJet04_fUniqueID[kMaxTrackJet04];   //[TrackJet04_]
   UInt_t          TrackJet04_fBits[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_PT[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_Eta[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_Phi[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_T[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_Mass[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_DeltaEta[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_DeltaPhi[kMaxTrackJet04];   //[TrackJet04_]
   UInt_t          TrackJet04_Flavor[kMaxTrackJet04];   //[TrackJet04_]
   UInt_t          TrackJet04_FlavorAlgo[kMaxTrackJet04];   //[TrackJet04_]
   UInt_t          TrackJet04_FlavorPhys[kMaxTrackJet04];   //[TrackJet04_]
   UInt_t          TrackJet04_BTag[kMaxTrackJet04];   //[TrackJet04_]
   UInt_t          TrackJet04_BTagAlgo[kMaxTrackJet04];   //[TrackJet04_]
   UInt_t          TrackJet04_BTagPhys[kMaxTrackJet04];   //[TrackJet04_]
   UInt_t          TrackJet04_TauTag[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_TauWeight[kMaxTrackJet04];   //[TrackJet04_]
   Int_t           TrackJet04_Charge[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_EhadOverEem[kMaxTrackJet04];   //[TrackJet04_]
   Int_t           TrackJet04_NCharged[kMaxTrackJet04];   //[TrackJet04_]
   Int_t           TrackJet04_NNeutrals[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_NeutralEnergyFraction[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_ChargedEnergyFraction[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_Beta[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_BetaStar[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_MeanSqDeltaR[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_PTD[kMaxTrackJet04];   //[TrackJet04_]
   Float_t         TrackJet04_FracPt[kMaxTrackJet04][5];   //[TrackJet04_]
   Float_t         TrackJet04_Tau[kMaxTrackJet04][5];   //[TrackJet04_]
   TLorentzVector  TrackJet04_SoftDroppedJet[kMaxTrackJet04];
   TLorentzVector  TrackJet04_SoftDroppedSubJet1[kMaxTrackJet04];
   TLorentzVector  TrackJet04_SoftDroppedSubJet2[kMaxTrackJet04];
   TLorentzVector  TrackJet04_TrimmedP4[5][kMaxTrackJet04];
   TLorentzVector  TrackJet04_PrunedP4[5][kMaxTrackJet04];
   TLorentzVector  TrackJet04_SoftDroppedP4[5][kMaxTrackJet04];
   Int_t           TrackJet04_NSubJetsTrimmed[kMaxTrackJet04];   //[TrackJet04_]
   Int_t           TrackJet04_NSubJetsPruned[kMaxTrackJet04];   //[TrackJet04_]
   Int_t           TrackJet04_NSubJetsSoftDropped[kMaxTrackJet04];   //[TrackJet04_]
   Double_t        TrackJet04_ExclYmerge23[kMaxTrackJet04];   //[TrackJet04_]
   Double_t        TrackJet04_ExclYmerge34[kMaxTrackJet04];   //[TrackJet04_]
   Double_t        TrackJet04_ExclYmerge45[kMaxTrackJet04];   //[TrackJet04_]
   Double_t        TrackJet04_ExclYmerge56[kMaxTrackJet04];   //[TrackJet04_]
   TRefArray       TrackJet04_Constituents[kMaxTrackJet04];
   TRefArray       TrackJet04_Particles[kMaxTrackJet04];
   TLorentzVector  TrackJet04_Area[kMaxTrackJet04];
   Int_t           TrackJet04_size;
   Int_t           TrackJet08_;
   UInt_t          TrackJet08_fUniqueID[kMaxTrackJet08];   //[TrackJet08_]
   UInt_t          TrackJet08_fBits[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_PT[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_Eta[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_Phi[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_T[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_Mass[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_DeltaEta[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_DeltaPhi[kMaxTrackJet08];   //[TrackJet08_]
   UInt_t          TrackJet08_Flavor[kMaxTrackJet08];   //[TrackJet08_]
   UInt_t          TrackJet08_FlavorAlgo[kMaxTrackJet08];   //[TrackJet08_]
   UInt_t          TrackJet08_FlavorPhys[kMaxTrackJet08];   //[TrackJet08_]
   UInt_t          TrackJet08_BTag[kMaxTrackJet08];   //[TrackJet08_]
   UInt_t          TrackJet08_BTagAlgo[kMaxTrackJet08];   //[TrackJet08_]
   UInt_t          TrackJet08_BTagPhys[kMaxTrackJet08];   //[TrackJet08_]
   UInt_t          TrackJet08_TauTag[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_TauWeight[kMaxTrackJet08];   //[TrackJet08_]
   Int_t           TrackJet08_Charge[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_EhadOverEem[kMaxTrackJet08];   //[TrackJet08_]
   Int_t           TrackJet08_NCharged[kMaxTrackJet08];   //[TrackJet08_]
   Int_t           TrackJet08_NNeutrals[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_NeutralEnergyFraction[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_ChargedEnergyFraction[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_Beta[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_BetaStar[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_MeanSqDeltaR[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_PTD[kMaxTrackJet08];   //[TrackJet08_]
   Float_t         TrackJet08_FracPt[kMaxTrackJet08][5];   //[TrackJet08_]
   Float_t         TrackJet08_Tau[kMaxTrackJet08][5];   //[TrackJet08_]
   TLorentzVector  TrackJet08_SoftDroppedJet[kMaxTrackJet08];
   TLorentzVector  TrackJet08_SoftDroppedSubJet1[kMaxTrackJet08];
   TLorentzVector  TrackJet08_SoftDroppedSubJet2[kMaxTrackJet08];
   TLorentzVector  TrackJet08_TrimmedP4[5][kMaxTrackJet08];
   TLorentzVector  TrackJet08_PrunedP4[5][kMaxTrackJet08];
   TLorentzVector  TrackJet08_SoftDroppedP4[5][kMaxTrackJet08];
   Int_t           TrackJet08_NSubJetsTrimmed[kMaxTrackJet08];   //[TrackJet08_]
   Int_t           TrackJet08_NSubJetsPruned[kMaxTrackJet08];   //[TrackJet08_]
   Int_t           TrackJet08_NSubJetsSoftDropped[kMaxTrackJet08];   //[TrackJet08_]
   Double_t        TrackJet08_ExclYmerge23[kMaxTrackJet08];   //[TrackJet08_]
   Double_t        TrackJet08_ExclYmerge34[kMaxTrackJet08];   //[TrackJet08_]
   Double_t        TrackJet08_ExclYmerge45[kMaxTrackJet08];   //[TrackJet08_]
   Double_t        TrackJet08_ExclYmerge56[kMaxTrackJet08];   //[TrackJet08_]
   TRefArray       TrackJet08_Constituents[kMaxTrackJet08];
   TRefArray       TrackJet08_Particles[kMaxTrackJet08];
   TLorentzVector  TrackJet08_Area[kMaxTrackJet08];
   Int_t           TrackJet08_size;
   Int_t           TrackJet15_;
   UInt_t          TrackJet15_fUniqueID[kMaxTrackJet15];   //[TrackJet15_]
   UInt_t          TrackJet15_fBits[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_PT[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_Eta[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_Phi[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_T[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_Mass[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_DeltaEta[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_DeltaPhi[kMaxTrackJet15];   //[TrackJet15_]
   UInt_t          TrackJet15_Flavor[kMaxTrackJet15];   //[TrackJet15_]
   UInt_t          TrackJet15_FlavorAlgo[kMaxTrackJet15];   //[TrackJet15_]
   UInt_t          TrackJet15_FlavorPhys[kMaxTrackJet15];   //[TrackJet15_]
   UInt_t          TrackJet15_BTag[kMaxTrackJet15];   //[TrackJet15_]
   UInt_t          TrackJet15_BTagAlgo[kMaxTrackJet15];   //[TrackJet15_]
   UInt_t          TrackJet15_BTagPhys[kMaxTrackJet15];   //[TrackJet15_]
   UInt_t          TrackJet15_TauTag[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_TauWeight[kMaxTrackJet15];   //[TrackJet15_]
   Int_t           TrackJet15_Charge[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_EhadOverEem[kMaxTrackJet15];   //[TrackJet15_]
   Int_t           TrackJet15_NCharged[kMaxTrackJet15];   //[TrackJet15_]
   Int_t           TrackJet15_NNeutrals[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_NeutralEnergyFraction[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_ChargedEnergyFraction[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_Beta[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_BetaStar[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_MeanSqDeltaR[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_PTD[kMaxTrackJet15];   //[TrackJet15_]
   Float_t         TrackJet15_FracPt[kMaxTrackJet15][5];   //[TrackJet15_]
   Float_t         TrackJet15_Tau[kMaxTrackJet15][5];   //[TrackJet15_]
   TLorentzVector  TrackJet15_SoftDroppedJet[kMaxTrackJet15];
   TLorentzVector  TrackJet15_SoftDroppedSubJet1[kMaxTrackJet15];
   TLorentzVector  TrackJet15_SoftDroppedSubJet2[kMaxTrackJet15];
   TLorentzVector  TrackJet15_TrimmedP4[5][kMaxTrackJet15];
   TLorentzVector  TrackJet15_PrunedP4[5][kMaxTrackJet15];
   TLorentzVector  TrackJet15_SoftDroppedP4[5][kMaxTrackJet15];
   Int_t           TrackJet15_NSubJetsTrimmed[kMaxTrackJet15];   //[TrackJet15_]
   Int_t           TrackJet15_NSubJetsPruned[kMaxTrackJet15];   //[TrackJet15_]
   Int_t           TrackJet15_NSubJetsSoftDropped[kMaxTrackJet15];   //[TrackJet15_]
   Double_t        TrackJet15_ExclYmerge23[kMaxTrackJet15];   //[TrackJet15_]
   Double_t        TrackJet15_ExclYmerge34[kMaxTrackJet15];   //[TrackJet15_]
   Double_t        TrackJet15_ExclYmerge45[kMaxTrackJet15];   //[TrackJet15_]
   Double_t        TrackJet15_ExclYmerge56[kMaxTrackJet15];   //[TrackJet15_]
   TRefArray       TrackJet15_Constituents[kMaxTrackJet15];
   TRefArray       TrackJet15_Particles[kMaxTrackJet15];
   TLorentzVector  TrackJet15_Area[kMaxTrackJet15];
   Int_t           TrackJet15_size;
   Int_t           MissingET_;
   UInt_t          MissingET_fUniqueID[kMaxMissingET];   //[MissingET_]
   UInt_t          MissingET_fBits[kMaxMissingET];   //[MissingET_]
   Float_t         MissingET_MET[kMaxMissingET];   //[MissingET_]
   Float_t         MissingET_Eta[kMaxMissingET];   //[MissingET_]
   Float_t         MissingET_Phi[kMaxMissingET];   //[MissingET_]
   Int_t           MissingET_size;
   Int_t           ScalarHT_;
   UInt_t          ScalarHT_fUniqueID[kMaxScalarHT];   //[ScalarHT_]
   UInt_t          ScalarHT_fBits[kMaxScalarHT];   //[ScalarHT_]
   Float_t         ScalarHT_HT[kMaxScalarHT];   //[ScalarHT_]
   Int_t           ScalarHT_size;

   // List of branches
   TBranch        *b_Event_;   //!
   TBranch        *b_Event_fUniqueID;   //!
   TBranch        *b_Event_fBits;   //!
   TBranch        *b_Event_Number;   //!
   TBranch        *b_Event_ReadTime;   //!
   TBranch        *b_Event_ProcTime;   //!
   TBranch        *b_Event_ProcessID;   //!
   TBranch        *b_Event_MPI;   //!
   TBranch        *b_Event_Weight;   //!
   TBranch        *b_Event_CrossSection;   //!
   TBranch        *b_Event_CrossSectionError;   //!
   TBranch        *b_Event_Scale;   //!
   TBranch        *b_Event_AlphaQED;   //!
   TBranch        *b_Event_AlphaQCD;   //!
   TBranch        *b_Event_ID1;   //!
   TBranch        *b_Event_ID2;   //!
   TBranch        *b_Event_X1;   //!
   TBranch        *b_Event_X2;   //!
   TBranch        *b_Event_ScalePDF;   //!
   TBranch        *b_Event_PDF1;   //!
   TBranch        *b_Event_PDF2;   //!
   TBranch        *b_Event_size;   //!
   TBranch        *b_Weight_;   //!
   TBranch        *b_Weight_fUniqueID;   //!
   TBranch        *b_Weight_fBits;   //!
   TBranch        *b_Weight_Weight;   //!
   TBranch        *b_Weight_size;   //!
   TBranch        *b_Particle_;   //!
   TBranch        *b_Particle_fUniqueID;   //!
   TBranch        *b_Particle_fBits;   //!
   TBranch        *b_Particle_PID;   //!
   TBranch        *b_Particle_Status;   //!
   TBranch        *b_Particle_IsPU;   //!
   TBranch        *b_Particle_M1;   //!
   TBranch        *b_Particle_M2;   //!
   TBranch        *b_Particle_D1;   //!
   TBranch        *b_Particle_D2;   //!
   TBranch        *b_Particle_Charge;   //!
   TBranch        *b_Particle_Mass;   //!
   TBranch        *b_Particle_E;   //!
   TBranch        *b_Particle_Px;   //!
   TBranch        *b_Particle_Py;   //!
   TBranch        *b_Particle_Pz;   //!
   TBranch        *b_Particle_P;   //!
   TBranch        *b_Particle_PT;   //!
   TBranch        *b_Particle_Eta;   //!
   TBranch        *b_Particle_Phi;   //!
   TBranch        *b_Particle_Rapidity;   //!
   TBranch        *b_Particle_T;   //!
   TBranch        *b_Particle_X;   //!
   TBranch        *b_Particle_Y;   //!
   TBranch        *b_Particle_Z;   //!
   TBranch        *b_Particle_size;   //!
   TBranch        *b_GenMissingET_;   //!
   TBranch        *b_GenMissingET_fUniqueID;   //!
   TBranch        *b_GenMissingET_fBits;   //!
   TBranch        *b_GenMissingET_MET;   //!
   TBranch        *b_GenMissingET_Eta;   //!
   TBranch        *b_GenMissingET_Phi;   //!
   TBranch        *b_GenMissingET_size;   //!
   TBranch        *b_EFlowTrack_;   //!
   TBranch        *b_EFlowTrack_fUniqueID;   //!
   TBranch        *b_EFlowTrack_fBits;   //!
   TBranch        *b_EFlowTrack_PID;   //!
   TBranch        *b_EFlowTrack_Charge;   //!
   TBranch        *b_EFlowTrack_P;   //!
   TBranch        *b_EFlowTrack_PT;   //!
   TBranch        *b_EFlowTrack_Eta;   //!
   TBranch        *b_EFlowTrack_Phi;   //!
   TBranch        *b_EFlowTrack_CtgTheta;   //!
   TBranch        *b_EFlowTrack_C;   //!
   TBranch        *b_EFlowTrack_Mass;   //!
   TBranch        *b_EFlowTrack_EtaOuter;   //!
   TBranch        *b_EFlowTrack_PhiOuter;   //!
   TBranch        *b_EFlowTrack_T;   //!
   TBranch        *b_EFlowTrack_X;   //!
   TBranch        *b_EFlowTrack_Y;   //!
   TBranch        *b_EFlowTrack_Z;   //!
   TBranch        *b_EFlowTrack_TOuter;   //!
   TBranch        *b_EFlowTrack_XOuter;   //!
   TBranch        *b_EFlowTrack_YOuter;   //!
   TBranch        *b_EFlowTrack_ZOuter;   //!
   TBranch        *b_EFlowTrack_Xd;   //!
   TBranch        *b_EFlowTrack_Yd;   //!
   TBranch        *b_EFlowTrack_Zd;   //!
   TBranch        *b_EFlowTrack_L;   //!
   TBranch        *b_EFlowTrack_D0;   //!
   TBranch        *b_EFlowTrack_DZ;   //!
   TBranch        *b_EFlowTrack_Nclusters;   //!
   TBranch        *b_EFlowTrack_dNdx;   //!
   TBranch        *b_EFlowTrack_ErrorP;   //!
   TBranch        *b_EFlowTrack_ErrorPT;   //!
   TBranch        *b_EFlowTrack_ErrorPhi;   //!
   TBranch        *b_EFlowTrack_ErrorCtgTheta;   //!
   TBranch        *b_EFlowTrack_ErrorT;   //!
   TBranch        *b_EFlowTrack_ErrorD0;   //!
   TBranch        *b_EFlowTrack_ErrorDZ;   //!
   TBranch        *b_EFlowTrack_ErrorC;   //!
   TBranch        *b_EFlowTrack_ErrorD0Phi;   //!
   TBranch        *b_EFlowTrack_ErrorD0C;   //!
   TBranch        *b_EFlowTrack_ErrorD0DZ;   //!
   TBranch        *b_EFlowTrack_ErrorD0CtgTheta;   //!
   TBranch        *b_EFlowTrack_ErrorPhiC;   //!
   TBranch        *b_EFlowTrack_ErrorPhiDZ;   //!
   TBranch        *b_EFlowTrack_ErrorPhiCtgTheta;   //!
   TBranch        *b_EFlowTrack_ErrorCDZ;   //!
   TBranch        *b_EFlowTrack_ErrorCCtgTheta;   //!
   TBranch        *b_EFlowTrack_ErrorDZCtgTheta;   //!
   TBranch        *b_EFlowTrack_Particle;   //!
   TBranch        *b_EFlowTrack_VertexIndex;   //!
   TBranch        *b_EFlowTrack_size;   //!
   TBranch        *b_EFlowPhoton_;   //!
   TBranch        *b_EFlowPhoton_fUniqueID;   //!
   TBranch        *b_EFlowPhoton_fBits;   //!
   TBranch        *b_EFlowPhoton_ET;   //!
   TBranch        *b_EFlowPhoton_Eta;   //!
   TBranch        *b_EFlowPhoton_Phi;   //!
   TBranch        *b_EFlowPhoton_E;   //!
   TBranch        *b_EFlowPhoton_T;   //!
   TBranch        *b_EFlowPhoton_NTimeHits;   //!
   TBranch        *b_EFlowPhoton_Eem;   //!
   TBranch        *b_EFlowPhoton_Ehad;   //!
   TBranch        *b_EFlowPhoton_Edges;   //!
   TBranch        *b_EFlowPhoton_Particles;   //!
   TBranch        *b_EFlowPhoton_size;   //!
   TBranch        *b_EFlowNeutralHadron_;   //!
   TBranch        *b_EFlowNeutralHadron_fUniqueID;   //!
   TBranch        *b_EFlowNeutralHadron_fBits;   //!
   TBranch        *b_EFlowNeutralHadron_ET;   //!
   TBranch        *b_EFlowNeutralHadron_Eta;   //!
   TBranch        *b_EFlowNeutralHadron_Phi;   //!
   TBranch        *b_EFlowNeutralHadron_E;   //!
   TBranch        *b_EFlowNeutralHadron_T;   //!
   TBranch        *b_EFlowNeutralHadron_NTimeHits;   //!
   TBranch        *b_EFlowNeutralHadron_Eem;   //!
   TBranch        *b_EFlowNeutralHadron_Ehad;   //!
   TBranch        *b_EFlowNeutralHadron_Edges;   //!
   TBranch        *b_EFlowNeutralHadron_Particles;   //!
   TBranch        *b_EFlowNeutralHadron_size;   //!
   TBranch        *b_Photon_;   //!
   TBranch        *b_Photon_fUniqueID;   //!
   TBranch        *b_Photon_fBits;   //!
   TBranch        *b_Photon_PT;   //!
   TBranch        *b_Photon_Eta;   //!
   TBranch        *b_Photon_Phi;   //!
   TBranch        *b_Photon_E;   //!
   TBranch        *b_Photon_T;   //!
   TBranch        *b_Photon_EhadOverEem;   //!
   TBranch        *b_Photon_Particles;   //!
   TBranch        *b_Photon_IsolationVar;   //!
   TBranch        *b_Photon_IsolationVarRhoCorr;   //!
   TBranch        *b_Photon_SumPtCharged;   //!
   TBranch        *b_Photon_SumPtNeutral;   //!
   TBranch        *b_Photon_SumPtChargedPU;   //!
   TBranch        *b_Photon_SumPt;   //!
   TBranch        *b_Photon_Status;   //!
   TBranch        *b_Photon_size;   //!
   TBranch        *b_Electron_;   //!
   TBranch        *b_Electron_fUniqueID;   //!
   TBranch        *b_Electron_fBits;   //!
   TBranch        *b_Electron_PT;   //!
   TBranch        *b_Electron_Eta;   //!
   TBranch        *b_Electron_Phi;   //!
   TBranch        *b_Electron_T;   //!
   TBranch        *b_Electron_Charge;   //!
   TBranch        *b_Electron_EhadOverEem;   //!
   TBranch        *b_Electron_Particle;   //!
   TBranch        *b_Electron_IsolationVar;   //!
   TBranch        *b_Electron_IsolationVarRhoCorr;   //!
   TBranch        *b_Electron_SumPtCharged;   //!
   TBranch        *b_Electron_SumPtNeutral;   //!
   TBranch        *b_Electron_SumPtChargedPU;   //!
   TBranch        *b_Electron_SumPt;   //!
   TBranch        *b_Electron_D0;   //!
   TBranch        *b_Electron_DZ;   //!
   TBranch        *b_Electron_ErrorD0;   //!
   TBranch        *b_Electron_ErrorDZ;   //!
   TBranch        *b_Electron_size;   //!
   TBranch        *b_Muon_;   //!
   TBranch        *b_Muon_fUniqueID;   //!
   TBranch        *b_Muon_fBits;   //!
   TBranch        *b_Muon_PT;   //!
   TBranch        *b_Muon_Eta;   //!
   TBranch        *b_Muon_Phi;   //!
   TBranch        *b_Muon_T;   //!
   TBranch        *b_Muon_Charge;   //!
   TBranch        *b_Muon_Particle;   //!
   TBranch        *b_Muon_IsolationVar;   //!
   TBranch        *b_Muon_IsolationVarRhoCorr;   //!
   TBranch        *b_Muon_SumPtCharged;   //!
   TBranch        *b_Muon_SumPtNeutral;   //!
   TBranch        *b_Muon_SumPtChargedPU;   //!
   TBranch        *b_Muon_SumPt;   //!
   TBranch        *b_Muon_D0;   //!
   TBranch        *b_Muon_DZ;   //!
   TBranch        *b_Muon_ErrorD0;   //!
   TBranch        *b_Muon_ErrorDZ;   //!
   TBranch        *b_Muon_size;   //!
   TBranch        *b_Jet_;   //!
   TBranch        *b_Jet_fUniqueID;   //!
   TBranch        *b_Jet_fBits;   //!
   TBranch        *b_Jet_PT;   //!
   TBranch        *b_Jet_Eta;   //!
   TBranch        *b_Jet_Phi;   //!
   TBranch        *b_Jet_T;   //!
   TBranch        *b_Jet_Mass;   //!
   TBranch        *b_Jet_DeltaEta;   //!
   TBranch        *b_Jet_DeltaPhi;   //!
   TBranch        *b_Jet_Flavor;   //!
   TBranch        *b_Jet_FlavorAlgo;   //!
   TBranch        *b_Jet_FlavorPhys;   //!
   TBranch        *b_Jet_BTag;   //!
   TBranch        *b_Jet_BTagAlgo;   //!
   TBranch        *b_Jet_BTagPhys;   //!
   TBranch        *b_Jet_TauTag;   //!
   TBranch        *b_Jet_TauWeight;   //!
   TBranch        *b_Jet_Charge;   //!
   TBranch        *b_Jet_EhadOverEem;   //!
   TBranch        *b_Jet_NCharged;   //!
   TBranch        *b_Jet_NNeutrals;   //!
   TBranch        *b_Jet_NeutralEnergyFraction;   //!
   TBranch        *b_Jet_ChargedEnergyFraction;   //!
   TBranch        *b_Jet_Beta;   //!
   TBranch        *b_Jet_BetaStar;   //!
   TBranch        *b_Jet_MeanSqDeltaR;   //!
   TBranch        *b_Jet_PTD;   //!
   TBranch        *b_Jet_FracPt;   //!
   TBranch        *b_Jet_Tau;   //!
   TBranch        *b_Jet_SoftDroppedJet;   //!
   TBranch        *b_Jet_SoftDroppedSubJet1;   //!
   TBranch        *b_Jet_SoftDroppedSubJet2;   //!
   TBranch        *b_Jet_TrimmedP4;   //!
   TBranch        *b_Jet_PrunedP4;   //!
   TBranch        *b_Jet_SoftDroppedP4;   //!
   TBranch        *b_Jet_NSubJetsTrimmed;   //!
   TBranch        *b_Jet_NSubJetsPruned;   //!
   TBranch        *b_Jet_NSubJetsSoftDropped;   //!
   TBranch        *b_Jet_ExclYmerge23;   //!
   TBranch        *b_Jet_ExclYmerge34;   //!
   TBranch        *b_Jet_ExclYmerge45;   //!
   TBranch        *b_Jet_ExclYmerge56;   //!
   TBranch        *b_Jet_Constituents;   //!
   TBranch        *b_Jet_Particles;   //!
   TBranch        *b_Jet_Area;   //!
   TBranch        *b_Jet_size;   //!
   TBranch        *b_PhotonNoIso_;   //!
   TBranch        *b_PhotonNoIso_fUniqueID;   //!
   TBranch        *b_PhotonNoIso_fBits;   //!
   TBranch        *b_PhotonNoIso_PT;   //!
   TBranch        *b_PhotonNoIso_Eta;   //!
   TBranch        *b_PhotonNoIso_Phi;   //!
   TBranch        *b_PhotonNoIso_E;   //!
   TBranch        *b_PhotonNoIso_T;   //!
   TBranch        *b_PhotonNoIso_EhadOverEem;   //!
   TBranch        *b_PhotonNoIso_Particles;   //!
   TBranch        *b_PhotonNoIso_IsolationVar;   //!
   TBranch        *b_PhotonNoIso_IsolationVarRhoCorr;   //!
   TBranch        *b_PhotonNoIso_SumPtCharged;   //!
   TBranch        *b_PhotonNoIso_SumPtNeutral;   //!
   TBranch        *b_PhotonNoIso_SumPtChargedPU;   //!
   TBranch        *b_PhotonNoIso_SumPt;   //!
   TBranch        *b_PhotonNoIso_Status;   //!
   TBranch        *b_PhotonNoIso_size;   //!
   TBranch        *b_ElectronNoIso_;   //!
   TBranch        *b_ElectronNoIso_fUniqueID;   //!
   TBranch        *b_ElectronNoIso_fBits;   //!
   TBranch        *b_ElectronNoIso_PT;   //!
   TBranch        *b_ElectronNoIso_Eta;   //!
   TBranch        *b_ElectronNoIso_Phi;   //!
   TBranch        *b_ElectronNoIso_T;   //!
   TBranch        *b_ElectronNoIso_Charge;   //!
   TBranch        *b_ElectronNoIso_EhadOverEem;   //!
   TBranch        *b_ElectronNoIso_Particle;   //!
   TBranch        *b_ElectronNoIso_IsolationVar;   //!
   TBranch        *b_ElectronNoIso_IsolationVarRhoCorr;   //!
   TBranch        *b_ElectronNoIso_SumPtCharged;   //!
   TBranch        *b_ElectronNoIso_SumPtNeutral;   //!
   TBranch        *b_ElectronNoIso_SumPtChargedPU;   //!
   TBranch        *b_ElectronNoIso_SumPt;   //!
   TBranch        *b_ElectronNoIso_D0;   //!
   TBranch        *b_ElectronNoIso_DZ;   //!
   TBranch        *b_ElectronNoIso_ErrorD0;   //!
   TBranch        *b_ElectronNoIso_ErrorDZ;   //!
   TBranch        *b_ElectronNoIso_size;   //!
   TBranch        *b_MuonNoIso_;   //!
   TBranch        *b_MuonNoIso_fUniqueID;   //!
   TBranch        *b_MuonNoIso_fBits;   //!
   TBranch        *b_MuonNoIso_PT;   //!
   TBranch        *b_MuonNoIso_Eta;   //!
   TBranch        *b_MuonNoIso_Phi;   //!
   TBranch        *b_MuonNoIso_T;   //!
   TBranch        *b_MuonNoIso_Charge;   //!
   TBranch        *b_MuonNoIso_Particle;   //!
   TBranch        *b_MuonNoIso_IsolationVar;   //!
   TBranch        *b_MuonNoIso_IsolationVarRhoCorr;   //!
   TBranch        *b_MuonNoIso_SumPtCharged;   //!
   TBranch        *b_MuonNoIso_SumPtNeutral;   //!
   TBranch        *b_MuonNoIso_SumPtChargedPU;   //!
   TBranch        *b_MuonNoIso_SumPt;   //!
   TBranch        *b_MuonNoIso_D0;   //!
   TBranch        *b_MuonNoIso_DZ;   //!
   TBranch        *b_MuonNoIso_ErrorD0;   //!
   TBranch        *b_MuonNoIso_ErrorDZ;   //!
   TBranch        *b_MuonNoIso_size;   //!
   TBranch        *b_JetNoIso_;   //!
   TBranch        *b_JetNoIso_fUniqueID;   //!
   TBranch        *b_JetNoIso_fBits;   //!
   TBranch        *b_JetNoIso_PT;   //!
   TBranch        *b_JetNoIso_Eta;   //!
   TBranch        *b_JetNoIso_Phi;   //!
   TBranch        *b_JetNoIso_T;   //!
   TBranch        *b_JetNoIso_Mass;   //!
   TBranch        *b_JetNoIso_DeltaEta;   //!
   TBranch        *b_JetNoIso_DeltaPhi;   //!
   TBranch        *b_JetNoIso_Flavor;   //!
   TBranch        *b_JetNoIso_FlavorAlgo;   //!
   TBranch        *b_JetNoIso_FlavorPhys;   //!
   TBranch        *b_JetNoIso_BTag;   //!
   TBranch        *b_JetNoIso_BTagAlgo;   //!
   TBranch        *b_JetNoIso_BTagPhys;   //!
   TBranch        *b_JetNoIso_TauTag;   //!
   TBranch        *b_JetNoIso_TauWeight;   //!
   TBranch        *b_JetNoIso_Charge;   //!
   TBranch        *b_JetNoIso_EhadOverEem;   //!
   TBranch        *b_JetNoIso_NCharged;   //!
   TBranch        *b_JetNoIso_NNeutrals;   //!
   TBranch        *b_JetNoIso_NeutralEnergyFraction;   //!
   TBranch        *b_JetNoIso_ChargedEnergyFraction;   //!
   TBranch        *b_JetNoIso_Beta;   //!
   TBranch        *b_JetNoIso_BetaStar;   //!
   TBranch        *b_JetNoIso_MeanSqDeltaR;   //!
   TBranch        *b_JetNoIso_PTD;   //!
   TBranch        *b_JetNoIso_FracPt;   //!
   TBranch        *b_JetNoIso_Tau;   //!
   TBranch        *b_JetNoIso_SoftDroppedJet;   //!
   TBranch        *b_JetNoIso_SoftDroppedSubJet1;   //!
   TBranch        *b_JetNoIso_SoftDroppedSubJet2;   //!
   TBranch        *b_JetNoIso_TrimmedP4;   //!
   TBranch        *b_JetNoIso_PrunedP4;   //!
   TBranch        *b_JetNoIso_SoftDroppedP4;   //!
   TBranch        *b_JetNoIso_NSubJetsTrimmed;   //!
   TBranch        *b_JetNoIso_NSubJetsPruned;   //!
   TBranch        *b_JetNoIso_NSubJetsSoftDropped;   //!
   TBranch        *b_JetNoIso_ExclYmerge23;   //!
   TBranch        *b_JetNoIso_ExclYmerge34;   //!
   TBranch        *b_JetNoIso_ExclYmerge45;   //!
   TBranch        *b_JetNoIso_ExclYmerge56;   //!
   TBranch        *b_JetNoIso_Constituents;   //!
   TBranch        *b_JetNoIso_Particles;   //!
   TBranch        *b_JetNoIso_Area;   //!
   TBranch        *b_JetNoIso_size;   //!
   TBranch        *b_GenJet02_;   //!
   TBranch        *b_GenJet02_fUniqueID;   //!
   TBranch        *b_GenJet02_fBits;   //!
   TBranch        *b_GenJet02_PT;   //!
   TBranch        *b_GenJet02_Eta;   //!
   TBranch        *b_GenJet02_Phi;   //!
   TBranch        *b_GenJet02_T;   //!
   TBranch        *b_GenJet02_Mass;   //!
   TBranch        *b_GenJet02_DeltaEta;   //!
   TBranch        *b_GenJet02_DeltaPhi;   //!
   TBranch        *b_GenJet02_Flavor;   //!
   TBranch        *b_GenJet02_FlavorAlgo;   //!
   TBranch        *b_GenJet02_FlavorPhys;   //!
   TBranch        *b_GenJet02_BTag;   //!
   TBranch        *b_GenJet02_BTagAlgo;   //!
   TBranch        *b_GenJet02_BTagPhys;   //!
   TBranch        *b_GenJet02_TauTag;   //!
   TBranch        *b_GenJet02_TauWeight;   //!
   TBranch        *b_GenJet02_Charge;   //!
   TBranch        *b_GenJet02_EhadOverEem;   //!
   TBranch        *b_GenJet02_NCharged;   //!
   TBranch        *b_GenJet02_NNeutrals;   //!
   TBranch        *b_GenJet02_NeutralEnergyFraction;   //!
   TBranch        *b_GenJet02_ChargedEnergyFraction;   //!
   TBranch        *b_GenJet02_Beta;   //!
   TBranch        *b_GenJet02_BetaStar;   //!
   TBranch        *b_GenJet02_MeanSqDeltaR;   //!
   TBranch        *b_GenJet02_PTD;   //!
   TBranch        *b_GenJet02_FracPt;   //!
   TBranch        *b_GenJet02_Tau;   //!
   TBranch        *b_GenJet02_SoftDroppedJet;   //!
   TBranch        *b_GenJet02_SoftDroppedSubJet1;   //!
   TBranch        *b_GenJet02_SoftDroppedSubJet2;   //!
   TBranch        *b_GenJet02_TrimmedP4;   //!
   TBranch        *b_GenJet02_PrunedP4;   //!
   TBranch        *b_GenJet02_SoftDroppedP4;   //!
   TBranch        *b_GenJet02_NSubJetsTrimmed;   //!
   TBranch        *b_GenJet02_NSubJetsPruned;   //!
   TBranch        *b_GenJet02_NSubJetsSoftDropped;   //!
   TBranch        *b_GenJet02_ExclYmerge23;   //!
   TBranch        *b_GenJet02_ExclYmerge34;   //!
   TBranch        *b_GenJet02_ExclYmerge45;   //!
   TBranch        *b_GenJet02_ExclYmerge56;   //!
   TBranch        *b_GenJet02_Constituents;   //!
   TBranch        *b_GenJet02_Particles;   //!
   TBranch        *b_GenJet02_Area;   //!
   TBranch        *b_GenJet02_size;   //!
   TBranch        *b_GenJet04_;   //!
   TBranch        *b_GenJet04_fUniqueID;   //!
   TBranch        *b_GenJet04_fBits;   //!
   TBranch        *b_GenJet04_PT;   //!
   TBranch        *b_GenJet04_Eta;   //!
   TBranch        *b_GenJet04_Phi;   //!
   TBranch        *b_GenJet04_T;   //!
   TBranch        *b_GenJet04_Mass;   //!
   TBranch        *b_GenJet04_DeltaEta;   //!
   TBranch        *b_GenJet04_DeltaPhi;   //!
   TBranch        *b_GenJet04_Flavor;   //!
   TBranch        *b_GenJet04_FlavorAlgo;   //!
   TBranch        *b_GenJet04_FlavorPhys;   //!
   TBranch        *b_GenJet04_BTag;   //!
   TBranch        *b_GenJet04_BTagAlgo;   //!
   TBranch        *b_GenJet04_BTagPhys;   //!
   TBranch        *b_GenJet04_TauTag;   //!
   TBranch        *b_GenJet04_TauWeight;   //!
   TBranch        *b_GenJet04_Charge;   //!
   TBranch        *b_GenJet04_EhadOverEem;   //!
   TBranch        *b_GenJet04_NCharged;   //!
   TBranch        *b_GenJet04_NNeutrals;   //!
   TBranch        *b_GenJet04_NeutralEnergyFraction;   //!
   TBranch        *b_GenJet04_ChargedEnergyFraction;   //!
   TBranch        *b_GenJet04_Beta;   //!
   TBranch        *b_GenJet04_BetaStar;   //!
   TBranch        *b_GenJet04_MeanSqDeltaR;   //!
   TBranch        *b_GenJet04_PTD;   //!
   TBranch        *b_GenJet04_FracPt;   //!
   TBranch        *b_GenJet04_Tau;   //!
   TBranch        *b_GenJet04_SoftDroppedJet;   //!
   TBranch        *b_GenJet04_SoftDroppedSubJet1;   //!
   TBranch        *b_GenJet04_SoftDroppedSubJet2;   //!
   TBranch        *b_GenJet04_TrimmedP4;   //!
   TBranch        *b_GenJet04_PrunedP4;   //!
   TBranch        *b_GenJet04_SoftDroppedP4;   //!
   TBranch        *b_GenJet04_NSubJetsTrimmed;   //!
   TBranch        *b_GenJet04_NSubJetsPruned;   //!
   TBranch        *b_GenJet04_NSubJetsSoftDropped;   //!
   TBranch        *b_GenJet04_ExclYmerge23;   //!
   TBranch        *b_GenJet04_ExclYmerge34;   //!
   TBranch        *b_GenJet04_ExclYmerge45;   //!
   TBranch        *b_GenJet04_ExclYmerge56;   //!
   TBranch        *b_GenJet04_Constituents;   //!
   TBranch        *b_GenJet04_Particles;   //!
   TBranch        *b_GenJet04_Area;   //!
   TBranch        *b_GenJet04_size;   //!
   TBranch        *b_GenJet08_;   //!
   TBranch        *b_GenJet08_fUniqueID;   //!
   TBranch        *b_GenJet08_fBits;   //!
   TBranch        *b_GenJet08_PT;   //!
   TBranch        *b_GenJet08_Eta;   //!
   TBranch        *b_GenJet08_Phi;   //!
   TBranch        *b_GenJet08_T;   //!
   TBranch        *b_GenJet08_Mass;   //!
   TBranch        *b_GenJet08_DeltaEta;   //!
   TBranch        *b_GenJet08_DeltaPhi;   //!
   TBranch        *b_GenJet08_Flavor;   //!
   TBranch        *b_GenJet08_FlavorAlgo;   //!
   TBranch        *b_GenJet08_FlavorPhys;   //!
   TBranch        *b_GenJet08_BTag;   //!
   TBranch        *b_GenJet08_BTagAlgo;   //!
   TBranch        *b_GenJet08_BTagPhys;   //!
   TBranch        *b_GenJet08_TauTag;   //!
   TBranch        *b_GenJet08_TauWeight;   //!
   TBranch        *b_GenJet08_Charge;   //!
   TBranch        *b_GenJet08_EhadOverEem;   //!
   TBranch        *b_GenJet08_NCharged;   //!
   TBranch        *b_GenJet08_NNeutrals;   //!
   TBranch        *b_GenJet08_NeutralEnergyFraction;   //!
   TBranch        *b_GenJet08_ChargedEnergyFraction;   //!
   TBranch        *b_GenJet08_Beta;   //!
   TBranch        *b_GenJet08_BetaStar;   //!
   TBranch        *b_GenJet08_MeanSqDeltaR;   //!
   TBranch        *b_GenJet08_PTD;   //!
   TBranch        *b_GenJet08_FracPt;   //!
   TBranch        *b_GenJet08_Tau;   //!
   TBranch        *b_GenJet08_SoftDroppedJet;   //!
   TBranch        *b_GenJet08_SoftDroppedSubJet1;   //!
   TBranch        *b_GenJet08_SoftDroppedSubJet2;   //!
   TBranch        *b_GenJet08_TrimmedP4;   //!
   TBranch        *b_GenJet08_PrunedP4;   //!
   TBranch        *b_GenJet08_SoftDroppedP4;   //!
   TBranch        *b_GenJet08_NSubJetsTrimmed;   //!
   TBranch        *b_GenJet08_NSubJetsPruned;   //!
   TBranch        *b_GenJet08_NSubJetsSoftDropped;   //!
   TBranch        *b_GenJet08_ExclYmerge23;   //!
   TBranch        *b_GenJet08_ExclYmerge34;   //!
   TBranch        *b_GenJet08_ExclYmerge45;   //!
   TBranch        *b_GenJet08_ExclYmerge56;   //!
   TBranch        *b_GenJet08_Constituents;   //!
   TBranch        *b_GenJet08_Particles;   //!
   TBranch        *b_GenJet08_Area;   //!
   TBranch        *b_GenJet08_size;   //!
   TBranch        *b_GenJet15_;   //!
   TBranch        *b_GenJet15_fUniqueID;   //!
   TBranch        *b_GenJet15_fBits;   //!
   TBranch        *b_GenJet15_PT;   //!
   TBranch        *b_GenJet15_Eta;   //!
   TBranch        *b_GenJet15_Phi;   //!
   TBranch        *b_GenJet15_T;   //!
   TBranch        *b_GenJet15_Mass;   //!
   TBranch        *b_GenJet15_DeltaEta;   //!
   TBranch        *b_GenJet15_DeltaPhi;   //!
   TBranch        *b_GenJet15_Flavor;   //!
   TBranch        *b_GenJet15_FlavorAlgo;   //!
   TBranch        *b_GenJet15_FlavorPhys;   //!
   TBranch        *b_GenJet15_BTag;   //!
   TBranch        *b_GenJet15_BTagAlgo;   //!
   TBranch        *b_GenJet15_BTagPhys;   //!
   TBranch        *b_GenJet15_TauTag;   //!
   TBranch        *b_GenJet15_TauWeight;   //!
   TBranch        *b_GenJet15_Charge;   //!
   TBranch        *b_GenJet15_EhadOverEem;   //!
   TBranch        *b_GenJet15_NCharged;   //!
   TBranch        *b_GenJet15_NNeutrals;   //!
   TBranch        *b_GenJet15_NeutralEnergyFraction;   //!
   TBranch        *b_GenJet15_ChargedEnergyFraction;   //!
   TBranch        *b_GenJet15_Beta;   //!
   TBranch        *b_GenJet15_BetaStar;   //!
   TBranch        *b_GenJet15_MeanSqDeltaR;   //!
   TBranch        *b_GenJet15_PTD;   //!
   TBranch        *b_GenJet15_FracPt;   //!
   TBranch        *b_GenJet15_Tau;   //!
   TBranch        *b_GenJet15_SoftDroppedJet;   //!
   TBranch        *b_GenJet15_SoftDroppedSubJet1;   //!
   TBranch        *b_GenJet15_SoftDroppedSubJet2;   //!
   TBranch        *b_GenJet15_TrimmedP4;   //!
   TBranch        *b_GenJet15_PrunedP4;   //!
   TBranch        *b_GenJet15_SoftDroppedP4;   //!
   TBranch        *b_GenJet15_NSubJetsTrimmed;   //!
   TBranch        *b_GenJet15_NSubJetsPruned;   //!
   TBranch        *b_GenJet15_NSubJetsSoftDropped;   //!
   TBranch        *b_GenJet15_ExclYmerge23;   //!
   TBranch        *b_GenJet15_ExclYmerge34;   //!
   TBranch        *b_GenJet15_ExclYmerge45;   //!
   TBranch        *b_GenJet15_ExclYmerge56;   //!
   TBranch        *b_GenJet15_Constituents;   //!
   TBranch        *b_GenJet15_Particles;   //!
   TBranch        *b_GenJet15_Area;   //!
   TBranch        *b_GenJet15_size;   //!
   TBranch        *b_ParticleFlowJet02_;   //!
   TBranch        *b_ParticleFlowJet02_fUniqueID;   //!
   TBranch        *b_ParticleFlowJet02_fBits;   //!
   TBranch        *b_ParticleFlowJet02_PT;   //!
   TBranch        *b_ParticleFlowJet02_Eta;   //!
   TBranch        *b_ParticleFlowJet02_Phi;   //!
   TBranch        *b_ParticleFlowJet02_T;   //!
   TBranch        *b_ParticleFlowJet02_Mass;   //!
   TBranch        *b_ParticleFlowJet02_DeltaEta;   //!
   TBranch        *b_ParticleFlowJet02_DeltaPhi;   //!
   TBranch        *b_ParticleFlowJet02_Flavor;   //!
   TBranch        *b_ParticleFlowJet02_FlavorAlgo;   //!
   TBranch        *b_ParticleFlowJet02_FlavorPhys;   //!
   TBranch        *b_ParticleFlowJet02_BTag;   //!
   TBranch        *b_ParticleFlowJet02_BTagAlgo;   //!
   TBranch        *b_ParticleFlowJet02_BTagPhys;   //!
   TBranch        *b_ParticleFlowJet02_TauTag;   //!
   TBranch        *b_ParticleFlowJet02_TauWeight;   //!
   TBranch        *b_ParticleFlowJet02_Charge;   //!
   TBranch        *b_ParticleFlowJet02_EhadOverEem;   //!
   TBranch        *b_ParticleFlowJet02_NCharged;   //!
   TBranch        *b_ParticleFlowJet02_NNeutrals;   //!
   TBranch        *b_ParticleFlowJet02_NeutralEnergyFraction;   //!
   TBranch        *b_ParticleFlowJet02_ChargedEnergyFraction;   //!
   TBranch        *b_ParticleFlowJet02_Beta;   //!
   TBranch        *b_ParticleFlowJet02_BetaStar;   //!
   TBranch        *b_ParticleFlowJet02_MeanSqDeltaR;   //!
   TBranch        *b_ParticleFlowJet02_PTD;   //!
   TBranch        *b_ParticleFlowJet02_FracPt;   //!
   TBranch        *b_ParticleFlowJet02_Tau;   //!
   TBranch        *b_ParticleFlowJet02_SoftDroppedJet;   //!
   TBranch        *b_ParticleFlowJet02_SoftDroppedSubJet1;   //!
   TBranch        *b_ParticleFlowJet02_SoftDroppedSubJet2;   //!
   TBranch        *b_ParticleFlowJet02_TrimmedP4;   //!
   TBranch        *b_ParticleFlowJet02_PrunedP4;   //!
   TBranch        *b_ParticleFlowJet02_SoftDroppedP4;   //!
   TBranch        *b_ParticleFlowJet02_NSubJetsTrimmed;   //!
   TBranch        *b_ParticleFlowJet02_NSubJetsPruned;   //!
   TBranch        *b_ParticleFlowJet02_NSubJetsSoftDropped;   //!
   TBranch        *b_ParticleFlowJet02_ExclYmerge23;   //!
   TBranch        *b_ParticleFlowJet02_ExclYmerge34;   //!
   TBranch        *b_ParticleFlowJet02_ExclYmerge45;   //!
   TBranch        *b_ParticleFlowJet02_ExclYmerge56;   //!
   TBranch        *b_ParticleFlowJet02_Constituents;   //!
   TBranch        *b_ParticleFlowJet02_Particles;   //!
   TBranch        *b_ParticleFlowJet02_Area;   //!
   TBranch        *b_ParticleFlowJet02_size;   //!
   TBranch        *b_ParticleFlowJet04_;   //!
   TBranch        *b_ParticleFlowJet04_fUniqueID;   //!
   TBranch        *b_ParticleFlowJet04_fBits;   //!
   TBranch        *b_ParticleFlowJet04_PT;   //!
   TBranch        *b_ParticleFlowJet04_Eta;   //!
   TBranch        *b_ParticleFlowJet04_Phi;   //!
   TBranch        *b_ParticleFlowJet04_T;   //!
   TBranch        *b_ParticleFlowJet04_Mass;   //!
   TBranch        *b_ParticleFlowJet04_DeltaEta;   //!
   TBranch        *b_ParticleFlowJet04_DeltaPhi;   //!
   TBranch        *b_ParticleFlowJet04_Flavor;   //!
   TBranch        *b_ParticleFlowJet04_FlavorAlgo;   //!
   TBranch        *b_ParticleFlowJet04_FlavorPhys;   //!
   TBranch        *b_ParticleFlowJet04_BTag;   //!
   TBranch        *b_ParticleFlowJet04_BTagAlgo;   //!
   TBranch        *b_ParticleFlowJet04_BTagPhys;   //!
   TBranch        *b_ParticleFlowJet04_TauTag;   //!
   TBranch        *b_ParticleFlowJet04_TauWeight;   //!
   TBranch        *b_ParticleFlowJet04_Charge;   //!
   TBranch        *b_ParticleFlowJet04_EhadOverEem;   //!
   TBranch        *b_ParticleFlowJet04_NCharged;   //!
   TBranch        *b_ParticleFlowJet04_NNeutrals;   //!
   TBranch        *b_ParticleFlowJet04_NeutralEnergyFraction;   //!
   TBranch        *b_ParticleFlowJet04_ChargedEnergyFraction;   //!
   TBranch        *b_ParticleFlowJet04_Beta;   //!
   TBranch        *b_ParticleFlowJet04_BetaStar;   //!
   TBranch        *b_ParticleFlowJet04_MeanSqDeltaR;   //!
   TBranch        *b_ParticleFlowJet04_PTD;   //!
   TBranch        *b_ParticleFlowJet04_FracPt;   //!
   TBranch        *b_ParticleFlowJet04_Tau;   //!
   TBranch        *b_ParticleFlowJet04_SoftDroppedJet;   //!
   TBranch        *b_ParticleFlowJet04_SoftDroppedSubJet1;   //!
   TBranch        *b_ParticleFlowJet04_SoftDroppedSubJet2;   //!
   TBranch        *b_ParticleFlowJet04_TrimmedP4;   //!
   TBranch        *b_ParticleFlowJet04_PrunedP4;   //!
   TBranch        *b_ParticleFlowJet04_SoftDroppedP4;   //!
   TBranch        *b_ParticleFlowJet04_NSubJetsTrimmed;   //!
   TBranch        *b_ParticleFlowJet04_NSubJetsPruned;   //!
   TBranch        *b_ParticleFlowJet04_NSubJetsSoftDropped;   //!
   TBranch        *b_ParticleFlowJet04_ExclYmerge23;   //!
   TBranch        *b_ParticleFlowJet04_ExclYmerge34;   //!
   TBranch        *b_ParticleFlowJet04_ExclYmerge45;   //!
   TBranch        *b_ParticleFlowJet04_ExclYmerge56;   //!
   TBranch        *b_ParticleFlowJet04_Constituents;   //!
   TBranch        *b_ParticleFlowJet04_Particles;   //!
   TBranch        *b_ParticleFlowJet04_Area;   //!
   TBranch        *b_ParticleFlowJet04_size;   //!
   TBranch        *b_ParticleFlowJet08_;   //!
   TBranch        *b_ParticleFlowJet08_fUniqueID;   //!
   TBranch        *b_ParticleFlowJet08_fBits;   //!
   TBranch        *b_ParticleFlowJet08_PT;   //!
   TBranch        *b_ParticleFlowJet08_Eta;   //!
   TBranch        *b_ParticleFlowJet08_Phi;   //!
   TBranch        *b_ParticleFlowJet08_T;   //!
   TBranch        *b_ParticleFlowJet08_Mass;   //!
   TBranch        *b_ParticleFlowJet08_DeltaEta;   //!
   TBranch        *b_ParticleFlowJet08_DeltaPhi;   //!
   TBranch        *b_ParticleFlowJet08_Flavor;   //!
   TBranch        *b_ParticleFlowJet08_FlavorAlgo;   //!
   TBranch        *b_ParticleFlowJet08_FlavorPhys;   //!
   TBranch        *b_ParticleFlowJet08_BTag;   //!
   TBranch        *b_ParticleFlowJet08_BTagAlgo;   //!
   TBranch        *b_ParticleFlowJet08_BTagPhys;   //!
   TBranch        *b_ParticleFlowJet08_TauTag;   //!
   TBranch        *b_ParticleFlowJet08_TauWeight;   //!
   TBranch        *b_ParticleFlowJet08_Charge;   //!
   TBranch        *b_ParticleFlowJet08_EhadOverEem;   //!
   TBranch        *b_ParticleFlowJet08_NCharged;   //!
   TBranch        *b_ParticleFlowJet08_NNeutrals;   //!
   TBranch        *b_ParticleFlowJet08_NeutralEnergyFraction;   //!
   TBranch        *b_ParticleFlowJet08_ChargedEnergyFraction;   //!
   TBranch        *b_ParticleFlowJet08_Beta;   //!
   TBranch        *b_ParticleFlowJet08_BetaStar;   //!
   TBranch        *b_ParticleFlowJet08_MeanSqDeltaR;   //!
   TBranch        *b_ParticleFlowJet08_PTD;   //!
   TBranch        *b_ParticleFlowJet08_FracPt;   //!
   TBranch        *b_ParticleFlowJet08_Tau;   //!
   TBranch        *b_ParticleFlowJet08_SoftDroppedJet;   //!
   TBranch        *b_ParticleFlowJet08_SoftDroppedSubJet1;   //!
   TBranch        *b_ParticleFlowJet08_SoftDroppedSubJet2;   //!
   TBranch        *b_ParticleFlowJet08_TrimmedP4;   //!
   TBranch        *b_ParticleFlowJet08_PrunedP4;   //!
   TBranch        *b_ParticleFlowJet08_SoftDroppedP4;   //!
   TBranch        *b_ParticleFlowJet08_NSubJetsTrimmed;   //!
   TBranch        *b_ParticleFlowJet08_NSubJetsPruned;   //!
   TBranch        *b_ParticleFlowJet08_NSubJetsSoftDropped;   //!
   TBranch        *b_ParticleFlowJet08_ExclYmerge23;   //!
   TBranch        *b_ParticleFlowJet08_ExclYmerge34;   //!
   TBranch        *b_ParticleFlowJet08_ExclYmerge45;   //!
   TBranch        *b_ParticleFlowJet08_ExclYmerge56;   //!
   TBranch        *b_ParticleFlowJet08_Constituents;   //!
   TBranch        *b_ParticleFlowJet08_Particles;   //!
   TBranch        *b_ParticleFlowJet08_Area;   //!
   TBranch        *b_ParticleFlowJet08_size;   //!
   TBranch        *b_ParticleFlowJet15_;   //!
   TBranch        *b_ParticleFlowJet15_fUniqueID;   //!
   TBranch        *b_ParticleFlowJet15_fBits;   //!
   TBranch        *b_ParticleFlowJet15_PT;   //!
   TBranch        *b_ParticleFlowJet15_Eta;   //!
   TBranch        *b_ParticleFlowJet15_Phi;   //!
   TBranch        *b_ParticleFlowJet15_T;   //!
   TBranch        *b_ParticleFlowJet15_Mass;   //!
   TBranch        *b_ParticleFlowJet15_DeltaEta;   //!
   TBranch        *b_ParticleFlowJet15_DeltaPhi;   //!
   TBranch        *b_ParticleFlowJet15_Flavor;   //!
   TBranch        *b_ParticleFlowJet15_FlavorAlgo;   //!
   TBranch        *b_ParticleFlowJet15_FlavorPhys;   //!
   TBranch        *b_ParticleFlowJet15_BTag;   //!
   TBranch        *b_ParticleFlowJet15_BTagAlgo;   //!
   TBranch        *b_ParticleFlowJet15_BTagPhys;   //!
   TBranch        *b_ParticleFlowJet15_TauTag;   //!
   TBranch        *b_ParticleFlowJet15_TauWeight;   //!
   TBranch        *b_ParticleFlowJet15_Charge;   //!
   TBranch        *b_ParticleFlowJet15_EhadOverEem;   //!
   TBranch        *b_ParticleFlowJet15_NCharged;   //!
   TBranch        *b_ParticleFlowJet15_NNeutrals;   //!
   TBranch        *b_ParticleFlowJet15_NeutralEnergyFraction;   //!
   TBranch        *b_ParticleFlowJet15_ChargedEnergyFraction;   //!
   TBranch        *b_ParticleFlowJet15_Beta;   //!
   TBranch        *b_ParticleFlowJet15_BetaStar;   //!
   TBranch        *b_ParticleFlowJet15_MeanSqDeltaR;   //!
   TBranch        *b_ParticleFlowJet15_PTD;   //!
   TBranch        *b_ParticleFlowJet15_FracPt;   //!
   TBranch        *b_ParticleFlowJet15_Tau;   //!
   TBranch        *b_ParticleFlowJet15_SoftDroppedJet;   //!
   TBranch        *b_ParticleFlowJet15_SoftDroppedSubJet1;   //!
   TBranch        *b_ParticleFlowJet15_SoftDroppedSubJet2;   //!
   TBranch        *b_ParticleFlowJet15_TrimmedP4;   //!
   TBranch        *b_ParticleFlowJet15_PrunedP4;   //!
   TBranch        *b_ParticleFlowJet15_SoftDroppedP4;   //!
   TBranch        *b_ParticleFlowJet15_NSubJetsTrimmed;   //!
   TBranch        *b_ParticleFlowJet15_NSubJetsPruned;   //!
   TBranch        *b_ParticleFlowJet15_NSubJetsSoftDropped;   //!
   TBranch        *b_ParticleFlowJet15_ExclYmerge23;   //!
   TBranch        *b_ParticleFlowJet15_ExclYmerge34;   //!
   TBranch        *b_ParticleFlowJet15_ExclYmerge45;   //!
   TBranch        *b_ParticleFlowJet15_ExclYmerge56;   //!
   TBranch        *b_ParticleFlowJet15_Constituents;   //!
   TBranch        *b_ParticleFlowJet15_Particles;   //!
   TBranch        *b_ParticleFlowJet15_Area;   //!
   TBranch        *b_ParticleFlowJet15_size;   //!
   TBranch        *b_CaloJet02_;   //!
   TBranch        *b_CaloJet02_fUniqueID;   //!
   TBranch        *b_CaloJet02_fBits;   //!
   TBranch        *b_CaloJet02_PT;   //!
   TBranch        *b_CaloJet02_Eta;   //!
   TBranch        *b_CaloJet02_Phi;   //!
   TBranch        *b_CaloJet02_T;   //!
   TBranch        *b_CaloJet02_Mass;   //!
   TBranch        *b_CaloJet02_DeltaEta;   //!
   TBranch        *b_CaloJet02_DeltaPhi;   //!
   TBranch        *b_CaloJet02_Flavor;   //!
   TBranch        *b_CaloJet02_FlavorAlgo;   //!
   TBranch        *b_CaloJet02_FlavorPhys;   //!
   TBranch        *b_CaloJet02_BTag;   //!
   TBranch        *b_CaloJet02_BTagAlgo;   //!
   TBranch        *b_CaloJet02_BTagPhys;   //!
   TBranch        *b_CaloJet02_TauTag;   //!
   TBranch        *b_CaloJet02_TauWeight;   //!
   TBranch        *b_CaloJet02_Charge;   //!
   TBranch        *b_CaloJet02_EhadOverEem;   //!
   TBranch        *b_CaloJet02_NCharged;   //!
   TBranch        *b_CaloJet02_NNeutrals;   //!
   TBranch        *b_CaloJet02_NeutralEnergyFraction;   //!
   TBranch        *b_CaloJet02_ChargedEnergyFraction;   //!
   TBranch        *b_CaloJet02_Beta;   //!
   TBranch        *b_CaloJet02_BetaStar;   //!
   TBranch        *b_CaloJet02_MeanSqDeltaR;   //!
   TBranch        *b_CaloJet02_PTD;   //!
   TBranch        *b_CaloJet02_FracPt;   //!
   TBranch        *b_CaloJet02_Tau;   //!
   TBranch        *b_CaloJet02_SoftDroppedJet;   //!
   TBranch        *b_CaloJet02_SoftDroppedSubJet1;   //!
   TBranch        *b_CaloJet02_SoftDroppedSubJet2;   //!
   TBranch        *b_CaloJet02_TrimmedP4;   //!
   TBranch        *b_CaloJet02_PrunedP4;   //!
   TBranch        *b_CaloJet02_SoftDroppedP4;   //!
   TBranch        *b_CaloJet02_NSubJetsTrimmed;   //!
   TBranch        *b_CaloJet02_NSubJetsPruned;   //!
   TBranch        *b_CaloJet02_NSubJetsSoftDropped;   //!
   TBranch        *b_CaloJet02_ExclYmerge23;   //!
   TBranch        *b_CaloJet02_ExclYmerge34;   //!
   TBranch        *b_CaloJet02_ExclYmerge45;   //!
   TBranch        *b_CaloJet02_ExclYmerge56;   //!
   TBranch        *b_CaloJet02_Constituents;   //!
   TBranch        *b_CaloJet02_Particles;   //!
   TBranch        *b_CaloJet02_Area;   //!
   TBranch        *b_CaloJet02_size;   //!
   TBranch        *b_CaloJet04_;   //!
   TBranch        *b_CaloJet04_fUniqueID;   //!
   TBranch        *b_CaloJet04_fBits;   //!
   TBranch        *b_CaloJet04_PT;   //!
   TBranch        *b_CaloJet04_Eta;   //!
   TBranch        *b_CaloJet04_Phi;   //!
   TBranch        *b_CaloJet04_T;   //!
   TBranch        *b_CaloJet04_Mass;   //!
   TBranch        *b_CaloJet04_DeltaEta;   //!
   TBranch        *b_CaloJet04_DeltaPhi;   //!
   TBranch        *b_CaloJet04_Flavor;   //!
   TBranch        *b_CaloJet04_FlavorAlgo;   //!
   TBranch        *b_CaloJet04_FlavorPhys;   //!
   TBranch        *b_CaloJet04_BTag;   //!
   TBranch        *b_CaloJet04_BTagAlgo;   //!
   TBranch        *b_CaloJet04_BTagPhys;   //!
   TBranch        *b_CaloJet04_TauTag;   //!
   TBranch        *b_CaloJet04_TauWeight;   //!
   TBranch        *b_CaloJet04_Charge;   //!
   TBranch        *b_CaloJet04_EhadOverEem;   //!
   TBranch        *b_CaloJet04_NCharged;   //!
   TBranch        *b_CaloJet04_NNeutrals;   //!
   TBranch        *b_CaloJet04_NeutralEnergyFraction;   //!
   TBranch        *b_CaloJet04_ChargedEnergyFraction;   //!
   TBranch        *b_CaloJet04_Beta;   //!
   TBranch        *b_CaloJet04_BetaStar;   //!
   TBranch        *b_CaloJet04_MeanSqDeltaR;   //!
   TBranch        *b_CaloJet04_PTD;   //!
   TBranch        *b_CaloJet04_FracPt;   //!
   TBranch        *b_CaloJet04_Tau;   //!
   TBranch        *b_CaloJet04_SoftDroppedJet;   //!
   TBranch        *b_CaloJet04_SoftDroppedSubJet1;   //!
   TBranch        *b_CaloJet04_SoftDroppedSubJet2;   //!
   TBranch        *b_CaloJet04_TrimmedP4;   //!
   TBranch        *b_CaloJet04_PrunedP4;   //!
   TBranch        *b_CaloJet04_SoftDroppedP4;   //!
   TBranch        *b_CaloJet04_NSubJetsTrimmed;   //!
   TBranch        *b_CaloJet04_NSubJetsPruned;   //!
   TBranch        *b_CaloJet04_NSubJetsSoftDropped;   //!
   TBranch        *b_CaloJet04_ExclYmerge23;   //!
   TBranch        *b_CaloJet04_ExclYmerge34;   //!
   TBranch        *b_CaloJet04_ExclYmerge45;   //!
   TBranch        *b_CaloJet04_ExclYmerge56;   //!
   TBranch        *b_CaloJet04_Constituents;   //!
   TBranch        *b_CaloJet04_Particles;   //!
   TBranch        *b_CaloJet04_Area;   //!
   TBranch        *b_CaloJet04_size;   //!
   TBranch        *b_CaloJet08_;   //!
   TBranch        *b_CaloJet08_fUniqueID;   //!
   TBranch        *b_CaloJet08_fBits;   //!
   TBranch        *b_CaloJet08_PT;   //!
   TBranch        *b_CaloJet08_Eta;   //!
   TBranch        *b_CaloJet08_Phi;   //!
   TBranch        *b_CaloJet08_T;   //!
   TBranch        *b_CaloJet08_Mass;   //!
   TBranch        *b_CaloJet08_DeltaEta;   //!
   TBranch        *b_CaloJet08_DeltaPhi;   //!
   TBranch        *b_CaloJet08_Flavor;   //!
   TBranch        *b_CaloJet08_FlavorAlgo;   //!
   TBranch        *b_CaloJet08_FlavorPhys;   //!
   TBranch        *b_CaloJet08_BTag;   //!
   TBranch        *b_CaloJet08_BTagAlgo;   //!
   TBranch        *b_CaloJet08_BTagPhys;   //!
   TBranch        *b_CaloJet08_TauTag;   //!
   TBranch        *b_CaloJet08_TauWeight;   //!
   TBranch        *b_CaloJet08_Charge;   //!
   TBranch        *b_CaloJet08_EhadOverEem;   //!
   TBranch        *b_CaloJet08_NCharged;   //!
   TBranch        *b_CaloJet08_NNeutrals;   //!
   TBranch        *b_CaloJet08_NeutralEnergyFraction;   //!
   TBranch        *b_CaloJet08_ChargedEnergyFraction;   //!
   TBranch        *b_CaloJet08_Beta;   //!
   TBranch        *b_CaloJet08_BetaStar;   //!
   TBranch        *b_CaloJet08_MeanSqDeltaR;   //!
   TBranch        *b_CaloJet08_PTD;   //!
   TBranch        *b_CaloJet08_FracPt;   //!
   TBranch        *b_CaloJet08_Tau;   //!
   TBranch        *b_CaloJet08_SoftDroppedJet;   //!
   TBranch        *b_CaloJet08_SoftDroppedSubJet1;   //!
   TBranch        *b_CaloJet08_SoftDroppedSubJet2;   //!
   TBranch        *b_CaloJet08_TrimmedP4;   //!
   TBranch        *b_CaloJet08_PrunedP4;   //!
   TBranch        *b_CaloJet08_SoftDroppedP4;   //!
   TBranch        *b_CaloJet08_NSubJetsTrimmed;   //!
   TBranch        *b_CaloJet08_NSubJetsPruned;   //!
   TBranch        *b_CaloJet08_NSubJetsSoftDropped;   //!
   TBranch        *b_CaloJet08_ExclYmerge23;   //!
   TBranch        *b_CaloJet08_ExclYmerge34;   //!
   TBranch        *b_CaloJet08_ExclYmerge45;   //!
   TBranch        *b_CaloJet08_ExclYmerge56;   //!
   TBranch        *b_CaloJet08_Constituents;   //!
   TBranch        *b_CaloJet08_Particles;   //!
   TBranch        *b_CaloJet08_Area;   //!
   TBranch        *b_CaloJet08_size;   //!
   TBranch        *b_CaloJet15_;   //!
   TBranch        *b_CaloJet15_fUniqueID;   //!
   TBranch        *b_CaloJet15_fBits;   //!
   TBranch        *b_CaloJet15_PT;   //!
   TBranch        *b_CaloJet15_Eta;   //!
   TBranch        *b_CaloJet15_Phi;   //!
   TBranch        *b_CaloJet15_T;   //!
   TBranch        *b_CaloJet15_Mass;   //!
   TBranch        *b_CaloJet15_DeltaEta;   //!
   TBranch        *b_CaloJet15_DeltaPhi;   //!
   TBranch        *b_CaloJet15_Flavor;   //!
   TBranch        *b_CaloJet15_FlavorAlgo;   //!
   TBranch        *b_CaloJet15_FlavorPhys;   //!
   TBranch        *b_CaloJet15_BTag;   //!
   TBranch        *b_CaloJet15_BTagAlgo;   //!
   TBranch        *b_CaloJet15_BTagPhys;   //!
   TBranch        *b_CaloJet15_TauTag;   //!
   TBranch        *b_CaloJet15_TauWeight;   //!
   TBranch        *b_CaloJet15_Charge;   //!
   TBranch        *b_CaloJet15_EhadOverEem;   //!
   TBranch        *b_CaloJet15_NCharged;   //!
   TBranch        *b_CaloJet15_NNeutrals;   //!
   TBranch        *b_CaloJet15_NeutralEnergyFraction;   //!
   TBranch        *b_CaloJet15_ChargedEnergyFraction;   //!
   TBranch        *b_CaloJet15_Beta;   //!
   TBranch        *b_CaloJet15_BetaStar;   //!
   TBranch        *b_CaloJet15_MeanSqDeltaR;   //!
   TBranch        *b_CaloJet15_PTD;   //!
   TBranch        *b_CaloJet15_FracPt;   //!
   TBranch        *b_CaloJet15_Tau;   //!
   TBranch        *b_CaloJet15_SoftDroppedJet;   //!
   TBranch        *b_CaloJet15_SoftDroppedSubJet1;   //!
   TBranch        *b_CaloJet15_SoftDroppedSubJet2;   //!
   TBranch        *b_CaloJet15_TrimmedP4;   //!
   TBranch        *b_CaloJet15_PrunedP4;   //!
   TBranch        *b_CaloJet15_SoftDroppedP4;   //!
   TBranch        *b_CaloJet15_NSubJetsTrimmed;   //!
   TBranch        *b_CaloJet15_NSubJetsPruned;   //!
   TBranch        *b_CaloJet15_NSubJetsSoftDropped;   //!
   TBranch        *b_CaloJet15_ExclYmerge23;   //!
   TBranch        *b_CaloJet15_ExclYmerge34;   //!
   TBranch        *b_CaloJet15_ExclYmerge45;   //!
   TBranch        *b_CaloJet15_ExclYmerge56;   //!
   TBranch        *b_CaloJet15_Constituents;   //!
   TBranch        *b_CaloJet15_Particles;   //!
   TBranch        *b_CaloJet15_Area;   //!
   TBranch        *b_CaloJet15_size;   //!
   TBranch        *b_TrackJet02_;   //!
   TBranch        *b_TrackJet02_fUniqueID;   //!
   TBranch        *b_TrackJet02_fBits;   //!
   TBranch        *b_TrackJet02_PT;   //!
   TBranch        *b_TrackJet02_Eta;   //!
   TBranch        *b_TrackJet02_Phi;   //!
   TBranch        *b_TrackJet02_T;   //!
   TBranch        *b_TrackJet02_Mass;   //!
   TBranch        *b_TrackJet02_DeltaEta;   //!
   TBranch        *b_TrackJet02_DeltaPhi;   //!
   TBranch        *b_TrackJet02_Flavor;   //!
   TBranch        *b_TrackJet02_FlavorAlgo;   //!
   TBranch        *b_TrackJet02_FlavorPhys;   //!
   TBranch        *b_TrackJet02_BTag;   //!
   TBranch        *b_TrackJet02_BTagAlgo;   //!
   TBranch        *b_TrackJet02_BTagPhys;   //!
   TBranch        *b_TrackJet02_TauTag;   //!
   TBranch        *b_TrackJet02_TauWeight;   //!
   TBranch        *b_TrackJet02_Charge;   //!
   TBranch        *b_TrackJet02_EhadOverEem;   //!
   TBranch        *b_TrackJet02_NCharged;   //!
   TBranch        *b_TrackJet02_NNeutrals;   //!
   TBranch        *b_TrackJet02_NeutralEnergyFraction;   //!
   TBranch        *b_TrackJet02_ChargedEnergyFraction;   //!
   TBranch        *b_TrackJet02_Beta;   //!
   TBranch        *b_TrackJet02_BetaStar;   //!
   TBranch        *b_TrackJet02_MeanSqDeltaR;   //!
   TBranch        *b_TrackJet02_PTD;   //!
   TBranch        *b_TrackJet02_FracPt;   //!
   TBranch        *b_TrackJet02_Tau;   //!
   TBranch        *b_TrackJet02_SoftDroppedJet;   //!
   TBranch        *b_TrackJet02_SoftDroppedSubJet1;   //!
   TBranch        *b_TrackJet02_SoftDroppedSubJet2;   //!
   TBranch        *b_TrackJet02_TrimmedP4;   //!
   TBranch        *b_TrackJet02_PrunedP4;   //!
   TBranch        *b_TrackJet02_SoftDroppedP4;   //!
   TBranch        *b_TrackJet02_NSubJetsTrimmed;   //!
   TBranch        *b_TrackJet02_NSubJetsPruned;   //!
   TBranch        *b_TrackJet02_NSubJetsSoftDropped;   //!
   TBranch        *b_TrackJet02_ExclYmerge23;   //!
   TBranch        *b_TrackJet02_ExclYmerge34;   //!
   TBranch        *b_TrackJet02_ExclYmerge45;   //!
   TBranch        *b_TrackJet02_ExclYmerge56;   //!
   TBranch        *b_TrackJet02_Constituents;   //!
   TBranch        *b_TrackJet02_Particles;   //!
   TBranch        *b_TrackJet02_Area;   //!
   TBranch        *b_TrackJet02_size;   //!
   TBranch        *b_TrackJet04_;   //!
   TBranch        *b_TrackJet04_fUniqueID;   //!
   TBranch        *b_TrackJet04_fBits;   //!
   TBranch        *b_TrackJet04_PT;   //!
   TBranch        *b_TrackJet04_Eta;   //!
   TBranch        *b_TrackJet04_Phi;   //!
   TBranch        *b_TrackJet04_T;   //!
   TBranch        *b_TrackJet04_Mass;   //!
   TBranch        *b_TrackJet04_DeltaEta;   //!
   TBranch        *b_TrackJet04_DeltaPhi;   //!
   TBranch        *b_TrackJet04_Flavor;   //!
   TBranch        *b_TrackJet04_FlavorAlgo;   //!
   TBranch        *b_TrackJet04_FlavorPhys;   //!
   TBranch        *b_TrackJet04_BTag;   //!
   TBranch        *b_TrackJet04_BTagAlgo;   //!
   TBranch        *b_TrackJet04_BTagPhys;   //!
   TBranch        *b_TrackJet04_TauTag;   //!
   TBranch        *b_TrackJet04_TauWeight;   //!
   TBranch        *b_TrackJet04_Charge;   //!
   TBranch        *b_TrackJet04_EhadOverEem;   //!
   TBranch        *b_TrackJet04_NCharged;   //!
   TBranch        *b_TrackJet04_NNeutrals;   //!
   TBranch        *b_TrackJet04_NeutralEnergyFraction;   //!
   TBranch        *b_TrackJet04_ChargedEnergyFraction;   //!
   TBranch        *b_TrackJet04_Beta;   //!
   TBranch        *b_TrackJet04_BetaStar;   //!
   TBranch        *b_TrackJet04_MeanSqDeltaR;   //!
   TBranch        *b_TrackJet04_PTD;   //!
   TBranch        *b_TrackJet04_FracPt;   //!
   TBranch        *b_TrackJet04_Tau;   //!
   TBranch        *b_TrackJet04_SoftDroppedJet;   //!
   TBranch        *b_TrackJet04_SoftDroppedSubJet1;   //!
   TBranch        *b_TrackJet04_SoftDroppedSubJet2;   //!
   TBranch        *b_TrackJet04_TrimmedP4;   //!
   TBranch        *b_TrackJet04_PrunedP4;   //!
   TBranch        *b_TrackJet04_SoftDroppedP4;   //!
   TBranch        *b_TrackJet04_NSubJetsTrimmed;   //!
   TBranch        *b_TrackJet04_NSubJetsPruned;   //!
   TBranch        *b_TrackJet04_NSubJetsSoftDropped;   //!
   TBranch        *b_TrackJet04_ExclYmerge23;   //!
   TBranch        *b_TrackJet04_ExclYmerge34;   //!
   TBranch        *b_TrackJet04_ExclYmerge45;   //!
   TBranch        *b_TrackJet04_ExclYmerge56;   //!
   TBranch        *b_TrackJet04_Constituents;   //!
   TBranch        *b_TrackJet04_Particles;   //!
   TBranch        *b_TrackJet04_Area;   //!
   TBranch        *b_TrackJet04_size;   //!
   TBranch        *b_TrackJet08_;   //!
   TBranch        *b_TrackJet08_fUniqueID;   //!
   TBranch        *b_TrackJet08_fBits;   //!
   TBranch        *b_TrackJet08_PT;   //!
   TBranch        *b_TrackJet08_Eta;   //!
   TBranch        *b_TrackJet08_Phi;   //!
   TBranch        *b_TrackJet08_T;   //!
   TBranch        *b_TrackJet08_Mass;   //!
   TBranch        *b_TrackJet08_DeltaEta;   //!
   TBranch        *b_TrackJet08_DeltaPhi;   //!
   TBranch        *b_TrackJet08_Flavor;   //!
   TBranch        *b_TrackJet08_FlavorAlgo;   //!
   TBranch        *b_TrackJet08_FlavorPhys;   //!
   TBranch        *b_TrackJet08_BTag;   //!
   TBranch        *b_TrackJet08_BTagAlgo;   //!
   TBranch        *b_TrackJet08_BTagPhys;   //!
   TBranch        *b_TrackJet08_TauTag;   //!
   TBranch        *b_TrackJet08_TauWeight;   //!
   TBranch        *b_TrackJet08_Charge;   //!
   TBranch        *b_TrackJet08_EhadOverEem;   //!
   TBranch        *b_TrackJet08_NCharged;   //!
   TBranch        *b_TrackJet08_NNeutrals;   //!
   TBranch        *b_TrackJet08_NeutralEnergyFraction;   //!
   TBranch        *b_TrackJet08_ChargedEnergyFraction;   //!
   TBranch        *b_TrackJet08_Beta;   //!
   TBranch        *b_TrackJet08_BetaStar;   //!
   TBranch        *b_TrackJet08_MeanSqDeltaR;   //!
   TBranch        *b_TrackJet08_PTD;   //!
   TBranch        *b_TrackJet08_FracPt;   //!
   TBranch        *b_TrackJet08_Tau;   //!
   TBranch        *b_TrackJet08_SoftDroppedJet;   //!
   TBranch        *b_TrackJet08_SoftDroppedSubJet1;   //!
   TBranch        *b_TrackJet08_SoftDroppedSubJet2;   //!
   TBranch        *b_TrackJet08_TrimmedP4;   //!
   TBranch        *b_TrackJet08_PrunedP4;   //!
   TBranch        *b_TrackJet08_SoftDroppedP4;   //!
   TBranch        *b_TrackJet08_NSubJetsTrimmed;   //!
   TBranch        *b_TrackJet08_NSubJetsPruned;   //!
   TBranch        *b_TrackJet08_NSubJetsSoftDropped;   //!
   TBranch        *b_TrackJet08_ExclYmerge23;   //!
   TBranch        *b_TrackJet08_ExclYmerge34;   //!
   TBranch        *b_TrackJet08_ExclYmerge45;   //!
   TBranch        *b_TrackJet08_ExclYmerge56;   //!
   TBranch        *b_TrackJet08_Constituents;   //!
   TBranch        *b_TrackJet08_Particles;   //!
   TBranch        *b_TrackJet08_Area;   //!
   TBranch        *b_TrackJet08_size;   //!
   TBranch        *b_TrackJet15_;   //!
   TBranch        *b_TrackJet15_fUniqueID;   //!
   TBranch        *b_TrackJet15_fBits;   //!
   TBranch        *b_TrackJet15_PT;   //!
   TBranch        *b_TrackJet15_Eta;   //!
   TBranch        *b_TrackJet15_Phi;   //!
   TBranch        *b_TrackJet15_T;   //!
   TBranch        *b_TrackJet15_Mass;   //!
   TBranch        *b_TrackJet15_DeltaEta;   //!
   TBranch        *b_TrackJet15_DeltaPhi;   //!
   TBranch        *b_TrackJet15_Flavor;   //!
   TBranch        *b_TrackJet15_FlavorAlgo;   //!
   TBranch        *b_TrackJet15_FlavorPhys;   //!
   TBranch        *b_TrackJet15_BTag;   //!
   TBranch        *b_TrackJet15_BTagAlgo;   //!
   TBranch        *b_TrackJet15_BTagPhys;   //!
   TBranch        *b_TrackJet15_TauTag;   //!
   TBranch        *b_TrackJet15_TauWeight;   //!
   TBranch        *b_TrackJet15_Charge;   //!
   TBranch        *b_TrackJet15_EhadOverEem;   //!
   TBranch        *b_TrackJet15_NCharged;   //!
   TBranch        *b_TrackJet15_NNeutrals;   //!
   TBranch        *b_TrackJet15_NeutralEnergyFraction;   //!
   TBranch        *b_TrackJet15_ChargedEnergyFraction;   //!
   TBranch        *b_TrackJet15_Beta;   //!
   TBranch        *b_TrackJet15_BetaStar;   //!
   TBranch        *b_TrackJet15_MeanSqDeltaR;   //!
   TBranch        *b_TrackJet15_PTD;   //!
   TBranch        *b_TrackJet15_FracPt;   //!
   TBranch        *b_TrackJet15_Tau;   //!
   TBranch        *b_TrackJet15_SoftDroppedJet;   //!
   TBranch        *b_TrackJet15_SoftDroppedSubJet1;   //!
   TBranch        *b_TrackJet15_SoftDroppedSubJet2;   //!
   TBranch        *b_TrackJet15_TrimmedP4;   //!
   TBranch        *b_TrackJet15_PrunedP4;   //!
   TBranch        *b_TrackJet15_SoftDroppedP4;   //!
   TBranch        *b_TrackJet15_NSubJetsTrimmed;   //!
   TBranch        *b_TrackJet15_NSubJetsPruned;   //!
   TBranch        *b_TrackJet15_NSubJetsSoftDropped;   //!
   TBranch        *b_TrackJet15_ExclYmerge23;   //!
   TBranch        *b_TrackJet15_ExclYmerge34;   //!
   TBranch        *b_TrackJet15_ExclYmerge45;   //!
   TBranch        *b_TrackJet15_ExclYmerge56;   //!
   TBranch        *b_TrackJet15_Constituents;   //!
   TBranch        *b_TrackJet15_Particles;   //!
   TBranch        *b_TrackJet15_Area;   //!
   TBranch        *b_TrackJet15_size;   //!
   TBranch        *b_MissingET_;   //!
   TBranch        *b_MissingET_fUniqueID;   //!
   TBranch        *b_MissingET_fBits;   //!
   TBranch        *b_MissingET_MET;   //!
   TBranch        *b_MissingET_Eta;   //!
   TBranch        *b_MissingET_Phi;   //!
   TBranch        *b_MissingET_size;   //!
   TBranch        *b_ScalarHT_;   //!
   TBranch        *b_ScalarHT_fUniqueID;   //!
   TBranch        *b_ScalarHT_fBits;   //!
   TBranch        *b_ScalarHT_HT;   //!
   TBranch        *b_ScalarHT_size;   //!

   Delphes(TTree *tree=0);
   virtual ~Delphes();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Delphes_cxx
Delphes::Delphes(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/work/project/physics/psriling/public/delphes_output_10.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/work/project/physics/psriling/public/delphes_output_10.root");
      }
      f->GetObject("Delphes",tree);

   }
   Init(tree);
}

Delphes::~Delphes()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Delphes::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Delphes::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Delphes::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event_, &b_Event_);
   fChain->SetBranchAddress("Event.fUniqueID", Event_fUniqueID, &b_Event_fUniqueID);
   fChain->SetBranchAddress("Event.fBits", Event_fBits, &b_Event_fBits);
   fChain->SetBranchAddress("Event.Number", Event_Number, &b_Event_Number);
   fChain->SetBranchAddress("Event.ReadTime", Event_ReadTime, &b_Event_ReadTime);
   fChain->SetBranchAddress("Event.ProcTime", Event_ProcTime, &b_Event_ProcTime);
   fChain->SetBranchAddress("Event.ProcessID", Event_ProcessID, &b_Event_ProcessID);
   fChain->SetBranchAddress("Event.MPI", Event_MPI, &b_Event_MPI);
   fChain->SetBranchAddress("Event.Weight", Event_Weight, &b_Event_Weight);
   fChain->SetBranchAddress("Event.CrossSection", Event_CrossSection, &b_Event_CrossSection);
   fChain->SetBranchAddress("Event.CrossSectionError", Event_CrossSectionError, &b_Event_CrossSectionError);
   fChain->SetBranchAddress("Event.Scale", Event_Scale, &b_Event_Scale);
   fChain->SetBranchAddress("Event.AlphaQED", Event_AlphaQED, &b_Event_AlphaQED);
   fChain->SetBranchAddress("Event.AlphaQCD", Event_AlphaQCD, &b_Event_AlphaQCD);
   fChain->SetBranchAddress("Event.ID1", Event_ID1, &b_Event_ID1);
   fChain->SetBranchAddress("Event.ID2", Event_ID2, &b_Event_ID2);
   fChain->SetBranchAddress("Event.X1", Event_X1, &b_Event_X1);
   fChain->SetBranchAddress("Event.X2", Event_X2, &b_Event_X2);
   fChain->SetBranchAddress("Event.ScalePDF", Event_ScalePDF, &b_Event_ScalePDF);
   fChain->SetBranchAddress("Event.PDF1", Event_PDF1, &b_Event_PDF1);
   fChain->SetBranchAddress("Event.PDF2", Event_PDF2, &b_Event_PDF2);
   fChain->SetBranchAddress("Event_size", &Event_size, &b_Event_size);
   fChain->SetBranchAddress("Weight", &Weight_, &b_Weight_);
   fChain->SetBranchAddress("Weight.fUniqueID", Weight_fUniqueID, &b_Weight_fUniqueID);
   fChain->SetBranchAddress("Weight.fBits", Weight_fBits, &b_Weight_fBits);
   fChain->SetBranchAddress("Weight.Weight", Weight_Weight, &b_Weight_Weight);
   fChain->SetBranchAddress("Weight_size", &Weight_size, &b_Weight_size);
   fChain->SetBranchAddress("Particle", &Particle_, &b_Particle_);
   fChain->SetBranchAddress("Particle.fUniqueID", Particle_fUniqueID, &b_Particle_fUniqueID);
   fChain->SetBranchAddress("Particle.fBits", Particle_fBits, &b_Particle_fBits);
   fChain->SetBranchAddress("Particle.PID", Particle_PID, &b_Particle_PID);
   fChain->SetBranchAddress("Particle.Status", Particle_Status, &b_Particle_Status);
   fChain->SetBranchAddress("Particle.IsPU", Particle_IsPU, &b_Particle_IsPU);
   fChain->SetBranchAddress("Particle.M1", Particle_M1, &b_Particle_M1);
   fChain->SetBranchAddress("Particle.M2", Particle_M2, &b_Particle_M2);
   fChain->SetBranchAddress("Particle.D1", Particle_D1, &b_Particle_D1);
   fChain->SetBranchAddress("Particle.D2", Particle_D2, &b_Particle_D2);
   fChain->SetBranchAddress("Particle.Charge", Particle_Charge, &b_Particle_Charge);
   fChain->SetBranchAddress("Particle.Mass", Particle_Mass, &b_Particle_Mass);
   fChain->SetBranchAddress("Particle.E", Particle_E, &b_Particle_E);
   fChain->SetBranchAddress("Particle.Px", Particle_Px, &b_Particle_Px);
   fChain->SetBranchAddress("Particle.Py", Particle_Py, &b_Particle_Py);
   fChain->SetBranchAddress("Particle.Pz", Particle_Pz, &b_Particle_Pz);
   fChain->SetBranchAddress("Particle.P", Particle_P, &b_Particle_P);
   fChain->SetBranchAddress("Particle.PT", Particle_PT, &b_Particle_PT);
   fChain->SetBranchAddress("Particle.Eta", Particle_Eta, &b_Particle_Eta);
   fChain->SetBranchAddress("Particle.Phi", Particle_Phi, &b_Particle_Phi);
   fChain->SetBranchAddress("Particle.Rapidity", Particle_Rapidity, &b_Particle_Rapidity);
   fChain->SetBranchAddress("Particle.T", Particle_T, &b_Particle_T);
   fChain->SetBranchAddress("Particle.X", Particle_X, &b_Particle_X);
   fChain->SetBranchAddress("Particle.Y", Particle_Y, &b_Particle_Y);
   fChain->SetBranchAddress("Particle.Z", Particle_Z, &b_Particle_Z);
   fChain->SetBranchAddress("Particle_size", &Particle_size, &b_Particle_size);
   fChain->SetBranchAddress("GenMissingET", &GenMissingET_, &b_GenMissingET_);
   fChain->SetBranchAddress("GenMissingET.fUniqueID", GenMissingET_fUniqueID, &b_GenMissingET_fUniqueID);
   fChain->SetBranchAddress("GenMissingET.fBits", GenMissingET_fBits, &b_GenMissingET_fBits);
   fChain->SetBranchAddress("GenMissingET.MET", GenMissingET_MET, &b_GenMissingET_MET);
   fChain->SetBranchAddress("GenMissingET.Eta", GenMissingET_Eta, &b_GenMissingET_Eta);
   fChain->SetBranchAddress("GenMissingET.Phi", GenMissingET_Phi, &b_GenMissingET_Phi);
   fChain->SetBranchAddress("GenMissingET_size", &GenMissingET_size, &b_GenMissingET_size);
   fChain->SetBranchAddress("EFlowTrack", &EFlowTrack_, &b_EFlowTrack_);
   fChain->SetBranchAddress("EFlowTrack.fUniqueID", EFlowTrack_fUniqueID, &b_EFlowTrack_fUniqueID);
   fChain->SetBranchAddress("EFlowTrack.fBits", EFlowTrack_fBits, &b_EFlowTrack_fBits);
   fChain->SetBranchAddress("EFlowTrack.PID", EFlowTrack_PID, &b_EFlowTrack_PID);
   fChain->SetBranchAddress("EFlowTrack.Charge", EFlowTrack_Charge, &b_EFlowTrack_Charge);
   fChain->SetBranchAddress("EFlowTrack.P", EFlowTrack_P, &b_EFlowTrack_P);
   fChain->SetBranchAddress("EFlowTrack.PT", EFlowTrack_PT, &b_EFlowTrack_PT);
   fChain->SetBranchAddress("EFlowTrack.Eta", EFlowTrack_Eta, &b_EFlowTrack_Eta);
   fChain->SetBranchAddress("EFlowTrack.Phi", EFlowTrack_Phi, &b_EFlowTrack_Phi);
   fChain->SetBranchAddress("EFlowTrack.CtgTheta", EFlowTrack_CtgTheta, &b_EFlowTrack_CtgTheta);
   fChain->SetBranchAddress("EFlowTrack.C", EFlowTrack_C, &b_EFlowTrack_C);
   fChain->SetBranchAddress("EFlowTrack.Mass", EFlowTrack_Mass, &b_EFlowTrack_Mass);
   fChain->SetBranchAddress("EFlowTrack.EtaOuter", EFlowTrack_EtaOuter, &b_EFlowTrack_EtaOuter);
   fChain->SetBranchAddress("EFlowTrack.PhiOuter", EFlowTrack_PhiOuter, &b_EFlowTrack_PhiOuter);
   fChain->SetBranchAddress("EFlowTrack.T", EFlowTrack_T, &b_EFlowTrack_T);
   fChain->SetBranchAddress("EFlowTrack.X", EFlowTrack_X, &b_EFlowTrack_X);
   fChain->SetBranchAddress("EFlowTrack.Y", EFlowTrack_Y, &b_EFlowTrack_Y);
   fChain->SetBranchAddress("EFlowTrack.Z", EFlowTrack_Z, &b_EFlowTrack_Z);
   fChain->SetBranchAddress("EFlowTrack.TOuter", EFlowTrack_TOuter, &b_EFlowTrack_TOuter);
   fChain->SetBranchAddress("EFlowTrack.XOuter", EFlowTrack_XOuter, &b_EFlowTrack_XOuter);
   fChain->SetBranchAddress("EFlowTrack.YOuter", EFlowTrack_YOuter, &b_EFlowTrack_YOuter);
   fChain->SetBranchAddress("EFlowTrack.ZOuter", EFlowTrack_ZOuter, &b_EFlowTrack_ZOuter);
   fChain->SetBranchAddress("EFlowTrack.Xd", EFlowTrack_Xd, &b_EFlowTrack_Xd);
   fChain->SetBranchAddress("EFlowTrack.Yd", EFlowTrack_Yd, &b_EFlowTrack_Yd);
   fChain->SetBranchAddress("EFlowTrack.Zd", EFlowTrack_Zd, &b_EFlowTrack_Zd);
   fChain->SetBranchAddress("EFlowTrack.L", EFlowTrack_L, &b_EFlowTrack_L);
   fChain->SetBranchAddress("EFlowTrack.D0", EFlowTrack_D0, &b_EFlowTrack_D0);
   fChain->SetBranchAddress("EFlowTrack.DZ", EFlowTrack_DZ, &b_EFlowTrack_DZ);
   fChain->SetBranchAddress("EFlowTrack.Nclusters", EFlowTrack_Nclusters, &b_EFlowTrack_Nclusters);
   fChain->SetBranchAddress("EFlowTrack.dNdx", EFlowTrack_dNdx, &b_EFlowTrack_dNdx);
   fChain->SetBranchAddress("EFlowTrack.ErrorP", EFlowTrack_ErrorP, &b_EFlowTrack_ErrorP);
   fChain->SetBranchAddress("EFlowTrack.ErrorPT", EFlowTrack_ErrorPT, &b_EFlowTrack_ErrorPT);
   fChain->SetBranchAddress("EFlowTrack.ErrorPhi", EFlowTrack_ErrorPhi, &b_EFlowTrack_ErrorPhi);
   fChain->SetBranchAddress("EFlowTrack.ErrorCtgTheta", EFlowTrack_ErrorCtgTheta, &b_EFlowTrack_ErrorCtgTheta);
   fChain->SetBranchAddress("EFlowTrack.ErrorT", EFlowTrack_ErrorT, &b_EFlowTrack_ErrorT);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0", EFlowTrack_ErrorD0, &b_EFlowTrack_ErrorD0);
   fChain->SetBranchAddress("EFlowTrack.ErrorDZ", EFlowTrack_ErrorDZ, &b_EFlowTrack_ErrorDZ);
   fChain->SetBranchAddress("EFlowTrack.ErrorC", EFlowTrack_ErrorC, &b_EFlowTrack_ErrorC);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0Phi", EFlowTrack_ErrorD0Phi, &b_EFlowTrack_ErrorD0Phi);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0C", EFlowTrack_ErrorD0C, &b_EFlowTrack_ErrorD0C);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0DZ", EFlowTrack_ErrorD0DZ, &b_EFlowTrack_ErrorD0DZ);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0CtgTheta", EFlowTrack_ErrorD0CtgTheta, &b_EFlowTrack_ErrorD0CtgTheta);
   fChain->SetBranchAddress("EFlowTrack.ErrorPhiC", EFlowTrack_ErrorPhiC, &b_EFlowTrack_ErrorPhiC);
   fChain->SetBranchAddress("EFlowTrack.ErrorPhiDZ", EFlowTrack_ErrorPhiDZ, &b_EFlowTrack_ErrorPhiDZ);
   fChain->SetBranchAddress("EFlowTrack.ErrorPhiCtgTheta", EFlowTrack_ErrorPhiCtgTheta, &b_EFlowTrack_ErrorPhiCtgTheta);
   fChain->SetBranchAddress("EFlowTrack.ErrorCDZ", EFlowTrack_ErrorCDZ, &b_EFlowTrack_ErrorCDZ);
   fChain->SetBranchAddress("EFlowTrack.ErrorCCtgTheta", EFlowTrack_ErrorCCtgTheta, &b_EFlowTrack_ErrorCCtgTheta);
   fChain->SetBranchAddress("EFlowTrack.ErrorDZCtgTheta", EFlowTrack_ErrorDZCtgTheta, &b_EFlowTrack_ErrorDZCtgTheta);
   fChain->SetBranchAddress("EFlowTrack.Particle", EFlowTrack_Particle, &b_EFlowTrack_Particle);
   fChain->SetBranchAddress("EFlowTrack.VertexIndex", EFlowTrack_VertexIndex, &b_EFlowTrack_VertexIndex);
   fChain->SetBranchAddress("EFlowTrack_size", &EFlowTrack_size, &b_EFlowTrack_size);
   fChain->SetBranchAddress("EFlowPhoton", &EFlowPhoton_, &b_EFlowPhoton_);
   fChain->SetBranchAddress("EFlowPhoton.fUniqueID", EFlowPhoton_fUniqueID, &b_EFlowPhoton_fUniqueID);
   fChain->SetBranchAddress("EFlowPhoton.fBits", EFlowPhoton_fBits, &b_EFlowPhoton_fBits);
   fChain->SetBranchAddress("EFlowPhoton.ET", EFlowPhoton_ET, &b_EFlowPhoton_ET);
   fChain->SetBranchAddress("EFlowPhoton.Eta", EFlowPhoton_Eta, &b_EFlowPhoton_Eta);
   fChain->SetBranchAddress("EFlowPhoton.Phi", EFlowPhoton_Phi, &b_EFlowPhoton_Phi);
   fChain->SetBranchAddress("EFlowPhoton.E", EFlowPhoton_E, &b_EFlowPhoton_E);
   fChain->SetBranchAddress("EFlowPhoton.T", EFlowPhoton_T, &b_EFlowPhoton_T);
   fChain->SetBranchAddress("EFlowPhoton.NTimeHits", EFlowPhoton_NTimeHits, &b_EFlowPhoton_NTimeHits);
   fChain->SetBranchAddress("EFlowPhoton.Eem", EFlowPhoton_Eem, &b_EFlowPhoton_Eem);
   fChain->SetBranchAddress("EFlowPhoton.Ehad", EFlowPhoton_Ehad, &b_EFlowPhoton_Ehad);
   fChain->SetBranchAddress("EFlowPhoton.Edges[4]", EFlowPhoton_Edges, &b_EFlowPhoton_Edges);
   fChain->SetBranchAddress("EFlowPhoton.Particles", EFlowPhoton_Particles, &b_EFlowPhoton_Particles);
   fChain->SetBranchAddress("EFlowPhoton_size", &EFlowPhoton_size, &b_EFlowPhoton_size);
   fChain->SetBranchAddress("EFlowNeutralHadron", &EFlowNeutralHadron_, &b_EFlowNeutralHadron_);
   fChain->SetBranchAddress("EFlowNeutralHadron.fUniqueID", EFlowNeutralHadron_fUniqueID, &b_EFlowNeutralHadron_fUniqueID);
   fChain->SetBranchAddress("EFlowNeutralHadron.fBits", EFlowNeutralHadron_fBits, &b_EFlowNeutralHadron_fBits);
   fChain->SetBranchAddress("EFlowNeutralHadron.ET", EFlowNeutralHadron_ET, &b_EFlowNeutralHadron_ET);
   fChain->SetBranchAddress("EFlowNeutralHadron.Eta", EFlowNeutralHadron_Eta, &b_EFlowNeutralHadron_Eta);
   fChain->SetBranchAddress("EFlowNeutralHadron.Phi", EFlowNeutralHadron_Phi, &b_EFlowNeutralHadron_Phi);
   fChain->SetBranchAddress("EFlowNeutralHadron.E", EFlowNeutralHadron_E, &b_EFlowNeutralHadron_E);
   fChain->SetBranchAddress("EFlowNeutralHadron.T", EFlowNeutralHadron_T, &b_EFlowNeutralHadron_T);
   fChain->SetBranchAddress("EFlowNeutralHadron.NTimeHits", EFlowNeutralHadron_NTimeHits, &b_EFlowNeutralHadron_NTimeHits);
   fChain->SetBranchAddress("EFlowNeutralHadron.Eem", EFlowNeutralHadron_Eem, &b_EFlowNeutralHadron_Eem);
   fChain->SetBranchAddress("EFlowNeutralHadron.Ehad", EFlowNeutralHadron_Ehad, &b_EFlowNeutralHadron_Ehad);
   fChain->SetBranchAddress("EFlowNeutralHadron.Edges[4]", EFlowNeutralHadron_Edges, &b_EFlowNeutralHadron_Edges);
   fChain->SetBranchAddress("EFlowNeutralHadron.Particles", EFlowNeutralHadron_Particles, &b_EFlowNeutralHadron_Particles);
   fChain->SetBranchAddress("EFlowNeutralHadron_size", &EFlowNeutralHadron_size, &b_EFlowNeutralHadron_size);
   fChain->SetBranchAddress("Photon", &Photon_, &b_Photon_);
   fChain->SetBranchAddress("Photon.fUniqueID", Photon_fUniqueID, &b_Photon_fUniqueID);
   fChain->SetBranchAddress("Photon.fBits", Photon_fBits, &b_Photon_fBits);
   fChain->SetBranchAddress("Photon.PT", Photon_PT, &b_Photon_PT);
   fChain->SetBranchAddress("Photon.Eta", Photon_Eta, &b_Photon_Eta);
   fChain->SetBranchAddress("Photon.Phi", Photon_Phi, &b_Photon_Phi);
   fChain->SetBranchAddress("Photon.E", Photon_E, &b_Photon_E);
   fChain->SetBranchAddress("Photon.T", Photon_T, &b_Photon_T);
   fChain->SetBranchAddress("Photon.EhadOverEem", Photon_EhadOverEem, &b_Photon_EhadOverEem);
   fChain->SetBranchAddress("Photon.Particles", Photon_Particles, &b_Photon_Particles);
   fChain->SetBranchAddress("Photon.IsolationVar", Photon_IsolationVar, &b_Photon_IsolationVar);
   fChain->SetBranchAddress("Photon.IsolationVarRhoCorr", Photon_IsolationVarRhoCorr, &b_Photon_IsolationVarRhoCorr);
   fChain->SetBranchAddress("Photon.SumPtCharged", Photon_SumPtCharged, &b_Photon_SumPtCharged);
   fChain->SetBranchAddress("Photon.SumPtNeutral", Photon_SumPtNeutral, &b_Photon_SumPtNeutral);
   fChain->SetBranchAddress("Photon.SumPtChargedPU", Photon_SumPtChargedPU, &b_Photon_SumPtChargedPU);
   fChain->SetBranchAddress("Photon.SumPt", Photon_SumPt, &b_Photon_SumPt);
   fChain->SetBranchAddress("Photon.Status", Photon_Status, &b_Photon_Status);
   fChain->SetBranchAddress("Photon_size", &Photon_size, &b_Photon_size);
   fChain->SetBranchAddress("Electron", &Electron_, &b_Electron_);
   fChain->SetBranchAddress("Electron.fUniqueID", Electron_fUniqueID, &b_Electron_fUniqueID);
   fChain->SetBranchAddress("Electron.fBits", Electron_fBits, &b_Electron_fBits);
   fChain->SetBranchAddress("Electron.PT", Electron_PT, &b_Electron_PT);
   fChain->SetBranchAddress("Electron.Eta", Electron_Eta, &b_Electron_Eta);
   fChain->SetBranchAddress("Electron.Phi", Electron_Phi, &b_Electron_Phi);
   fChain->SetBranchAddress("Electron.T", Electron_T, &b_Electron_T);
   fChain->SetBranchAddress("Electron.Charge", Electron_Charge, &b_Electron_Charge);
   fChain->SetBranchAddress("Electron.EhadOverEem", Electron_EhadOverEem, &b_Electron_EhadOverEem);
   fChain->SetBranchAddress("Electron.Particle", Electron_Particle, &b_Electron_Particle);
   fChain->SetBranchAddress("Electron.IsolationVar", Electron_IsolationVar, &b_Electron_IsolationVar);
   fChain->SetBranchAddress("Electron.IsolationVarRhoCorr", Electron_IsolationVarRhoCorr, &b_Electron_IsolationVarRhoCorr);
   fChain->SetBranchAddress("Electron.SumPtCharged", Electron_SumPtCharged, &b_Electron_SumPtCharged);
   fChain->SetBranchAddress("Electron.SumPtNeutral", Electron_SumPtNeutral, &b_Electron_SumPtNeutral);
   fChain->SetBranchAddress("Electron.SumPtChargedPU", Electron_SumPtChargedPU, &b_Electron_SumPtChargedPU);
   fChain->SetBranchAddress("Electron.SumPt", Electron_SumPt, &b_Electron_SumPt);
   fChain->SetBranchAddress("Electron.D0", Electron_D0, &b_Electron_D0);
   fChain->SetBranchAddress("Electron.DZ", Electron_DZ, &b_Electron_DZ);
   fChain->SetBranchAddress("Electron.ErrorD0", Electron_ErrorD0, &b_Electron_ErrorD0);
   fChain->SetBranchAddress("Electron.ErrorDZ", Electron_ErrorDZ, &b_Electron_ErrorDZ);
   fChain->SetBranchAddress("Electron_size", &Electron_size, &b_Electron_size);
   fChain->SetBranchAddress("Muon", &Muon_, &b_Muon_);
   fChain->SetBranchAddress("Muon.fUniqueID", Muon_fUniqueID, &b_Muon_fUniqueID);
   fChain->SetBranchAddress("Muon.fBits", Muon_fBits, &b_Muon_fBits);
   fChain->SetBranchAddress("Muon.PT", Muon_PT, &b_Muon_PT);
   fChain->SetBranchAddress("Muon.Eta", Muon_Eta, &b_Muon_Eta);
   fChain->SetBranchAddress("Muon.Phi", Muon_Phi, &b_Muon_Phi);
   fChain->SetBranchAddress("Muon.T", Muon_T, &b_Muon_T);
   fChain->SetBranchAddress("Muon.Charge", Muon_Charge, &b_Muon_Charge);
   fChain->SetBranchAddress("Muon.Particle", Muon_Particle, &b_Muon_Particle);
   fChain->SetBranchAddress("Muon.IsolationVar", Muon_IsolationVar, &b_Muon_IsolationVar);
   fChain->SetBranchAddress("Muon.IsolationVarRhoCorr", Muon_IsolationVarRhoCorr, &b_Muon_IsolationVarRhoCorr);
   fChain->SetBranchAddress("Muon.SumPtCharged", Muon_SumPtCharged, &b_Muon_SumPtCharged);
   fChain->SetBranchAddress("Muon.SumPtNeutral", Muon_SumPtNeutral, &b_Muon_SumPtNeutral);
   fChain->SetBranchAddress("Muon.SumPtChargedPU", Muon_SumPtChargedPU, &b_Muon_SumPtChargedPU);
   fChain->SetBranchAddress("Muon.SumPt", Muon_SumPt, &b_Muon_SumPt);
   fChain->SetBranchAddress("Muon.D0", Muon_D0, &b_Muon_D0);
   fChain->SetBranchAddress("Muon.DZ", Muon_DZ, &b_Muon_DZ);
   fChain->SetBranchAddress("Muon.ErrorD0", Muon_ErrorD0, &b_Muon_ErrorD0);
   fChain->SetBranchAddress("Muon.ErrorDZ", Muon_ErrorDZ, &b_Muon_ErrorDZ);
   fChain->SetBranchAddress("Muon_size", &Muon_size, &b_Muon_size);
   fChain->SetBranchAddress("Jet", &Jet_, &b_Jet_);
   fChain->SetBranchAddress("Jet.fUniqueID", Jet_fUniqueID, &b_Jet_fUniqueID);
   fChain->SetBranchAddress("Jet.fBits", Jet_fBits, &b_Jet_fBits);
   fChain->SetBranchAddress("Jet.PT", Jet_PT, &b_Jet_PT);
   fChain->SetBranchAddress("Jet.Eta", Jet_Eta, &b_Jet_Eta);
   fChain->SetBranchAddress("Jet.Phi", Jet_Phi, &b_Jet_Phi);
   fChain->SetBranchAddress("Jet.T", Jet_T, &b_Jet_T);
   fChain->SetBranchAddress("Jet.Mass", Jet_Mass, &b_Jet_Mass);
   fChain->SetBranchAddress("Jet.DeltaEta", Jet_DeltaEta, &b_Jet_DeltaEta);
   fChain->SetBranchAddress("Jet.DeltaPhi", Jet_DeltaPhi, &b_Jet_DeltaPhi);
   fChain->SetBranchAddress("Jet.Flavor", Jet_Flavor, &b_Jet_Flavor);
   fChain->SetBranchAddress("Jet.FlavorAlgo", Jet_FlavorAlgo, &b_Jet_FlavorAlgo);
   fChain->SetBranchAddress("Jet.FlavorPhys", Jet_FlavorPhys, &b_Jet_FlavorPhys);
   fChain->SetBranchAddress("Jet.BTag", Jet_BTag, &b_Jet_BTag);
   fChain->SetBranchAddress("Jet.BTagAlgo", Jet_BTagAlgo, &b_Jet_BTagAlgo);
   fChain->SetBranchAddress("Jet.BTagPhys", Jet_BTagPhys, &b_Jet_BTagPhys);
   fChain->SetBranchAddress("Jet.TauTag", Jet_TauTag, &b_Jet_TauTag);
   fChain->SetBranchAddress("Jet.TauWeight", Jet_TauWeight, &b_Jet_TauWeight);
   fChain->SetBranchAddress("Jet.Charge", Jet_Charge, &b_Jet_Charge);
   fChain->SetBranchAddress("Jet.EhadOverEem", Jet_EhadOverEem, &b_Jet_EhadOverEem);
   fChain->SetBranchAddress("Jet.NCharged", Jet_NCharged, &b_Jet_NCharged);
   fChain->SetBranchAddress("Jet.NNeutrals", Jet_NNeutrals, &b_Jet_NNeutrals);
   fChain->SetBranchAddress("Jet.NeutralEnergyFraction", Jet_NeutralEnergyFraction, &b_Jet_NeutralEnergyFraction);
   fChain->SetBranchAddress("Jet.ChargedEnergyFraction", Jet_ChargedEnergyFraction, &b_Jet_ChargedEnergyFraction);
   fChain->SetBranchAddress("Jet.Beta", Jet_Beta, &b_Jet_Beta);
   fChain->SetBranchAddress("Jet.BetaStar", Jet_BetaStar, &b_Jet_BetaStar);
   fChain->SetBranchAddress("Jet.MeanSqDeltaR", Jet_MeanSqDeltaR, &b_Jet_MeanSqDeltaR);
   fChain->SetBranchAddress("Jet.PTD", Jet_PTD, &b_Jet_PTD);
   fChain->SetBranchAddress("Jet.FracPt[5]", Jet_FracPt, &b_Jet_FracPt);
   fChain->SetBranchAddress("Jet.Tau[5]", Jet_Tau, &b_Jet_Tau);
   fChain->SetBranchAddress("Jet.SoftDroppedJet", Jet_SoftDroppedJet, &b_Jet_SoftDroppedJet);
   fChain->SetBranchAddress("Jet.SoftDroppedSubJet1", Jet_SoftDroppedSubJet1, &b_Jet_SoftDroppedSubJet1);
   fChain->SetBranchAddress("Jet.SoftDroppedSubJet2", Jet_SoftDroppedSubJet2, &b_Jet_SoftDroppedSubJet2);
   fChain->SetBranchAddress("Jet.TrimmedP4[5]", Jet_TrimmedP4, &b_Jet_TrimmedP4);
   fChain->SetBranchAddress("Jet.PrunedP4[5]", Jet_PrunedP4, &b_Jet_PrunedP4);
   fChain->SetBranchAddress("Jet.SoftDroppedP4[5]", Jet_SoftDroppedP4, &b_Jet_SoftDroppedP4);
   fChain->SetBranchAddress("Jet.NSubJetsTrimmed", Jet_NSubJetsTrimmed, &b_Jet_NSubJetsTrimmed);
   fChain->SetBranchAddress("Jet.NSubJetsPruned", Jet_NSubJetsPruned, &b_Jet_NSubJetsPruned);
   fChain->SetBranchAddress("Jet.NSubJetsSoftDropped", Jet_NSubJetsSoftDropped, &b_Jet_NSubJetsSoftDropped);
   fChain->SetBranchAddress("Jet.ExclYmerge23", Jet_ExclYmerge23, &b_Jet_ExclYmerge23);
   fChain->SetBranchAddress("Jet.ExclYmerge34", Jet_ExclYmerge34, &b_Jet_ExclYmerge34);
   fChain->SetBranchAddress("Jet.ExclYmerge45", Jet_ExclYmerge45, &b_Jet_ExclYmerge45);
   fChain->SetBranchAddress("Jet.ExclYmerge56", Jet_ExclYmerge56, &b_Jet_ExclYmerge56);
   fChain->SetBranchAddress("Jet.Constituents", Jet_Constituents, &b_Jet_Constituents);
   fChain->SetBranchAddress("Jet.Particles", Jet_Particles, &b_Jet_Particles);
   fChain->SetBranchAddress("Jet.Area", Jet_Area, &b_Jet_Area);
   fChain->SetBranchAddress("Jet_size", &Jet_size, &b_Jet_size);
   fChain->SetBranchAddress("PhotonNoIso", &PhotonNoIso_, &b_PhotonNoIso_);
   fChain->SetBranchAddress("PhotonNoIso.fUniqueID", PhotonNoIso_fUniqueID, &b_PhotonNoIso_fUniqueID);
   fChain->SetBranchAddress("PhotonNoIso.fBits", PhotonNoIso_fBits, &b_PhotonNoIso_fBits);
   fChain->SetBranchAddress("PhotonNoIso.PT", PhotonNoIso_PT, &b_PhotonNoIso_PT);
   fChain->SetBranchAddress("PhotonNoIso.Eta", PhotonNoIso_Eta, &b_PhotonNoIso_Eta);
   fChain->SetBranchAddress("PhotonNoIso.Phi", PhotonNoIso_Phi, &b_PhotonNoIso_Phi);
   fChain->SetBranchAddress("PhotonNoIso.E", PhotonNoIso_E, &b_PhotonNoIso_E);
   fChain->SetBranchAddress("PhotonNoIso.T", PhotonNoIso_T, &b_PhotonNoIso_T);
   fChain->SetBranchAddress("PhotonNoIso.EhadOverEem", PhotonNoIso_EhadOverEem, &b_PhotonNoIso_EhadOverEem);
   fChain->SetBranchAddress("PhotonNoIso.Particles", PhotonNoIso_Particles, &b_PhotonNoIso_Particles);
   fChain->SetBranchAddress("PhotonNoIso.IsolationVar", PhotonNoIso_IsolationVar, &b_PhotonNoIso_IsolationVar);
   fChain->SetBranchAddress("PhotonNoIso.IsolationVarRhoCorr", PhotonNoIso_IsolationVarRhoCorr, &b_PhotonNoIso_IsolationVarRhoCorr);
   fChain->SetBranchAddress("PhotonNoIso.SumPtCharged", PhotonNoIso_SumPtCharged, &b_PhotonNoIso_SumPtCharged);
   fChain->SetBranchAddress("PhotonNoIso.SumPtNeutral", PhotonNoIso_SumPtNeutral, &b_PhotonNoIso_SumPtNeutral);
   fChain->SetBranchAddress("PhotonNoIso.SumPtChargedPU", PhotonNoIso_SumPtChargedPU, &b_PhotonNoIso_SumPtChargedPU);
   fChain->SetBranchAddress("PhotonNoIso.SumPt", PhotonNoIso_SumPt, &b_PhotonNoIso_SumPt);
   fChain->SetBranchAddress("PhotonNoIso.Status", PhotonNoIso_Status, &b_PhotonNoIso_Status);
   fChain->SetBranchAddress("PhotonNoIso_size", &PhotonNoIso_size, &b_PhotonNoIso_size);
   fChain->SetBranchAddress("ElectronNoIso", &ElectronNoIso_, &b_ElectronNoIso_);
   fChain->SetBranchAddress("ElectronNoIso.fUniqueID", ElectronNoIso_fUniqueID, &b_ElectronNoIso_fUniqueID);
   fChain->SetBranchAddress("ElectronNoIso.fBits", ElectronNoIso_fBits, &b_ElectronNoIso_fBits);
   fChain->SetBranchAddress("ElectronNoIso.PT", ElectronNoIso_PT, &b_ElectronNoIso_PT);
   fChain->SetBranchAddress("ElectronNoIso.Eta", ElectronNoIso_Eta, &b_ElectronNoIso_Eta);
   fChain->SetBranchAddress("ElectronNoIso.Phi", ElectronNoIso_Phi, &b_ElectronNoIso_Phi);
   fChain->SetBranchAddress("ElectronNoIso.T", ElectronNoIso_T, &b_ElectronNoIso_T);
   fChain->SetBranchAddress("ElectronNoIso.Charge", ElectronNoIso_Charge, &b_ElectronNoIso_Charge);
   fChain->SetBranchAddress("ElectronNoIso.EhadOverEem", ElectronNoIso_EhadOverEem, &b_ElectronNoIso_EhadOverEem);
   fChain->SetBranchAddress("ElectronNoIso.Particle", ElectronNoIso_Particle, &b_ElectronNoIso_Particle);
   fChain->SetBranchAddress("ElectronNoIso.IsolationVar", ElectronNoIso_IsolationVar, &b_ElectronNoIso_IsolationVar);
   fChain->SetBranchAddress("ElectronNoIso.IsolationVarRhoCorr", ElectronNoIso_IsolationVarRhoCorr, &b_ElectronNoIso_IsolationVarRhoCorr);
   fChain->SetBranchAddress("ElectronNoIso.SumPtCharged", ElectronNoIso_SumPtCharged, &b_ElectronNoIso_SumPtCharged);
   fChain->SetBranchAddress("ElectronNoIso.SumPtNeutral", ElectronNoIso_SumPtNeutral, &b_ElectronNoIso_SumPtNeutral);
   fChain->SetBranchAddress("ElectronNoIso.SumPtChargedPU", ElectronNoIso_SumPtChargedPU, &b_ElectronNoIso_SumPtChargedPU);
   fChain->SetBranchAddress("ElectronNoIso.SumPt", ElectronNoIso_SumPt, &b_ElectronNoIso_SumPt);
   fChain->SetBranchAddress("ElectronNoIso.D0", ElectronNoIso_D0, &b_ElectronNoIso_D0);
   fChain->SetBranchAddress("ElectronNoIso.DZ", ElectronNoIso_DZ, &b_ElectronNoIso_DZ);
   fChain->SetBranchAddress("ElectronNoIso.ErrorD0", ElectronNoIso_ErrorD0, &b_ElectronNoIso_ErrorD0);
   fChain->SetBranchAddress("ElectronNoIso.ErrorDZ", ElectronNoIso_ErrorDZ, &b_ElectronNoIso_ErrorDZ);
   fChain->SetBranchAddress("ElectronNoIso_size", &ElectronNoIso_size, &b_ElectronNoIso_size);
   fChain->SetBranchAddress("MuonNoIso", &MuonNoIso_, &b_MuonNoIso_);
   fChain->SetBranchAddress("MuonNoIso.fUniqueID", MuonNoIso_fUniqueID, &b_MuonNoIso_fUniqueID);
   fChain->SetBranchAddress("MuonNoIso.fBits", MuonNoIso_fBits, &b_MuonNoIso_fBits);
   fChain->SetBranchAddress("MuonNoIso.PT", MuonNoIso_PT, &b_MuonNoIso_PT);
   fChain->SetBranchAddress("MuonNoIso.Eta", MuonNoIso_Eta, &b_MuonNoIso_Eta);
   fChain->SetBranchAddress("MuonNoIso.Phi", MuonNoIso_Phi, &b_MuonNoIso_Phi);
   fChain->SetBranchAddress("MuonNoIso.T", MuonNoIso_T, &b_MuonNoIso_T);
   fChain->SetBranchAddress("MuonNoIso.Charge", MuonNoIso_Charge, &b_MuonNoIso_Charge);
   fChain->SetBranchAddress("MuonNoIso.Particle", MuonNoIso_Particle, &b_MuonNoIso_Particle);
   fChain->SetBranchAddress("MuonNoIso.IsolationVar", MuonNoIso_IsolationVar, &b_MuonNoIso_IsolationVar);
   fChain->SetBranchAddress("MuonNoIso.IsolationVarRhoCorr", MuonNoIso_IsolationVarRhoCorr, &b_MuonNoIso_IsolationVarRhoCorr);
   fChain->SetBranchAddress("MuonNoIso.SumPtCharged", MuonNoIso_SumPtCharged, &b_MuonNoIso_SumPtCharged);
   fChain->SetBranchAddress("MuonNoIso.SumPtNeutral", MuonNoIso_SumPtNeutral, &b_MuonNoIso_SumPtNeutral);
   fChain->SetBranchAddress("MuonNoIso.SumPtChargedPU", MuonNoIso_SumPtChargedPU, &b_MuonNoIso_SumPtChargedPU);
   fChain->SetBranchAddress("MuonNoIso.SumPt", MuonNoIso_SumPt, &b_MuonNoIso_SumPt);
   fChain->SetBranchAddress("MuonNoIso.D0", MuonNoIso_D0, &b_MuonNoIso_D0);
   fChain->SetBranchAddress("MuonNoIso.DZ", MuonNoIso_DZ, &b_MuonNoIso_DZ);
   fChain->SetBranchAddress("MuonNoIso.ErrorD0", MuonNoIso_ErrorD0, &b_MuonNoIso_ErrorD0);
   fChain->SetBranchAddress("MuonNoIso.ErrorDZ", MuonNoIso_ErrorDZ, &b_MuonNoIso_ErrorDZ);
   fChain->SetBranchAddress("MuonNoIso_size", &MuonNoIso_size, &b_MuonNoIso_size);
   fChain->SetBranchAddress("JetNoIso", &JetNoIso_, &b_JetNoIso_);
   fChain->SetBranchAddress("JetNoIso.fUniqueID", JetNoIso_fUniqueID, &b_JetNoIso_fUniqueID);
   fChain->SetBranchAddress("JetNoIso.fBits", JetNoIso_fBits, &b_JetNoIso_fBits);
   fChain->SetBranchAddress("JetNoIso.PT", JetNoIso_PT, &b_JetNoIso_PT);
   fChain->SetBranchAddress("JetNoIso.Eta", JetNoIso_Eta, &b_JetNoIso_Eta);
   fChain->SetBranchAddress("JetNoIso.Phi", JetNoIso_Phi, &b_JetNoIso_Phi);
   fChain->SetBranchAddress("JetNoIso.T", JetNoIso_T, &b_JetNoIso_T);
   fChain->SetBranchAddress("JetNoIso.Mass", JetNoIso_Mass, &b_JetNoIso_Mass);
   fChain->SetBranchAddress("JetNoIso.DeltaEta", JetNoIso_DeltaEta, &b_JetNoIso_DeltaEta);
   fChain->SetBranchAddress("JetNoIso.DeltaPhi", JetNoIso_DeltaPhi, &b_JetNoIso_DeltaPhi);
   fChain->SetBranchAddress("JetNoIso.Flavor", JetNoIso_Flavor, &b_JetNoIso_Flavor);
   fChain->SetBranchAddress("JetNoIso.FlavorAlgo", JetNoIso_FlavorAlgo, &b_JetNoIso_FlavorAlgo);
   fChain->SetBranchAddress("JetNoIso.FlavorPhys", JetNoIso_FlavorPhys, &b_JetNoIso_FlavorPhys);
   fChain->SetBranchAddress("JetNoIso.BTag", JetNoIso_BTag, &b_JetNoIso_BTag);
   fChain->SetBranchAddress("JetNoIso.BTagAlgo", JetNoIso_BTagAlgo, &b_JetNoIso_BTagAlgo);
   fChain->SetBranchAddress("JetNoIso.BTagPhys", JetNoIso_BTagPhys, &b_JetNoIso_BTagPhys);
   fChain->SetBranchAddress("JetNoIso.TauTag", JetNoIso_TauTag, &b_JetNoIso_TauTag);
   fChain->SetBranchAddress("JetNoIso.TauWeight", JetNoIso_TauWeight, &b_JetNoIso_TauWeight);
   fChain->SetBranchAddress("JetNoIso.Charge", JetNoIso_Charge, &b_JetNoIso_Charge);
   fChain->SetBranchAddress("JetNoIso.EhadOverEem", JetNoIso_EhadOverEem, &b_JetNoIso_EhadOverEem);
   fChain->SetBranchAddress("JetNoIso.NCharged", JetNoIso_NCharged, &b_JetNoIso_NCharged);
   fChain->SetBranchAddress("JetNoIso.NNeutrals", JetNoIso_NNeutrals, &b_JetNoIso_NNeutrals);
   fChain->SetBranchAddress("JetNoIso.NeutralEnergyFraction", JetNoIso_NeutralEnergyFraction, &b_JetNoIso_NeutralEnergyFraction);
   fChain->SetBranchAddress("JetNoIso.ChargedEnergyFraction", JetNoIso_ChargedEnergyFraction, &b_JetNoIso_ChargedEnergyFraction);
   fChain->SetBranchAddress("JetNoIso.Beta", JetNoIso_Beta, &b_JetNoIso_Beta);
   fChain->SetBranchAddress("JetNoIso.BetaStar", JetNoIso_BetaStar, &b_JetNoIso_BetaStar);
   fChain->SetBranchAddress("JetNoIso.MeanSqDeltaR", JetNoIso_MeanSqDeltaR, &b_JetNoIso_MeanSqDeltaR);
   fChain->SetBranchAddress("JetNoIso.PTD", JetNoIso_PTD, &b_JetNoIso_PTD);
   fChain->SetBranchAddress("JetNoIso.FracPt[5]", JetNoIso_FracPt, &b_JetNoIso_FracPt);
   fChain->SetBranchAddress("JetNoIso.Tau[5]", JetNoIso_Tau, &b_JetNoIso_Tau);
   fChain->SetBranchAddress("JetNoIso.SoftDroppedJet", JetNoIso_SoftDroppedJet, &b_JetNoIso_SoftDroppedJet);
   fChain->SetBranchAddress("JetNoIso.SoftDroppedSubJet1", JetNoIso_SoftDroppedSubJet1, &b_JetNoIso_SoftDroppedSubJet1);
   fChain->SetBranchAddress("JetNoIso.SoftDroppedSubJet2", JetNoIso_SoftDroppedSubJet2, &b_JetNoIso_SoftDroppedSubJet2);
   fChain->SetBranchAddress("JetNoIso.TrimmedP4[5]", JetNoIso_TrimmedP4, &b_JetNoIso_TrimmedP4);
   fChain->SetBranchAddress("JetNoIso.PrunedP4[5]", JetNoIso_PrunedP4, &b_JetNoIso_PrunedP4);
   fChain->SetBranchAddress("JetNoIso.SoftDroppedP4[5]", JetNoIso_SoftDroppedP4, &b_JetNoIso_SoftDroppedP4);
   fChain->SetBranchAddress("JetNoIso.NSubJetsTrimmed", JetNoIso_NSubJetsTrimmed, &b_JetNoIso_NSubJetsTrimmed);
   fChain->SetBranchAddress("JetNoIso.NSubJetsPruned", JetNoIso_NSubJetsPruned, &b_JetNoIso_NSubJetsPruned);
   fChain->SetBranchAddress("JetNoIso.NSubJetsSoftDropped", JetNoIso_NSubJetsSoftDropped, &b_JetNoIso_NSubJetsSoftDropped);
   fChain->SetBranchAddress("JetNoIso.ExclYmerge23", JetNoIso_ExclYmerge23, &b_JetNoIso_ExclYmerge23);
   fChain->SetBranchAddress("JetNoIso.ExclYmerge34", JetNoIso_ExclYmerge34, &b_JetNoIso_ExclYmerge34);
   fChain->SetBranchAddress("JetNoIso.ExclYmerge45", JetNoIso_ExclYmerge45, &b_JetNoIso_ExclYmerge45);
   fChain->SetBranchAddress("JetNoIso.ExclYmerge56", JetNoIso_ExclYmerge56, &b_JetNoIso_ExclYmerge56);
   fChain->SetBranchAddress("JetNoIso.Constituents", JetNoIso_Constituents, &b_JetNoIso_Constituents);
   fChain->SetBranchAddress("JetNoIso.Particles", JetNoIso_Particles, &b_JetNoIso_Particles);
   fChain->SetBranchAddress("JetNoIso.Area", JetNoIso_Area, &b_JetNoIso_Area);
   fChain->SetBranchAddress("JetNoIso_size", &JetNoIso_size, &b_JetNoIso_size);
   fChain->SetBranchAddress("GenJet02", &GenJet02_, &b_GenJet02_);
   fChain->SetBranchAddress("GenJet02.fUniqueID", GenJet02_fUniqueID, &b_GenJet02_fUniqueID);
   fChain->SetBranchAddress("GenJet02.fBits", GenJet02_fBits, &b_GenJet02_fBits);
   fChain->SetBranchAddress("GenJet02.PT", GenJet02_PT, &b_GenJet02_PT);
   fChain->SetBranchAddress("GenJet02.Eta", GenJet02_Eta, &b_GenJet02_Eta);
   fChain->SetBranchAddress("GenJet02.Phi", GenJet02_Phi, &b_GenJet02_Phi);
   fChain->SetBranchAddress("GenJet02.T", GenJet02_T, &b_GenJet02_T);
   fChain->SetBranchAddress("GenJet02.Mass", GenJet02_Mass, &b_GenJet02_Mass);
   fChain->SetBranchAddress("GenJet02.DeltaEta", GenJet02_DeltaEta, &b_GenJet02_DeltaEta);
   fChain->SetBranchAddress("GenJet02.DeltaPhi", GenJet02_DeltaPhi, &b_GenJet02_DeltaPhi);
   fChain->SetBranchAddress("GenJet02.Flavor", GenJet02_Flavor, &b_GenJet02_Flavor);
   fChain->SetBranchAddress("GenJet02.FlavorAlgo", GenJet02_FlavorAlgo, &b_GenJet02_FlavorAlgo);
   fChain->SetBranchAddress("GenJet02.FlavorPhys", GenJet02_FlavorPhys, &b_GenJet02_FlavorPhys);
   fChain->SetBranchAddress("GenJet02.BTag", GenJet02_BTag, &b_GenJet02_BTag);
   fChain->SetBranchAddress("GenJet02.BTagAlgo", GenJet02_BTagAlgo, &b_GenJet02_BTagAlgo);
   fChain->SetBranchAddress("GenJet02.BTagPhys", GenJet02_BTagPhys, &b_GenJet02_BTagPhys);
   fChain->SetBranchAddress("GenJet02.TauTag", GenJet02_TauTag, &b_GenJet02_TauTag);
   fChain->SetBranchAddress("GenJet02.TauWeight", GenJet02_TauWeight, &b_GenJet02_TauWeight);
   fChain->SetBranchAddress("GenJet02.Charge", GenJet02_Charge, &b_GenJet02_Charge);
   fChain->SetBranchAddress("GenJet02.EhadOverEem", GenJet02_EhadOverEem, &b_GenJet02_EhadOverEem);
   fChain->SetBranchAddress("GenJet02.NCharged", GenJet02_NCharged, &b_GenJet02_NCharged);
   fChain->SetBranchAddress("GenJet02.NNeutrals", GenJet02_NNeutrals, &b_GenJet02_NNeutrals);
   fChain->SetBranchAddress("GenJet02.NeutralEnergyFraction", GenJet02_NeutralEnergyFraction, &b_GenJet02_NeutralEnergyFraction);
   fChain->SetBranchAddress("GenJet02.ChargedEnergyFraction", GenJet02_ChargedEnergyFraction, &b_GenJet02_ChargedEnergyFraction);
   fChain->SetBranchAddress("GenJet02.Beta", GenJet02_Beta, &b_GenJet02_Beta);
   fChain->SetBranchAddress("GenJet02.BetaStar", GenJet02_BetaStar, &b_GenJet02_BetaStar);
   fChain->SetBranchAddress("GenJet02.MeanSqDeltaR", GenJet02_MeanSqDeltaR, &b_GenJet02_MeanSqDeltaR);
   fChain->SetBranchAddress("GenJet02.PTD", GenJet02_PTD, &b_GenJet02_PTD);
   fChain->SetBranchAddress("GenJet02.FracPt[5]", GenJet02_FracPt, &b_GenJet02_FracPt);
   fChain->SetBranchAddress("GenJet02.Tau[5]", GenJet02_Tau, &b_GenJet02_Tau);
   fChain->SetBranchAddress("GenJet02.SoftDroppedJet", GenJet02_SoftDroppedJet, &b_GenJet02_SoftDroppedJet);
   fChain->SetBranchAddress("GenJet02.SoftDroppedSubJet1", GenJet02_SoftDroppedSubJet1, &b_GenJet02_SoftDroppedSubJet1);
   fChain->SetBranchAddress("GenJet02.SoftDroppedSubJet2", GenJet02_SoftDroppedSubJet2, &b_GenJet02_SoftDroppedSubJet2);
   fChain->SetBranchAddress("GenJet02.TrimmedP4[5]", GenJet02_TrimmedP4, &b_GenJet02_TrimmedP4);
   fChain->SetBranchAddress("GenJet02.PrunedP4[5]", GenJet02_PrunedP4, &b_GenJet02_PrunedP4);
   fChain->SetBranchAddress("GenJet02.SoftDroppedP4[5]", GenJet02_SoftDroppedP4, &b_GenJet02_SoftDroppedP4);
   fChain->SetBranchAddress("GenJet02.NSubJetsTrimmed", GenJet02_NSubJetsTrimmed, &b_GenJet02_NSubJetsTrimmed);
   fChain->SetBranchAddress("GenJet02.NSubJetsPruned", GenJet02_NSubJetsPruned, &b_GenJet02_NSubJetsPruned);
   fChain->SetBranchAddress("GenJet02.NSubJetsSoftDropped", GenJet02_NSubJetsSoftDropped, &b_GenJet02_NSubJetsSoftDropped);
   fChain->SetBranchAddress("GenJet02.ExclYmerge23", GenJet02_ExclYmerge23, &b_GenJet02_ExclYmerge23);
   fChain->SetBranchAddress("GenJet02.ExclYmerge34", GenJet02_ExclYmerge34, &b_GenJet02_ExclYmerge34);
   fChain->SetBranchAddress("GenJet02.ExclYmerge45", GenJet02_ExclYmerge45, &b_GenJet02_ExclYmerge45);
   fChain->SetBranchAddress("GenJet02.ExclYmerge56", GenJet02_ExclYmerge56, &b_GenJet02_ExclYmerge56);
   fChain->SetBranchAddress("GenJet02.Constituents", GenJet02_Constituents, &b_GenJet02_Constituents);
   fChain->SetBranchAddress("GenJet02.Particles", GenJet02_Particles, &b_GenJet02_Particles);
   fChain->SetBranchAddress("GenJet02.Area", GenJet02_Area, &b_GenJet02_Area);
   fChain->SetBranchAddress("GenJet02_size", &GenJet02_size, &b_GenJet02_size);
   fChain->SetBranchAddress("GenJet04", &GenJet04_, &b_GenJet04_);
   fChain->SetBranchAddress("GenJet04.fUniqueID", GenJet04_fUniqueID, &b_GenJet04_fUniqueID);
   fChain->SetBranchAddress("GenJet04.fBits", GenJet04_fBits, &b_GenJet04_fBits);
   fChain->SetBranchAddress("GenJet04.PT", GenJet04_PT, &b_GenJet04_PT);
   fChain->SetBranchAddress("GenJet04.Eta", GenJet04_Eta, &b_GenJet04_Eta);
   fChain->SetBranchAddress("GenJet04.Phi", GenJet04_Phi, &b_GenJet04_Phi);
   fChain->SetBranchAddress("GenJet04.T", GenJet04_T, &b_GenJet04_T);
   fChain->SetBranchAddress("GenJet04.Mass", GenJet04_Mass, &b_GenJet04_Mass);
   fChain->SetBranchAddress("GenJet04.DeltaEta", GenJet04_DeltaEta, &b_GenJet04_DeltaEta);
   fChain->SetBranchAddress("GenJet04.DeltaPhi", GenJet04_DeltaPhi, &b_GenJet04_DeltaPhi);
   fChain->SetBranchAddress("GenJet04.Flavor", GenJet04_Flavor, &b_GenJet04_Flavor);
   fChain->SetBranchAddress("GenJet04.FlavorAlgo", GenJet04_FlavorAlgo, &b_GenJet04_FlavorAlgo);
   fChain->SetBranchAddress("GenJet04.FlavorPhys", GenJet04_FlavorPhys, &b_GenJet04_FlavorPhys);
   fChain->SetBranchAddress("GenJet04.BTag", GenJet04_BTag, &b_GenJet04_BTag);
   fChain->SetBranchAddress("GenJet04.BTagAlgo", GenJet04_BTagAlgo, &b_GenJet04_BTagAlgo);
   fChain->SetBranchAddress("GenJet04.BTagPhys", GenJet04_BTagPhys, &b_GenJet04_BTagPhys);
   fChain->SetBranchAddress("GenJet04.TauTag", GenJet04_TauTag, &b_GenJet04_TauTag);
   fChain->SetBranchAddress("GenJet04.TauWeight", GenJet04_TauWeight, &b_GenJet04_TauWeight);
   fChain->SetBranchAddress("GenJet04.Charge", GenJet04_Charge, &b_GenJet04_Charge);
   fChain->SetBranchAddress("GenJet04.EhadOverEem", GenJet04_EhadOverEem, &b_GenJet04_EhadOverEem);
   fChain->SetBranchAddress("GenJet04.NCharged", GenJet04_NCharged, &b_GenJet04_NCharged);
   fChain->SetBranchAddress("GenJet04.NNeutrals", GenJet04_NNeutrals, &b_GenJet04_NNeutrals);
   fChain->SetBranchAddress("GenJet04.NeutralEnergyFraction", GenJet04_NeutralEnergyFraction, &b_GenJet04_NeutralEnergyFraction);
   fChain->SetBranchAddress("GenJet04.ChargedEnergyFraction", GenJet04_ChargedEnergyFraction, &b_GenJet04_ChargedEnergyFraction);
   fChain->SetBranchAddress("GenJet04.Beta", GenJet04_Beta, &b_GenJet04_Beta);
   fChain->SetBranchAddress("GenJet04.BetaStar", GenJet04_BetaStar, &b_GenJet04_BetaStar);
   fChain->SetBranchAddress("GenJet04.MeanSqDeltaR", GenJet04_MeanSqDeltaR, &b_GenJet04_MeanSqDeltaR);
   fChain->SetBranchAddress("GenJet04.PTD", GenJet04_PTD, &b_GenJet04_PTD);
   fChain->SetBranchAddress("GenJet04.FracPt[5]", GenJet04_FracPt, &b_GenJet04_FracPt);
   fChain->SetBranchAddress("GenJet04.Tau[5]", GenJet04_Tau, &b_GenJet04_Tau);
   fChain->SetBranchAddress("GenJet04.SoftDroppedJet", GenJet04_SoftDroppedJet, &b_GenJet04_SoftDroppedJet);
   fChain->SetBranchAddress("GenJet04.SoftDroppedSubJet1", GenJet04_SoftDroppedSubJet1, &b_GenJet04_SoftDroppedSubJet1);
   fChain->SetBranchAddress("GenJet04.SoftDroppedSubJet2", GenJet04_SoftDroppedSubJet2, &b_GenJet04_SoftDroppedSubJet2);
   fChain->SetBranchAddress("GenJet04.TrimmedP4[5]", GenJet04_TrimmedP4, &b_GenJet04_TrimmedP4);
   fChain->SetBranchAddress("GenJet04.PrunedP4[5]", GenJet04_PrunedP4, &b_GenJet04_PrunedP4);
   fChain->SetBranchAddress("GenJet04.SoftDroppedP4[5]", GenJet04_SoftDroppedP4, &b_GenJet04_SoftDroppedP4);
   fChain->SetBranchAddress("GenJet04.NSubJetsTrimmed", GenJet04_NSubJetsTrimmed, &b_GenJet04_NSubJetsTrimmed);
   fChain->SetBranchAddress("GenJet04.NSubJetsPruned", GenJet04_NSubJetsPruned, &b_GenJet04_NSubJetsPruned);
   fChain->SetBranchAddress("GenJet04.NSubJetsSoftDropped", GenJet04_NSubJetsSoftDropped, &b_GenJet04_NSubJetsSoftDropped);
   fChain->SetBranchAddress("GenJet04.ExclYmerge23", GenJet04_ExclYmerge23, &b_GenJet04_ExclYmerge23);
   fChain->SetBranchAddress("GenJet04.ExclYmerge34", GenJet04_ExclYmerge34, &b_GenJet04_ExclYmerge34);
   fChain->SetBranchAddress("GenJet04.ExclYmerge45", GenJet04_ExclYmerge45, &b_GenJet04_ExclYmerge45);
   fChain->SetBranchAddress("GenJet04.ExclYmerge56", GenJet04_ExclYmerge56, &b_GenJet04_ExclYmerge56);
   fChain->SetBranchAddress("GenJet04.Constituents", GenJet04_Constituents, &b_GenJet04_Constituents);
   fChain->SetBranchAddress("GenJet04.Particles", GenJet04_Particles, &b_GenJet04_Particles);
   fChain->SetBranchAddress("GenJet04.Area", GenJet04_Area, &b_GenJet04_Area);
   fChain->SetBranchAddress("GenJet04_size", &GenJet04_size, &b_GenJet04_size);
   fChain->SetBranchAddress("GenJet08", &GenJet08_, &b_GenJet08_);
   fChain->SetBranchAddress("GenJet08.fUniqueID", GenJet08_fUniqueID, &b_GenJet08_fUniqueID);
   fChain->SetBranchAddress("GenJet08.fBits", GenJet08_fBits, &b_GenJet08_fBits);
   fChain->SetBranchAddress("GenJet08.PT", GenJet08_PT, &b_GenJet08_PT);
   fChain->SetBranchAddress("GenJet08.Eta", GenJet08_Eta, &b_GenJet08_Eta);
   fChain->SetBranchAddress("GenJet08.Phi", GenJet08_Phi, &b_GenJet08_Phi);
   fChain->SetBranchAddress("GenJet08.T", GenJet08_T, &b_GenJet08_T);
   fChain->SetBranchAddress("GenJet08.Mass", GenJet08_Mass, &b_GenJet08_Mass);
   fChain->SetBranchAddress("GenJet08.DeltaEta", GenJet08_DeltaEta, &b_GenJet08_DeltaEta);
   fChain->SetBranchAddress("GenJet08.DeltaPhi", GenJet08_DeltaPhi, &b_GenJet08_DeltaPhi);
   fChain->SetBranchAddress("GenJet08.Flavor", GenJet08_Flavor, &b_GenJet08_Flavor);
   fChain->SetBranchAddress("GenJet08.FlavorAlgo", GenJet08_FlavorAlgo, &b_GenJet08_FlavorAlgo);
   fChain->SetBranchAddress("GenJet08.FlavorPhys", GenJet08_FlavorPhys, &b_GenJet08_FlavorPhys);
   fChain->SetBranchAddress("GenJet08.BTag", GenJet08_BTag, &b_GenJet08_BTag);
   fChain->SetBranchAddress("GenJet08.BTagAlgo", GenJet08_BTagAlgo, &b_GenJet08_BTagAlgo);
   fChain->SetBranchAddress("GenJet08.BTagPhys", GenJet08_BTagPhys, &b_GenJet08_BTagPhys);
   fChain->SetBranchAddress("GenJet08.TauTag", GenJet08_TauTag, &b_GenJet08_TauTag);
   fChain->SetBranchAddress("GenJet08.TauWeight", GenJet08_TauWeight, &b_GenJet08_TauWeight);
   fChain->SetBranchAddress("GenJet08.Charge", GenJet08_Charge, &b_GenJet08_Charge);
   fChain->SetBranchAddress("GenJet08.EhadOverEem", GenJet08_EhadOverEem, &b_GenJet08_EhadOverEem);
   fChain->SetBranchAddress("GenJet08.NCharged", GenJet08_NCharged, &b_GenJet08_NCharged);
   fChain->SetBranchAddress("GenJet08.NNeutrals", GenJet08_NNeutrals, &b_GenJet08_NNeutrals);
   fChain->SetBranchAddress("GenJet08.NeutralEnergyFraction", GenJet08_NeutralEnergyFraction, &b_GenJet08_NeutralEnergyFraction);
   fChain->SetBranchAddress("GenJet08.ChargedEnergyFraction", GenJet08_ChargedEnergyFraction, &b_GenJet08_ChargedEnergyFraction);
   fChain->SetBranchAddress("GenJet08.Beta", GenJet08_Beta, &b_GenJet08_Beta);
   fChain->SetBranchAddress("GenJet08.BetaStar", GenJet08_BetaStar, &b_GenJet08_BetaStar);
   fChain->SetBranchAddress("GenJet08.MeanSqDeltaR", GenJet08_MeanSqDeltaR, &b_GenJet08_MeanSqDeltaR);
   fChain->SetBranchAddress("GenJet08.PTD", GenJet08_PTD, &b_GenJet08_PTD);
   fChain->SetBranchAddress("GenJet08.FracPt[5]", GenJet08_FracPt, &b_GenJet08_FracPt);
   fChain->SetBranchAddress("GenJet08.Tau[5]", GenJet08_Tau, &b_GenJet08_Tau);
   fChain->SetBranchAddress("GenJet08.SoftDroppedJet", GenJet08_SoftDroppedJet, &b_GenJet08_SoftDroppedJet);
   fChain->SetBranchAddress("GenJet08.SoftDroppedSubJet1", GenJet08_SoftDroppedSubJet1, &b_GenJet08_SoftDroppedSubJet1);
   fChain->SetBranchAddress("GenJet08.SoftDroppedSubJet2", GenJet08_SoftDroppedSubJet2, &b_GenJet08_SoftDroppedSubJet2);
   fChain->SetBranchAddress("GenJet08.TrimmedP4[5]", GenJet08_TrimmedP4, &b_GenJet08_TrimmedP4);
   fChain->SetBranchAddress("GenJet08.PrunedP4[5]", GenJet08_PrunedP4, &b_GenJet08_PrunedP4);
   fChain->SetBranchAddress("GenJet08.SoftDroppedP4[5]", GenJet08_SoftDroppedP4, &b_GenJet08_SoftDroppedP4);
   fChain->SetBranchAddress("GenJet08.NSubJetsTrimmed", GenJet08_NSubJetsTrimmed, &b_GenJet08_NSubJetsTrimmed);
   fChain->SetBranchAddress("GenJet08.NSubJetsPruned", GenJet08_NSubJetsPruned, &b_GenJet08_NSubJetsPruned);
   fChain->SetBranchAddress("GenJet08.NSubJetsSoftDropped", GenJet08_NSubJetsSoftDropped, &b_GenJet08_NSubJetsSoftDropped);
   fChain->SetBranchAddress("GenJet08.ExclYmerge23", GenJet08_ExclYmerge23, &b_GenJet08_ExclYmerge23);
   fChain->SetBranchAddress("GenJet08.ExclYmerge34", GenJet08_ExclYmerge34, &b_GenJet08_ExclYmerge34);
   fChain->SetBranchAddress("GenJet08.ExclYmerge45", GenJet08_ExclYmerge45, &b_GenJet08_ExclYmerge45);
   fChain->SetBranchAddress("GenJet08.ExclYmerge56", GenJet08_ExclYmerge56, &b_GenJet08_ExclYmerge56);
   fChain->SetBranchAddress("GenJet08.Constituents", GenJet08_Constituents, &b_GenJet08_Constituents);
   fChain->SetBranchAddress("GenJet08.Particles", GenJet08_Particles, &b_GenJet08_Particles);
   fChain->SetBranchAddress("GenJet08.Area", GenJet08_Area, &b_GenJet08_Area);
   fChain->SetBranchAddress("GenJet08_size", &GenJet08_size, &b_GenJet08_size);
   fChain->SetBranchAddress("GenJet15", &GenJet15_, &b_GenJet15_);
   fChain->SetBranchAddress("GenJet15.fUniqueID", GenJet15_fUniqueID, &b_GenJet15_fUniqueID);
   fChain->SetBranchAddress("GenJet15.fBits", GenJet15_fBits, &b_GenJet15_fBits);
   fChain->SetBranchAddress("GenJet15.PT", GenJet15_PT, &b_GenJet15_PT);
   fChain->SetBranchAddress("GenJet15.Eta", GenJet15_Eta, &b_GenJet15_Eta);
   fChain->SetBranchAddress("GenJet15.Phi", GenJet15_Phi, &b_GenJet15_Phi);
   fChain->SetBranchAddress("GenJet15.T", GenJet15_T, &b_GenJet15_T);
   fChain->SetBranchAddress("GenJet15.Mass", GenJet15_Mass, &b_GenJet15_Mass);
   fChain->SetBranchAddress("GenJet15.DeltaEta", GenJet15_DeltaEta, &b_GenJet15_DeltaEta);
   fChain->SetBranchAddress("GenJet15.DeltaPhi", GenJet15_DeltaPhi, &b_GenJet15_DeltaPhi);
   fChain->SetBranchAddress("GenJet15.Flavor", GenJet15_Flavor, &b_GenJet15_Flavor);
   fChain->SetBranchAddress("GenJet15.FlavorAlgo", GenJet15_FlavorAlgo, &b_GenJet15_FlavorAlgo);
   fChain->SetBranchAddress("GenJet15.FlavorPhys", GenJet15_FlavorPhys, &b_GenJet15_FlavorPhys);
   fChain->SetBranchAddress("GenJet15.BTag", GenJet15_BTag, &b_GenJet15_BTag);
   fChain->SetBranchAddress("GenJet15.BTagAlgo", GenJet15_BTagAlgo, &b_GenJet15_BTagAlgo);
   fChain->SetBranchAddress("GenJet15.BTagPhys", GenJet15_BTagPhys, &b_GenJet15_BTagPhys);
   fChain->SetBranchAddress("GenJet15.TauTag", GenJet15_TauTag, &b_GenJet15_TauTag);
   fChain->SetBranchAddress("GenJet15.TauWeight", GenJet15_TauWeight, &b_GenJet15_TauWeight);
   fChain->SetBranchAddress("GenJet15.Charge", GenJet15_Charge, &b_GenJet15_Charge);
   fChain->SetBranchAddress("GenJet15.EhadOverEem", GenJet15_EhadOverEem, &b_GenJet15_EhadOverEem);
   fChain->SetBranchAddress("GenJet15.NCharged", GenJet15_NCharged, &b_GenJet15_NCharged);
   fChain->SetBranchAddress("GenJet15.NNeutrals", GenJet15_NNeutrals, &b_GenJet15_NNeutrals);
   fChain->SetBranchAddress("GenJet15.NeutralEnergyFraction", GenJet15_NeutralEnergyFraction, &b_GenJet15_NeutralEnergyFraction);
   fChain->SetBranchAddress("GenJet15.ChargedEnergyFraction", GenJet15_ChargedEnergyFraction, &b_GenJet15_ChargedEnergyFraction);
   fChain->SetBranchAddress("GenJet15.Beta", GenJet15_Beta, &b_GenJet15_Beta);
   fChain->SetBranchAddress("GenJet15.BetaStar", GenJet15_BetaStar, &b_GenJet15_BetaStar);
   fChain->SetBranchAddress("GenJet15.MeanSqDeltaR", GenJet15_MeanSqDeltaR, &b_GenJet15_MeanSqDeltaR);
   fChain->SetBranchAddress("GenJet15.PTD", GenJet15_PTD, &b_GenJet15_PTD);
   fChain->SetBranchAddress("GenJet15.FracPt[5]", GenJet15_FracPt, &b_GenJet15_FracPt);
   fChain->SetBranchAddress("GenJet15.Tau[5]", GenJet15_Tau, &b_GenJet15_Tau);
   fChain->SetBranchAddress("GenJet15.SoftDroppedJet", GenJet15_SoftDroppedJet, &b_GenJet15_SoftDroppedJet);
   fChain->SetBranchAddress("GenJet15.SoftDroppedSubJet1", GenJet15_SoftDroppedSubJet1, &b_GenJet15_SoftDroppedSubJet1);
   fChain->SetBranchAddress("GenJet15.SoftDroppedSubJet2", GenJet15_SoftDroppedSubJet2, &b_GenJet15_SoftDroppedSubJet2);
   fChain->SetBranchAddress("GenJet15.TrimmedP4[5]", GenJet15_TrimmedP4, &b_GenJet15_TrimmedP4);
   fChain->SetBranchAddress("GenJet15.PrunedP4[5]", GenJet15_PrunedP4, &b_GenJet15_PrunedP4);
   fChain->SetBranchAddress("GenJet15.SoftDroppedP4[5]", GenJet15_SoftDroppedP4, &b_GenJet15_SoftDroppedP4);
   fChain->SetBranchAddress("GenJet15.NSubJetsTrimmed", GenJet15_NSubJetsTrimmed, &b_GenJet15_NSubJetsTrimmed);
   fChain->SetBranchAddress("GenJet15.NSubJetsPruned", GenJet15_NSubJetsPruned, &b_GenJet15_NSubJetsPruned);
   fChain->SetBranchAddress("GenJet15.NSubJetsSoftDropped", GenJet15_NSubJetsSoftDropped, &b_GenJet15_NSubJetsSoftDropped);
   fChain->SetBranchAddress("GenJet15.ExclYmerge23", GenJet15_ExclYmerge23, &b_GenJet15_ExclYmerge23);
   fChain->SetBranchAddress("GenJet15.ExclYmerge34", GenJet15_ExclYmerge34, &b_GenJet15_ExclYmerge34);
   fChain->SetBranchAddress("GenJet15.ExclYmerge45", GenJet15_ExclYmerge45, &b_GenJet15_ExclYmerge45);
   fChain->SetBranchAddress("GenJet15.ExclYmerge56", GenJet15_ExclYmerge56, &b_GenJet15_ExclYmerge56);
   fChain->SetBranchAddress("GenJet15.Constituents", GenJet15_Constituents, &b_GenJet15_Constituents);
   fChain->SetBranchAddress("GenJet15.Particles", GenJet15_Particles, &b_GenJet15_Particles);
   fChain->SetBranchAddress("GenJet15.Area", GenJet15_Area, &b_GenJet15_Area);
   fChain->SetBranchAddress("GenJet15_size", &GenJet15_size, &b_GenJet15_size);
   fChain->SetBranchAddress("ParticleFlowJet02", &ParticleFlowJet02_, &b_ParticleFlowJet02_);
   fChain->SetBranchAddress("ParticleFlowJet02.fUniqueID", ParticleFlowJet02_fUniqueID, &b_ParticleFlowJet02_fUniqueID);
   fChain->SetBranchAddress("ParticleFlowJet02.fBits", ParticleFlowJet02_fBits, &b_ParticleFlowJet02_fBits);
   fChain->SetBranchAddress("ParticleFlowJet02.PT", ParticleFlowJet02_PT, &b_ParticleFlowJet02_PT);
   fChain->SetBranchAddress("ParticleFlowJet02.Eta", ParticleFlowJet02_Eta, &b_ParticleFlowJet02_Eta);
   fChain->SetBranchAddress("ParticleFlowJet02.Phi", ParticleFlowJet02_Phi, &b_ParticleFlowJet02_Phi);
   fChain->SetBranchAddress("ParticleFlowJet02.T", ParticleFlowJet02_T, &b_ParticleFlowJet02_T);
   fChain->SetBranchAddress("ParticleFlowJet02.Mass", ParticleFlowJet02_Mass, &b_ParticleFlowJet02_Mass);
   fChain->SetBranchAddress("ParticleFlowJet02.DeltaEta", ParticleFlowJet02_DeltaEta, &b_ParticleFlowJet02_DeltaEta);
   fChain->SetBranchAddress("ParticleFlowJet02.DeltaPhi", ParticleFlowJet02_DeltaPhi, &b_ParticleFlowJet02_DeltaPhi);
   fChain->SetBranchAddress("ParticleFlowJet02.Flavor", ParticleFlowJet02_Flavor, &b_ParticleFlowJet02_Flavor);
   fChain->SetBranchAddress("ParticleFlowJet02.FlavorAlgo", ParticleFlowJet02_FlavorAlgo, &b_ParticleFlowJet02_FlavorAlgo);
   fChain->SetBranchAddress("ParticleFlowJet02.FlavorPhys", ParticleFlowJet02_FlavorPhys, &b_ParticleFlowJet02_FlavorPhys);
   fChain->SetBranchAddress("ParticleFlowJet02.BTag", ParticleFlowJet02_BTag, &b_ParticleFlowJet02_BTag);
   fChain->SetBranchAddress("ParticleFlowJet02.BTagAlgo", ParticleFlowJet02_BTagAlgo, &b_ParticleFlowJet02_BTagAlgo);
   fChain->SetBranchAddress("ParticleFlowJet02.BTagPhys", ParticleFlowJet02_BTagPhys, &b_ParticleFlowJet02_BTagPhys);
   fChain->SetBranchAddress("ParticleFlowJet02.TauTag", ParticleFlowJet02_TauTag, &b_ParticleFlowJet02_TauTag);
   fChain->SetBranchAddress("ParticleFlowJet02.TauWeight", ParticleFlowJet02_TauWeight, &b_ParticleFlowJet02_TauWeight);
   fChain->SetBranchAddress("ParticleFlowJet02.Charge", ParticleFlowJet02_Charge, &b_ParticleFlowJet02_Charge);
   fChain->SetBranchAddress("ParticleFlowJet02.EhadOverEem", ParticleFlowJet02_EhadOverEem, &b_ParticleFlowJet02_EhadOverEem);
   fChain->SetBranchAddress("ParticleFlowJet02.NCharged", ParticleFlowJet02_NCharged, &b_ParticleFlowJet02_NCharged);
   fChain->SetBranchAddress("ParticleFlowJet02.NNeutrals", ParticleFlowJet02_NNeutrals, &b_ParticleFlowJet02_NNeutrals);
   fChain->SetBranchAddress("ParticleFlowJet02.NeutralEnergyFraction", ParticleFlowJet02_NeutralEnergyFraction, &b_ParticleFlowJet02_NeutralEnergyFraction);
   fChain->SetBranchAddress("ParticleFlowJet02.ChargedEnergyFraction", ParticleFlowJet02_ChargedEnergyFraction, &b_ParticleFlowJet02_ChargedEnergyFraction);
   fChain->SetBranchAddress("ParticleFlowJet02.Beta", ParticleFlowJet02_Beta, &b_ParticleFlowJet02_Beta);
   fChain->SetBranchAddress("ParticleFlowJet02.BetaStar", ParticleFlowJet02_BetaStar, &b_ParticleFlowJet02_BetaStar);
   fChain->SetBranchAddress("ParticleFlowJet02.MeanSqDeltaR", ParticleFlowJet02_MeanSqDeltaR, &b_ParticleFlowJet02_MeanSqDeltaR);
   fChain->SetBranchAddress("ParticleFlowJet02.PTD", ParticleFlowJet02_PTD, &b_ParticleFlowJet02_PTD);
   fChain->SetBranchAddress("ParticleFlowJet02.FracPt[5]", ParticleFlowJet02_FracPt, &b_ParticleFlowJet02_FracPt);
   fChain->SetBranchAddress("ParticleFlowJet02.Tau[5]", ParticleFlowJet02_Tau, &b_ParticleFlowJet02_Tau);
   fChain->SetBranchAddress("ParticleFlowJet02.SoftDroppedJet", ParticleFlowJet02_SoftDroppedJet, &b_ParticleFlowJet02_SoftDroppedJet);
   fChain->SetBranchAddress("ParticleFlowJet02.SoftDroppedSubJet1", ParticleFlowJet02_SoftDroppedSubJet1, &b_ParticleFlowJet02_SoftDroppedSubJet1);
   fChain->SetBranchAddress("ParticleFlowJet02.SoftDroppedSubJet2", ParticleFlowJet02_SoftDroppedSubJet2, &b_ParticleFlowJet02_SoftDroppedSubJet2);
   fChain->SetBranchAddress("ParticleFlowJet02.TrimmedP4[5]", ParticleFlowJet02_TrimmedP4, &b_ParticleFlowJet02_TrimmedP4);
   fChain->SetBranchAddress("ParticleFlowJet02.PrunedP4[5]", ParticleFlowJet02_PrunedP4, &b_ParticleFlowJet02_PrunedP4);
   fChain->SetBranchAddress("ParticleFlowJet02.SoftDroppedP4[5]", ParticleFlowJet02_SoftDroppedP4, &b_ParticleFlowJet02_SoftDroppedP4);
   fChain->SetBranchAddress("ParticleFlowJet02.NSubJetsTrimmed", ParticleFlowJet02_NSubJetsTrimmed, &b_ParticleFlowJet02_NSubJetsTrimmed);
   fChain->SetBranchAddress("ParticleFlowJet02.NSubJetsPruned", ParticleFlowJet02_NSubJetsPruned, &b_ParticleFlowJet02_NSubJetsPruned);
   fChain->SetBranchAddress("ParticleFlowJet02.NSubJetsSoftDropped", ParticleFlowJet02_NSubJetsSoftDropped, &b_ParticleFlowJet02_NSubJetsSoftDropped);
   fChain->SetBranchAddress("ParticleFlowJet02.ExclYmerge23", ParticleFlowJet02_ExclYmerge23, &b_ParticleFlowJet02_ExclYmerge23);
   fChain->SetBranchAddress("ParticleFlowJet02.ExclYmerge34", ParticleFlowJet02_ExclYmerge34, &b_ParticleFlowJet02_ExclYmerge34);
   fChain->SetBranchAddress("ParticleFlowJet02.ExclYmerge45", ParticleFlowJet02_ExclYmerge45, &b_ParticleFlowJet02_ExclYmerge45);
   fChain->SetBranchAddress("ParticleFlowJet02.ExclYmerge56", ParticleFlowJet02_ExclYmerge56, &b_ParticleFlowJet02_ExclYmerge56);
   fChain->SetBranchAddress("ParticleFlowJet02.Constituents", ParticleFlowJet02_Constituents, &b_ParticleFlowJet02_Constituents);
   fChain->SetBranchAddress("ParticleFlowJet02.Particles", ParticleFlowJet02_Particles, &b_ParticleFlowJet02_Particles);
   fChain->SetBranchAddress("ParticleFlowJet02.Area", ParticleFlowJet02_Area, &b_ParticleFlowJet02_Area);
   fChain->SetBranchAddress("ParticleFlowJet02_size", &ParticleFlowJet02_size, &b_ParticleFlowJet02_size);
   fChain->SetBranchAddress("ParticleFlowJet04", &ParticleFlowJet04_, &b_ParticleFlowJet04_);
   fChain->SetBranchAddress("ParticleFlowJet04.fUniqueID", ParticleFlowJet04_fUniqueID, &b_ParticleFlowJet04_fUniqueID);
   fChain->SetBranchAddress("ParticleFlowJet04.fBits", ParticleFlowJet04_fBits, &b_ParticleFlowJet04_fBits);
   fChain->SetBranchAddress("ParticleFlowJet04.PT", ParticleFlowJet04_PT, &b_ParticleFlowJet04_PT);
   fChain->SetBranchAddress("ParticleFlowJet04.Eta", ParticleFlowJet04_Eta, &b_ParticleFlowJet04_Eta);
   fChain->SetBranchAddress("ParticleFlowJet04.Phi", ParticleFlowJet04_Phi, &b_ParticleFlowJet04_Phi);
   fChain->SetBranchAddress("ParticleFlowJet04.T", ParticleFlowJet04_T, &b_ParticleFlowJet04_T);
   fChain->SetBranchAddress("ParticleFlowJet04.Mass", ParticleFlowJet04_Mass, &b_ParticleFlowJet04_Mass);
   fChain->SetBranchAddress("ParticleFlowJet04.DeltaEta", ParticleFlowJet04_DeltaEta, &b_ParticleFlowJet04_DeltaEta);
   fChain->SetBranchAddress("ParticleFlowJet04.DeltaPhi", ParticleFlowJet04_DeltaPhi, &b_ParticleFlowJet04_DeltaPhi);
   fChain->SetBranchAddress("ParticleFlowJet04.Flavor", ParticleFlowJet04_Flavor, &b_ParticleFlowJet04_Flavor);
   fChain->SetBranchAddress("ParticleFlowJet04.FlavorAlgo", ParticleFlowJet04_FlavorAlgo, &b_ParticleFlowJet04_FlavorAlgo);
   fChain->SetBranchAddress("ParticleFlowJet04.FlavorPhys", ParticleFlowJet04_FlavorPhys, &b_ParticleFlowJet04_FlavorPhys);
   fChain->SetBranchAddress("ParticleFlowJet04.BTag", ParticleFlowJet04_BTag, &b_ParticleFlowJet04_BTag);
   fChain->SetBranchAddress("ParticleFlowJet04.BTagAlgo", ParticleFlowJet04_BTagAlgo, &b_ParticleFlowJet04_BTagAlgo);
   fChain->SetBranchAddress("ParticleFlowJet04.BTagPhys", ParticleFlowJet04_BTagPhys, &b_ParticleFlowJet04_BTagPhys);
   fChain->SetBranchAddress("ParticleFlowJet04.TauTag", ParticleFlowJet04_TauTag, &b_ParticleFlowJet04_TauTag);
   fChain->SetBranchAddress("ParticleFlowJet04.TauWeight", ParticleFlowJet04_TauWeight, &b_ParticleFlowJet04_TauWeight);
   fChain->SetBranchAddress("ParticleFlowJet04.Charge", ParticleFlowJet04_Charge, &b_ParticleFlowJet04_Charge);
   fChain->SetBranchAddress("ParticleFlowJet04.EhadOverEem", ParticleFlowJet04_EhadOverEem, &b_ParticleFlowJet04_EhadOverEem);
   fChain->SetBranchAddress("ParticleFlowJet04.NCharged", ParticleFlowJet04_NCharged, &b_ParticleFlowJet04_NCharged);
   fChain->SetBranchAddress("ParticleFlowJet04.NNeutrals", ParticleFlowJet04_NNeutrals, &b_ParticleFlowJet04_NNeutrals);
   fChain->SetBranchAddress("ParticleFlowJet04.NeutralEnergyFraction", ParticleFlowJet04_NeutralEnergyFraction, &b_ParticleFlowJet04_NeutralEnergyFraction);
   fChain->SetBranchAddress("ParticleFlowJet04.ChargedEnergyFraction", ParticleFlowJet04_ChargedEnergyFraction, &b_ParticleFlowJet04_ChargedEnergyFraction);
   fChain->SetBranchAddress("ParticleFlowJet04.Beta", ParticleFlowJet04_Beta, &b_ParticleFlowJet04_Beta);
   fChain->SetBranchAddress("ParticleFlowJet04.BetaStar", ParticleFlowJet04_BetaStar, &b_ParticleFlowJet04_BetaStar);
   fChain->SetBranchAddress("ParticleFlowJet04.MeanSqDeltaR", ParticleFlowJet04_MeanSqDeltaR, &b_ParticleFlowJet04_MeanSqDeltaR);
   fChain->SetBranchAddress("ParticleFlowJet04.PTD", ParticleFlowJet04_PTD, &b_ParticleFlowJet04_PTD);
   fChain->SetBranchAddress("ParticleFlowJet04.FracPt[5]", ParticleFlowJet04_FracPt, &b_ParticleFlowJet04_FracPt);
   fChain->SetBranchAddress("ParticleFlowJet04.Tau[5]", ParticleFlowJet04_Tau, &b_ParticleFlowJet04_Tau);
   fChain->SetBranchAddress("ParticleFlowJet04.SoftDroppedJet", ParticleFlowJet04_SoftDroppedJet, &b_ParticleFlowJet04_SoftDroppedJet);
   fChain->SetBranchAddress("ParticleFlowJet04.SoftDroppedSubJet1", ParticleFlowJet04_SoftDroppedSubJet1, &b_ParticleFlowJet04_SoftDroppedSubJet1);
   fChain->SetBranchAddress("ParticleFlowJet04.SoftDroppedSubJet2", ParticleFlowJet04_SoftDroppedSubJet2, &b_ParticleFlowJet04_SoftDroppedSubJet2);
   fChain->SetBranchAddress("ParticleFlowJet04.TrimmedP4[5]", ParticleFlowJet04_TrimmedP4, &b_ParticleFlowJet04_TrimmedP4);
   fChain->SetBranchAddress("ParticleFlowJet04.PrunedP4[5]", ParticleFlowJet04_PrunedP4, &b_ParticleFlowJet04_PrunedP4);
   fChain->SetBranchAddress("ParticleFlowJet04.SoftDroppedP4[5]", ParticleFlowJet04_SoftDroppedP4, &b_ParticleFlowJet04_SoftDroppedP4);
   fChain->SetBranchAddress("ParticleFlowJet04.NSubJetsTrimmed", ParticleFlowJet04_NSubJetsTrimmed, &b_ParticleFlowJet04_NSubJetsTrimmed);
   fChain->SetBranchAddress("ParticleFlowJet04.NSubJetsPruned", ParticleFlowJet04_NSubJetsPruned, &b_ParticleFlowJet04_NSubJetsPruned);
   fChain->SetBranchAddress("ParticleFlowJet04.NSubJetsSoftDropped", ParticleFlowJet04_NSubJetsSoftDropped, &b_ParticleFlowJet04_NSubJetsSoftDropped);
   fChain->SetBranchAddress("ParticleFlowJet04.ExclYmerge23", ParticleFlowJet04_ExclYmerge23, &b_ParticleFlowJet04_ExclYmerge23);
   fChain->SetBranchAddress("ParticleFlowJet04.ExclYmerge34", ParticleFlowJet04_ExclYmerge34, &b_ParticleFlowJet04_ExclYmerge34);
   fChain->SetBranchAddress("ParticleFlowJet04.ExclYmerge45", ParticleFlowJet04_ExclYmerge45, &b_ParticleFlowJet04_ExclYmerge45);
   fChain->SetBranchAddress("ParticleFlowJet04.ExclYmerge56", ParticleFlowJet04_ExclYmerge56, &b_ParticleFlowJet04_ExclYmerge56);
   fChain->SetBranchAddress("ParticleFlowJet04.Constituents", ParticleFlowJet04_Constituents, &b_ParticleFlowJet04_Constituents);
   fChain->SetBranchAddress("ParticleFlowJet04.Particles", ParticleFlowJet04_Particles, &b_ParticleFlowJet04_Particles);
   fChain->SetBranchAddress("ParticleFlowJet04.Area", ParticleFlowJet04_Area, &b_ParticleFlowJet04_Area);
   fChain->SetBranchAddress("ParticleFlowJet04_size", &ParticleFlowJet04_size, &b_ParticleFlowJet04_size);
   fChain->SetBranchAddress("ParticleFlowJet08", &ParticleFlowJet08_, &b_ParticleFlowJet08_);
   fChain->SetBranchAddress("ParticleFlowJet08.fUniqueID", ParticleFlowJet08_fUniqueID, &b_ParticleFlowJet08_fUniqueID);
   fChain->SetBranchAddress("ParticleFlowJet08.fBits", ParticleFlowJet08_fBits, &b_ParticleFlowJet08_fBits);
   fChain->SetBranchAddress("ParticleFlowJet08.PT", ParticleFlowJet08_PT, &b_ParticleFlowJet08_PT);
   fChain->SetBranchAddress("ParticleFlowJet08.Eta", ParticleFlowJet08_Eta, &b_ParticleFlowJet08_Eta);
   fChain->SetBranchAddress("ParticleFlowJet08.Phi", ParticleFlowJet08_Phi, &b_ParticleFlowJet08_Phi);
   fChain->SetBranchAddress("ParticleFlowJet08.T", ParticleFlowJet08_T, &b_ParticleFlowJet08_T);
   fChain->SetBranchAddress("ParticleFlowJet08.Mass", ParticleFlowJet08_Mass, &b_ParticleFlowJet08_Mass);
   fChain->SetBranchAddress("ParticleFlowJet08.DeltaEta", ParticleFlowJet08_DeltaEta, &b_ParticleFlowJet08_DeltaEta);
   fChain->SetBranchAddress("ParticleFlowJet08.DeltaPhi", ParticleFlowJet08_DeltaPhi, &b_ParticleFlowJet08_DeltaPhi);
   fChain->SetBranchAddress("ParticleFlowJet08.Flavor", ParticleFlowJet08_Flavor, &b_ParticleFlowJet08_Flavor);
   fChain->SetBranchAddress("ParticleFlowJet08.FlavorAlgo", ParticleFlowJet08_FlavorAlgo, &b_ParticleFlowJet08_FlavorAlgo);
   fChain->SetBranchAddress("ParticleFlowJet08.FlavorPhys", ParticleFlowJet08_FlavorPhys, &b_ParticleFlowJet08_FlavorPhys);
   fChain->SetBranchAddress("ParticleFlowJet08.BTag", ParticleFlowJet08_BTag, &b_ParticleFlowJet08_BTag);
   fChain->SetBranchAddress("ParticleFlowJet08.BTagAlgo", ParticleFlowJet08_BTagAlgo, &b_ParticleFlowJet08_BTagAlgo);
   fChain->SetBranchAddress("ParticleFlowJet08.BTagPhys", ParticleFlowJet08_BTagPhys, &b_ParticleFlowJet08_BTagPhys);
   fChain->SetBranchAddress("ParticleFlowJet08.TauTag", ParticleFlowJet08_TauTag, &b_ParticleFlowJet08_TauTag);
   fChain->SetBranchAddress("ParticleFlowJet08.TauWeight", ParticleFlowJet08_TauWeight, &b_ParticleFlowJet08_TauWeight);
   fChain->SetBranchAddress("ParticleFlowJet08.Charge", ParticleFlowJet08_Charge, &b_ParticleFlowJet08_Charge);
   fChain->SetBranchAddress("ParticleFlowJet08.EhadOverEem", ParticleFlowJet08_EhadOverEem, &b_ParticleFlowJet08_EhadOverEem);
   fChain->SetBranchAddress("ParticleFlowJet08.NCharged", ParticleFlowJet08_NCharged, &b_ParticleFlowJet08_NCharged);
   fChain->SetBranchAddress("ParticleFlowJet08.NNeutrals", ParticleFlowJet08_NNeutrals, &b_ParticleFlowJet08_NNeutrals);
   fChain->SetBranchAddress("ParticleFlowJet08.NeutralEnergyFraction", ParticleFlowJet08_NeutralEnergyFraction, &b_ParticleFlowJet08_NeutralEnergyFraction);
   fChain->SetBranchAddress("ParticleFlowJet08.ChargedEnergyFraction", ParticleFlowJet08_ChargedEnergyFraction, &b_ParticleFlowJet08_ChargedEnergyFraction);
   fChain->SetBranchAddress("ParticleFlowJet08.Beta", ParticleFlowJet08_Beta, &b_ParticleFlowJet08_Beta);
   fChain->SetBranchAddress("ParticleFlowJet08.BetaStar", ParticleFlowJet08_BetaStar, &b_ParticleFlowJet08_BetaStar);
   fChain->SetBranchAddress("ParticleFlowJet08.MeanSqDeltaR", ParticleFlowJet08_MeanSqDeltaR, &b_ParticleFlowJet08_MeanSqDeltaR);
   fChain->SetBranchAddress("ParticleFlowJet08.PTD", ParticleFlowJet08_PTD, &b_ParticleFlowJet08_PTD);
   fChain->SetBranchAddress("ParticleFlowJet08.FracPt[5]", ParticleFlowJet08_FracPt, &b_ParticleFlowJet08_FracPt);
   fChain->SetBranchAddress("ParticleFlowJet08.Tau[5]", ParticleFlowJet08_Tau, &b_ParticleFlowJet08_Tau);
   fChain->SetBranchAddress("ParticleFlowJet08.SoftDroppedJet", ParticleFlowJet08_SoftDroppedJet, &b_ParticleFlowJet08_SoftDroppedJet);
   fChain->SetBranchAddress("ParticleFlowJet08.SoftDroppedSubJet1", ParticleFlowJet08_SoftDroppedSubJet1, &b_ParticleFlowJet08_SoftDroppedSubJet1);
   fChain->SetBranchAddress("ParticleFlowJet08.SoftDroppedSubJet2", ParticleFlowJet08_SoftDroppedSubJet2, &b_ParticleFlowJet08_SoftDroppedSubJet2);
   fChain->SetBranchAddress("ParticleFlowJet08.TrimmedP4[5]", ParticleFlowJet08_TrimmedP4, &b_ParticleFlowJet08_TrimmedP4);
   fChain->SetBranchAddress("ParticleFlowJet08.PrunedP4[5]", ParticleFlowJet08_PrunedP4, &b_ParticleFlowJet08_PrunedP4);
   fChain->SetBranchAddress("ParticleFlowJet08.SoftDroppedP4[5]", ParticleFlowJet08_SoftDroppedP4, &b_ParticleFlowJet08_SoftDroppedP4);
   fChain->SetBranchAddress("ParticleFlowJet08.NSubJetsTrimmed", ParticleFlowJet08_NSubJetsTrimmed, &b_ParticleFlowJet08_NSubJetsTrimmed);
   fChain->SetBranchAddress("ParticleFlowJet08.NSubJetsPruned", ParticleFlowJet08_NSubJetsPruned, &b_ParticleFlowJet08_NSubJetsPruned);
   fChain->SetBranchAddress("ParticleFlowJet08.NSubJetsSoftDropped", ParticleFlowJet08_NSubJetsSoftDropped, &b_ParticleFlowJet08_NSubJetsSoftDropped);
   fChain->SetBranchAddress("ParticleFlowJet08.ExclYmerge23", ParticleFlowJet08_ExclYmerge23, &b_ParticleFlowJet08_ExclYmerge23);
   fChain->SetBranchAddress("ParticleFlowJet08.ExclYmerge34", ParticleFlowJet08_ExclYmerge34, &b_ParticleFlowJet08_ExclYmerge34);
   fChain->SetBranchAddress("ParticleFlowJet08.ExclYmerge45", ParticleFlowJet08_ExclYmerge45, &b_ParticleFlowJet08_ExclYmerge45);
   fChain->SetBranchAddress("ParticleFlowJet08.ExclYmerge56", ParticleFlowJet08_ExclYmerge56, &b_ParticleFlowJet08_ExclYmerge56);
   fChain->SetBranchAddress("ParticleFlowJet08.Constituents", ParticleFlowJet08_Constituents, &b_ParticleFlowJet08_Constituents);
   fChain->SetBranchAddress("ParticleFlowJet08.Particles", ParticleFlowJet08_Particles, &b_ParticleFlowJet08_Particles);
   fChain->SetBranchAddress("ParticleFlowJet08.Area", ParticleFlowJet08_Area, &b_ParticleFlowJet08_Area);
   fChain->SetBranchAddress("ParticleFlowJet08_size", &ParticleFlowJet08_size, &b_ParticleFlowJet08_size);
   fChain->SetBranchAddress("ParticleFlowJet15", &ParticleFlowJet15_, &b_ParticleFlowJet15_);
   fChain->SetBranchAddress("ParticleFlowJet15.fUniqueID", ParticleFlowJet15_fUniqueID, &b_ParticleFlowJet15_fUniqueID);
   fChain->SetBranchAddress("ParticleFlowJet15.fBits", ParticleFlowJet15_fBits, &b_ParticleFlowJet15_fBits);
   fChain->SetBranchAddress("ParticleFlowJet15.PT", ParticleFlowJet15_PT, &b_ParticleFlowJet15_PT);
   fChain->SetBranchAddress("ParticleFlowJet15.Eta", ParticleFlowJet15_Eta, &b_ParticleFlowJet15_Eta);
   fChain->SetBranchAddress("ParticleFlowJet15.Phi", ParticleFlowJet15_Phi, &b_ParticleFlowJet15_Phi);
   fChain->SetBranchAddress("ParticleFlowJet15.T", ParticleFlowJet15_T, &b_ParticleFlowJet15_T);
   fChain->SetBranchAddress("ParticleFlowJet15.Mass", ParticleFlowJet15_Mass, &b_ParticleFlowJet15_Mass);
   fChain->SetBranchAddress("ParticleFlowJet15.DeltaEta", ParticleFlowJet15_DeltaEta, &b_ParticleFlowJet15_DeltaEta);
   fChain->SetBranchAddress("ParticleFlowJet15.DeltaPhi", ParticleFlowJet15_DeltaPhi, &b_ParticleFlowJet15_DeltaPhi);
   fChain->SetBranchAddress("ParticleFlowJet15.Flavor", ParticleFlowJet15_Flavor, &b_ParticleFlowJet15_Flavor);
   fChain->SetBranchAddress("ParticleFlowJet15.FlavorAlgo", ParticleFlowJet15_FlavorAlgo, &b_ParticleFlowJet15_FlavorAlgo);
   fChain->SetBranchAddress("ParticleFlowJet15.FlavorPhys", ParticleFlowJet15_FlavorPhys, &b_ParticleFlowJet15_FlavorPhys);
   fChain->SetBranchAddress("ParticleFlowJet15.BTag", ParticleFlowJet15_BTag, &b_ParticleFlowJet15_BTag);
   fChain->SetBranchAddress("ParticleFlowJet15.BTagAlgo", ParticleFlowJet15_BTagAlgo, &b_ParticleFlowJet15_BTagAlgo);
   fChain->SetBranchAddress("ParticleFlowJet15.BTagPhys", ParticleFlowJet15_BTagPhys, &b_ParticleFlowJet15_BTagPhys);
   fChain->SetBranchAddress("ParticleFlowJet15.TauTag", ParticleFlowJet15_TauTag, &b_ParticleFlowJet15_TauTag);
   fChain->SetBranchAddress("ParticleFlowJet15.TauWeight", ParticleFlowJet15_TauWeight, &b_ParticleFlowJet15_TauWeight);
   fChain->SetBranchAddress("ParticleFlowJet15.Charge", ParticleFlowJet15_Charge, &b_ParticleFlowJet15_Charge);
   fChain->SetBranchAddress("ParticleFlowJet15.EhadOverEem", ParticleFlowJet15_EhadOverEem, &b_ParticleFlowJet15_EhadOverEem);
   fChain->SetBranchAddress("ParticleFlowJet15.NCharged", ParticleFlowJet15_NCharged, &b_ParticleFlowJet15_NCharged);
   fChain->SetBranchAddress("ParticleFlowJet15.NNeutrals", ParticleFlowJet15_NNeutrals, &b_ParticleFlowJet15_NNeutrals);
   fChain->SetBranchAddress("ParticleFlowJet15.NeutralEnergyFraction", ParticleFlowJet15_NeutralEnergyFraction, &b_ParticleFlowJet15_NeutralEnergyFraction);
   fChain->SetBranchAddress("ParticleFlowJet15.ChargedEnergyFraction", ParticleFlowJet15_ChargedEnergyFraction, &b_ParticleFlowJet15_ChargedEnergyFraction);
   fChain->SetBranchAddress("ParticleFlowJet15.Beta", ParticleFlowJet15_Beta, &b_ParticleFlowJet15_Beta);
   fChain->SetBranchAddress("ParticleFlowJet15.BetaStar", ParticleFlowJet15_BetaStar, &b_ParticleFlowJet15_BetaStar);
   fChain->SetBranchAddress("ParticleFlowJet15.MeanSqDeltaR", ParticleFlowJet15_MeanSqDeltaR, &b_ParticleFlowJet15_MeanSqDeltaR);
   fChain->SetBranchAddress("ParticleFlowJet15.PTD", ParticleFlowJet15_PTD, &b_ParticleFlowJet15_PTD);
   fChain->SetBranchAddress("ParticleFlowJet15.FracPt[5]", ParticleFlowJet15_FracPt, &b_ParticleFlowJet15_FracPt);
   fChain->SetBranchAddress("ParticleFlowJet15.Tau[5]", ParticleFlowJet15_Tau, &b_ParticleFlowJet15_Tau);
   fChain->SetBranchAddress("ParticleFlowJet15.SoftDroppedJet", ParticleFlowJet15_SoftDroppedJet, &b_ParticleFlowJet15_SoftDroppedJet);
   fChain->SetBranchAddress("ParticleFlowJet15.SoftDroppedSubJet1", ParticleFlowJet15_SoftDroppedSubJet1, &b_ParticleFlowJet15_SoftDroppedSubJet1);
   fChain->SetBranchAddress("ParticleFlowJet15.SoftDroppedSubJet2", ParticleFlowJet15_SoftDroppedSubJet2, &b_ParticleFlowJet15_SoftDroppedSubJet2);
   fChain->SetBranchAddress("ParticleFlowJet15.TrimmedP4[5]", ParticleFlowJet15_TrimmedP4, &b_ParticleFlowJet15_TrimmedP4);
   fChain->SetBranchAddress("ParticleFlowJet15.PrunedP4[5]", ParticleFlowJet15_PrunedP4, &b_ParticleFlowJet15_PrunedP4);
   fChain->SetBranchAddress("ParticleFlowJet15.SoftDroppedP4[5]", ParticleFlowJet15_SoftDroppedP4, &b_ParticleFlowJet15_SoftDroppedP4);
   fChain->SetBranchAddress("ParticleFlowJet15.NSubJetsTrimmed", ParticleFlowJet15_NSubJetsTrimmed, &b_ParticleFlowJet15_NSubJetsTrimmed);
   fChain->SetBranchAddress("ParticleFlowJet15.NSubJetsPruned", ParticleFlowJet15_NSubJetsPruned, &b_ParticleFlowJet15_NSubJetsPruned);
   fChain->SetBranchAddress("ParticleFlowJet15.NSubJetsSoftDropped", ParticleFlowJet15_NSubJetsSoftDropped, &b_ParticleFlowJet15_NSubJetsSoftDropped);
   fChain->SetBranchAddress("ParticleFlowJet15.ExclYmerge23", ParticleFlowJet15_ExclYmerge23, &b_ParticleFlowJet15_ExclYmerge23);
   fChain->SetBranchAddress("ParticleFlowJet15.ExclYmerge34", ParticleFlowJet15_ExclYmerge34, &b_ParticleFlowJet15_ExclYmerge34);
   fChain->SetBranchAddress("ParticleFlowJet15.ExclYmerge45", ParticleFlowJet15_ExclYmerge45, &b_ParticleFlowJet15_ExclYmerge45);
   fChain->SetBranchAddress("ParticleFlowJet15.ExclYmerge56", ParticleFlowJet15_ExclYmerge56, &b_ParticleFlowJet15_ExclYmerge56);
   fChain->SetBranchAddress("ParticleFlowJet15.Constituents", ParticleFlowJet15_Constituents, &b_ParticleFlowJet15_Constituents);
   fChain->SetBranchAddress("ParticleFlowJet15.Particles", ParticleFlowJet15_Particles, &b_ParticleFlowJet15_Particles);
   fChain->SetBranchAddress("ParticleFlowJet15.Area", ParticleFlowJet15_Area, &b_ParticleFlowJet15_Area);
   fChain->SetBranchAddress("ParticleFlowJet15_size", &ParticleFlowJet15_size, &b_ParticleFlowJet15_size);
   fChain->SetBranchAddress("CaloJet02", &CaloJet02_, &b_CaloJet02_);
   fChain->SetBranchAddress("CaloJet02.fUniqueID", CaloJet02_fUniqueID, &b_CaloJet02_fUniqueID);
   fChain->SetBranchAddress("CaloJet02.fBits", CaloJet02_fBits, &b_CaloJet02_fBits);
   fChain->SetBranchAddress("CaloJet02.PT", CaloJet02_PT, &b_CaloJet02_PT);
   fChain->SetBranchAddress("CaloJet02.Eta", CaloJet02_Eta, &b_CaloJet02_Eta);
   fChain->SetBranchAddress("CaloJet02.Phi", CaloJet02_Phi, &b_CaloJet02_Phi);
   fChain->SetBranchAddress("CaloJet02.T", CaloJet02_T, &b_CaloJet02_T);
   fChain->SetBranchAddress("CaloJet02.Mass", CaloJet02_Mass, &b_CaloJet02_Mass);
   fChain->SetBranchAddress("CaloJet02.DeltaEta", CaloJet02_DeltaEta, &b_CaloJet02_DeltaEta);
   fChain->SetBranchAddress("CaloJet02.DeltaPhi", CaloJet02_DeltaPhi, &b_CaloJet02_DeltaPhi);
   fChain->SetBranchAddress("CaloJet02.Flavor", CaloJet02_Flavor, &b_CaloJet02_Flavor);
   fChain->SetBranchAddress("CaloJet02.FlavorAlgo", CaloJet02_FlavorAlgo, &b_CaloJet02_FlavorAlgo);
   fChain->SetBranchAddress("CaloJet02.FlavorPhys", CaloJet02_FlavorPhys, &b_CaloJet02_FlavorPhys);
   fChain->SetBranchAddress("CaloJet02.BTag", CaloJet02_BTag, &b_CaloJet02_BTag);
   fChain->SetBranchAddress("CaloJet02.BTagAlgo", CaloJet02_BTagAlgo, &b_CaloJet02_BTagAlgo);
   fChain->SetBranchAddress("CaloJet02.BTagPhys", CaloJet02_BTagPhys, &b_CaloJet02_BTagPhys);
   fChain->SetBranchAddress("CaloJet02.TauTag", CaloJet02_TauTag, &b_CaloJet02_TauTag);
   fChain->SetBranchAddress("CaloJet02.TauWeight", CaloJet02_TauWeight, &b_CaloJet02_TauWeight);
   fChain->SetBranchAddress("CaloJet02.Charge", CaloJet02_Charge, &b_CaloJet02_Charge);
   fChain->SetBranchAddress("CaloJet02.EhadOverEem", CaloJet02_EhadOverEem, &b_CaloJet02_EhadOverEem);
   fChain->SetBranchAddress("CaloJet02.NCharged", CaloJet02_NCharged, &b_CaloJet02_NCharged);
   fChain->SetBranchAddress("CaloJet02.NNeutrals", CaloJet02_NNeutrals, &b_CaloJet02_NNeutrals);
   fChain->SetBranchAddress("CaloJet02.NeutralEnergyFraction", CaloJet02_NeutralEnergyFraction, &b_CaloJet02_NeutralEnergyFraction);
   fChain->SetBranchAddress("CaloJet02.ChargedEnergyFraction", CaloJet02_ChargedEnergyFraction, &b_CaloJet02_ChargedEnergyFraction);
   fChain->SetBranchAddress("CaloJet02.Beta", CaloJet02_Beta, &b_CaloJet02_Beta);
   fChain->SetBranchAddress("CaloJet02.BetaStar", CaloJet02_BetaStar, &b_CaloJet02_BetaStar);
   fChain->SetBranchAddress("CaloJet02.MeanSqDeltaR", CaloJet02_MeanSqDeltaR, &b_CaloJet02_MeanSqDeltaR);
   fChain->SetBranchAddress("CaloJet02.PTD", CaloJet02_PTD, &b_CaloJet02_PTD);
   fChain->SetBranchAddress("CaloJet02.FracPt[5]", CaloJet02_FracPt, &b_CaloJet02_FracPt);
   fChain->SetBranchAddress("CaloJet02.Tau[5]", CaloJet02_Tau, &b_CaloJet02_Tau);
   fChain->SetBranchAddress("CaloJet02.SoftDroppedJet", CaloJet02_SoftDroppedJet, &b_CaloJet02_SoftDroppedJet);
   fChain->SetBranchAddress("CaloJet02.SoftDroppedSubJet1", CaloJet02_SoftDroppedSubJet1, &b_CaloJet02_SoftDroppedSubJet1);
   fChain->SetBranchAddress("CaloJet02.SoftDroppedSubJet2", CaloJet02_SoftDroppedSubJet2, &b_CaloJet02_SoftDroppedSubJet2);
   fChain->SetBranchAddress("CaloJet02.TrimmedP4[5]", CaloJet02_TrimmedP4, &b_CaloJet02_TrimmedP4);
   fChain->SetBranchAddress("CaloJet02.PrunedP4[5]", CaloJet02_PrunedP4, &b_CaloJet02_PrunedP4);
   fChain->SetBranchAddress("CaloJet02.SoftDroppedP4[5]", CaloJet02_SoftDroppedP4, &b_CaloJet02_SoftDroppedP4);
   fChain->SetBranchAddress("CaloJet02.NSubJetsTrimmed", CaloJet02_NSubJetsTrimmed, &b_CaloJet02_NSubJetsTrimmed);
   fChain->SetBranchAddress("CaloJet02.NSubJetsPruned", CaloJet02_NSubJetsPruned, &b_CaloJet02_NSubJetsPruned);
   fChain->SetBranchAddress("CaloJet02.NSubJetsSoftDropped", CaloJet02_NSubJetsSoftDropped, &b_CaloJet02_NSubJetsSoftDropped);
   fChain->SetBranchAddress("CaloJet02.ExclYmerge23", CaloJet02_ExclYmerge23, &b_CaloJet02_ExclYmerge23);
   fChain->SetBranchAddress("CaloJet02.ExclYmerge34", CaloJet02_ExclYmerge34, &b_CaloJet02_ExclYmerge34);
   fChain->SetBranchAddress("CaloJet02.ExclYmerge45", CaloJet02_ExclYmerge45, &b_CaloJet02_ExclYmerge45);
   fChain->SetBranchAddress("CaloJet02.ExclYmerge56", CaloJet02_ExclYmerge56, &b_CaloJet02_ExclYmerge56);
   fChain->SetBranchAddress("CaloJet02.Constituents", CaloJet02_Constituents, &b_CaloJet02_Constituents);
   fChain->SetBranchAddress("CaloJet02.Particles", CaloJet02_Particles, &b_CaloJet02_Particles);
   fChain->SetBranchAddress("CaloJet02.Area", CaloJet02_Area, &b_CaloJet02_Area);
   fChain->SetBranchAddress("CaloJet02_size", &CaloJet02_size, &b_CaloJet02_size);
   fChain->SetBranchAddress("CaloJet04", &CaloJet04_, &b_CaloJet04_);
   fChain->SetBranchAddress("CaloJet04.fUniqueID", CaloJet04_fUniqueID, &b_CaloJet04_fUniqueID);
   fChain->SetBranchAddress("CaloJet04.fBits", CaloJet04_fBits, &b_CaloJet04_fBits);
   fChain->SetBranchAddress("CaloJet04.PT", CaloJet04_PT, &b_CaloJet04_PT);
   fChain->SetBranchAddress("CaloJet04.Eta", CaloJet04_Eta, &b_CaloJet04_Eta);
   fChain->SetBranchAddress("CaloJet04.Phi", CaloJet04_Phi, &b_CaloJet04_Phi);
   fChain->SetBranchAddress("CaloJet04.T", CaloJet04_T, &b_CaloJet04_T);
   fChain->SetBranchAddress("CaloJet04.Mass", CaloJet04_Mass, &b_CaloJet04_Mass);
   fChain->SetBranchAddress("CaloJet04.DeltaEta", CaloJet04_DeltaEta, &b_CaloJet04_DeltaEta);
   fChain->SetBranchAddress("CaloJet04.DeltaPhi", CaloJet04_DeltaPhi, &b_CaloJet04_DeltaPhi);
   fChain->SetBranchAddress("CaloJet04.Flavor", CaloJet04_Flavor, &b_CaloJet04_Flavor);
   fChain->SetBranchAddress("CaloJet04.FlavorAlgo", CaloJet04_FlavorAlgo, &b_CaloJet04_FlavorAlgo);
   fChain->SetBranchAddress("CaloJet04.FlavorPhys", CaloJet04_FlavorPhys, &b_CaloJet04_FlavorPhys);
   fChain->SetBranchAddress("CaloJet04.BTag", CaloJet04_BTag, &b_CaloJet04_BTag);
   fChain->SetBranchAddress("CaloJet04.BTagAlgo", CaloJet04_BTagAlgo, &b_CaloJet04_BTagAlgo);
   fChain->SetBranchAddress("CaloJet04.BTagPhys", CaloJet04_BTagPhys, &b_CaloJet04_BTagPhys);
   fChain->SetBranchAddress("CaloJet04.TauTag", CaloJet04_TauTag, &b_CaloJet04_TauTag);
   fChain->SetBranchAddress("CaloJet04.TauWeight", CaloJet04_TauWeight, &b_CaloJet04_TauWeight);
   fChain->SetBranchAddress("CaloJet04.Charge", CaloJet04_Charge, &b_CaloJet04_Charge);
   fChain->SetBranchAddress("CaloJet04.EhadOverEem", CaloJet04_EhadOverEem, &b_CaloJet04_EhadOverEem);
   fChain->SetBranchAddress("CaloJet04.NCharged", CaloJet04_NCharged, &b_CaloJet04_NCharged);
   fChain->SetBranchAddress("CaloJet04.NNeutrals", CaloJet04_NNeutrals, &b_CaloJet04_NNeutrals);
   fChain->SetBranchAddress("CaloJet04.NeutralEnergyFraction", CaloJet04_NeutralEnergyFraction, &b_CaloJet04_NeutralEnergyFraction);
   fChain->SetBranchAddress("CaloJet04.ChargedEnergyFraction", CaloJet04_ChargedEnergyFraction, &b_CaloJet04_ChargedEnergyFraction);
   fChain->SetBranchAddress("CaloJet04.Beta", CaloJet04_Beta, &b_CaloJet04_Beta);
   fChain->SetBranchAddress("CaloJet04.BetaStar", CaloJet04_BetaStar, &b_CaloJet04_BetaStar);
   fChain->SetBranchAddress("CaloJet04.MeanSqDeltaR", CaloJet04_MeanSqDeltaR, &b_CaloJet04_MeanSqDeltaR);
   fChain->SetBranchAddress("CaloJet04.PTD", CaloJet04_PTD, &b_CaloJet04_PTD);
   fChain->SetBranchAddress("CaloJet04.FracPt[5]", CaloJet04_FracPt, &b_CaloJet04_FracPt);
   fChain->SetBranchAddress("CaloJet04.Tau[5]", CaloJet04_Tau, &b_CaloJet04_Tau);
   fChain->SetBranchAddress("CaloJet04.SoftDroppedJet", CaloJet04_SoftDroppedJet, &b_CaloJet04_SoftDroppedJet);
   fChain->SetBranchAddress("CaloJet04.SoftDroppedSubJet1", CaloJet04_SoftDroppedSubJet1, &b_CaloJet04_SoftDroppedSubJet1);
   fChain->SetBranchAddress("CaloJet04.SoftDroppedSubJet2", CaloJet04_SoftDroppedSubJet2, &b_CaloJet04_SoftDroppedSubJet2);
   fChain->SetBranchAddress("CaloJet04.TrimmedP4[5]", CaloJet04_TrimmedP4, &b_CaloJet04_TrimmedP4);
   fChain->SetBranchAddress("CaloJet04.PrunedP4[5]", CaloJet04_PrunedP4, &b_CaloJet04_PrunedP4);
   fChain->SetBranchAddress("CaloJet04.SoftDroppedP4[5]", CaloJet04_SoftDroppedP4, &b_CaloJet04_SoftDroppedP4);
   fChain->SetBranchAddress("CaloJet04.NSubJetsTrimmed", CaloJet04_NSubJetsTrimmed, &b_CaloJet04_NSubJetsTrimmed);
   fChain->SetBranchAddress("CaloJet04.NSubJetsPruned", CaloJet04_NSubJetsPruned, &b_CaloJet04_NSubJetsPruned);
   fChain->SetBranchAddress("CaloJet04.NSubJetsSoftDropped", CaloJet04_NSubJetsSoftDropped, &b_CaloJet04_NSubJetsSoftDropped);
   fChain->SetBranchAddress("CaloJet04.ExclYmerge23", CaloJet04_ExclYmerge23, &b_CaloJet04_ExclYmerge23);
   fChain->SetBranchAddress("CaloJet04.ExclYmerge34", CaloJet04_ExclYmerge34, &b_CaloJet04_ExclYmerge34);
   fChain->SetBranchAddress("CaloJet04.ExclYmerge45", CaloJet04_ExclYmerge45, &b_CaloJet04_ExclYmerge45);
   fChain->SetBranchAddress("CaloJet04.ExclYmerge56", CaloJet04_ExclYmerge56, &b_CaloJet04_ExclYmerge56);
   fChain->SetBranchAddress("CaloJet04.Constituents", CaloJet04_Constituents, &b_CaloJet04_Constituents);
   fChain->SetBranchAddress("CaloJet04.Particles", CaloJet04_Particles, &b_CaloJet04_Particles);
   fChain->SetBranchAddress("CaloJet04.Area", CaloJet04_Area, &b_CaloJet04_Area);
   fChain->SetBranchAddress("CaloJet04_size", &CaloJet04_size, &b_CaloJet04_size);
   fChain->SetBranchAddress("CaloJet08", &CaloJet08_, &b_CaloJet08_);
   fChain->SetBranchAddress("CaloJet08.fUniqueID", CaloJet08_fUniqueID, &b_CaloJet08_fUniqueID);
   fChain->SetBranchAddress("CaloJet08.fBits", CaloJet08_fBits, &b_CaloJet08_fBits);
   fChain->SetBranchAddress("CaloJet08.PT", CaloJet08_PT, &b_CaloJet08_PT);
   fChain->SetBranchAddress("CaloJet08.Eta", CaloJet08_Eta, &b_CaloJet08_Eta);
   fChain->SetBranchAddress("CaloJet08.Phi", CaloJet08_Phi, &b_CaloJet08_Phi);
   fChain->SetBranchAddress("CaloJet08.T", CaloJet08_T, &b_CaloJet08_T);
   fChain->SetBranchAddress("CaloJet08.Mass", CaloJet08_Mass, &b_CaloJet08_Mass);
   fChain->SetBranchAddress("CaloJet08.DeltaEta", CaloJet08_DeltaEta, &b_CaloJet08_DeltaEta);
   fChain->SetBranchAddress("CaloJet08.DeltaPhi", CaloJet08_DeltaPhi, &b_CaloJet08_DeltaPhi);
   fChain->SetBranchAddress("CaloJet08.Flavor", CaloJet08_Flavor, &b_CaloJet08_Flavor);
   fChain->SetBranchAddress("CaloJet08.FlavorAlgo", CaloJet08_FlavorAlgo, &b_CaloJet08_FlavorAlgo);
   fChain->SetBranchAddress("CaloJet08.FlavorPhys", CaloJet08_FlavorPhys, &b_CaloJet08_FlavorPhys);
   fChain->SetBranchAddress("CaloJet08.BTag", CaloJet08_BTag, &b_CaloJet08_BTag);
   fChain->SetBranchAddress("CaloJet08.BTagAlgo", CaloJet08_BTagAlgo, &b_CaloJet08_BTagAlgo);
   fChain->SetBranchAddress("CaloJet08.BTagPhys", CaloJet08_BTagPhys, &b_CaloJet08_BTagPhys);
   fChain->SetBranchAddress("CaloJet08.TauTag", CaloJet08_TauTag, &b_CaloJet08_TauTag);
   fChain->SetBranchAddress("CaloJet08.TauWeight", CaloJet08_TauWeight, &b_CaloJet08_TauWeight);
   fChain->SetBranchAddress("CaloJet08.Charge", CaloJet08_Charge, &b_CaloJet08_Charge);
   fChain->SetBranchAddress("CaloJet08.EhadOverEem", CaloJet08_EhadOverEem, &b_CaloJet08_EhadOverEem);
   fChain->SetBranchAddress("CaloJet08.NCharged", CaloJet08_NCharged, &b_CaloJet08_NCharged);
   fChain->SetBranchAddress("CaloJet08.NNeutrals", CaloJet08_NNeutrals, &b_CaloJet08_NNeutrals);
   fChain->SetBranchAddress("CaloJet08.NeutralEnergyFraction", CaloJet08_NeutralEnergyFraction, &b_CaloJet08_NeutralEnergyFraction);
   fChain->SetBranchAddress("CaloJet08.ChargedEnergyFraction", CaloJet08_ChargedEnergyFraction, &b_CaloJet08_ChargedEnergyFraction);
   fChain->SetBranchAddress("CaloJet08.Beta", CaloJet08_Beta, &b_CaloJet08_Beta);
   fChain->SetBranchAddress("CaloJet08.BetaStar", CaloJet08_BetaStar, &b_CaloJet08_BetaStar);
   fChain->SetBranchAddress("CaloJet08.MeanSqDeltaR", CaloJet08_MeanSqDeltaR, &b_CaloJet08_MeanSqDeltaR);
   fChain->SetBranchAddress("CaloJet08.PTD", CaloJet08_PTD, &b_CaloJet08_PTD);
   fChain->SetBranchAddress("CaloJet08.FracPt[5]", CaloJet08_FracPt, &b_CaloJet08_FracPt);
   fChain->SetBranchAddress("CaloJet08.Tau[5]", CaloJet08_Tau, &b_CaloJet08_Tau);
   fChain->SetBranchAddress("CaloJet08.SoftDroppedJet", CaloJet08_SoftDroppedJet, &b_CaloJet08_SoftDroppedJet);
   fChain->SetBranchAddress("CaloJet08.SoftDroppedSubJet1", CaloJet08_SoftDroppedSubJet1, &b_CaloJet08_SoftDroppedSubJet1);
   fChain->SetBranchAddress("CaloJet08.SoftDroppedSubJet2", CaloJet08_SoftDroppedSubJet2, &b_CaloJet08_SoftDroppedSubJet2);
   fChain->SetBranchAddress("CaloJet08.TrimmedP4[5]", CaloJet08_TrimmedP4, &b_CaloJet08_TrimmedP4);
   fChain->SetBranchAddress("CaloJet08.PrunedP4[5]", CaloJet08_PrunedP4, &b_CaloJet08_PrunedP4);
   fChain->SetBranchAddress("CaloJet08.SoftDroppedP4[5]", CaloJet08_SoftDroppedP4, &b_CaloJet08_SoftDroppedP4);
   fChain->SetBranchAddress("CaloJet08.NSubJetsTrimmed", CaloJet08_NSubJetsTrimmed, &b_CaloJet08_NSubJetsTrimmed);
   fChain->SetBranchAddress("CaloJet08.NSubJetsPruned", CaloJet08_NSubJetsPruned, &b_CaloJet08_NSubJetsPruned);
   fChain->SetBranchAddress("CaloJet08.NSubJetsSoftDropped", CaloJet08_NSubJetsSoftDropped, &b_CaloJet08_NSubJetsSoftDropped);
   fChain->SetBranchAddress("CaloJet08.ExclYmerge23", CaloJet08_ExclYmerge23, &b_CaloJet08_ExclYmerge23);
   fChain->SetBranchAddress("CaloJet08.ExclYmerge34", CaloJet08_ExclYmerge34, &b_CaloJet08_ExclYmerge34);
   fChain->SetBranchAddress("CaloJet08.ExclYmerge45", CaloJet08_ExclYmerge45, &b_CaloJet08_ExclYmerge45);
   fChain->SetBranchAddress("CaloJet08.ExclYmerge56", CaloJet08_ExclYmerge56, &b_CaloJet08_ExclYmerge56);
   fChain->SetBranchAddress("CaloJet08.Constituents", CaloJet08_Constituents, &b_CaloJet08_Constituents);
   fChain->SetBranchAddress("CaloJet08.Particles", CaloJet08_Particles, &b_CaloJet08_Particles);
   fChain->SetBranchAddress("CaloJet08.Area", CaloJet08_Area, &b_CaloJet08_Area);
   fChain->SetBranchAddress("CaloJet08_size", &CaloJet08_size, &b_CaloJet08_size);
   fChain->SetBranchAddress("CaloJet15", &CaloJet15_, &b_CaloJet15_);
   fChain->SetBranchAddress("CaloJet15.fUniqueID", CaloJet15_fUniqueID, &b_CaloJet15_fUniqueID);
   fChain->SetBranchAddress("CaloJet15.fBits", CaloJet15_fBits, &b_CaloJet15_fBits);
   fChain->SetBranchAddress("CaloJet15.PT", CaloJet15_PT, &b_CaloJet15_PT);
   fChain->SetBranchAddress("CaloJet15.Eta", CaloJet15_Eta, &b_CaloJet15_Eta);
   fChain->SetBranchAddress("CaloJet15.Phi", CaloJet15_Phi, &b_CaloJet15_Phi);
   fChain->SetBranchAddress("CaloJet15.T", CaloJet15_T, &b_CaloJet15_T);
   fChain->SetBranchAddress("CaloJet15.Mass", CaloJet15_Mass, &b_CaloJet15_Mass);
   fChain->SetBranchAddress("CaloJet15.DeltaEta", CaloJet15_DeltaEta, &b_CaloJet15_DeltaEta);
   fChain->SetBranchAddress("CaloJet15.DeltaPhi", CaloJet15_DeltaPhi, &b_CaloJet15_DeltaPhi);
   fChain->SetBranchAddress("CaloJet15.Flavor", CaloJet15_Flavor, &b_CaloJet15_Flavor);
   fChain->SetBranchAddress("CaloJet15.FlavorAlgo", CaloJet15_FlavorAlgo, &b_CaloJet15_FlavorAlgo);
   fChain->SetBranchAddress("CaloJet15.FlavorPhys", CaloJet15_FlavorPhys, &b_CaloJet15_FlavorPhys);
   fChain->SetBranchAddress("CaloJet15.BTag", CaloJet15_BTag, &b_CaloJet15_BTag);
   fChain->SetBranchAddress("CaloJet15.BTagAlgo", CaloJet15_BTagAlgo, &b_CaloJet15_BTagAlgo);
   fChain->SetBranchAddress("CaloJet15.BTagPhys", CaloJet15_BTagPhys, &b_CaloJet15_BTagPhys);
   fChain->SetBranchAddress("CaloJet15.TauTag", CaloJet15_TauTag, &b_CaloJet15_TauTag);
   fChain->SetBranchAddress("CaloJet15.TauWeight", CaloJet15_TauWeight, &b_CaloJet15_TauWeight);
   fChain->SetBranchAddress("CaloJet15.Charge", CaloJet15_Charge, &b_CaloJet15_Charge);
   fChain->SetBranchAddress("CaloJet15.EhadOverEem", CaloJet15_EhadOverEem, &b_CaloJet15_EhadOverEem);
   fChain->SetBranchAddress("CaloJet15.NCharged", CaloJet15_NCharged, &b_CaloJet15_NCharged);
   fChain->SetBranchAddress("CaloJet15.NNeutrals", CaloJet15_NNeutrals, &b_CaloJet15_NNeutrals);
   fChain->SetBranchAddress("CaloJet15.NeutralEnergyFraction", CaloJet15_NeutralEnergyFraction, &b_CaloJet15_NeutralEnergyFraction);
   fChain->SetBranchAddress("CaloJet15.ChargedEnergyFraction", CaloJet15_ChargedEnergyFraction, &b_CaloJet15_ChargedEnergyFraction);
   fChain->SetBranchAddress("CaloJet15.Beta", CaloJet15_Beta, &b_CaloJet15_Beta);
   fChain->SetBranchAddress("CaloJet15.BetaStar", CaloJet15_BetaStar, &b_CaloJet15_BetaStar);
   fChain->SetBranchAddress("CaloJet15.MeanSqDeltaR", CaloJet15_MeanSqDeltaR, &b_CaloJet15_MeanSqDeltaR);
   fChain->SetBranchAddress("CaloJet15.PTD", CaloJet15_PTD, &b_CaloJet15_PTD);
   fChain->SetBranchAddress("CaloJet15.FracPt[5]", CaloJet15_FracPt, &b_CaloJet15_FracPt);
   fChain->SetBranchAddress("CaloJet15.Tau[5]", CaloJet15_Tau, &b_CaloJet15_Tau);
   fChain->SetBranchAddress("CaloJet15.SoftDroppedJet", CaloJet15_SoftDroppedJet, &b_CaloJet15_SoftDroppedJet);
   fChain->SetBranchAddress("CaloJet15.SoftDroppedSubJet1", CaloJet15_SoftDroppedSubJet1, &b_CaloJet15_SoftDroppedSubJet1);
   fChain->SetBranchAddress("CaloJet15.SoftDroppedSubJet2", CaloJet15_SoftDroppedSubJet2, &b_CaloJet15_SoftDroppedSubJet2);
   fChain->SetBranchAddress("CaloJet15.TrimmedP4[5]", CaloJet15_TrimmedP4, &b_CaloJet15_TrimmedP4);
   fChain->SetBranchAddress("CaloJet15.PrunedP4[5]", CaloJet15_PrunedP4, &b_CaloJet15_PrunedP4);
   fChain->SetBranchAddress("CaloJet15.SoftDroppedP4[5]", CaloJet15_SoftDroppedP4, &b_CaloJet15_SoftDroppedP4);
   fChain->SetBranchAddress("CaloJet15.NSubJetsTrimmed", CaloJet15_NSubJetsTrimmed, &b_CaloJet15_NSubJetsTrimmed);
   fChain->SetBranchAddress("CaloJet15.NSubJetsPruned", CaloJet15_NSubJetsPruned, &b_CaloJet15_NSubJetsPruned);
   fChain->SetBranchAddress("CaloJet15.NSubJetsSoftDropped", CaloJet15_NSubJetsSoftDropped, &b_CaloJet15_NSubJetsSoftDropped);
   fChain->SetBranchAddress("CaloJet15.ExclYmerge23", CaloJet15_ExclYmerge23, &b_CaloJet15_ExclYmerge23);
   fChain->SetBranchAddress("CaloJet15.ExclYmerge34", CaloJet15_ExclYmerge34, &b_CaloJet15_ExclYmerge34);
   fChain->SetBranchAddress("CaloJet15.ExclYmerge45", CaloJet15_ExclYmerge45, &b_CaloJet15_ExclYmerge45);
   fChain->SetBranchAddress("CaloJet15.ExclYmerge56", CaloJet15_ExclYmerge56, &b_CaloJet15_ExclYmerge56);
   fChain->SetBranchAddress("CaloJet15.Constituents", CaloJet15_Constituents, &b_CaloJet15_Constituents);
   fChain->SetBranchAddress("CaloJet15.Particles", CaloJet15_Particles, &b_CaloJet15_Particles);
   fChain->SetBranchAddress("CaloJet15.Area", CaloJet15_Area, &b_CaloJet15_Area);
   fChain->SetBranchAddress("CaloJet15_size", &CaloJet15_size, &b_CaloJet15_size);
   fChain->SetBranchAddress("TrackJet02", &TrackJet02_, &b_TrackJet02_);
   fChain->SetBranchAddress("TrackJet02.fUniqueID", TrackJet02_fUniqueID, &b_TrackJet02_fUniqueID);
   fChain->SetBranchAddress("TrackJet02.fBits", TrackJet02_fBits, &b_TrackJet02_fBits);
   fChain->SetBranchAddress("TrackJet02.PT", TrackJet02_PT, &b_TrackJet02_PT);
   fChain->SetBranchAddress("TrackJet02.Eta", TrackJet02_Eta, &b_TrackJet02_Eta);
   fChain->SetBranchAddress("TrackJet02.Phi", TrackJet02_Phi, &b_TrackJet02_Phi);
   fChain->SetBranchAddress("TrackJet02.T", TrackJet02_T, &b_TrackJet02_T);
   fChain->SetBranchAddress("TrackJet02.Mass", TrackJet02_Mass, &b_TrackJet02_Mass);
   fChain->SetBranchAddress("TrackJet02.DeltaEta", TrackJet02_DeltaEta, &b_TrackJet02_DeltaEta);
   fChain->SetBranchAddress("TrackJet02.DeltaPhi", TrackJet02_DeltaPhi, &b_TrackJet02_DeltaPhi);
   fChain->SetBranchAddress("TrackJet02.Flavor", TrackJet02_Flavor, &b_TrackJet02_Flavor);
   fChain->SetBranchAddress("TrackJet02.FlavorAlgo", TrackJet02_FlavorAlgo, &b_TrackJet02_FlavorAlgo);
   fChain->SetBranchAddress("TrackJet02.FlavorPhys", TrackJet02_FlavorPhys, &b_TrackJet02_FlavorPhys);
   fChain->SetBranchAddress("TrackJet02.BTag", TrackJet02_BTag, &b_TrackJet02_BTag);
   fChain->SetBranchAddress("TrackJet02.BTagAlgo", TrackJet02_BTagAlgo, &b_TrackJet02_BTagAlgo);
   fChain->SetBranchAddress("TrackJet02.BTagPhys", TrackJet02_BTagPhys, &b_TrackJet02_BTagPhys);
   fChain->SetBranchAddress("TrackJet02.TauTag", TrackJet02_TauTag, &b_TrackJet02_TauTag);
   fChain->SetBranchAddress("TrackJet02.TauWeight", TrackJet02_TauWeight, &b_TrackJet02_TauWeight);
   fChain->SetBranchAddress("TrackJet02.Charge", TrackJet02_Charge, &b_TrackJet02_Charge);
   fChain->SetBranchAddress("TrackJet02.EhadOverEem", TrackJet02_EhadOverEem, &b_TrackJet02_EhadOverEem);
   fChain->SetBranchAddress("TrackJet02.NCharged", TrackJet02_NCharged, &b_TrackJet02_NCharged);
   fChain->SetBranchAddress("TrackJet02.NNeutrals", TrackJet02_NNeutrals, &b_TrackJet02_NNeutrals);
   fChain->SetBranchAddress("TrackJet02.NeutralEnergyFraction", TrackJet02_NeutralEnergyFraction, &b_TrackJet02_NeutralEnergyFraction);
   fChain->SetBranchAddress("TrackJet02.ChargedEnergyFraction", TrackJet02_ChargedEnergyFraction, &b_TrackJet02_ChargedEnergyFraction);
   fChain->SetBranchAddress("TrackJet02.Beta", TrackJet02_Beta, &b_TrackJet02_Beta);
   fChain->SetBranchAddress("TrackJet02.BetaStar", TrackJet02_BetaStar, &b_TrackJet02_BetaStar);
   fChain->SetBranchAddress("TrackJet02.MeanSqDeltaR", TrackJet02_MeanSqDeltaR, &b_TrackJet02_MeanSqDeltaR);
   fChain->SetBranchAddress("TrackJet02.PTD", TrackJet02_PTD, &b_TrackJet02_PTD);
   fChain->SetBranchAddress("TrackJet02.FracPt[5]", TrackJet02_FracPt, &b_TrackJet02_FracPt);
   fChain->SetBranchAddress("TrackJet02.Tau[5]", TrackJet02_Tau, &b_TrackJet02_Tau);
   fChain->SetBranchAddress("TrackJet02.SoftDroppedJet", TrackJet02_SoftDroppedJet, &b_TrackJet02_SoftDroppedJet);
   fChain->SetBranchAddress("TrackJet02.SoftDroppedSubJet1", TrackJet02_SoftDroppedSubJet1, &b_TrackJet02_SoftDroppedSubJet1);
   fChain->SetBranchAddress("TrackJet02.SoftDroppedSubJet2", TrackJet02_SoftDroppedSubJet2, &b_TrackJet02_SoftDroppedSubJet2);
   fChain->SetBranchAddress("TrackJet02.TrimmedP4[5]", TrackJet02_TrimmedP4, &b_TrackJet02_TrimmedP4);
   fChain->SetBranchAddress("TrackJet02.PrunedP4[5]", TrackJet02_PrunedP4, &b_TrackJet02_PrunedP4);
   fChain->SetBranchAddress("TrackJet02.SoftDroppedP4[5]", TrackJet02_SoftDroppedP4, &b_TrackJet02_SoftDroppedP4);
   fChain->SetBranchAddress("TrackJet02.NSubJetsTrimmed", TrackJet02_NSubJetsTrimmed, &b_TrackJet02_NSubJetsTrimmed);
   fChain->SetBranchAddress("TrackJet02.NSubJetsPruned", TrackJet02_NSubJetsPruned, &b_TrackJet02_NSubJetsPruned);
   fChain->SetBranchAddress("TrackJet02.NSubJetsSoftDropped", TrackJet02_NSubJetsSoftDropped, &b_TrackJet02_NSubJetsSoftDropped);
   fChain->SetBranchAddress("TrackJet02.ExclYmerge23", TrackJet02_ExclYmerge23, &b_TrackJet02_ExclYmerge23);
   fChain->SetBranchAddress("TrackJet02.ExclYmerge34", TrackJet02_ExclYmerge34, &b_TrackJet02_ExclYmerge34);
   fChain->SetBranchAddress("TrackJet02.ExclYmerge45", TrackJet02_ExclYmerge45, &b_TrackJet02_ExclYmerge45);
   fChain->SetBranchAddress("TrackJet02.ExclYmerge56", TrackJet02_ExclYmerge56, &b_TrackJet02_ExclYmerge56);
   fChain->SetBranchAddress("TrackJet02.Constituents", TrackJet02_Constituents, &b_TrackJet02_Constituents);
   fChain->SetBranchAddress("TrackJet02.Particles", TrackJet02_Particles, &b_TrackJet02_Particles);
   fChain->SetBranchAddress("TrackJet02.Area", TrackJet02_Area, &b_TrackJet02_Area);
   fChain->SetBranchAddress("TrackJet02_size", &TrackJet02_size, &b_TrackJet02_size);
   fChain->SetBranchAddress("TrackJet04", &TrackJet04_, &b_TrackJet04_);
   fChain->SetBranchAddress("TrackJet04.fUniqueID", TrackJet04_fUniqueID, &b_TrackJet04_fUniqueID);
   fChain->SetBranchAddress("TrackJet04.fBits", TrackJet04_fBits, &b_TrackJet04_fBits);
   fChain->SetBranchAddress("TrackJet04.PT", TrackJet04_PT, &b_TrackJet04_PT);
   fChain->SetBranchAddress("TrackJet04.Eta", TrackJet04_Eta, &b_TrackJet04_Eta);
   fChain->SetBranchAddress("TrackJet04.Phi", TrackJet04_Phi, &b_TrackJet04_Phi);
   fChain->SetBranchAddress("TrackJet04.T", TrackJet04_T, &b_TrackJet04_T);
   fChain->SetBranchAddress("TrackJet04.Mass", TrackJet04_Mass, &b_TrackJet04_Mass);
   fChain->SetBranchAddress("TrackJet04.DeltaEta", TrackJet04_DeltaEta, &b_TrackJet04_DeltaEta);
   fChain->SetBranchAddress("TrackJet04.DeltaPhi", TrackJet04_DeltaPhi, &b_TrackJet04_DeltaPhi);
   fChain->SetBranchAddress("TrackJet04.Flavor", TrackJet04_Flavor, &b_TrackJet04_Flavor);
   fChain->SetBranchAddress("TrackJet04.FlavorAlgo", TrackJet04_FlavorAlgo, &b_TrackJet04_FlavorAlgo);
   fChain->SetBranchAddress("TrackJet04.FlavorPhys", TrackJet04_FlavorPhys, &b_TrackJet04_FlavorPhys);
   fChain->SetBranchAddress("TrackJet04.BTag", TrackJet04_BTag, &b_TrackJet04_BTag);
   fChain->SetBranchAddress("TrackJet04.BTagAlgo", TrackJet04_BTagAlgo, &b_TrackJet04_BTagAlgo);
   fChain->SetBranchAddress("TrackJet04.BTagPhys", TrackJet04_BTagPhys, &b_TrackJet04_BTagPhys);
   fChain->SetBranchAddress("TrackJet04.TauTag", TrackJet04_TauTag, &b_TrackJet04_TauTag);
   fChain->SetBranchAddress("TrackJet04.TauWeight", TrackJet04_TauWeight, &b_TrackJet04_TauWeight);
   fChain->SetBranchAddress("TrackJet04.Charge", TrackJet04_Charge, &b_TrackJet04_Charge);
   fChain->SetBranchAddress("TrackJet04.EhadOverEem", TrackJet04_EhadOverEem, &b_TrackJet04_EhadOverEem);
   fChain->SetBranchAddress("TrackJet04.NCharged", TrackJet04_NCharged, &b_TrackJet04_NCharged);
   fChain->SetBranchAddress("TrackJet04.NNeutrals", TrackJet04_NNeutrals, &b_TrackJet04_NNeutrals);
   fChain->SetBranchAddress("TrackJet04.NeutralEnergyFraction", TrackJet04_NeutralEnergyFraction, &b_TrackJet04_NeutralEnergyFraction);
   fChain->SetBranchAddress("TrackJet04.ChargedEnergyFraction", TrackJet04_ChargedEnergyFraction, &b_TrackJet04_ChargedEnergyFraction);
   fChain->SetBranchAddress("TrackJet04.Beta", TrackJet04_Beta, &b_TrackJet04_Beta);
   fChain->SetBranchAddress("TrackJet04.BetaStar", TrackJet04_BetaStar, &b_TrackJet04_BetaStar);
   fChain->SetBranchAddress("TrackJet04.MeanSqDeltaR", TrackJet04_MeanSqDeltaR, &b_TrackJet04_MeanSqDeltaR);
   fChain->SetBranchAddress("TrackJet04.PTD", TrackJet04_PTD, &b_TrackJet04_PTD);
   fChain->SetBranchAddress("TrackJet04.FracPt[5]", TrackJet04_FracPt, &b_TrackJet04_FracPt);
   fChain->SetBranchAddress("TrackJet04.Tau[5]", TrackJet04_Tau, &b_TrackJet04_Tau);
   fChain->SetBranchAddress("TrackJet04.SoftDroppedJet", TrackJet04_SoftDroppedJet, &b_TrackJet04_SoftDroppedJet);
   fChain->SetBranchAddress("TrackJet04.SoftDroppedSubJet1", TrackJet04_SoftDroppedSubJet1, &b_TrackJet04_SoftDroppedSubJet1);
   fChain->SetBranchAddress("TrackJet04.SoftDroppedSubJet2", TrackJet04_SoftDroppedSubJet2, &b_TrackJet04_SoftDroppedSubJet2);
   fChain->SetBranchAddress("TrackJet04.TrimmedP4[5]", TrackJet04_TrimmedP4, &b_TrackJet04_TrimmedP4);
   fChain->SetBranchAddress("TrackJet04.PrunedP4[5]", TrackJet04_PrunedP4, &b_TrackJet04_PrunedP4);
   fChain->SetBranchAddress("TrackJet04.SoftDroppedP4[5]", TrackJet04_SoftDroppedP4, &b_TrackJet04_SoftDroppedP4);
   fChain->SetBranchAddress("TrackJet04.NSubJetsTrimmed", TrackJet04_NSubJetsTrimmed, &b_TrackJet04_NSubJetsTrimmed);
   fChain->SetBranchAddress("TrackJet04.NSubJetsPruned", TrackJet04_NSubJetsPruned, &b_TrackJet04_NSubJetsPruned);
   fChain->SetBranchAddress("TrackJet04.NSubJetsSoftDropped", TrackJet04_NSubJetsSoftDropped, &b_TrackJet04_NSubJetsSoftDropped);
   fChain->SetBranchAddress("TrackJet04.ExclYmerge23", TrackJet04_ExclYmerge23, &b_TrackJet04_ExclYmerge23);
   fChain->SetBranchAddress("TrackJet04.ExclYmerge34", TrackJet04_ExclYmerge34, &b_TrackJet04_ExclYmerge34);
   fChain->SetBranchAddress("TrackJet04.ExclYmerge45", TrackJet04_ExclYmerge45, &b_TrackJet04_ExclYmerge45);
   fChain->SetBranchAddress("TrackJet04.ExclYmerge56", TrackJet04_ExclYmerge56, &b_TrackJet04_ExclYmerge56);
   fChain->SetBranchAddress("TrackJet04.Constituents", TrackJet04_Constituents, &b_TrackJet04_Constituents);
   fChain->SetBranchAddress("TrackJet04.Particles", TrackJet04_Particles, &b_TrackJet04_Particles);
   fChain->SetBranchAddress("TrackJet04.Area", TrackJet04_Area, &b_TrackJet04_Area);
   fChain->SetBranchAddress("TrackJet04_size", &TrackJet04_size, &b_TrackJet04_size);
   fChain->SetBranchAddress("TrackJet08", &TrackJet08_, &b_TrackJet08_);
   fChain->SetBranchAddress("TrackJet08.fUniqueID", TrackJet08_fUniqueID, &b_TrackJet08_fUniqueID);
   fChain->SetBranchAddress("TrackJet08.fBits", TrackJet08_fBits, &b_TrackJet08_fBits);
   fChain->SetBranchAddress("TrackJet08.PT", TrackJet08_PT, &b_TrackJet08_PT);
   fChain->SetBranchAddress("TrackJet08.Eta", TrackJet08_Eta, &b_TrackJet08_Eta);
   fChain->SetBranchAddress("TrackJet08.Phi", TrackJet08_Phi, &b_TrackJet08_Phi);
   fChain->SetBranchAddress("TrackJet08.T", TrackJet08_T, &b_TrackJet08_T);
   fChain->SetBranchAddress("TrackJet08.Mass", TrackJet08_Mass, &b_TrackJet08_Mass);
   fChain->SetBranchAddress("TrackJet08.DeltaEta", TrackJet08_DeltaEta, &b_TrackJet08_DeltaEta);
   fChain->SetBranchAddress("TrackJet08.DeltaPhi", TrackJet08_DeltaPhi, &b_TrackJet08_DeltaPhi);
   fChain->SetBranchAddress("TrackJet08.Flavor", TrackJet08_Flavor, &b_TrackJet08_Flavor);
   fChain->SetBranchAddress("TrackJet08.FlavorAlgo", TrackJet08_FlavorAlgo, &b_TrackJet08_FlavorAlgo);
   fChain->SetBranchAddress("TrackJet08.FlavorPhys", TrackJet08_FlavorPhys, &b_TrackJet08_FlavorPhys);
   fChain->SetBranchAddress("TrackJet08.BTag", TrackJet08_BTag, &b_TrackJet08_BTag);
   fChain->SetBranchAddress("TrackJet08.BTagAlgo", TrackJet08_BTagAlgo, &b_TrackJet08_BTagAlgo);
   fChain->SetBranchAddress("TrackJet08.BTagPhys", TrackJet08_BTagPhys, &b_TrackJet08_BTagPhys);
   fChain->SetBranchAddress("TrackJet08.TauTag", TrackJet08_TauTag, &b_TrackJet08_TauTag);
   fChain->SetBranchAddress("TrackJet08.TauWeight", TrackJet08_TauWeight, &b_TrackJet08_TauWeight);
   fChain->SetBranchAddress("TrackJet08.Charge", TrackJet08_Charge, &b_TrackJet08_Charge);
   fChain->SetBranchAddress("TrackJet08.EhadOverEem", TrackJet08_EhadOverEem, &b_TrackJet08_EhadOverEem);
   fChain->SetBranchAddress("TrackJet08.NCharged", TrackJet08_NCharged, &b_TrackJet08_NCharged);
   fChain->SetBranchAddress("TrackJet08.NNeutrals", TrackJet08_NNeutrals, &b_TrackJet08_NNeutrals);
   fChain->SetBranchAddress("TrackJet08.NeutralEnergyFraction", TrackJet08_NeutralEnergyFraction, &b_TrackJet08_NeutralEnergyFraction);
   fChain->SetBranchAddress("TrackJet08.ChargedEnergyFraction", TrackJet08_ChargedEnergyFraction, &b_TrackJet08_ChargedEnergyFraction);
   fChain->SetBranchAddress("TrackJet08.Beta", TrackJet08_Beta, &b_TrackJet08_Beta);
   fChain->SetBranchAddress("TrackJet08.BetaStar", TrackJet08_BetaStar, &b_TrackJet08_BetaStar);
   fChain->SetBranchAddress("TrackJet08.MeanSqDeltaR", TrackJet08_MeanSqDeltaR, &b_TrackJet08_MeanSqDeltaR);
   fChain->SetBranchAddress("TrackJet08.PTD", TrackJet08_PTD, &b_TrackJet08_PTD);
   fChain->SetBranchAddress("TrackJet08.FracPt[5]", TrackJet08_FracPt, &b_TrackJet08_FracPt);
   fChain->SetBranchAddress("TrackJet08.Tau[5]", TrackJet08_Tau, &b_TrackJet08_Tau);
   fChain->SetBranchAddress("TrackJet08.SoftDroppedJet", TrackJet08_SoftDroppedJet, &b_TrackJet08_SoftDroppedJet);
   fChain->SetBranchAddress("TrackJet08.SoftDroppedSubJet1", TrackJet08_SoftDroppedSubJet1, &b_TrackJet08_SoftDroppedSubJet1);
   fChain->SetBranchAddress("TrackJet08.SoftDroppedSubJet2", TrackJet08_SoftDroppedSubJet2, &b_TrackJet08_SoftDroppedSubJet2);
   fChain->SetBranchAddress("TrackJet08.TrimmedP4[5]", TrackJet08_TrimmedP4, &b_TrackJet08_TrimmedP4);
   fChain->SetBranchAddress("TrackJet08.PrunedP4[5]", TrackJet08_PrunedP4, &b_TrackJet08_PrunedP4);
   fChain->SetBranchAddress("TrackJet08.SoftDroppedP4[5]", TrackJet08_SoftDroppedP4, &b_TrackJet08_SoftDroppedP4);
   fChain->SetBranchAddress("TrackJet08.NSubJetsTrimmed", TrackJet08_NSubJetsTrimmed, &b_TrackJet08_NSubJetsTrimmed);
   fChain->SetBranchAddress("TrackJet08.NSubJetsPruned", TrackJet08_NSubJetsPruned, &b_TrackJet08_NSubJetsPruned);
   fChain->SetBranchAddress("TrackJet08.NSubJetsSoftDropped", TrackJet08_NSubJetsSoftDropped, &b_TrackJet08_NSubJetsSoftDropped);
   fChain->SetBranchAddress("TrackJet08.ExclYmerge23", TrackJet08_ExclYmerge23, &b_TrackJet08_ExclYmerge23);
   fChain->SetBranchAddress("TrackJet08.ExclYmerge34", TrackJet08_ExclYmerge34, &b_TrackJet08_ExclYmerge34);
   fChain->SetBranchAddress("TrackJet08.ExclYmerge45", TrackJet08_ExclYmerge45, &b_TrackJet08_ExclYmerge45);
   fChain->SetBranchAddress("TrackJet08.ExclYmerge56", TrackJet08_ExclYmerge56, &b_TrackJet08_ExclYmerge56);
   fChain->SetBranchAddress("TrackJet08.Constituents", TrackJet08_Constituents, &b_TrackJet08_Constituents);
   fChain->SetBranchAddress("TrackJet08.Particles", TrackJet08_Particles, &b_TrackJet08_Particles);
   fChain->SetBranchAddress("TrackJet08.Area", TrackJet08_Area, &b_TrackJet08_Area);
   fChain->SetBranchAddress("TrackJet08_size", &TrackJet08_size, &b_TrackJet08_size);
   fChain->SetBranchAddress("TrackJet15", &TrackJet15_, &b_TrackJet15_);
   fChain->SetBranchAddress("TrackJet15.fUniqueID", TrackJet15_fUniqueID, &b_TrackJet15_fUniqueID);
   fChain->SetBranchAddress("TrackJet15.fBits", TrackJet15_fBits, &b_TrackJet15_fBits);
   fChain->SetBranchAddress("TrackJet15.PT", TrackJet15_PT, &b_TrackJet15_PT);
   fChain->SetBranchAddress("TrackJet15.Eta", TrackJet15_Eta, &b_TrackJet15_Eta);
   fChain->SetBranchAddress("TrackJet15.Phi", TrackJet15_Phi, &b_TrackJet15_Phi);
   fChain->SetBranchAddress("TrackJet15.T", TrackJet15_T, &b_TrackJet15_T);
   fChain->SetBranchAddress("TrackJet15.Mass", TrackJet15_Mass, &b_TrackJet15_Mass);
   fChain->SetBranchAddress("TrackJet15.DeltaEta", TrackJet15_DeltaEta, &b_TrackJet15_DeltaEta);
   fChain->SetBranchAddress("TrackJet15.DeltaPhi", TrackJet15_DeltaPhi, &b_TrackJet15_DeltaPhi);
   fChain->SetBranchAddress("TrackJet15.Flavor", TrackJet15_Flavor, &b_TrackJet15_Flavor);
   fChain->SetBranchAddress("TrackJet15.FlavorAlgo", TrackJet15_FlavorAlgo, &b_TrackJet15_FlavorAlgo);
   fChain->SetBranchAddress("TrackJet15.FlavorPhys", TrackJet15_FlavorPhys, &b_TrackJet15_FlavorPhys);
   fChain->SetBranchAddress("TrackJet15.BTag", TrackJet15_BTag, &b_TrackJet15_BTag);
   fChain->SetBranchAddress("TrackJet15.BTagAlgo", TrackJet15_BTagAlgo, &b_TrackJet15_BTagAlgo);
   fChain->SetBranchAddress("TrackJet15.BTagPhys", TrackJet15_BTagPhys, &b_TrackJet15_BTagPhys);
   fChain->SetBranchAddress("TrackJet15.TauTag", TrackJet15_TauTag, &b_TrackJet15_TauTag);
   fChain->SetBranchAddress("TrackJet15.TauWeight", TrackJet15_TauWeight, &b_TrackJet15_TauWeight);
   fChain->SetBranchAddress("TrackJet15.Charge", TrackJet15_Charge, &b_TrackJet15_Charge);
   fChain->SetBranchAddress("TrackJet15.EhadOverEem", TrackJet15_EhadOverEem, &b_TrackJet15_EhadOverEem);
   fChain->SetBranchAddress("TrackJet15.NCharged", TrackJet15_NCharged, &b_TrackJet15_NCharged);
   fChain->SetBranchAddress("TrackJet15.NNeutrals", TrackJet15_NNeutrals, &b_TrackJet15_NNeutrals);
   fChain->SetBranchAddress("TrackJet15.NeutralEnergyFraction", TrackJet15_NeutralEnergyFraction, &b_TrackJet15_NeutralEnergyFraction);
   fChain->SetBranchAddress("TrackJet15.ChargedEnergyFraction", TrackJet15_ChargedEnergyFraction, &b_TrackJet15_ChargedEnergyFraction);
   fChain->SetBranchAddress("TrackJet15.Beta", TrackJet15_Beta, &b_TrackJet15_Beta);
   fChain->SetBranchAddress("TrackJet15.BetaStar", TrackJet15_BetaStar, &b_TrackJet15_BetaStar);
   fChain->SetBranchAddress("TrackJet15.MeanSqDeltaR", TrackJet15_MeanSqDeltaR, &b_TrackJet15_MeanSqDeltaR);
   fChain->SetBranchAddress("TrackJet15.PTD", TrackJet15_PTD, &b_TrackJet15_PTD);
   fChain->SetBranchAddress("TrackJet15.FracPt[5]", TrackJet15_FracPt, &b_TrackJet15_FracPt);
   fChain->SetBranchAddress("TrackJet15.Tau[5]", TrackJet15_Tau, &b_TrackJet15_Tau);
   fChain->SetBranchAddress("TrackJet15.SoftDroppedJet", TrackJet15_SoftDroppedJet, &b_TrackJet15_SoftDroppedJet);
   fChain->SetBranchAddress("TrackJet15.SoftDroppedSubJet1", TrackJet15_SoftDroppedSubJet1, &b_TrackJet15_SoftDroppedSubJet1);
   fChain->SetBranchAddress("TrackJet15.SoftDroppedSubJet2", TrackJet15_SoftDroppedSubJet2, &b_TrackJet15_SoftDroppedSubJet2);
   fChain->SetBranchAddress("TrackJet15.TrimmedP4[5]", TrackJet15_TrimmedP4, &b_TrackJet15_TrimmedP4);
   fChain->SetBranchAddress("TrackJet15.PrunedP4[5]", TrackJet15_PrunedP4, &b_TrackJet15_PrunedP4);
   fChain->SetBranchAddress("TrackJet15.SoftDroppedP4[5]", TrackJet15_SoftDroppedP4, &b_TrackJet15_SoftDroppedP4);
   fChain->SetBranchAddress("TrackJet15.NSubJetsTrimmed", TrackJet15_NSubJetsTrimmed, &b_TrackJet15_NSubJetsTrimmed);
   fChain->SetBranchAddress("TrackJet15.NSubJetsPruned", TrackJet15_NSubJetsPruned, &b_TrackJet15_NSubJetsPruned);
   fChain->SetBranchAddress("TrackJet15.NSubJetsSoftDropped", TrackJet15_NSubJetsSoftDropped, &b_TrackJet15_NSubJetsSoftDropped);
   fChain->SetBranchAddress("TrackJet15.ExclYmerge23", TrackJet15_ExclYmerge23, &b_TrackJet15_ExclYmerge23);
   fChain->SetBranchAddress("TrackJet15.ExclYmerge34", TrackJet15_ExclYmerge34, &b_TrackJet15_ExclYmerge34);
   fChain->SetBranchAddress("TrackJet15.ExclYmerge45", TrackJet15_ExclYmerge45, &b_TrackJet15_ExclYmerge45);
   fChain->SetBranchAddress("TrackJet15.ExclYmerge56", TrackJet15_ExclYmerge56, &b_TrackJet15_ExclYmerge56);
   fChain->SetBranchAddress("TrackJet15.Constituents", TrackJet15_Constituents, &b_TrackJet15_Constituents);
   fChain->SetBranchAddress("TrackJet15.Particles", TrackJet15_Particles, &b_TrackJet15_Particles);
   fChain->SetBranchAddress("TrackJet15.Area", TrackJet15_Area, &b_TrackJet15_Area);
   fChain->SetBranchAddress("TrackJet15_size", &TrackJet15_size, &b_TrackJet15_size);
   fChain->SetBranchAddress("MissingET", &MissingET_, &b_MissingET_);
   fChain->SetBranchAddress("MissingET.fUniqueID", MissingET_fUniqueID, &b_MissingET_fUniqueID);
   fChain->SetBranchAddress("MissingET.fBits", MissingET_fBits, &b_MissingET_fBits);
   fChain->SetBranchAddress("MissingET.MET", MissingET_MET, &b_MissingET_MET);
   fChain->SetBranchAddress("MissingET.Eta", MissingET_Eta, &b_MissingET_Eta);
   fChain->SetBranchAddress("MissingET.Phi", MissingET_Phi, &b_MissingET_Phi);
   fChain->SetBranchAddress("MissingET_size", &MissingET_size, &b_MissingET_size);
   fChain->SetBranchAddress("ScalarHT", &ScalarHT_, &b_ScalarHT_);
   fChain->SetBranchAddress("ScalarHT.fUniqueID", ScalarHT_fUniqueID, &b_ScalarHT_fUniqueID);
   fChain->SetBranchAddress("ScalarHT.fBits", ScalarHT_fBits, &b_ScalarHT_fBits);
   fChain->SetBranchAddress("ScalarHT.HT", ScalarHT_HT, &b_ScalarHT_HT);
   fChain->SetBranchAddress("ScalarHT_size", &ScalarHT_size, &b_ScalarHT_size);
   Notify();
}

Bool_t Delphes::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Delphes::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Delphes::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Delphes_cxx
