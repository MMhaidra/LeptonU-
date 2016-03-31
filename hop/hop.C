#include<iostream>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TH1F.h"

using namespace std;



void AddHop (){
  cout << "Now we add Hop ! " << endl;
  //TFile * signal_file = new TFile("LeptonU-data-electrons.root", "UPDATE");
  TFile *signal_file = new TFile ("/exp/LHCb/amhis/LeptonU-/tuples/montecarlo/spring16/mars/mc-15154001.root", "UPDATE");
 
  TTree * signalcopy = (TTree*)signal_file->Get("Tuple_Bu2LLK_eeLine2/DecayTree");

  cout << "Number of events in the small tree : " << signalcopy->GetEntries() << endl;
  

ULong64_t eventNumber;


Double_t Lambdab_ENDVERTEX_X ;
Double_t Lambdab_ENDVERTEX_Y ;
Double_t Lambdab_ENDVERTEX_Z ;

Double_t Lambdab_OWNPV_X ;
Double_t Lambdab_OWNPV_Y ;
Double_t Lambdab_OWNPV_Z ;

Double_t Lambdab_FD_OWNPV ;

Double_t Lambdastar_PX ;
Double_t Lambdastar_PY ;
Double_t Lambdastar_PZ ;
Double_t Lambdastar_P ;
Double_t Lambdastar_PE ;

Double_t Jpsi_PX ;
Double_t Jpsi_PY ;
Double_t Jpsi_PZ ;
Double_t Jpsi_P ;

Double_t L1_PX ;
Double_t L1_PY ;
Double_t L1_PZ ;
Double_t L1_P ;

Double_t L2_PX ;
Double_t L2_PY ;
Double_t L2_PZ ;
Double_t L2_P ;



Double_t Proton_PX ;
Double_t Proton_PY ;
Double_t Proton_PZ ;
Double_t Proton_P ;

Double_t Proton_PE ;


Double_t Kaon_PX ;
Double_t Kaon_PY ;
Double_t Kaon_PZ ;
Double_t Kaon_P ;
Double_t Kaon_PE ;


signalcopy->SetBranchAddress("Lambdab_ENDVERTEX_X" ,&Lambdab_ENDVERTEX_X);
signalcopy->SetBranchAddress("Lambdab_ENDVERTEX_Y" ,&Lambdab_ENDVERTEX_Y);
signalcopy->SetBranchAddress("Lambdab_ENDVERTEX_Z" ,&Lambdab_ENDVERTEX_Z);

signalcopy->SetBranchAddress("Lambdab_OWNPV_X" ,&Lambdab_OWNPV_X);
signalcopy->SetBranchAddress("Lambdab_OWNPV_Y" ,&Lambdab_OWNPV_Y);
signalcopy->SetBranchAddress("Lambdab_OWNPV_Z" ,&Lambdab_OWNPV_Z);

signalcopy->SetBranchAddress("Lambdab_FD_OWNPV" ,&Lambdab_FD_OWNPV);

signalcopy->SetBranchAddress("Lambdastar_PX" ,&Lambdastar_PX) ;
signalcopy->SetBranchAddress("Lambdastar_PY" ,&Lambdastar_PY) ;
signalcopy->SetBranchAddress("Lambdastar_PZ" ,&Lambdastar_PZ) ;
signalcopy->SetBranchAddress("Lambdastar_P" ,&Lambdastar_P) ;
signalcopy->SetBranchAddress("Lambdastar_PE" ,&Lambdastar_PE) ;

signalcopy->SetBranchAddress("Jpsi_PX" ,&Jpsi_PX) ;
signalcopy->SetBranchAddress("Jpsi_PY" ,&Jpsi_PY) ;
signalcopy->SetBranchAddress("Jpsi_PZ" ,&Jpsi_PZ) ;
signalcopy->SetBranchAddress("Jpsi_P" ,&Jpsi_P) ;


signalcopy->SetBranchAddress("L1_PX" ,&L1_PX) ;
signalcopy->SetBranchAddress("L1_PY" ,&L1_PY) ;
signalcopy->SetBranchAddress("L1_PZ" ,&L1_PZ) ;
signalcopy->SetBranchAddress("L1_P" ,&L1_P) ;

signalcopy->SetBranchAddress("L2_PX" ,&L2_PX) ;
signalcopy->SetBranchAddress("L2_PY" ,&L2_PY) ;
signalcopy->SetBranchAddress("L2_PZ" ,&L2_PZ) ;
signalcopy->SetBranchAddress("L2_P" ,&L2_P) ;


signalcopy->SetBranchAddress("Kaon_PX" ,&Kaon_PX) ;
signalcopy->SetBranchAddress("Kaon_PY" ,&Kaon_PY) ;
signalcopy->SetBranchAddress("Kaon_PZ" ,&Kaon_PZ) ;
signalcopy->SetBranchAddress("Kaon_P" , &Kaon_P) ;
signalcopy->SetBranchAddress("Kaon_PE" ,&Kaon_PE) ;



signalcopy->SetBranchAddress("Proton_PX" ,&Proton_PX) ;
signalcopy->SetBranchAddress("Proton_PY" ,&Proton_PY) ;
signalcopy->SetBranchAddress("Proton_PZ" ,&Proton_PZ) ;
signalcopy->SetBranchAddress("Proton_P" ,&Proton_P) ;
signalcopy->SetBranchAddress("Proton_PE" ,&Proton_PE) ;






Int_t signumberOfEntries = signalcopy->GetEntries();

const Double_t PDG_e_M = 0.510998910 ;
Float_t Jpsi_SF ;
TLorentzVector KstarMom ;
TLorentzVector ScaledL1Mom ;
TLorentzVector ScaledL2Mom ;
TLorentzVector ScaledJpsiMom ;
TLorentzVector ScaledBMom ;


Float_t val =0.0;
Float_t valMee =0.0;
Float_t HopVal =0.0;

TBranch *newbranch1 = signalcopy->Branch("HOP_Lambdab_MM", &val, "HOP_Lambdab_MM");
TBranch *newbranch2 = signalcopy->Branch("HOP", &HopVal, "HOP") ;
TBranch *newbranch3 = signalcopy->Branch("HOP_Jpsi_MM", &valMee, "HOP_Jpsi_MM");
//TBranch *newbranch4 = signalcopy->Branch("phi_InvMass",&phi_InvMass ,"phi_InvMass");

for (Int_t loopie=0; loopie < signumberOfEntries; ++loopie){
  signalcopy->GetEntry(loopie);


  // compute the pT of the K* wrt to the B line of flight

  Double_t cosThetaKstar = ((Lambdab_ENDVERTEX_X-Lambdab_OWNPV_X)*Lambdastar_PX+(Lambdab_ENDVERTEX_Y-Lambdab_OWNPV_Y)*Lambdastar_PY+(Lambdab_ENDVERTEX_Z-Lambdab_OWNPV_Z)*Lambdastar_PZ)/(Lambdastar_P*Lambdab_FD_OWNPV) ;
  Double_t pTKstar = Lambdastar_P*TMath::Sqrt(1.-cosThetaKstar*cosThetaKstar) ;
  Double_t cosThetaY = ((Lambdab_ENDVERTEX_X-Lambdab_OWNPV_X)*Jpsi_PX+(Lambdab_ENDVERTEX_Y-Lambdab_OWNPV_Y)*Jpsi_PY+ (Lambdab_ENDVERTEX_Z-Lambdab_OWNPV_Z)*Jpsi_PZ)/(Jpsi_P*Lambdab_FD_OWNPV);
  Double_t pTY = Jpsi_P*TMath::Sqrt(1.-cosThetaY*cosThetaY) ;

  Jpsi_SF = pTKstar/pTY;

  HopVal = Jpsi_SF ;

  ScaledL1Mom.SetXYZM( Jpsi_SF*L1_PX, Jpsi_SF*L1_PY, Jpsi_SF*L1_PZ, PDG_e_M );
  ScaledL2Mom.SetXYZM( Jpsi_SF*L2_PX, Jpsi_SF*L2_PY, Jpsi_SF*L2_PZ, PDG_e_M );
  KstarMom.SetPxPyPzE( Lambdastar_PX, Lambdastar_PY, Lambdastar_PZ, Lambdastar_PE );
  ScaledBMom = ScaledL1Mom + ScaledL2Mom + KstarMom  ;
  ScaledJpsiMom = ScaledL1Mom + ScaledL2Mom ;
  val = ScaledBMom.M();
  valMee = ScaledJpsiMom.M();
  newbranch1->Fill();
  newbranch2->Fill();
  newbranch3->Fill();
  //newbranch4->Fill();
}
signal_file->Write();







  return;


}

