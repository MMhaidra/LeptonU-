// this is a little script to reduce the data for MVA input
// and count the luminosity just because we can 
#include "TFile.h"
#include "TCut.h"
#include "TTree.h"
#include <iostream>


void reduce(){




  bool Jpsi = true; 
  bool Photon = false;
  bool Penguin = false;
  TFile * f = new TFile ("/exp/LHCb/amhis/LeptonU/tuples/data/LeptonU-total-electrons-11122015.root");
  TTree * tree = (TTree*)f->Get("TupleFromDataKst_Bu2LLK_eeLine2/DecayTree");
  
  cout << "Total number of entries : "  << tree->GetEntries () << endl;
  
  
  // Mass windows
  TCut UpperSideBand = "B_M0123_Subst3_pi2p";
  TCut JpsiW = "Jpsi_MM > 2200 && Jpsi_MM < 3400 "; 
  TCut PhotonW  = "Jpsi_MM < 20"; 
  TCut LambdastarW = "B_M23_Subst3_pi2p< 2000";
  
  //PID cuts
  TCut ElectronID = "L1_PIDe>1.2&& L1_PIDe>1.2";
  TCut ProtonID = "Pion_ProbNNp > 0.1";					
  TCut KaonID = "Kaon_ProbNNk> 0.1";
 

  //Apply Trigger cuts
  //TCut L0 = 
  //TCut signal_selection = "J_psi_1S_MM > 2800. && J_psi_1S_MM < 3400. &&fabs(Lambda_1520_0_MM-1520)<150&&fabs(Lambda_b0_MM-5620)<1500&& Kplus_ProbNNk>0.4&&p~minus_ProbNNp>0.4&&eplus_PIDe>1.2&&eminus_PIDe>1.2";
  
  TCut ApplyThis = LambdastarW && ElectronID && ProtonID && KaonID;
  
  if (Jpsi == true) {ApplyThis = ApplyThis&&JpsiW;}
  if (Photon  == true){ApplyThis = ApplyThis &&PhotonW;}
  
  //Apply this to see something around the J/psi
  // root [10] tree->Draw("B_M0123_Subst3_pi2p", "Pion_ProbNNp>0.4&& Kaon_ProbNNk >0.4 && B_VtxChi2_0123 < 9 && B_IPChi2_0123<9")
  
  cout <<"For now we use the following cuts : " << ApplyThis << endl;
  TFile * signal_file = new TFile("/exp/LHCb/amhis/LeptonU/tuples/data/forfit/real-data-electrons-jpsi.root", "recreate");
  TTree * signal_copy = tree->CopyTree( ApplyThis);
  signal_copy->Write();
  
  
  
  
  
  TTree* lumi =  (TTree*)f->Get("GetIntegratedLuminosity/LumiTuple");
  cout << lumi->GetEntries()  << endl;
  double IntegratedLuminosity = 0;
  lumi->SetBranchAddress("IntegratedLuminosity", &IntegratedLuminosity);
  double total =0 ;
  for (int i = 0 ; i < lumi->GetEntries(); i++){
    lumi->GetEntry(i);
    
    total = total + IntegratedLuminosity;
    
  }
  cout << "-----------------" << endl;
  cout << "Total luminosity " << total/100. << " fb-1 "<<  endl;
  cout << "-----------------" << endl;
  
  return;
  
}
