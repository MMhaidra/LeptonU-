
#include "TFile.h"
#include "TCut.h"
#include "TTree.h"
#include <iostream>

// What do we do : this is a little script to reduce the data for MVA input
// and count the luminosity just because we can. 

void FirstReduction(){
  bool Jpsi = true; 
  bool Photon = false;
  bool Penguin = false;
  TFile * f = new TFile ("/exp/LHCb/amhis/LeptonU-/tuples/data/LeptonU-total-electrons-23022016.root");
  TTree * tree = (TTree*)f->Get("TupleFromDataKst_Bu2LLK_eeLine2/DecayTree");
  
  cout << "Total number of entries : "  << tree->GetEntries () << endl;
  
  
  // Mass windows
  TCut UpperSideBand = "B_M0123_Subst3_pi2p > 5700";
  TCut JpsiW = "Jpsi_MM > 2200 && Jpsi_MM < 3400 "; 
  TCut PhotonW  = "Jpsi_MM < 20"; 
  TCut LambdastarW = "B_M23_Subst3_pi2p< 2000";
  
  //PID cuts
  TCut ElectronID = "L1_PIDe>1.2&& L1_PIDe>1.2";
  TCut ProtonID = "Pion_ProbNNp > 0.1";				
  TCut KaonID = "Kaon_ProbNNk> 0.1";
 

  //Apply Trigger cuts
   
  //  TCut ApplyThis = LambdastarW && ElectronID && ProtonID && KaonID && UpperSideBand;
   TCut ApplyThis =  LambdastarW && ElectronID && ProtonID && KaonID ;

  
  if (Jpsi == true) {ApplyThis = ApplyThis&&JpsiW;}
  if (Photon  == true){ApplyThis = ApplyThis &&PhotonW;}
  //Apply this to see something around the J/psi
  // root [10] tree->Draw("B_M0123_Subst3_pi2p", "Pion_ProbNNp>0.4&& Kaon_ProbNNk >0.4 && B_VtxChi2_0123 < 9 && B_IPChi2_0123<9")
  cout <<"For now we use the following cuts : " << ApplyThis << endl;
    TFile * signal_file = new TFile("/exp/LHCb/amhis/LeptonU-/tuples/data/forfit/real-data-electrons-jpsi-23022016.root", "recreate");
  // TFile * signal_file = new TFile("/exp/LHCb/amhis/LeptonU-/tuples/data/forfit/real-data-electrons-jpsi-uppersideband.root", "recreate");
  

  TTree * signal_copy = tree->CopyTree( ApplyThis);
  signal_copy->Write();
  cout << "We are done with the first reduction :-) " << endl;
  return;
}


void SecondReduction(){
    
  TFile * f = new TFile ("/exp/LHCb/amhis/LeptonU-/tuples/data/forfit/real-data-electrons-jpsi-23022016.root");
  TTree * tree = (TTree*)f->Get("DecayTree");
  

  TCut TighterCuts = "Pion_ProbNNp>0.4&& Kaon_ProbNNk >0.4 && B_VtxChi2_0123 < 9 && B_IPChi2_0123<9";
  
  TFile * Tighterfile = new TFile ("/exp/LHCb/amhis/LeptonU-/tuples/data/forfit/real-data-electrons-jpsi-tight-23022016.root", "recreate"); 
  TTree * Tightertree = tree->CopyTree(TighterCuts);



  Tightertree->Write();
  cout << "Now we also applied this : " << TighterCuts  << endl;
  cout << "We are done with the second reduction yippie ! " << endl;


  return;
}
