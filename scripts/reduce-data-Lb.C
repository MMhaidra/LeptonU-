
#include "TFile.h"
#include "TCut.h"
#include "TTree.h"
#include <iostream>

// What do we do : this is a little script to reduce the data for MVA input
// and count the luminosity just because we can. 

void FirstReduction(){
  bool Jpsi = false; 
  bool Photon = true;
  bool Penguin = false;
  //  TString inputfile = ""
  TFile * f = new TFile ("/exp/LHCb/amhis/LeptonU-/tuples/data/LeptonU-data-electrons-24032016.root");
  TTree * tree = (TTree*)f->Get("TupleFromData_Bu2LLK_eeLine2/DecayTree");
  
  cout << "Total number of entries : "  << tree->GetEntries () << endl;
  
  
  // Mass windows
  TCut UpperSideBand = "Lambdab> 5700";
  TCut JpsiW = "Jpsi_MM > 2200 && Jpsi_MM < 3400 "; 
  TCut PhotonW  = "Jpsi_MM < 20"; 
  TCut LambdastarW = "Lambdastar_MM< 2000";
  
  //PID cuts
  TCut ElectronID = "L1_PIDe>1.2&& L1_PIDe>1.2";
  TCut ProtonID = "Proton_ProbNNp > 0.1";				
  TCut KaonID = "Kaon_ProbNNk> 0.1";
 

  //Apply Trigger cuts

  //L0
  TCut L0Electron=  "(Lambdab_L0ElectronDecision_TOS) ";
  TCut L0Hadron = " (Lambdab_L0HadronDecision_TOS && !Lambdab_L0ElectronDecision_TOS) ";
  TCut L0TIS = "(Lambdab_L0Global_TIS && !(Lambdab_L0ElectronDecision_TOS || Lambdab_L0HadronDecision_TOS))";
  TCut L0 = L0Electron || L0Hadron || L0TIS;
  
  //HLt1
  TCut Hlt1 = "(Lambdab_Hlt1TrackAllL0Decision_TOS)";
  
  
  //Hlt2 
  TCut Hlt2 = "(Lambdab_Hlt2Topo2BodyBBDTDecision_TOS ||  Lambdab_Hlt2Topo4BodyBBDTDecision_TOS || Lambdab_Hlt2TopoE2BodyBBDTDecision_TOS  ||   Lambdab_Hlt2TopoE4BodyBBDTDecision_TOS )";
  
  
  TCut  Trigger = L0 && Hlt1 && Hlt2 ;
  
  // cout <<   LambdastarW  << "  " << ElectronID << " " << ProtonID << " " << KaonID << "  " << L0Electron  << "  " << L0Hadron << " " << L0TIS <<  " "  << Hlt1 << " " << Hlt2 << endl;
    TCut ApplyThis =  LambdastarW && ElectronID && ProtonID && KaonID && L0 ;

    
  
  if (Jpsi == true) {ApplyThis = ApplyThis&&JpsiW;}
  if (Photon  == true){ApplyThis = ApplyThis &&PhotonW;}
 


  //Apply this to see something around the J/psi
  // root [10] tree->Draw("B_M0123_Subst3_pi2p", "Pion_ProbNNp>0.4&& Kaon_ProbNNk >0.4 && B_VtxChi2_0123 < 9 && B_IPChi2_0123<9")
  cout <<"For now we use the following cuts : " << ApplyThis << endl;
  TFile * signal_file = new TFile("/exp/LHCb/amhis/LeptonU-/tuples/data/forfit/real-data-electrons-jpsi-24032016.root", "recreate");
  // TFile * signal_file = new TFile("/exp/LHCb/amhis/LeptonU-/tuples/data/forfit/real-data-electrons-jpsi-uppersideband.root", "recreate");
  
  
  TTree * signal_copy = tree->CopyTree( ApplyThis);
  signal_copy->Write();
  cout << "We are done with the first reduction :-)  and now we have : " << signal_copy -> GetEntries()  << " ie : a reduction of : "  <<    ((float)signal_copy -> GetEntries()  / (float)tree -> GetEntries() )* 100 << " % "<< endl;
  return;
}


void SecondReduction(){
    
  TFile * f = new TFile ("/exp/LHCb/amhis/LeptonU-/tuples/data/forfit/real-data-electrons-photon-24032016.root");
  TTree * tree = (TTree*)f->Get("DecayTree");
  
  cout << "There are  " << tree->GetEntries()   << " entries in the tree after the first reduction " << endl;
  //  TCut TighterCuts = "Pion_ProbNNp>0.4&& Kaon_ProbNNk >0.4 && B_VtxChi2_0123 < 9 && B_IPChi2_0123<9";
  TCut TighterCuts = "Proton_ProbNNp>0.1&& Kaon_ProbNNk >0.1 && Lambdab_ENDVERTEX_CHI2< 9 && Lambdab_LOKI_DTF_CHI2NDOF < 9 ";
  
  TFile * Tighterfile = new TFile ("/exp/LHCb/amhis/LeptonU-/tuples/data/forfit/real-data-electrons-photon-tight-24032016.root", "recreate"); 
  //  TCut TighterCuts = "Pion_ProbNNp>0.4&& Kaon_ProbNNk >0.4 && B_VtxChi2_0123 < 9 && B_IPChi2_0123<9";
  TTree * Tightertree = tree->CopyTree(TighterCuts);



  Tightertree->Write();
  cout << "Now we also applied this : " << TighterCuts  << endl;
  cout << "We are done with the second reduction yippie ! " << endl;


  return;
}
