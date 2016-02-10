
#include "TFile.h"
#include "TCut.h"
#include "TTree.h"
#include <iostream>


void efficiencies(){
  bool Jpsi = true; 
  bool Photon = false;
  bool Penguin = false;
  TFile * f = new TFile ("/exp/LHCb/amhis/LeptonU-/tuples/montecarlo/spring16/mc-15454101-leptonU.root");
  TTree * tree = (TTree*)f->Get("Tuple_Bu2LLK_eeLine2/DecayTree");
  
  cout << "Total number of entries : "  << tree->GetEntries () << endl;
  
 
  return;

}
