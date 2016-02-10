#include "TFile.h"
#include "TCut.h"
#include "TTree.h"
#include <iostream>

void tellMeTheLumi (){


  TFile * f = new TFile ("/exp/LHCb/amhis/LeptonU-/tuples/data/LeptonU-total-electrons-11122015.root");
  TTree* lumi =  (TTree*)f->Get("GetIntegratedLuminosity/LumiTuple");
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
