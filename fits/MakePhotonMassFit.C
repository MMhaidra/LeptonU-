
// April 2016
// What do I do : 
// The aim of this code is to fit  Lb --> pK gamma
// Voila



#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

//standard libs
#include<iostream>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"




//RooFit librairies




using namespace std;


using namespace RooFit ;

#include "RooBinning.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCBShape.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooPlot.h"



void ApplyMoreCutsBeforeFit(){
  //------------------------------------------------
    TFile * f = new TFile ("/exp/LHCb/amhis/LeptonU-/tuples/data/forfit/real-data-electrons-photon-30032016-reduced.root");
  //TFile * f = new TFile ("/exp/LHCb/amhis/LeptonU-/tuples/data/forfit/real-data-electrons-jpsi-30032016-reduced.root");
  TTree * signaltree = (TTree*)f->Get("DecayTree");
  TCut PID = "Kaon_MC12TuneV3_ProbNNk > 0.1 && Proton_MC12TuneV3_ProbNNp>0.1&&L2_PIDe>1.2&&L1_PIDe>1.2";
  
  
  // Cuts taken from  pk mumu analysis
  //  TCut HarderPID= "Proton_MC12TuneV3_ProbNNp >0.2 &&Proton_MC12TuneV3_ProbNNk<0.8 && Proton_MC12TuneV3_ProbNNpi<0.7 && Kaon_MC12TuneV3_ProbNNk>0.2 && Kaon_MC12TuneV3_ProbNNp < 0.8 &&L2_PIDe>1.2&&L1_PIDe>1.2" ;
  

  TCut HarderPID= "Proton_ProbNNp >0.2 &&Proton_ProbNNk<0.8 && Proton_ProbNNpi<0.7 && Kaon_ProbNNk>0.2 && Kaon_ProbNNp < 0.8 &&L2_PIDe>1.2&&L1_PIDe>1.2" ;
  
  TCut Kinematics = "Proton_PT >500 && Lambdab_DTF_PV_chi2<10  && Lambdab_PT>2000 ";

  TCut Lstar = "Lambdastar_M<1550 && Lambdastar_M>1450";
 
  //  TCut HopCut = "Lambdab_MM> (4000 + 344 +  100*log(Lambdab_FDCHI2_OWNPV)) ";

  TCut HopCut = "Lambdab_MM> (4000 + 344 +  75*log(Lambdab_FDCHI2_OWNPV)) ";

  TCut RemovePhi = "Lambdab_M01_Subst0_p2K>1050";
  //  TCut RemoveKst = "Lambdab_M01_Subst0_p2pi< 850 && Lambdab_M01_Subst0_p2pi>920";
  TCut RemoveKst = "Lambdab_M01_Subst0_p2pi< 850 ";
  
  //  TCut ApplyTheseCuts = HarderPID && Kinematics&& HopCut && RemovePhi && RemoveKst;
  TCut ApplyTheseCuts = HarderPID && Kinematics ;
  

  TFile * localFileForFit = new TFile("localFileForFit.root", "RECREATE"); 
  cout <<" These are the cuts that we applied : " <<  ApplyTheseCuts << endl;
  TTree * smallerTreeForFit = signaltree->CopyTree(ApplyTheseCuts);
  smallerTreeForFit->Write();
  f->Close();
  return;
}
  



void MakeQuickPlots (){

  TFile * localFileForFit = new TFile("localFileForFit-photon.root");
  TTree * smallerTreeForFit = (TTree*)localFileForFit->Get("DecayTree");
  

  Double_t Lambdab_M =0;
  smallerTreeForFit->SetBranchAddress("Lambdab_M", &Lambdab_M);

  Double_t Lambdab_LOKI_MASS_JpsiConstr =0;
  smallerTreeForFit->SetBranchAddress("Lambdab_LOKI_MASS_JpsiConstr", &Lambdab_LOKI_MASS_JpsiConstr);
 

  Double_t Lambdab_M01 =0; 
  smallerTreeForFit->SetBranchAddress("Lambdab_M01", &Lambdab_M01);
    
  Double_t Lambdab_M01_Subst0_p2pi=0;
    smallerTreeForFit->SetBranchAddress("Lambdab_M01_Subst0_p2pi", &Lambdab_M01_Subst0_p2pi);

  Double_t Lambdab_M01_Subst0_p2K =0;
  smallerTreeForFit->SetBranchAddress("Lambdab_M01_Subst0_p2K", &Lambdab_M01_Subst0_p2K);
    

  Double_t Lambdab_M0123_Subst0_p2K =0;
  smallerTreeForFit->SetBranchAddress("Lambdab_M0123_Subst0_p2K", &Lambdab_M0123_Subst0_p2K);

  TH1F * LbMHisto = new TH1F ("LbMHisto", "LbMHisto", 100, 4300, 7000);
  TH1F * LbMHisto_red= new TH1F ("LbMHisto_red", "LbMHisto_red", 100, 4300, 7000);



 TH1F * LbMConstHisto = new TH1F ("LbMConstHisto", "LbMConstHisto", 100, 4300, 7000);
  TH1F * LbMConstHisto_red = new TH1F ("LbMConstHisto_red", "LbMConstHisto_red", 100, 4300, 7000);
  TH1F * LbMConstHisto_redphi = new TH1F ("LbMConstHisto_redphi", "LbMConstHisto_redphi", 100, 4300, 7000);

  TH2F * LbMConstVspK = new TH2F("LbMConstVspK", "LbMconstVspK", 100, 4300, 7000, 100, 1400, 2100);
  TH2F * LbMConstVsKK = new TH2F("LbMConstVsKK", "LbMconstVsKK", 100, 4300, 7000, 100, 1000, 1500);
  TH2F * LbMConstVsKpi = new TH2F("LbMConstVsKpi", "LbMconstVsKK", 100, 4300, 7000, 100, 600, 1200);
  
  
  
  
  for (int i = 0; i < smallerTreeForFit->GetEntries(); i++){
    smallerTreeForFit->GetEntry(i);
    LbMHisto->Fill(Lambdab_M);

    
    LbMConstHisto->Fill(Lambdab_LOKI_MASS_JpsiConstr);
    




    // LbMConstVspK->Fill(Lambdab_LOKI_MASS_JpsiConstr, Lambdab_M01);
    // LbMConstVsKK->Fill(Lambdab_LOKI_MASS_JpsiConstr, Lambdab_M01_Subst0_p2K);
    // LbMConstVsKpi->Fill(Lambdab_LOKI_MASS_JpsiConstr, Lambdab_M01_Subst0_p2pi);
    
    LbMConstVspK->Fill(Lambdab_M, Lambdab_M01);
    LbMConstVsKK->Fill(Lambdab_M, Lambdab_M01_Subst0_p2K);
    LbMConstVsKpi->Fill(Lambdab_M, Lambdab_M01_Subst0_p2pi);
    

    if (Lambdab_M01_Subst0_p2K > 1030) { 
      LbMHisto_red->Fill(Lambdab_M);
      LbMConstHisto_red->Fill(Lambdab_LOKI_MASS_JpsiConstr);
    }
    
  }

  TCanvas * c0 = new TCanvas ("c0", "c0", 600,400);
  c0->Divide(2,1);
  c0->cd(1);

  LbMHisto->Draw();
  LbMHisto_red->Draw("SAME");
  LbMHisto_red->SetLineColor(2);
  c0->cd(2);
  LbMConstHisto->Draw();
  LbMConstHisto_red->Draw("SAME");
  LbMConstHisto_red->SetLineColor(2);



  TCanvas *c1 = new TCanvas("c1", "c1", 700,400);
  c1->Divide(3,1);
  c1->cd(1);
  LbMConstVspK->Draw("ZCOL");
  c1->cd(2);
  LbMConstVsKK->Draw("ZCOL");
  c1->cd(3);
  LbMConstVsKpi->Draw("ZCOL");
  
  return;
}


void MakePhotonMassFit (){

  TFile * localFileForFit = new TFile("localFileForFit.root");
  TTree * smallerTreeForFit = (TTree*)localFileForFit->Get("DecayTree");
  //declare the variables to be fitted

  //RooRealVar mass("Lambdab_M", "#gamma p K ", 4600., 7000., "MeV/c^{2}" );
  RooRealVar mass("Lambdab_LOKI_MASS_JpsiConstr", "p K J/#psi ", 4600., 7000., "MeV/c^{2}" );
  RooRealVar pKMass ("Lambdastar", "Lambdastar", 1450,1550, "MeV/c^{2}");
  
  //RooDataSet *data  = RooDataSet::read("data.txt",  RooArgList(mass));
  RooDataSet *data_electrons  =  new RooDataSet("data_electrons", "data_electrons", smallerTreeForFit,  RooArgList(mass, pKMass));
  RooDataSet *data_muons  =  new RooDataSet("data_muons", "data_muons", smallerTreeForFit,  RooArgList(mass));
  
  //signal variables
  
  RooRealVar meanMass_S ("meanMass_S","Mean of the mass ",5679.5 ,5000,6200);
  RooRealVar SigmaMassCore_S("sigmaMassCore_S","Sigma Core of the mass ",50,10,300);
  RooRealVar SigmaMassTail_S("sigmaMassTail_S","Sigma Core of the mass ",80,10,400);
  
  RooRealVar Alpha_S ("Alpha_S", "Alpha_S",0,0,100);
  RooRealVar nCB_S("nCB_S","nCB_S",3, 0, 100);
  //background variables
  RooRealVar slope_Mass_B("slope_Mass_B", "slope_Mass_B", 0, -100, 100);
  //number of signal and background events to be fitted
  RooRealVar N_Sig ("N_Sig","# of signal events",100,0,30000);
  RooRealVar N_Bkg("N_Bkg","# of background events",100,0,30000);
  //-----------------------------------------------------
  //signal pdf
 RooGaussian GaussMass_pdf_S ("GaussMass_pdf_S","GaussMass S", mass , meanMass_S , SigmaMassTail_S  );
 RooCBShape CB_pdf_S ("CB_pdf_S", "CB_pdf_S", mass , meanMass_S , SigmaMassCore_S, Alpha_S, nCB_S  );
 RooRealVar frac ("frac", "frac", 0.5, 0,1);
 RooAddPdf signal_tot ( "signal_tot", "signal_tot", RooArgList(CB_pdf_S, GaussMass_pdf_S), RooArgSet(frac));
 //background pdf
 //RooPolynomial flatComb_pdf_B("flatComb_pdf_B","flatComb",bMass_var, RooArgList(slope_Mass_B) ) ;
 RooChebychev flatComb_pdf_B("flatComb_pdf_B","flatComb",mass, RooArgList(slope_Mass_B) ) ;
 RooExponential exp_pdf_B("exp_pdf_B", "exp_pdf_B", mass, slope_Mass_B);
 //total pdf
 
 
 // RooAddPdf totalpdf("totalpdf","totalpdf", RooArgList(GaussMass_pdf_S,exp_pdf_B), RooArgList(N_Sig,N_Bkg) );
 RooAddPdf totalpdf("totalpdf","totalpdf", RooArgList(CB_pdf_S, exp_pdf_B), RooArgList(N_Sig,N_Bkg) );
 //RooAddPdf totalpdf("totalpdf","totalpdf", RooArgList(signal_tot, exp_pdf_B), RooArgList(N_Sig,N_Bkg) );
 
 //-------------------------------------------------
 //fitting
 // totalpdf.fitTo(*data_electrons , Extended());
 

 // RooBinning massBins(50);
 ///
 totalpdf.fitTo(*data_electrons, "e");
 //plotting;
 TCanvas * mass_canvas = new TCanvas("mass", "mass", 600,600);
 mass_canvas->cd();
 
 // RooPlot* mass_plot = mass.frame(Title(""));
 RooPlot* mass_plot = mass.frame();

  mass_plot->GetYaxis()->SetTitleOffset (1.5);
  data_electrons->plotOn(mass_plot ,Binning(50));

  
  //  data_electrons->plotOn(mass_plot);
 
 // mass.setRange("signal_region", 5500., 5700.); //We define here the signal region.
 // RooAbsReal* igx_sig = exp_pdf_B.createIntegral(mass,NormSet(mass),Range("signal_region")) ;
 
 // cout << "----------" << endl;
 // cout << "Fraction of background events in the signal region =   " << igx_sig->getVal() << endl ;
 // cout << "Number of signal events in the signal region " << igx_sig->getVal()*N_Bkg.getVal() << endl;
 // cout << "-------" << endl;
 

 totalpdf.plotOn(mass_plot);

 //totalpdf.plotOn(mass_plot, Components(GaussMass_pdf_S), LineColor(801));
 //totalpdf.plotOn(mass_plot, Components(signal_tot), LineColor(613));
 totalpdf.plotOn(mass_plot,Components(exp_pdf_B),  DrawOption("F"), FillColor(921));
 // totalpdf.plotOn(mass_plot,Components(exp_pdf_B),LineStyle(kDashed), LineColor(kBlue));
 totalpdf.plotOn(mass_plot, Components(CB_pdf_S), LineColor(882));

 totalpdf.plotOn(mass_plot);
   data_electrons->plotOn(mass_plot ,Binning(50));

 mass_plot->Draw();
 mass_plot->SetTitle("");
 
 
 return ;
};
