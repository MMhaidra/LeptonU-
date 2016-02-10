
//standard libs
#include<iostream>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TH1F.h"
//RooFit librairies




using namespace std;
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCBShape.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooPlot.h"

using namespace RooFit ;


void extractR(){
  //------------------------------------------------
  TFile * f = new TFile ("/exp/LHCb/amhis/LeptonU-/tuples/data/forfit/real-data-electrons-jpsi-tight.root");
  TTree * signaltree = (TTree*)f->Get("DecayTree");
  
  //declare the variables to be fitted
  
  RooRealVar mass("B_M0123_Subst3_pi2p", "J/#psi (e^{+}e^{-}) p K ", 4800., 6500., "MeV/c^{2}" );
  
  //RooRealVar mass("Lambda_b0_LOKI_MASS_JpsiConstr", "J/#psi (e^{+}e^{-}) p K ", 5000., 6200., "MeV/c^{2}" );
  
  //RooDataSet *data  = RooDataSet::read("data.txt",  RooArgList(mass));
  RooDataSet *data_electrons  =  new RooDataSet("data_electrons", "data_electrons", signaltree,  RooArgList(mass));
  RooDataSet *data_muons  =  new RooDataSet("data_muons", "data_muons", signaltree,  RooArgList(mass));
  
  //signal variables
  RooRealVar meanMass_S ("meanMass_S","Mean of the mass ",5679.5 ,5000,6200);
  RooRealVar SigmaMassCore_S("sigmaMassCore_S","Sigma Core of the mass ",50,10,300);
  RooRealVar SigmaMassTail_S("sigmaMassTail_S","Sigma Core of the mass ",80,10,400);
  
  RooRealVar Alpha_S ("Alpha_S", "Alpha_S",2,0,100);
  RooRealVar nCB_S("nCB_S","nCB_S",0, 0, 10);
  //background variables
  RooRealVar slope_Mass_B("slope_Mass_B", "slope_Mass_B", 0, -100, 100);
  //number of signal and background events to be fitted
  RooRealVar N_Sig ("N_Sig","# of signal events",100,0,3000);
  RooRealVar N_Bkg("N_Bkg","# of background events",100,0,3000);
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
 
 
 RooAddPdf totalpdf("totalpdf","totalpdf", RooArgList(GaussMass_pdf_S,exp_pdf_B), RooArgList(N_Sig,N_Bkg) );
 //RooAddPdf totalpdf("totalpdf","totalpdf", RooArgList(CB_pdf_S, exp_pdf_B), RooArgList(N_Sig,N_Bkg) );
 //RooAddPdf totalpdf("totalpdf","totalpdf", RooArgList(signal_tot, exp_pdf_B), RooArgList(N_Sig,N_Bkg) );
 
 //-------------------------------------------------
 //fitting
 totalpdf.fitTo(*data_electrons,Extended());
 //plotting
 TCanvas * mass_canvas = new TCanvas("mass", "mass", 600,600);
 mass_canvas->cd();
 
 RooPlot* mass_plot = mass.frame(Title(""));
 mass_plot->GetYaxis()->SetTitleOffset (1.5);
 data_electrons->plotOn(mass_plot,Binning(50));
 
 mass.setRange("signal_region", 5500., 5700.); //We define here the signal region.
 RooAbsReal* igx_sig = exp_pdf_B.createIntegral(mass,NormSet(mass),Range("signal_region")) ;
 
 cout << "----------" << endl;
 cout << "Fraction of background events in the signal region =   " << igx_sig->getVal() << endl ;
 cout << "Number of signal events in the signal region " << igx_sig->getVal()*N_Bkg.getVal() << endl;
 cout << "-------" << endl;
 

 totalpdf.plotOn(mass_plot);
 totalpdf.plotOn(mass_plot, Components(CB_pdf_S), LineColor(801));
 //totalpdf.plotOn(mass_plot, Components(GaussMass_pdf_S), LineColor(801));
 //totalpdf.plotOn(mass_plot, Components(signal_tot), LineColor(613));
 totalpdf.plotOn(mass_plot,Components(exp_pdf_B),LineStyle(kDashed), LineColor(kBlue));
 totalpdf.plotOn(mass_plot);
 
 mass_plot->Draw();
 mass_plot->SetTitle("");
 
 
 return ;
};
