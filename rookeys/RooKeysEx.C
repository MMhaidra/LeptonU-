
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "RooPlot.h"
using namespace RooFit ;

//This is an example of how to use RooKeys PDF

void RooKeysEx()
{

  cout << "--------------------------"<< endl;
  cout << "Start the RooKeys Example "<< endl;
  cout << "--------------------------"<< endl;

  // This is the variable that we will use ex: it could be the invariant mass  Lb- >pK gamma
  RooRealVar x("K#pi#mu#mu","K#pi#mu#mu",1800,1940, "MeV/c2") ;



   // This is some constant number that we use to generate the toy data
  RooPolynomial p("p","p",x,RooArgList(RooConst(0.01),RooConst(-0.01),RooConst(0.0004))) ;
 // These are some PDF that we use to generate the toy day
  RooExponential exp("exp", "exp", x, RooConst(-0.01));
  RooGaussian signal("signal", "signal",  x, RooConst(1866), RooConst(8));
  RooGaussian G1("G1", "G1",x,RooConst(1850),RooConst(10) );
  RooGaussian G2("G2", "G2",x,RooConst(1840),RooConst(15) );


  RooAddPdf SumG1andG2("SumG1andG2", "SumG1andG2",RooArgList(G1,G2), RooArgList(RooConst(0.6) ) );
  //RooAddPdf TotalPDF ("TotalPDF", "TotalPDF", RooArgList(signal, SumG1andG2),RooArgList(RooConst(0.6)));
  RooAddPdf TotalPDF("TotalPDF", "TotalPDF", RooArgList(signal, SumG1andG2, exp),RooArgList(RooConst(0.4), RooConst(0.2) ));

  // Prepare the today data
  //RooDataSet* data1 = p.generate(x,10000) ;
  //RooDataSet* data1 = G.generate(x,5000) ;
  //RooDataSet* data1 = SumG1andG2.generate(x,500) ;

  // Here we generate some events in the x variable using the PDFs that we defined before
  RooDataSet* data1= SumG1andG2.generate(x,1000);

  RooDataSet* data2 = TotalPDF.generate(x,10000) ; // generate signal and background


  cout << "--------------------------"<< endl;
  cout << " Done generating the data " << endl;
  cout << "--------------------------"<< endl;


  // is mirrored over the boundaries to minimize edge effects in distribution
  // that do not fall to zero towards the edges
  RooKeysPdf kest1("kest1","kest1",x,*data1,RooKeysPdf::MirrorBoth) ;
  // An adaptive kernel estimation pdf on the same data without mirroring option
  // for comparison
  RooKeysPdf kest2("kest2","kest2",x,*data1,RooKeysPdf::NoMirror) ;
  // Adaptive kernel estimation pdf with increased bandwidth scale factor
  // (promotes smoothness over detail preservation)
  RooKeysPdf kest3("kest3","kest3",x,*data1,RooKeysPdf::MirrorBoth,2) ;



 //===========================================================================
  // Here are the PDF and variable that we will use in the fit :
  //===========================================================================
  RooRealVar signal_mean("signal_mean","signal_mean", 1850, 1840, 1900);
  RooRealVar signal_width("signal_width", "signal_width", 10, 0, 20);
  //----
  RooRealVar comb_slope ("comb_slope", "comb_slope", 0, -1., 1);
  //----
  RooRealVar Nsig("Nsig", "Nsig", 100, 0, 10000);
  RooRealVar Nbkg_comb("Nbkg_comb", "Nbkg_comb", 100, 0, 10000);
  RooRealVar Nbkg_misID("Nbkg_misID", "Nbkg_misID", 100, 0, 10000);
  //----
  RooGaussian signal_fitted ("signal_fitted", "signal_fitted", x, signal_mean, signal_width);
  //  RooKeysPdf background_fitted ("background_fitted", "background_fitted",  x, *data1, RooKeysPdf::MirrorBoth );
  //----
  RooExponential comb_PDF ("comb_PDF", "comb_PDF", x, comb_slope);
  //---
  //RooAddPdf TotalPDF_fitted ("TotalPDF_fitted", "TotalPDF_fitted", RooArgList(signal_fitted, background_fitted), RooArgList(Nsig,Nbkg));
  RooAddPdf TotalPDF_fitted ("TotalPDF_fitted", "TotalPDF_fitted", RooArgList(signal_fitted, kest1, comb_PDF  ), RooArgList(Nsig,Nbkg_misID, Nbkg_comb));
  TotalPDF_fitted -> fitTo(*data2);
  cout << "----------------------------"<< endl;
  cout << " Done doing the fit to data " << endl;
  cout << "----------------------------"<< endl;




  // Here we are going to plot stuff


  RooPlot* frame0 = x.frame(Title("Adaptive kernel estimation pdf with and w/o mirroring"),Bins(20)) ;

  //  G1.plotOn(frame0, LineColor(1)) ;
  //G2.plotOn(frame0, LineColor(2)) ;
  SumG1andG2.plotOn(frame0) ;





  // Plot kernel estimation pdfs with and without mirroring over data
  RooPlot* frame = x.frame(Title("Adaptive kernel estimation pdf with and w/o mirroring"),Bins(20)) ;
  data1->plotOn(frame) ;
  kest1.plotOn(frame, LineColor(875)) ;
  //kest2.plotOn(frame, LineColor(429)) ;
  //kest3.plotOn(frame, LineColor(625)) ;

  // Plot kernel estimation pdfs with regular and increased bandwidth
  RooPlot* frame2 = x.frame(Title("Adaptive kernel estimation pdf with regular, increased bandwidth")) ;
  kest1.plotOn(frame2, LineColor(875)) ;
  kest3.plotOn(frame2, LineColor(807)) ;



  RooPlot* frame3 = x.frame(Title(""));
   data2->plotOn(frame3);

   TotalPDF_fitted.plotOn(frame3, Components(kest1), LineColor(875));
   TotalPDF_fitted.plotOn(frame3, Components(signal_fitted), LineColor(429));
   TotalPDF_fitted.plotOn(frame3, Components(comb_PDF), LineColor(801));
   TotalPDF_fitted.plotOn(frame3);

  TCanvas* c = new TCanvas("c","c",600,600) ;
  c->Divide(2,2) ;
  c->cd(1) ;
  frame0->GetYaxis()->SetTitleOffset(1.4) ; frame0->Draw() ;
  c->cd(2) ;
  frame->GetYaxis()->SetTitleOffset(1.4) ; frame->Draw() ;
  c->cd(3) ;
  frame2->GetYaxis()->SetTitleOffset(1.8) ; frame2->Draw() ;
  c->cd(4) ;
  frame3->GetYaxis()->SetTitleOffset(1.8) ; frame3->Draw() ;
  c->SaveAs("Keys.pdf");

  TCanvas * c1 = new TCanvas("c1","c1", 400,400);
  //c1->cd();
  frame3->GetYaxis()->SetTitleOffset(1.8) ; frame3->Draw() ;
  frame3->SetTitle("");
  c1->SaveAs("KeysResult.pdf");

  cout<< "Fin"<< endl;








}
