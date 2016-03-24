#include "RooTDalitz/RooDalitzAmplitude.h"
#include "TList.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "TSystem.h"
#include "RooDataSet.h"
#include "TFile.h"
#include "RooFitResult.h"
#include "TDatime.h"
#include "RooRealConstant.h"
#include "RooMinuit.h"
#include "RooAddition.h"
#include "TLorentzVector.h"
#include "TBranch.h"
//#include "RooTDalitz/RooSUMVar.h"

using namespace RooFit;

double xmkp;
double xcostheta;
double xcostheta1;
double xcostheta2;
double xphi1;
double xphi2;
double xcosthetaB;
double xcosthetaZ;
double xcosthetaPsi;
double xphiZ;
double xphiPsi;
double xphiMu;
double xcosthetap;
double xalphaMu;
double xmjpsip;
 
inline bool exists(const std::string& name) {
  return ( access( name.c_str(), F_OK ) != -1 );
//  struct stat buffer;   
//  return (stat (name.c_str(), &buffer) == 0); 
}

void helicityJpsiLam(TLorentzVector *Pmu1, TLorentzVector *Pmu2, TLorentzVector *PJ,
                     TLorentzVector *Pproton, TLorentzVector *Pkaon)
{  
  TLorentzVector PJpsi = *Pmu1 + *Pmu2;
  TLorentzVector PLam = *Pproton+ *Pkaon;
  TLorentzVector PB = PJpsi+PLam;
  TLorentzVector PJpsi_MassCon = *PJ;
//  PJpsi_MassCon.SetVectM(PJpsi.Vect(), 3096.916);
  xmkp = PLam.Mag()/1000.;
  xmjpsip = (PJpsi_MassCon + *Pproton).Mag()/1000.;
  TVector3 PBunit = PB.Vect().Unit();
  TVector3 Pbeam(0.0,0.0,1.0);
  //  TVector3  nvec = PBunit.Cross(TVector3(0,0,1));
  
  TVector3 boosttoparent = -PB.BoostVector();
  
  //  std::cout << "before " <<std::endl;
  //  PLam.Print();
  Pmu1->Boost(boosttoparent);
  Pmu2->Boost(boosttoparent);
  Pproton->Boost(boosttoparent);
  Pkaon->Boost(boosttoparent);
  PJpsi.Boost(boosttoparent);
  PLam.Boost(boosttoparent);
  
  TVector3 LUnit = PLam.Vect().Unit();
  TVector3 psiUnit = PJpsi.Vect().Unit();

  xcostheta = PBunit.Dot(LUnit);
//  phi = 0;
  //Lam rest frame
  boosttoparent = -PLam.BoostVector();
  Pkaon->Boost(boosttoparent);

  TVector3 Unit = Pkaon->Vect().Unit();
  xcostheta1 = -(psiUnit).Dot(Unit);
  
  TVector3 aboost = -(PB.Vect() - ((PB.Vect()).Dot(LUnit))*LUnit).Unit();
  double cosphi = aboost.Dot(Unit);
  double sinphi = (((-psiUnit).Cross(aboost)).Unit()).Dot(Unit);
  xphi1 = atan2(sinphi,cosphi);
  
  //Jpsi rest frame
  boosttoparent = -PJpsi.BoostVector();
  Pmu2->Boost(boosttoparent);
  Unit = Pmu2->Vect().Unit();
  xcostheta2 = (psiUnit).Dot(Unit);
  
  cosphi = aboost.Dot(Unit);
  sinphi = (((psiUnit).Cross(aboost)).Unit()).Dot(Unit);
  xphi2 = atan2(sinphi,cosphi);
     
}

void helicityZK(TLorentzVector *Pmu1, TLorentzVector *Pmu2,
                     TLorentzVector *Pproton, TLorentzVector *Pkaon)
{
  TLorentzVector PJpsi = *Pmu1 + *Pmu2;
  TLorentzVector PLam = *Pproton+ *Pkaon;
  TLorentzVector PB = PJpsi+PLam;
  TLorentzVector PZ = PJpsi + *Pproton;
  TVector3 PBunit = PB.Vect().Unit();
  TVector3 Pbeam(0.0,0.0,1.0);
  //  TVector3  nvec = PBunit.Cross(TVector3(0,0,1));
  
  TVector3 boosttoparent = -PB.BoostVector();
  
  //  std::cout << "before " <<std::endl;
  //  PLam.Print();
  Pmu1->Boost(boosttoparent);
  Pmu2->Boost(boosttoparent);
  Pproton->Boost(boosttoparent);
  Pkaon->Boost(boosttoparent);
  PJpsi.Boost(boosttoparent);
  PLam.Boost(boosttoparent);
  PZ.Boost(boosttoparent);
  
  TVector3 LUnit = PLam.Vect().Unit();
  TVector3 psiUnit = PJpsi.Vect().Unit();
  TVector3 ZUnit = PZ.Vect().Unit();

  xcosthetaB = PBunit.Dot(ZUnit);
  TVector3 aboost = (PLam.Vect()-(PLam.Vect().Dot(PBunit))*PBunit).Unit();
  double cosphi = aboost.Dot(ZUnit);
  double sinphi = ((PBunit.Cross(aboost)).Unit()).Dot(ZUnit);
  xphiZ = atan2(sinphi,cosphi);
  
  //Z rest frame
  boosttoparent = -PZ.BoostVector();
  PJpsi.Boost(boosttoparent);
  Pmu1->Boost(boosttoparent);
  Pmu2->Boost(boosttoparent);
  Pproton->Boost(boosttoparent);
  Pkaon->Boost(boosttoparent);
  
  TVector3 Unit = PJpsi.Vect().Unit();
  xcosthetaZ = -(Pkaon->Vect().Unit()).Dot(Unit);
  
  aboost = -(PB.Vect() - ((PB.Vect()).Dot(ZUnit))*ZUnit).Unit();
  cosphi = aboost.Dot(Unit);
  sinphi = (((-Pkaon->Vect().Unit()).Cross(aboost)).Unit()).Dot(Unit);
  xphiPsi = atan2(sinphi,cosphi);
  
  //Jpsi rest frame
  psiUnit = PJpsi.Vect().Unit();
  boosttoparent = -PJpsi.BoostVector();
  Pmu2->Boost(boosttoparent);
  Unit = Pmu2->Vect().Unit();
  xcosthetaPsi = (psiUnit).Dot(Unit);
  aboost = -(-Pkaon->Vect() + ((Pkaon->Vect()).Dot(psiUnit))*(psiUnit)).Unit();
  cosphi = aboost.Dot(Unit);
  sinphi = (((psiUnit).Cross(aboost)).Unit()).Dot(Unit);
  xphiMu = atan2(sinphi,cosphi);
}


void helicityTwoFrame(TLorentzVector *Pmu1, TLorentzVector *Pmu2,
                     TLorentzVector *Pproton, TLorentzVector *Pkaon)
{
  TLorentzVector PJpsi = *Pmu1 + *Pmu2;
  TLorentzVector PLam = *Pproton+ *Pkaon;
  TLorentzVector PB = PJpsi+PLam;
  TLorentzVector PZ = PJpsi + *Pproton;
  TVector3 PBunit = PB.Vect().Unit();
  TVector3 Pbeam(0.0,0.0,1.0);
  //  TVector3  nvec = PBunit.Cross(TVector3(0,0,1));
  
  TVector3 boosttoparent = -PB.BoostVector();
  
  //  std::cout << "before " <<std::endl;
  //  PLam.Print();
  Pmu1->Boost(boosttoparent);
  Pmu2->Boost(boosttoparent);
  Pproton->Boost(boosttoparent);
  Pkaon->Boost(boosttoparent);
  PJpsi.Boost(boosttoparent);
  PLam.Boost(boosttoparent);
  PZ.Boost(boosttoparent);
  
   
  TLorentzVector PKaon_p = *Pkaon;
  TLorentzVector PJpsi_p = PJpsi;
  
  boosttoparent = -Pproton->BoostVector();
  PKaon_p.Boost(boosttoparent);
  PJpsi_p.Boost(boosttoparent);
  /*
  boosttoparent = - PLam.BoostVector();
  TLorentzVector PKaon_p = *Pkaon;
  TLorentzVector Pproton_p = *Pproton;
  PKaon_p.Boost(boosttoparent);
  Pproton_p.Boost(boosttoparent);
  boosttoparent = - Pproton_p.BoostVector();
  PKaon_p.Boost(boosttoparent);
  
  boosttoparent = -PZ.BoostVector();
  TLorentzVector PJpsi_p = PJpsi;
  Pproton_p = *Pproton;
  PJpsi_p.Boost(boosttoparent);
  Pproton_p.Boost(boosttoparent);
  boosttoparent = - Pproton_p.BoostVector();
  PJpsi_p.Boost(boosttoparent);
   */ 
  xcosthetap = (PKaon_p.Vect().Unit()).Dot(PJpsi_p.Vect().Unit());
    
/*  
  TLorentzVector PKaon_Z = *Pkaon;
  TLorentzVector PJpsi_Z = PJpsi;
  TLorentzVector PKaon_L = *Pkaon;
  TLorentzVector PJpsi_L = PJpsi;
  boosttoparent = -PZ.BoostVector();
  PKaon_Z.Boost(boosttoparent);
  PJpsi_Z.Boost(boosttoparent);  
  boosttoparent = -PLam.BoostVector();
  PKaon_L.Boost(boosttoparent);
  PJpsi_L.Boost(boosttoparent);  
  
  TVector3 x0_Z = (-PKaon_Z.Vect()+(PKaon_Z.Vect().Dot(PJpsi_Z.Vect().Unit()))*(PJpsi_Z.Vect().Unit())).Unit();
  TVector3 x0_L = (-PJpsi_L.Vect()+(PJpsi_L.Vect().Dot(PKaon_L.Vect().Unit()))*(PKaon_L.Vect().Unit())).Unit();  
  TVector3 z0_Z = -PJpsi_p.Vect().Unit();
  std::cout << (z0_Z.Cross(x0_Z)).Unit().Dot(x0_L) << std::endl;
  */
  
  
  //Jpsi rest frame
  boosttoparent = -PJpsi.BoostVector();
  Pmu2->Boost(boosttoparent);
  Pproton->Boost(boosttoparent);
  PLam.Boost(boosttoparent);
  TVector3 z3 = Pmu2->Vect().Unit();
  TVector3 x3_Z = -(-Pproton->Vect()+(Pproton->Vect().Dot(Pmu2->Vect().Unit()))*Pmu2->Vect().Unit()).Unit();    
  TVector3 x3_L = -(-PLam.Vect()+(PLam.Vect().Dot(Pmu2->Vect().Unit()))*Pmu2->Vect().Unit()).Unit();      
  
  double sinphi = ((z3.Cross(x3_Z)).Unit()).Dot(x3_L);
  double cosphi = x3_Z.Dot(x3_L);
  
  xalphaMu = atan2(sinphi,cosphi);
  
}

//set how many SL of Lambda to float
void SetLPar(RooArgList* argli, int maxind, bool fixfirst=false)
{
    int  spin = ((RooAbsReal&)(*argli)[argli->getSize()-2]).getVal();
//    std::cout << spin << std::endl;
  int maxi = maxind-1;
  if(spin==1&&maxi>3) maxi = 3;
//    std::cout << spin << std::endl;
//  if(spin==1&&ind>3)  return;
  int mini = 0;
  if(fixfirst) mini = 1;
  for(int ind=mini; ind<=maxi; ++ind) {
    RooRealVar *var = (RooRealVar*)(argli->at(2*ind));
    var->setError(0.1); var->setConstant(0);
    RooRealVar *var1 = (RooRealVar*)(argli->at(2*ind+1));
    var1->setError(0.1); var1->setConstant(0);
  }
  for(int ind=maxi+1; ind<6; ++ind) {
    RooRealVar *var = (RooRealVar*)(argli->at(2*ind));
    var->setVal(0.); var->setConstant(1);
    RooRealVar *var1 = (RooRealVar*)(argli->at(2*ind+1));
    var1->setVal(0.); var1->setConstant(1);
  }
}

void FloatPar(RooArgList* argli, int ind)
{
  if(ind>5) return;
    int  spin = ((RooAbsReal&)(*argli)[argli->getSize()-2]).getVal();
//    std::cout << spin << std::endl;
  if(spin==1&&ind>3)  return;
  RooRealVar *var = (RooRealVar*)(argli->at(2*ind));
  var->setError(0.1); var->setConstant(0);
  RooRealVar *var1 = (RooRealVar*)(argli->at(2*ind+1));
  var1->setError(0.1); var1->setConstant(0);
}

void ResetPar(RooArgList* argli)
{
  int  J = ((RooAbsReal&)(*argli)[argli->getSize()-2]).getVal();
  int ind(0);
  for(int S = abs(J-2); S<=J+2; S+=2) {
    for(int L = S - 1; L<= S + 1; L+=2) {       
      RooRealVar *var = (RooRealVar*)(argli->at(2*ind));
      var->setVal(var->getVal()*sqrt(2./((double)S+1.)));
      RooRealVar *var1 = (RooRealVar*)(argli->at(2*ind+1));
      var1->setVal(var1->getVal()*sqrt(2./((double)S+1.)));
      ind++;
   }
  }
}


void FloatZPar(RooArgList* argli)
{
  int  spin = ((RooAbsReal&)(*argli)[argli->getSize()-2]).getVal();
//    std::cout << spin << std::endl;
//  if(spin==1&&ind>3)  return;
  for(int ind=0; ind<=3; ++ind) {
    RooRealVar *var = (RooRealVar*)(argli->at(2*ind));
    var->setError(0.1); var->setConstant(0);
    RooRealVar *var1 = (RooRealVar*)(argli->at(2*ind+1));
    var1->setError(0.1); var1->setConstant(0);
  }
  if(abs(spin)<=1) {
    int ind=3;
    RooRealVar *var = (RooRealVar*)(argli->at(2*ind));
    var->setVal(0.); var->setConstant(1);
    RooRealVar *var1 = (RooRealVar*)(argli->at(2*ind+1));
    var1->setVal(0.); var1->setConstant(1);
  }
}

int main(int argc, char **argv)
{
 TDatime BeginTime;
 std::cout << "Time(begin)  " << BeginTime.GetYear() << "." << BeginTime.GetMonth() << "." << BeginTime.GetDay() << "    " << BeginTime.GetHour() << ":" << BeginTime.GetMinute() << ":" << BeginTime.GetSecond() << std::endl; 

  
//#ifdef __CINT__
//  gROOT->ProcessLineSync(".x RooDalitzAmplitude.cxx+");
//#endif
  if( argc < 2){
     fprintf(stderr,"usage: Fit Z2[J] Z2[P] Z[J] Z[P]\n");
     exit(3);
   }
  int ZS = 5;//atoi( argv[3] );
  int  ZP = 1;//atoi( argv[4] );
  std::cout << "Fit Z J P " << ZS << " " << ZP << std::endl;  

  int NoCPU(1);

  char *inputname = argv[1];
  int ZS2 = 3;//atoi( argv[1] );
  int  ZP2 = -1;//atoi( argv[2] );  
  
  std::cout << "Fit Z2 J P " << ZS2 << " " << ZP2 << std::endl;    
  gSystem->Load("RooIpatia2_cxx.so");
          
          
  RooRealVar mkp("mkp","m(K^{-}p)",1.4,2.6);
  RooRealVar mjpsip("mjpsip","",4.,5.2);
//  cosk:cosmu:chi
  RooRealVar cosTheta_L("cosTheta_L","cosTheta_L",-1,1);
  RooRealVar cosTheta_Jpsi("cosTheta_Jpsi","cosTheta_Jpsi",-1,1);
  RooRealVar cosTheta_Lb("cosTheta_Lb","cosTheta_Lb",-1,1);
  RooRealVar Z_cosTheta_Lb("Z_cosTheta_Lb","Z_cosTheta_Lb",-1,1);
  RooRealVar Z_cosTheta_Z("Z_cosTheta_Z","Z_cosTheta_Z",-1,1);
  RooRealVar Z_cosTheta_Jpsi("Z_cosTheta_Jpsi","Z_cosTheta_Jpsi",-1,1);
  RooRealVar cosTheta_p("cosTheta_p","cosTheta_p",-1,1);
  RooRealVar phiK("phiK","phiK",-TMath::Pi(), TMath::Pi());
  RooRealVar phiMu("phiMu","phiMu",-TMath::Pi(), TMath::Pi());
  RooRealVar Z_phiZ("Z_phiZ","Z_phiZ",-TMath::Pi(), TMath::Pi());
  RooRealVar Z_phiJpsi("Z_phiJpsi","Z_phiJpsi",-TMath::Pi(), TMath::Pi());
  RooRealVar Z_phiMu("Z_phiMu","Z_phiMu",-TMath::Pi(), TMath::Pi());
  RooRealVar alpha_Mu("alpha_Mu","alpha_Mu",-TMath::Pi(), TMath::Pi());
  RooArgList *obs = new RooArgList(mkp,cosTheta_Lb,cosTheta_L,cosTheta_Jpsi,phiK,phiMu);  
  obs->add(Z_cosTheta_Lb);
  obs->add(Z_cosTheta_Z);
  obs->add(Z_cosTheta_Jpsi);
  obs->add(Z_phiZ);
  obs->add(Z_phiJpsi);
  obs->add(Z_phiMu);
  obs->add(cosTheta_p);
  obs->add(alpha_Mu);
  obs->add(mjpsip);
  
  RooRealVar sw("sw","sw",0);
  RooArgList *obs1 = new RooArgList(mkp,cosTheta_Lb,cosTheta_L,cosTheta_Jpsi,phiK,phiMu);  
  obs1->add(Z_cosTheta_Lb);
  obs1->add(Z_cosTheta_Z);
  obs1->add(Z_cosTheta_Jpsi);
  obs1->add(Z_phiZ);
  obs1->add(Z_phiJpsi);
  obs1->add(Z_phiMu);
  obs1->add(cosTheta_p);
  obs1->add(alpha_Mu);
  obs1->add(mjpsip);

  obs1->add(sw);
  TFile *fdata = new TFile("sPlot.root");
  RooDataSet *datars1 = (RooDataSet*)fdata->Get("datars");
  RooDataSet *datars = (RooDataSet*)datars1->reduce("mkp<2.522584&&mjpsip>4.0351880460&&mjpsip<5.125823");
  RooRealVar *index = new RooRealVar("index","index",0);
  RooDataSet *IND = new RooDataSet("IND","", RooArgSet(*index));
  

  double nev = datars->numEntries();
  std::cout << "sum " << nev << std::endl;
  for(int i=0; i<nev; ++i) {
     *index = i;
//     if(i%1000==0) std::cout << "i " << i << std::endl;
     IND->add(RooArgSet(*index));

  }
  datars->merge(IND); 
//  datars->Print("V");  
  RooDataSet * data_fit = new  RooDataSet(TString(datars->GetName())+TString("new"),datars->GetTitle(),datars,*datars->get(),0,"nsig_sw");
//  data_fit->Print("V");  
  
  obs->add(*index);

  TFile *fmc = new TFile("mcsw.root");//mcswbin.root");
  TTree *tree = (TTree*)fmc->Get("h1");
  RooDataSet * mcsw = new RooDataSet("mcsw","",tree, *obs1);
  RooDataSet *INDMC = new RooDataSet("INDMC","", RooArgSet(*index));
  nev = mcsw->numEntries();

  for(int i=0; i<nev; ++i) {
     *index = i;
//     if(i%1000==0) std::cout << "i " << i << std::endl;
     INDMC->add(RooArgSet(*index));
  }
  mcsw->merge(INDMC); 
//  mcsw->Print("V");  
  RooDataSet * mcnorm = new  RooDataSet(TString(mcsw->GetName())+TString("new"),mcsw->GetTitle(),mcsw,*mcsw->get(),0,"sw");
//  mcnorm->Print("V");  
  std::cout.precision(20);
  std::cout << "# mc" << nev << " sum weighted " <<  std::endl;
  
  double mcr = nev/(double)datars->numEntries();
  int cpu1 = (double) NoCPU - (double)NoCPU/(1.+1./sqrt(mcr));
  int cpu2 = NoCPU-cpu1;
  std::cout << "CPU DIV " << cpu1 << " " << cpu2 << std::endl;
  
  datars->Delete();
  datars1->Delete();
  mcsw->Delete();
//  exit(0);  
  double w(0.);
  double w2(0.);
  for( int i = 0; i < data_fit->numEntries(); ++i ){
    data_fit->get(i);
    double nsw = data_fit->weight();
    w += nsw;
    w2 += pow(nsw,2);
  }
  double swf2 = w/w2;
  double swf  = sqrt(swf2);

  double Wmc = mcnorm->sumEntries();
  double Wda = data_fit->sumEntries();
  RooRealVar Weimc("Weimc","Wei of mc norm", Wmc);
  RooRealVar Weida("Weida","Wei of data", Wda);
  RooRealVar SWF2("SWF2","Scale ", swf2);
  
  std::cout << "w data " << w << " " << Wda << std::endl;
  std::cout << "alpha to factor the -2lnL" << swf2 << std::endl;  
  TList *listZ = new TList();
  //===Lambda(Z) 1/2-
  RooRealVar m0_Z("m0_Z","m0",4.45,4.4,4.5);
  RooRealVar width_Z("width_Z","width",0.04,0,0.1);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_Bp_Z("a_Bp_Z","",0.5,-5,5);
  RooRealVar p_Bp_Z("p_Bp_Z","",0.5,-5,5);
  RooRealVar a_L0_Z("a_L0_Z","",0.03);
  RooRealVar p_L0_Z("p_L0_Z","",0.04);
  RooRealVar a_L1_Z("a_L1_Z","",0.02);
  RooRealVar p_L1_Z("p_L1_Z","",0.01);
  RooRealVar a_L2_Z("a_L2_Z","",0.02);
  RooRealVar p_L2_Z("p_L2_Z","",0.01);

  RooArgList* L_Z = new RooArgList(a_Bp_Z,p_Bp_Z,
				      a_L0_Z,p_L0_Z,
				      a_L1_Z,p_L1_Z,
				      a_L2_Z,p_L2_Z,"L_Z");
  L_Z->add(m0_Z);
  L_Z->add(width_Z);
  L_Z->add(RooRealConstant::value(ZS));//2 x spin
  L_Z->add(RooRealConstant::value(ZP)); //parity
  
  listZ->Add(L_Z);
  
  //===Lambda(Z) 1/2-
  RooRealVar m0_Z2("m0_Z2","m0",4.35,4.2,4.4);
  RooRealVar width_Z2("width_Z2","width",0.2,0,0.5);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_Bp_Z2("a_Bp_Z2","",0.5,-5,5);
  RooRealVar p_Bp_Z2("p_Bp_Z2","",0.5,-5,5);
  RooRealVar a_L0_Z2("a_L0_Z2","",0.5);
  RooRealVar p_L0_Z2("p_L0_Z2","",0.5);
  RooRealVar a_L1_Z2("a_L1_Z2","",0.5);
  RooRealVar p_L1_Z2("p_L1_Z2","",0.5);
  RooRealVar a_L2_Z2("a_L2_Z2","",0);
  RooRealVar p_L2_Z2("p_L2_Z2","",0);

  RooArgList* L_Z2 = new RooArgList(a_Bp_Z2,p_Bp_Z2,
				      a_L0_Z2,p_L0_Z2,
				      a_L1_Z2,p_L1_Z2,
				      a_L2_Z2,p_L2_Z2,"L_Z2");
  L_Z2->add(m0_Z2);
  L_Z2->add(width_Z2);
  L_Z2->add(RooRealConstant::value(ZS2));//2 x spin
  L_Z2->add(RooRealConstant::value(ZP2)); //parity
  
  listZ->Add(L_Z2);

  TList *list = new TList();
  //===Lambda(1520) 3/2-
  RooRealVar m0_1520("m0_1520","m0",1.5195);
  RooRealVar width_1520("width_1520","width",0.0156);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_p1p1_1520("a_p1p1_1520","",1.);
  RooRealVar p_p1p1_1520("p_p1p1_1520","",0.);
  RooRealVar a_p100_1520("a_p100_1520","",0);
  RooRealVar p_p100_1520("p_p100_1520","",0);
  RooRealVar a_m100_1520("a_m100_1520","",0);
  RooRealVar p_m100_1520("p_m100_1520","",0);
  RooRealVar a_m1m1_1520("a_m1m1_1520","",0);
  RooRealVar p_m1m1_1520("p_m1m1_1520","",0);
  RooRealVar a_p3p1_1520("a_p3p1_1520","",0);
  RooRealVar p_p3p1_1520("p_p3p1_1520","",0);
  RooRealVar a_m3m1_1520("a_m3m1_1520","",0);
  RooRealVar p_m3m1_1520("p_m3m1_1520","",0);

  RooArgList* L_1520 = new RooArgList(a_p1p1_1520,p_p1p1_1520,
				      a_p100_1520,p_p100_1520,
				      a_m100_1520,p_m100_1520,
				      a_m1m1_1520,p_m1m1_1520,"L_1520");
  L_1520->add(a_p3p1_1520);
  L_1520->add(p_p3p1_1520);
  L_1520->add(a_m3m1_1520);
  L_1520->add(p_m3m1_1520);
  L_1520->add(m0_1520);
  L_1520->add(width_1520);
  L_1520->add(RooRealConstant::value(3.));//2 x spin
  L_1520->add(RooRealConstant::value(-1)); //parity
  


    //===Lambda(1600) 1/2+
  RooRealVar m0_1600("m0_1600","m0",1.6);
  RooRealVar width_1600("width_1600","width",0.15);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_p1p1_1600("a_p1p1_1600","",40.39/100,-300,300);
  RooRealVar p_p1p1_1600("p_p1p1_1600","",-1.018/100,-300,300);
  RooRealVar a_p100_1600("a_p100_1600","",0);
  RooRealVar p_p100_1600("p_p100_1600","",0);
  RooRealVar a_m100_1600("a_m100_1600","",0);
  RooRealVar p_m100_1600("p_m100_1600","",0);
  RooRealVar a_m1m1_1600("a_m1m1_1600","",0);
  RooRealVar p_m1m1_1600("p_m1m1_1600","",0);
  RooRealVar a_p3p1_1600("a_p3p1_1600","",0);
  RooRealVar p_p3p1_1600("p_p3p1_1600","",0);
  RooRealVar a_m3m1_1600("a_m3m1_1600","",0);
  RooRealVar p_m3m1_1600("p_m3m1_1600","",0);

  RooArgList* L_1600 = new RooArgList(a_p1p1_1600,p_p1p1_1600,
				      a_p100_1600,p_p100_1600,
				      a_m100_1600,p_m100_1600,
				      a_m1m1_1600,p_m1m1_1600,"L_1600");
  L_1600->add(a_p3p1_1600);
  L_1600->add(p_p3p1_1600);
  L_1600->add(a_m3m1_1600);
  L_1600->add(p_m3m1_1600);
  L_1600->add(m0_1600);
  L_1600->add(width_1600);
  L_1600->add(RooRealConstant::value(1.));//2 x spin
  L_1600->add(RooRealConstant::value(1)); //parity


  //===Lambda(1670) 1/2-
  RooRealVar m0_1670("m0_1670","m0",1.67);
  RooRealVar width_1670("width_1670","width",0.035);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_p1p1_1670("a_p1p1_1670","",-0.5976/100,-300,300);
  RooRealVar p_p1p1_1670("p_p1p1_1670","",1.218/100,-300,300);
  RooRealVar a_p100_1670("a_p100_1670","",0);
  RooRealVar p_p100_1670("p_p100_1670","",0);
  RooRealVar a_m100_1670("a_m100_1670","",0);
  RooRealVar p_m100_1670("p_m100_1670","",0);
  RooRealVar a_m1m1_1670("a_m1m1_1670","",0);
  RooRealVar p_m1m1_1670("p_m1m1_1670","",0);
  RooRealVar a_p3p1_1670("a_p3p1_1670","",0);
  RooRealVar p_p3p1_1670("p_p3p1_1670","",0);
  RooRealVar a_m3m1_1670("a_m3m1_1670","",0);
  RooRealVar p_m3m1_1670("p_m3m1_1670","",0);

  RooArgList* L_1670 = new RooArgList(a_p1p1_1670,p_p1p1_1670,
				      a_p100_1670,p_p100_1670,
				      a_m100_1670,p_m100_1670,
				      a_m1m1_1670,p_m1m1_1670,"L_1670");
  L_1670->add(a_p3p1_1670);
  L_1670->add(p_p3p1_1670);
  L_1670->add(a_m3m1_1670);
  L_1670->add(p_m3m1_1670);
  L_1670->add(m0_1670);
  L_1670->add(width_1670);
  L_1670->add(RooRealConstant::value(1.));//2 x spin
  L_1670->add(RooRealConstant::value(-1)); //parity



  
  //===Lambda(1690) 3/2-
  RooRealVar m0_1690("m0_1690","m0",1.715);
  RooRealVar width_1690("width_1690","width",0.06);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_p1p1_1690("a_p1p1_1690","",42.5/100,-300,300);
  RooRealVar p_p1p1_1690("p_p1p1_1690","",2.471/100,-300,300);
  RooRealVar a_p100_1690("a_p100_1690","",0);
  RooRealVar p_p100_1690("p_p100_1690","",0);
  RooRealVar a_m100_1690("a_m100_1690","",0);
  RooRealVar p_m100_1690("p_m100_1690","",0);
  RooRealVar a_m1m1_1690("a_m1m1_1690","",0);
  RooRealVar p_m1m1_1690("p_m1m1_1690","",0);
  RooRealVar a_p3p1_1690("a_p3p1_1690","",0);
  RooRealVar p_p3p1_1690("p_p3p1_1690","",0);
  RooRealVar a_m3m1_1690("a_m3m1_1690","",0);
  RooRealVar p_m3m1_1690("p_m3m1_1690","",0);

  RooArgList* L_1690 = new RooArgList(a_p1p1_1690,p_p1p1_1690,
				      a_p100_1690,p_p100_1690,
				      a_m100_1690,p_m100_1690,
				      a_m1m1_1690,p_m1m1_1690,"L_1690");
  L_1690->add(a_p3p1_1690);
  L_1690->add(p_p3p1_1690);
  L_1690->add(a_m3m1_1690);
  L_1690->add(p_m3m1_1690);
  L_1690->add(m0_1690);
  L_1690->add(width_1690);
  L_1690->add(RooRealConstant::value(3.));//2 x spin
  L_1690->add(RooRealConstant::value(-1)); //parity
  



  //===Lambda(1800) 1/2-
  RooRealVar m0_1800("m0_1800","m0",1.8);
  RooRealVar width_1800("width_1800","width",0.3);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_p1p1_1800("a_p1p1_1800","",0.00,-300,300);
  RooRealVar p_p1p1_1800("p_p1p1_1800","",0.00,-300,300);
  RooRealVar a_p100_1800("a_p100_1800","",0);
  RooRealVar p_p100_1800("p_p100_1800","",0);
  RooRealVar a_m100_1800("a_m100_1800","",0);
  RooRealVar p_m100_1800("p_m100_1800","",0);
  RooRealVar a_m1m1_1800("a_m1m1_1800","",0);
  RooRealVar p_m1m1_1800("p_m1m1_1800","",0);
  RooRealVar a_p3p1_1800("a_p3p1_1800","",0);
  RooRealVar p_p3p1_1800("p_p3p1_1800","",0);
  RooRealVar a_m3m1_1800("a_m3m1_1800","",0);
  RooRealVar p_m3m1_1800("p_m3m1_1800","",0);

  RooArgList* L_1800 = new RooArgList(a_p1p1_1800,p_p1p1_1800,
				      a_p100_1800,p_p100_1800,
				      a_m100_1800,p_m100_1800,
				      a_m1m1_1800,p_m1m1_1800,"L_1800");
  L_1800->add(a_p3p1_1800);
  L_1800->add(p_p3p1_1800);
  L_1800->add(a_m3m1_1800);
  L_1800->add(p_m3m1_1800);
  L_1800->add(m0_1800);
  L_1800->add(width_1800);
  L_1800->add(RooRealConstant::value(1.));//2 x spin
  L_1800->add(RooRealConstant::value(-1)); //parity


  //===Lambda(1810) 1/2+
  RooRealVar m0_1810("m0_1810","m0",1.81);
  RooRealVar width_1810("width_1810","width",0.15);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_p1p1_1810("a_p1p1_1810","",-3.752/100,-300,300);
  RooRealVar p_p1p1_1810("p_p1p1_1810","",-10.28/100,-300,300);
  RooRealVar a_p100_1810("a_p100_1810","",0);
  RooRealVar p_p100_1810("p_p100_1810","",0);
  RooRealVar a_m100_1810("a_m100_1810","",0);
  RooRealVar p_m100_1810("p_m100_1810","",0);
  RooRealVar a_m1m1_1810("a_m1m1_1810","",0);
  RooRealVar p_m1m1_1810("p_m1m1_1810","",0);
  RooRealVar a_p3p1_1810("a_p3p1_1810","",0);
  RooRealVar p_p3p1_1810("p_p3p1_1810","",0);
  RooRealVar a_m3m1_1810("a_m3m1_1810","",0);
  RooRealVar p_m3m1_1810("p_m3m1_1810","",0);

  RooArgList* L_1810 = new RooArgList(a_p1p1_1810,p_p1p1_1810,
				      a_p100_1810,p_p100_1810,
				      a_m100_1810,p_m100_1810,
				      a_m1m1_1810,p_m1m1_1810,"L_1810");
  L_1810->add(a_p3p1_1810);
  L_1810->add(p_p3p1_1810);
  L_1810->add(a_m3m1_1810);
  L_1810->add(p_m3m1_1810);
  L_1810->add(m0_1810);
  L_1810->add(width_1810);
  L_1810->add(RooRealConstant::value(1.));//2 x spin
  L_1810->add(RooRealConstant::value(1)); //parity



  //===Lambda(1820) 5/2+
  RooRealVar m0_1820("m0_1820","m0",1.82);
  RooRealVar width_1820("width_1820","width",0.08);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_p1p1_1820("a_p1p1_1820","",0.2,-300,300);
  RooRealVar p_p1p1_1820("p_p1p1_1820","",-10.28/100,-300,300);
  RooRealVar a_p100_1820("a_p100_1820","",0);
  RooRealVar p_p100_1820("p_p100_1820","",0);
  RooRealVar a_m100_1820("a_m100_1820","",0);
  RooRealVar p_m100_1820("p_m100_1820","",0);
  RooRealVar a_m1m1_1820("a_m1m1_1820","",0);
  RooRealVar p_m1m1_1820("p_m1m1_1820","",0);
  RooRealVar a_p3p1_1820("a_p3p1_1820","",0);
  RooRealVar p_p3p1_1820("p_p3p1_1820","",0);
  RooRealVar a_m3m1_1820("a_m3m1_1820","",0);
  RooRealVar p_m3m1_1820("p_m3m1_1820","",0);

  RooArgList* L_1820 = new RooArgList(a_p1p1_1820,p_p1p1_1820,
				      a_p100_1820,p_p100_1820,
				      a_m100_1820,p_m100_1820,
				      a_m1m1_1820,p_m1m1_1820,"L_1820");
  L_1820->add(a_p3p1_1820);
  L_1820->add(p_p3p1_1820);
  L_1820->add(a_m3m1_1820);
  L_1820->add(p_m3m1_1820);
  L_1820->add(m0_1820);
  L_1820->add(width_1820);
  L_1820->add(RooRealConstant::value(5.));//2 x spin
  L_1820->add(RooRealConstant::value(1)); //parity
  
//  list->Add(L_1820);


  //===Lambda(1830) 5/2-
  RooRealVar m0_1830("m0_1830","m0",1.83);
  RooRealVar width_1830("width_1830","width",0.095);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_p1p1_1830("a_p1p1_1830","",0.2,-300,300);
  RooRealVar p_p1p1_1830("p_p1p1_1830","",-10.28/100,-300,300);
  RooRealVar a_p100_1830("a_p100_1830","",0);
  RooRealVar p_p100_1830("p_p100_1830","",0);
  RooRealVar a_m100_1830("a_m100_1830","",0);
  RooRealVar p_m100_1830("p_m100_1830","",0);
  RooRealVar a_m1m1_1830("a_m1m1_1830","",0);
  RooRealVar p_m1m1_1830("p_m1m1_1830","",0);
  RooRealVar a_p3p1_1830("a_p3p1_1830","",0);
  RooRealVar p_p3p1_1830("p_p3p1_1830","",0);
  RooRealVar a_m3m1_1830("a_m3m1_1830","",0);
  RooRealVar p_m3m1_1830("p_m3m1_1830","",0);

  RooArgList* L_1830 = new RooArgList(a_p1p1_1830,p_p1p1_1830,
				      a_p100_1830,p_p100_1830,
				      a_m100_1830,p_m100_1830,
				      a_m1m1_1830,p_m1m1_1830,"L_1830");
  L_1830->add(a_p3p1_1830);
  L_1830->add(p_p3p1_1830);
  L_1830->add(a_m3m1_1830);
  L_1830->add(p_m3m1_1830);
  L_1830->add(m0_1830);
  L_1830->add(width_1830);
  L_1830->add(RooRealConstant::value(5.));//2 x spin
  L_1830->add(RooRealConstant::value(-1)); //parity
  
//  list->Add(L_1830);

  //===Lambda(1890) 3/2+
  RooRealVar m0_1890("m0_1890","m0",1.89);
  RooRealVar width_1890("width_1890","width",0.1);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_p1p1_1890("a_p1p1_1890","",0.2,-300,300);
  RooRealVar p_p1p1_1890("p_p1p1_1890","",-10.28/100,-300,300);
  RooRealVar a_p100_1890("a_p100_1890","",0);
  RooRealVar p_p100_1890("p_p100_1890","",0);
  RooRealVar a_m100_1890("a_m100_1890","",0);
  RooRealVar p_m100_1890("p_m100_1890","",0);
  RooRealVar a_m1m1_1890("a_m1m1_1890","",0);
  RooRealVar p_m1m1_1890("p_m1m1_1890","",0);
  RooRealVar a_p3p1_1890("a_p3p1_1890","",0);
  RooRealVar p_p3p1_1890("p_p3p1_1890","",0);
  RooRealVar a_m3m1_1890("a_m3m1_1890","",0);
  RooRealVar p_m3m1_1890("p_m3m1_1890","",0);

  RooArgList* L_1890 = new RooArgList(a_p1p1_1890,p_p1p1_1890,
				      a_p100_1890,p_p100_1890,
				      a_m100_1890,p_m100_1890,
				      a_m1m1_1890,p_m1m1_1890,"L_1890");
  L_1890->add(a_p3p1_1890);
  L_1890->add(p_p3p1_1890);
  L_1890->add(a_m3m1_1890);
  L_1890->add(p_m3m1_1890);
  L_1890->add(m0_1890);
  L_1890->add(width_1890);
  L_1890->add(RooRealConstant::value(3.));//2 x spin
  L_1890->add(RooRealConstant::value(1)); //parity
  
//  list->Add(L_1890);

  //2110 5/2+
  RooRealVar m0_2110("m0_2110","m0",2.11);
  RooRealVar width_2110("width_2110","width",0.2);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_p1p1_2110("a_p1p1_2110","",1.,-300,300);
  RooRealVar p_p1p1_2110("p_p1p1_2110","",0.,-300,300);
  RooRealVar a_p100_2110("a_p100_2110","",0.);
  RooRealVar p_p100_2110("p_p100_2110","",0.);
  RooRealVar a_m100_2110("a_m100_2110","",0.);
  RooRealVar p_m100_2110("p_m100_2110","",0.);
  RooRealVar a_m1m1_2110("a_m1m1_2110","",0.);
  RooRealVar p_m1m1_2110("p_m1m1_2110","",0.);
  RooRealVar a_p3p1_2110("a_p3p1_2110","",0.);
  RooRealVar p_p3p1_2110("p_p3p1_2110","",0.);
  RooRealVar a_m3m1_2110("a_m3m1_2110","",0.);
  RooRealVar p_m3m1_2110("p_m3m1_2110","",0.);

  RooArgList* L_2110 = new RooArgList(a_p1p1_2110,p_p1p1_2110,
				      a_p100_2110,p_p100_2110,
				      a_m100_2110,p_m100_2110,
				      a_m1m1_2110,p_m1m1_2110,"L_2110");
  L_2110->add(a_p3p1_2110);
  L_2110->add(p_p3p1_2110);
  L_2110->add(a_m3m1_2110);
  L_2110->add(p_m3m1_2110);
  L_2110->add(m0_2110);
  L_2110->add(width_2110);
  L_2110->add(RooRealConstant::value(5.));//2 x spin
  L_2110->add(RooRealConstant::value(1)); //parity

  //2100 7/2-
  RooRealVar m0_2100("m0_2100","m0",2.10);
  RooRealVar width_2100("width_2100","width",0.2);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_p1p1_2100("a_p1p1_2100","",1.,-300,300);
  RooRealVar p_p1p1_2100("p_p1p1_2100","",0.,-300,300);
  RooRealVar a_p100_2100("a_p100_2100","",0.);
  RooRealVar p_p100_2100("p_p100_2100","",0.);
  RooRealVar a_m100_2100("a_m100_2100","",0.);
  RooRealVar p_m100_2100("p_m100_2100","",0.);
  RooRealVar a_m1m1_2100("a_m1m1_2100","",0.);
  RooRealVar p_m1m1_2100("p_m1m1_2100","",0.);
  RooRealVar a_p3p1_2100("a_p3p1_2100","",0.);
  RooRealVar p_p3p1_2100("p_p3p1_2100","",0.);
  RooRealVar a_m3m1_2100("a_m3m1_2100","",0.);
  RooRealVar p_m3m1_2100("p_m3m1_2100","",0.);

  RooArgList* L_2100 = new RooArgList(a_p1p1_2100,p_p1p1_2100,
				      a_p100_2100,p_p100_2100,
				      a_m100_2100,p_m100_2100,
				      a_m1m1_2100,p_m1m1_2100,"L_2100");
  L_2100->add(a_p3p1_2100);
  L_2100->add(p_p3p1_2100);
  L_2100->add(a_m3m1_2100);
  L_2100->add(p_m3m1_2100);
  L_2100->add(m0_2100);
  L_2100->add(width_2100);
  L_2100->add(RooRealConstant::value(7.));//2 x spin
  L_2100->add(RooRealConstant::value(-1)); //parity


  //===Lambda(1405) 1/2-
  RooRealVar m0_1405("m0_1405","m0",1.4051);
  RooRealVar width_1405("width_1405","width",0.0505);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_p1p1_1405("a_p1p1_1405","",0.001,-300,300);
  RooRealVar p_p1p1_1405("p_p1p1_1405","",0.,-300,300);
  RooRealVar a_p100_1405("a_p100_1405","",0);
  RooRealVar p_p100_1405("p_p100_1405","",0);
  RooRealVar a_m100_1405("a_m100_1405","",0);
  RooRealVar p_m100_1405("p_m100_1405","",0);
  RooRealVar a_m1m1_1405("a_m1m1_1405","",0);
  RooRealVar p_m1m1_1405("p_m1m1_1405","",0);
  RooRealVar a_p3p1_1405("a_p3p1_1405","",0);
  RooRealVar p_p3p1_1405("p_p3p1_1405","",0);
  RooRealVar a_m3m1_1405("a_m3m1_1405","",0);
  RooRealVar p_m3m1_1405("p_m3m1_1405","",0);

  RooArgList* L_1405 = new RooArgList(a_p1p1_1405,p_p1p1_1405,
				      a_p100_1405,p_p100_1405,
				      a_m100_1405,p_m100_1405,
				      a_m1m1_1405,p_m1m1_1405,"L_1405");
  L_1405->add(a_p3p1_1405);
  L_1405->add(p_p3p1_1405);
  L_1405->add(a_m3m1_1405);
  L_1405->add(p_m3m1_1405);
  L_1405->add(m0_1405);
  L_1405->add(width_1405);
  L_1405->add(RooRealConstant::value(1.));//2 x spin
  L_1405->add(RooRealConstant::value(-1)); //parity


  //===Lambda(2350) 9/2+
  RooRealVar m0_2350("m0_2350","m0",2.35);
  RooRealVar width_2350("width_2350","width",0.15);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_p1p1_2350("a_p1p1_2350","",0.001,-10,10);
  RooRealVar p_p1p1_2350("p_p1p1_2350","",0.,-10,10);
  RooRealVar a_p100_2350("a_p100_2350","",0);
  RooRealVar p_p100_2350("p_p100_2350","",0);
  RooRealVar a_m100_2350("a_m100_2350","",0);
  RooRealVar p_m100_2350("p_m100_2350","",0);
  RooRealVar a_m1m1_2350("a_m1m1_2350","",0);
  RooRealVar p_m1m1_2350("p_m1m1_2350","",0);
  RooRealVar a_p3p1_2350("a_p3p1_2350","",0);
  RooRealVar p_p3p1_2350("p_p3p1_2350","",0);
  RooRealVar a_m3m1_2350("a_m3m1_2350","",0);
  RooRealVar p_m3m1_2350("p_m3m1_2350","",0);

  RooArgList* L_2350 = new RooArgList(a_p1p1_2350,p_p1p1_2350,
				      a_p100_2350,p_p100_2350,
				      a_m100_2350,p_m100_2350,
				      a_m1m1_2350,p_m1m1_2350,"L_2350");
  L_2350->add(a_p3p1_2350);
  L_2350->add(p_p3p1_2350);
  L_2350->add(a_m3m1_2350);
  L_2350->add(p_m3m1_2350);
  L_2350->add(m0_2350);
  L_2350->add(width_2350);
  L_2350->add(RooRealConstant::value(9.));//2 x spin
  L_2350->add(RooRealConstant::value(1)); //parity
  
//  list->Add(L_2350);

  RooRealVar m0_2585("m0_2585","m0",2.585);
  RooRealVar width_2585("width_2585","width",0.3);
  // a &phase [1/2,1] [1/2,0] [-1/2, 0] [-1/2, -1] [3/2,1] [-3/2, -1]
  RooRealVar a_p1p1_2585("a_p1p1_2585","",1.,-300,300);
  RooRealVar p_p1p1_2585("p_p1p1_2585","",0.,-300,300);
  RooRealVar a_p100_2585("a_p100_2585","",0.);
  RooRealVar p_p100_2585("p_p100_2585","",0.);
  RooRealVar a_m100_2585("a_m100_2585","",0.);
  RooRealVar p_m100_2585("p_m100_2585","",0.);
  RooRealVar a_m1m1_2585("a_m1m1_2585","",0.);
  RooRealVar p_m1m1_2585("p_m1m1_2585","",0.);
  RooRealVar a_p3p1_2585("a_p3p1_2585","",0.);
  RooRealVar p_p3p1_2585("p_p3p1_2585","",0.);
  RooRealVar a_m3m1_2585("a_m3m1_2585","",0.);
  RooRealVar p_m3m1_2585("p_m3m1_2585","",0.);

  RooArgList* L_2585 = new RooArgList(a_p1p1_2585,p_p1p1_2585,
				      a_p100_2585,p_p100_2585,
				      a_m100_2585,p_m100_2585,
				      a_m1m1_2585,p_m1m1_2585,"L_2585");
  L_2585->add(a_p3p1_2585);
  L_2585->add(p_p3p1_2585);
  L_2585->add(a_m3m1_2585);
  L_2585->add(p_m3m1_2585);
  L_2585->add(m0_2585);
  L_2585->add(width_2585);
  L_2585->add(RooRealConstant::value(5.));//2 x spin
  L_2585->add(RooRealConstant::value(-1)); //parity
  
//  list->Add(L_2585);  
  list->Add(L_1405);  
  list->Add(L_1520);  
  list->Add(L_1600);  
  list->Add(L_1670);  
  list->Add(L_1690);  
  list->Add(L_1800);  
  list->Add(L_1810);  
  list->Add(L_1820);  
  list->Add(L_1830);
  list->Add(L_1890);
  list->Add(L_2100);    
  list->Add(L_2110);    
//  list->Add(L_2350);    
  
//  return ;
  RooDalitzAmplitude *sig = new RooDalitzAmplitude("sig","", *obs, list, listZ, "MCFlatAcc.root",*data_fit);
  //sig->genToy("Test");
  RooArgSet* setdlz = sig->getParameters(*data_fit);
//  RooArgSet* setdlz1 = signorm->getParameters(*data_fit);  

   
//   setdlz->readFromFile("result/fitz-m1.func");   
  char name[10];
  if(ZP>0) {
   sprintf(name,"p%i",ZS);
  } else {
   sprintf(name,"m%i",ZS);  
  }
  
  char nameZ2[10];
  if(ZP2>0) {
   sprintf(nameZ2,"p%i",ZS2);
  } else {
   sprintf(nameZ2,"m%i",ZS2);  
  }  
  //setdlz->readFromFile(TString::Format("result/fitz-%s.func",name).Data());//"fitz.func");   
//  setdlz->readFromFile(TString::Format("fit2zall-3half-%s.func",name).Data());


  setdlz->readFromFile("fit2zall-m3-p5-res5.func");  
  
  //setdlz->readFromFile("func/fit2zall-m3-p5-res8.func");  
//  if(exists(TString::Format("func/fit2zall-%s-%s-res8.func",nameZ2,name).Data())) {  
//    setdlz->readFromFile(TString::Format("funcold/fit2zall-%s-%s-res8.func",nameZ2,name).Data());    
/*    a_Bp_Z2.setVal(a_Bp_Z2.getVal()/10.);
    p_Bp_Z2.setVal(p_Bp_Z2.getVal()/10.);    
    a_L0_Z2.setVal(a_L0_Z2.getVal()/10.);
    p_L0_Z2.setVal(p_L0_Z2.getVal()/10.);    
    a_L1_Z2.setVal(a_L1_Z2.getVal()/10.);    
    p_L1_Z2.setVal(p_L1_Z2.getVal()/10.);    
    a_Bp_Z.setVal(a_Bp_Z.getVal()/10.);
    p_Bp_Z.setVal(p_Bp_Z.getVal()/10.);    
    a_L0_Z.setVal(a_L0_Z.getVal()/10.);
    p_L0_Z.setVal(p_L0_Z.getVal()/10.);    
    a_L1_Z.setVal(a_L1_Z.getVal()/10.);    
    p_L1_Z.setVal(p_L1_Z.getVal()/10.);   
    */
//  }           
  
  SetLPar(L_1405,3);      
  SetLPar(L_1520,5,true);  
  SetLPar(L_1600,3);  
  SetLPar(L_1670,3);  
  SetLPar(L_1690,5);  
  SetLPar(L_1800,5);  
  SetLPar(L_1810,3);  
  SetLPar(L_1820,1);  
  SetLPar(L_1830,1);
  SetLPar(L_1890,3);
  SetLPar(L_2110,1);    
  SetLPar(L_2100,1);      
//  SetLPar(L_2350,3);  
  FloatZPar(L_Z);
  FloatZPar(L_Z2);
  
//  setdlz->readFromFile("fit2zall-1half-p1.func");
  setdlz->Print("V"); 
  int nres = list->GetSize();//+listZ->GetSize();
  double Dsum[20];
  sig->getInt(Dsum);
  char resname[20];
  double sum(0);
  std::cout << "======Fit Fraction======" << std::endl;
  for(int i=0; i<nres; ++i) {
     printf("%10s %6.2f\n", (list->At(i))->GetName(),Dsum[i]*100.0 );
     sum += Dsum[i];
  }
  for(int i=0; i<listZ->GetSize(); ++i) {
    printf("%10s %6.2f\n", (listZ->At(i))->GetName(),Dsum[i+nres]*100.0 );
    sum += Dsum[i+nres];    
  }
  printf("Total  %6.2f\n", sum*100.0 );  

//  sig->genToy("Test");           
//  RooFitResult* fitRes = sig->fitTo(*data_fit, Hesse(0), Verbose(1), NumCPU(NoCPU), Strategy(1), Timer(1), Minimizer("Minuit","Migrad"),  Save(), SumW2Error(0));

  TFile *flhcbmc =new TFile(TString(inputname)+".root","update");
  TTree *fChain = (TTree*)flhcbmc->Get("DecayTreeTuple/DecayTree");
  
   Int_t           B0_TRUEID;
   Double_t        B0_TRUEP_E;
   Double_t        B0_TRUEP_X;
   Double_t        B0_TRUEP_Y;
   Double_t        B0_TRUEP_Z;
   Double_t        B0_TRUEPT;
   Double_t        B0_TRUEORIGINVERTEX_X;
   Double_t        B0_TRUEORIGINVERTEX_Y;
   Double_t        B0_TRUEORIGINVERTEX_Z;
   Double_t        B0_TRUEENDVERTEX_X;
   Double_t        B0_TRUEENDVERTEX_Y;
   Double_t        B0_TRUEENDVERTEX_Z;
   Bool_t          B0_TRUEISSTABLE;
   Double_t        B0_TRUETAU;
   Int_t           JPs_TRUEID;
   Double_t        JPs_TRUEP_E;
   Double_t        JPs_TRUEP_X;
   Double_t        JPs_TRUEP_Y;
   Double_t        JPs_TRUEP_Z;
   Double_t        JPs_TRUEPT;
   Double_t        JPs_TRUEORIGINVERTEX_X;
   Double_t        JPs_TRUEORIGINVERTEX_Y;
   Double_t        JPs_TRUEORIGINVERTEX_Z;
   Double_t        JPs_TRUEENDVERTEX_X;
   Double_t        JPs_TRUEENDVERTEX_Y;
   Double_t        JPs_TRUEENDVERTEX_Z;
   Bool_t          JPs_TRUEISSTABLE;
   Double_t        JPs_TRUETAU;
   Int_t           E1_TRUEID;
   Double_t        E1_TRUEP_E;
   Double_t        E1_TRUEP_X;
   Double_t        E1_TRUEP_Y;
   Double_t        E1_TRUEP_Z;
   Double_t        E1_TRUEPT;
   Double_t        E1_TRUEORIGINVERTEX_X;
   Double_t        E1_TRUEORIGINVERTEX_Y;
   Double_t        E1_TRUEORIGINVERTEX_Z;
   Double_t        E1_TRUEENDVERTEX_X;
   Double_t        E1_TRUEENDVERTEX_Y;
   Double_t        E1_TRUEENDVERTEX_Z;
   Bool_t          E1_TRUEISSTABLE;
   Double_t        E1_TRUETAU;
   Int_t           E2_TRUEID;
   Double_t        E2_TRUEP_E;
   Double_t        E2_TRUEP_X;
   Double_t        E2_TRUEP_Y;
   Double_t        E2_TRUEP_Z;
   Double_t        E2_TRUEPT;
   Double_t        E2_TRUEORIGINVERTEX_X;
   Double_t        E2_TRUEORIGINVERTEX_Y;
   Double_t        E2_TRUEORIGINVERTEX_Z;
   Double_t        E2_TRUEENDVERTEX_X;
   Double_t        E2_TRUEENDVERTEX_Y;
   Double_t        E2_TRUEENDVERTEX_Z;
   Bool_t          E2_TRUEISSTABLE;
   Double_t        E2_TRUETAU;
   Int_t           Kst_TRUEID;
   Double_t        Kst_TRUEP_E;
   Double_t        Kst_TRUEP_X;
   Double_t        Kst_TRUEP_Y;
   Double_t        Kst_TRUEP_Z;
   Double_t        Kst_TRUEPT;
   Double_t        Kst_TRUEORIGINVERTEX_X;
   Double_t        Kst_TRUEORIGINVERTEX_Y;
   Double_t        Kst_TRUEORIGINVERTEX_Z;
   Double_t        Kst_TRUEENDVERTEX_X;
   Double_t        Kst_TRUEENDVERTEX_Y;
   Double_t        Kst_TRUEENDVERTEX_Z;
   Bool_t          Kst_TRUEISSTABLE;
   Double_t        Kst_TRUETAU;
   Int_t           K_TRUEID;
   Double_t        K_TRUEP_E;
   Double_t        K_TRUEP_X;
   Double_t        K_TRUEP_Y;
   Double_t        K_TRUEP_Z;
   Double_t        K_TRUEPT;
   Double_t        K_TRUEORIGINVERTEX_X;
   Double_t        K_TRUEORIGINVERTEX_Y;
   Double_t        K_TRUEORIGINVERTEX_Z;
   Double_t        K_TRUEENDVERTEX_X;
   Double_t        K_TRUEENDVERTEX_Y;
   Double_t        K_TRUEENDVERTEX_Z;
   Bool_t          K_TRUEISSTABLE;
   Double_t        K_TRUETAU;
   Int_t           Pi_TRUEID;
   Double_t        Pi_TRUEP_E;
   Double_t        Pi_TRUEP_X;
   Double_t        Pi_TRUEP_Y;
   Double_t        Pi_TRUEP_Z;
   Double_t        Pi_TRUEPT;
   Double_t        Pi_TRUEORIGINVERTEX_X;
   Double_t        Pi_TRUEORIGINVERTEX_Y;
   Double_t        Pi_TRUEORIGINVERTEX_Z;
   Double_t        Pi_TRUEENDVERTEX_X;
   Double_t        Pi_TRUEENDVERTEX_Y;
   Double_t        Pi_TRUEENDVERTEX_Z;
   Bool_t          Pi_TRUEISSTABLE;
   Double_t        Pi_TRUETAU;
   TBranch        *b_B0_TRUEID;   //!
   TBranch        *b_B0_TRUEP_E;   //!
   TBranch        *b_B0_TRUEP_X;   //!
   TBranch        *b_B0_TRUEP_Y;   //!
   TBranch        *b_B0_TRUEP_Z;   //!
   TBranch        *b_B0_TRUEPT;   //!
   TBranch        *b_B0_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_B0_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_B0_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_B0_TRUEENDVERTEX_X;   //!
   TBranch        *b_B0_TRUEENDVERTEX_Y;   //!
   TBranch        *b_B0_TRUEENDVERTEX_Z;   //!
   TBranch        *b_B0_TRUEISSTABLE;   //!
   TBranch        *b_B0_TRUETAU;   //!
   TBranch        *b_JPs_TRUEID;   //!
   TBranch        *b_JPs_TRUEP_E;   //!
   TBranch        *b_JPs_TRUEP_X;   //!
   TBranch        *b_JPs_TRUEP_Y;   //!
   TBranch        *b_JPs_TRUEP_Z;   //!
   TBranch        *b_JPs_TRUEPT;   //!
   TBranch        *b_JPs_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_JPs_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_JPs_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_JPs_TRUEENDVERTEX_X;   //!
   TBranch        *b_JPs_TRUEENDVERTEX_Y;   //!
   TBranch        *b_JPs_TRUEENDVERTEX_Z;   //!
   TBranch        *b_JPs_TRUEISSTABLE;   //!
   TBranch        *b_JPs_TRUETAU;   //!
   TBranch        *b_E1_TRUEID;   //!
   TBranch        *b_E1_TRUEP_E;   //!
   TBranch        *b_E1_TRUEP_X;   //!
   TBranch        *b_E1_TRUEP_Y;   //!
   TBranch        *b_E1_TRUEP_Z;   //!
   TBranch        *b_E1_TRUEPT;   //!
   TBranch        *b_E1_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_E1_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_E1_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_E1_TRUEENDVERTEX_X;   //!
   TBranch        *b_E1_TRUEENDVERTEX_Y;   //!
   TBranch        *b_E1_TRUEENDVERTEX_Z;   //!
   TBranch        *b_E1_TRUEISSTABLE;   //!
   TBranch        *b_E1_TRUETAU;   //!
   TBranch        *b_E2_TRUEID;   //!
   TBranch        *b_E2_TRUEP_E;   //!
   TBranch        *b_E2_TRUEP_X;   //!
   TBranch        *b_E2_TRUEP_Y;   //!
   TBranch        *b_E2_TRUEP_Z;   //!
   TBranch        *b_E2_TRUEPT;   //!
   TBranch        *b_E2_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_E2_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_E2_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_E2_TRUEENDVERTEX_X;   //!
   TBranch        *b_E2_TRUEENDVERTEX_Y;   //!
   TBranch        *b_E2_TRUEENDVERTEX_Z;   //!
   TBranch        *b_E2_TRUEISSTABLE;   //!
   TBranch        *b_E2_TRUETAU;   //!
   TBranch        *b_Kst_TRUEID;   //!
   TBranch        *b_Kst_TRUEP_E;   //!
   TBranch        *b_Kst_TRUEP_X;   //!
   TBranch        *b_Kst_TRUEP_Y;   //!
   TBranch        *b_Kst_TRUEP_Z;   //!
   TBranch        *b_Kst_TRUEPT;   //!
   TBranch        *b_Kst_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_Kst_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_Kst_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_Kst_TRUEENDVERTEX_X;   //!
   TBranch        *b_Kst_TRUEENDVERTEX_Y;   //!
   TBranch        *b_Kst_TRUEENDVERTEX_Z;   //!
   TBranch        *b_Kst_TRUEISSTABLE;   //!
   TBranch        *b_Kst_TRUETAU;   //!
   TBranch        *b_K_TRUEID;   //!
   TBranch        *b_K_TRUEP_E;   //!
   TBranch        *b_K_TRUEP_X;   //!
   TBranch        *b_K_TRUEP_Y;   //!
   TBranch        *b_K_TRUEP_Z;   //!
   TBranch        *b_K_TRUEPT;   //!
   TBranch        *b_K_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_K_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_K_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_K_TRUEENDVERTEX_X;   //!
   TBranch        *b_K_TRUEENDVERTEX_Y;   //!
   TBranch        *b_K_TRUEENDVERTEX_Z;   //!
   TBranch        *b_K_TRUEISSTABLE;   //!
   TBranch        *b_K_TRUETAU;   //!
   TBranch        *b_Pi_TRUEID;   //!
   TBranch        *b_Pi_TRUEP_E;   //!
   TBranch        *b_Pi_TRUEP_X;   //!
   TBranch        *b_Pi_TRUEP_Y;   //!
   TBranch        *b_Pi_TRUEP_Z;   //!
   TBranch        *b_Pi_TRUEPT;   //!
   TBranch        *b_Pi_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_Pi_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_Pi_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_Pi_TRUEENDVERTEX_X;   //!
   TBranch        *b_Pi_TRUEENDVERTEX_Y;   //!
   TBranch        *b_Pi_TRUEENDVERTEX_Z;   //!
   TBranch        *b_Pi_TRUEISSTABLE;   //!
   TBranch        *b_Pi_TRUETAU;   //!
   fChain->SetBranchAddress("B0_TRUEID", &B0_TRUEID, &b_B0_TRUEID);
   fChain->SetBranchAddress("B0_TRUEP_E", &B0_TRUEP_E, &b_B0_TRUEP_E);
   fChain->SetBranchAddress("B0_TRUEP_X", &B0_TRUEP_X, &b_B0_TRUEP_X);
   fChain->SetBranchAddress("B0_TRUEP_Y", &B0_TRUEP_Y, &b_B0_TRUEP_Y);
   fChain->SetBranchAddress("B0_TRUEP_Z", &B0_TRUEP_Z, &b_B0_TRUEP_Z);
   fChain->SetBranchAddress("B0_TRUEPT", &B0_TRUEPT, &b_B0_TRUEPT);
   fChain->SetBranchAddress("B0_TRUEORIGINVERTEX_X", &B0_TRUEORIGINVERTEX_X, &b_B0_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("B0_TRUEORIGINVERTEX_Y", &B0_TRUEORIGINVERTEX_Y, &b_B0_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("B0_TRUEORIGINVERTEX_Z", &B0_TRUEORIGINVERTEX_Z, &b_B0_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("B0_TRUEENDVERTEX_X", &B0_TRUEENDVERTEX_X, &b_B0_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("B0_TRUEENDVERTEX_Y", &B0_TRUEENDVERTEX_Y, &b_B0_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("B0_TRUEENDVERTEX_Z", &B0_TRUEENDVERTEX_Z, &b_B0_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("B0_TRUEISSTABLE", &B0_TRUEISSTABLE, &b_B0_TRUEISSTABLE);
   fChain->SetBranchAddress("B0_TRUETAU", &B0_TRUETAU, &b_B0_TRUETAU);
   fChain->SetBranchAddress("JPs_TRUEID", &JPs_TRUEID, &b_JPs_TRUEID);
   fChain->SetBranchAddress("JPs_TRUEP_E", &JPs_TRUEP_E, &b_JPs_TRUEP_E);
   fChain->SetBranchAddress("JPs_TRUEP_X", &JPs_TRUEP_X, &b_JPs_TRUEP_X);
   fChain->SetBranchAddress("JPs_TRUEP_Y", &JPs_TRUEP_Y, &b_JPs_TRUEP_Y);
   fChain->SetBranchAddress("JPs_TRUEP_Z", &JPs_TRUEP_Z, &b_JPs_TRUEP_Z);
   fChain->SetBranchAddress("JPs_TRUEPT", &JPs_TRUEPT, &b_JPs_TRUEPT);
   fChain->SetBranchAddress("JPs_TRUEORIGINVERTEX_X", &JPs_TRUEORIGINVERTEX_X, &b_JPs_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("JPs_TRUEORIGINVERTEX_Y", &JPs_TRUEORIGINVERTEX_Y, &b_JPs_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("JPs_TRUEORIGINVERTEX_Z", &JPs_TRUEORIGINVERTEX_Z, &b_JPs_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("JPs_TRUEENDVERTEX_X", &JPs_TRUEENDVERTEX_X, &b_JPs_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("JPs_TRUEENDVERTEX_Y", &JPs_TRUEENDVERTEX_Y, &b_JPs_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("JPs_TRUEENDVERTEX_Z", &JPs_TRUEENDVERTEX_Z, &b_JPs_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("JPs_TRUEISSTABLE", &JPs_TRUEISSTABLE, &b_JPs_TRUEISSTABLE);
   fChain->SetBranchAddress("JPs_TRUETAU", &JPs_TRUETAU, &b_JPs_TRUETAU);
   fChain->SetBranchAddress("E1_TRUEID", &E1_TRUEID, &b_E1_TRUEID);
   fChain->SetBranchAddress("E1_TRUEP_E", &E1_TRUEP_E, &b_E1_TRUEP_E);
   fChain->SetBranchAddress("E1_TRUEP_X", &E1_TRUEP_X, &b_E1_TRUEP_X);
   fChain->SetBranchAddress("E1_TRUEP_Y", &E1_TRUEP_Y, &b_E1_TRUEP_Y);
   fChain->SetBranchAddress("E1_TRUEP_Z", &E1_TRUEP_Z, &b_E1_TRUEP_Z);
   fChain->SetBranchAddress("E1_TRUEPT", &E1_TRUEPT, &b_E1_TRUEPT);
   fChain->SetBranchAddress("E1_TRUEORIGINVERTEX_X", &E1_TRUEORIGINVERTEX_X, &b_E1_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("E1_TRUEORIGINVERTEX_Y", &E1_TRUEORIGINVERTEX_Y, &b_E1_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("E1_TRUEORIGINVERTEX_Z", &E1_TRUEORIGINVERTEX_Z, &b_E1_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("E1_TRUEENDVERTEX_X", &E1_TRUEENDVERTEX_X, &b_E1_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("E1_TRUEENDVERTEX_Y", &E1_TRUEENDVERTEX_Y, &b_E1_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("E1_TRUEENDVERTEX_Z", &E1_TRUEENDVERTEX_Z, &b_E1_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("E1_TRUEISSTABLE", &E1_TRUEISSTABLE, &b_E1_TRUEISSTABLE);
   fChain->SetBranchAddress("E1_TRUETAU", &E1_TRUETAU, &b_E1_TRUETAU);
   fChain->SetBranchAddress("E2_TRUEID", &E2_TRUEID, &b_E2_TRUEID);
   fChain->SetBranchAddress("E2_TRUEP_E", &E2_TRUEP_E, &b_E2_TRUEP_E);
   fChain->SetBranchAddress("E2_TRUEP_X", &E2_TRUEP_X, &b_E2_TRUEP_X);
   fChain->SetBranchAddress("E2_TRUEP_Y", &E2_TRUEP_Y, &b_E2_TRUEP_Y);
   fChain->SetBranchAddress("E2_TRUEP_Z", &E2_TRUEP_Z, &b_E2_TRUEP_Z);
   fChain->SetBranchAddress("E2_TRUEPT", &E2_TRUEPT, &b_E2_TRUEPT);
   fChain->SetBranchAddress("E2_TRUEORIGINVERTEX_X", &E2_TRUEORIGINVERTEX_X, &b_E2_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("E2_TRUEORIGINVERTEX_Y", &E2_TRUEORIGINVERTEX_Y, &b_E2_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("E2_TRUEORIGINVERTEX_Z", &E2_TRUEORIGINVERTEX_Z, &b_E2_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("E2_TRUEENDVERTEX_X", &E2_TRUEENDVERTEX_X, &b_E2_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("E2_TRUEENDVERTEX_Y", &E2_TRUEENDVERTEX_Y, &b_E2_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("E2_TRUEENDVERTEX_Z", &E2_TRUEENDVERTEX_Z, &b_E2_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("E2_TRUEISSTABLE", &E2_TRUEISSTABLE, &b_E2_TRUEISSTABLE);
   fChain->SetBranchAddress("E2_TRUETAU", &E2_TRUETAU, &b_E2_TRUETAU);
   fChain->SetBranchAddress("Kst_TRUEID", &Kst_TRUEID, &b_Kst_TRUEID);
   fChain->SetBranchAddress("Kst_TRUEP_E", &Kst_TRUEP_E, &b_Kst_TRUEP_E);
   fChain->SetBranchAddress("Kst_TRUEP_X", &Kst_TRUEP_X, &b_Kst_TRUEP_X);
   fChain->SetBranchAddress("Kst_TRUEP_Y", &Kst_TRUEP_Y, &b_Kst_TRUEP_Y);
   fChain->SetBranchAddress("Kst_TRUEP_Z", &Kst_TRUEP_Z, &b_Kst_TRUEP_Z);
   fChain->SetBranchAddress("Kst_TRUEPT", &Kst_TRUEPT, &b_Kst_TRUEPT);
   fChain->SetBranchAddress("Kst_TRUEORIGINVERTEX_X", &Kst_TRUEORIGINVERTEX_X, &b_Kst_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("Kst_TRUEORIGINVERTEX_Y", &Kst_TRUEORIGINVERTEX_Y, &b_Kst_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("Kst_TRUEORIGINVERTEX_Z", &Kst_TRUEORIGINVERTEX_Z, &b_Kst_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("Kst_TRUEENDVERTEX_X", &Kst_TRUEENDVERTEX_X, &b_Kst_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("Kst_TRUEENDVERTEX_Y", &Kst_TRUEENDVERTEX_Y, &b_Kst_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("Kst_TRUEENDVERTEX_Z", &Kst_TRUEENDVERTEX_Z, &b_Kst_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("Kst_TRUEISSTABLE", &Kst_TRUEISSTABLE, &b_Kst_TRUEISSTABLE);
   fChain->SetBranchAddress("Kst_TRUETAU", &Kst_TRUETAU, &b_Kst_TRUETAU);
   fChain->SetBranchAddress("K_TRUEID", &K_TRUEID, &b_K_TRUEID);
   fChain->SetBranchAddress("K_TRUEP_E", &K_TRUEP_E, &b_K_TRUEP_E);
   fChain->SetBranchAddress("K_TRUEP_X", &K_TRUEP_X, &b_K_TRUEP_X);
   fChain->SetBranchAddress("K_TRUEP_Y", &K_TRUEP_Y, &b_K_TRUEP_Y);
   fChain->SetBranchAddress("K_TRUEP_Z", &K_TRUEP_Z, &b_K_TRUEP_Z);
   fChain->SetBranchAddress("K_TRUEPT", &K_TRUEPT, &b_K_TRUEPT);
   fChain->SetBranchAddress("K_TRUEORIGINVERTEX_X", &K_TRUEORIGINVERTEX_X, &b_K_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("K_TRUEORIGINVERTEX_Y", &K_TRUEORIGINVERTEX_Y, &b_K_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("K_TRUEORIGINVERTEX_Z", &K_TRUEORIGINVERTEX_Z, &b_K_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("K_TRUEENDVERTEX_X", &K_TRUEENDVERTEX_X, &b_K_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("K_TRUEENDVERTEX_Y", &K_TRUEENDVERTEX_Y, &b_K_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("K_TRUEENDVERTEX_Z", &K_TRUEENDVERTEX_Z, &b_K_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("K_TRUEISSTABLE", &K_TRUEISSTABLE, &b_K_TRUEISSTABLE);
   fChain->SetBranchAddress("K_TRUETAU", &K_TRUETAU, &b_K_TRUETAU);
   fChain->SetBranchAddress("Pi_TRUEID", &Pi_TRUEID, &b_Pi_TRUEID);
   fChain->SetBranchAddress("Pi_TRUEP_E", &Pi_TRUEP_E, &b_Pi_TRUEP_E);
   fChain->SetBranchAddress("Pi_TRUEP_X", &Pi_TRUEP_X, &b_Pi_TRUEP_X);
   fChain->SetBranchAddress("Pi_TRUEP_Y", &Pi_TRUEP_Y, &b_Pi_TRUEP_Y);
   fChain->SetBranchAddress("Pi_TRUEP_Z", &Pi_TRUEP_Z, &b_Pi_TRUEP_Z);
   fChain->SetBranchAddress("Pi_TRUEPT", &Pi_TRUEPT, &b_Pi_TRUEPT);
   fChain->SetBranchAddress("Pi_TRUEORIGINVERTEX_X", &Pi_TRUEORIGINVERTEX_X, &b_Pi_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("Pi_TRUEORIGINVERTEX_Y", &Pi_TRUEORIGINVERTEX_Y, &b_Pi_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("Pi_TRUEORIGINVERTEX_Z", &Pi_TRUEORIGINVERTEX_Z, &b_Pi_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("Pi_TRUEENDVERTEX_X", &Pi_TRUEENDVERTEX_X, &b_Pi_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("Pi_TRUEENDVERTEX_Y", &Pi_TRUEENDVERTEX_Y, &b_Pi_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("Pi_TRUEENDVERTEX_Z", &Pi_TRUEENDVERTEX_Z, &b_Pi_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("Pi_TRUEISSTABLE", &Pi_TRUEISSTABLE, &b_Pi_TRUEISSTABLE);
   fChain->SetBranchAddress("Pi_TRUETAU", &Pi_TRUETAU, &b_Pi_TRUETAU);
  Double_t wdp;
  TBranch *newBranch1 = fChain->Branch("mkp", &xmkp, "mkp/D");
  TBranch *newBranch2 = fChain->Branch("mjpsip", &xmjpsip, "mjpsip/D");
  TBranch *newBranch = fChain->Branch("wdp", &wdp, "wdp/D");
  Long64_t numevt = fChain->GetEntries();
 // numevt = 100;
  for(Long64_t i=0; i<numevt; ++i) {
    fChain->GetEntry(i);
    wdp = -1.0;
    xmkp = -1.0;
    xmjpsip = -1.0;
    if(abs(B0_TRUEID)==5122)  {
      double  muplusPX, muplusPY, muplusPZ, muplusPE;
      double muminusPX, muminusPY, muminusPZ, muminusPE;
      double protonPX, protonPY, protonPZ, protonPE;
      double kaonPX, kaonPY, kaonPZ, kaonPE(-1);
      
      TLorentzVector p[4];
      Int_t id[4];
      p[0] = TLorentzVector(E1_TRUEP_X,E1_TRUEP_Y,E1_TRUEP_Z,E1_TRUEP_E); id[0] = E1_TRUEID;
      p[1] = TLorentzVector(E2_TRUEP_X,E2_TRUEP_Y,E2_TRUEP_Z,E2_TRUEP_E); id[1] = E2_TRUEID;      
      p[2] = TLorentzVector(K_TRUEP_X,K_TRUEP_Y,K_TRUEP_Z,K_TRUEP_E); id[2] = K_TRUEID;      
      p[3] = TLorentzVector(Pi_TRUEP_X,Pi_TRUEP_Y,Pi_TRUEP_Z,Pi_TRUEP_E); id[3] = Pi_TRUEID;      
      Int_t matchid[4] ={-1,-1,-1,-1};
      Int_t match(0);
      for(int ip=0; ip<4; ++ip) {
         if(abs(id[ip])==11&&B0_TRUEID*id[ip]<0) {
            if(matchid[0]<0) match++;         
            matchid[0] = ip;
         } else if(abs(id[ip])==11&&B0_TRUEID*id[ip]>0) {
            if(matchid[1]<0) match++;         
            matchid[1] = ip;
         } else if(abs(id[ip])==321) {   
            if(matchid[2]<0) match++;         
            matchid[2] = ip;
         }else if(abs(id[ip])==2212) {                    
            if(matchid[3]<0) match++;         
            matchid[3] = ip;
         }
      }
      if(match==4) {
	muplusPE = p[matchid[0]].E();
	muplusPX = p[matchid[0]].Px();
	muplusPY = p[matchid[0]].Py();
	muplusPZ = p[matchid[0]].Pz();
	muminusPE = p[matchid[1]].E();
	muminusPX = p[matchid[1]].Px();
	muminusPY = p[matchid[1]].Py();
	muminusPZ = p[matchid[1]].Pz();
	kaonPE = p[matchid[2]].E();
	kaonPX = p[matchid[2]].Px();
	kaonPY = p[matchid[2]].Py();
	kaonPZ = p[matchid[2]].Pz();
	protonPE = p[matchid[3]].E();
	protonPX = p[matchid[3]].Px();
	protonPY = p[matchid[3]].Py();
	protonPZ = p[matchid[3]].Pz();
      
      //      if(kaonPE>0)
      if(kaonPE>0) {
      int ID = B0_TRUEID/abs(B0_TRUEID);
      
      TLorentzVector *Pmu1 = new TLorentzVector( muplusPX, muplusPY, muplusPZ, muplusPE);
      TLorentzVector *Pmu2 = new TLorentzVector( muminusPX, muminusPY, muminusPZ, muminusPE);      
      TLorentzVector *Pproton = new TLorentzVector( protonPX, protonPY, protonPZ, protonPE);
      TLorentzVector *Pkaon = new TLorentzVector( kaonPX, kaonPY, kaonPZ, kaonPE);
      TLorentzVector *PJ = new TLorentzVector(muplusPX+muminusPX,muplusPY+muminusPY,muplusPZ+muminusPZ,muplusPE+muminusPE);
      if(JPs_TRUEID==443) {
        PJ =  new TLorentzVector(JPs_TRUEP_X, JPs_TRUEP_Y,JPs_TRUEP_Z, JPs_TRUEP_E);      
      } else {
        PJ = new TLorentzVector(B0_TRUEP_X-kaonPX-protonPX, B0_TRUEP_Y-kaonPY-protonPY, B0_TRUEP_Z-kaonPZ-protonPZ, B0_TRUEP_E-kaonPE-protonPE);
      }
      helicityJpsiLam(Pmu1, Pmu2, PJ, Pproton, Pkaon);
      Pmu1 = new TLorentzVector( muplusPX, muplusPY, muplusPZ, muplusPE);
      Pmu2 = new TLorentzVector( muminusPX, muminusPY, muminusPZ, muminusPE);      
      Pproton = new TLorentzVector( protonPX, protonPY, protonPZ, protonPE);
      Pkaon = new TLorentzVector( kaonPX, kaonPY, kaonPZ, kaonPE);
        
      helicityZK(Pmu1, Pmu2, Pproton, Pkaon);
      Pmu1 = new TLorentzVector( muplusPX, muplusPY, muplusPZ, muplusPE);
      Pmu2 = new TLorentzVector( muminusPX, muminusPY, muminusPZ, muminusPE);      
      Pproton = new TLorentzVector( protonPX, protonPY, protonPZ, protonPE);
      Pkaon = new TLorentzVector( kaonPX, kaonPY, kaonPZ, kaonPE);        
        helicityTwoFrame(Pmu1, Pmu2, Pproton, Pkaon); 
        Double_t Vdlz[15];
        Vdlz[0] = xmkp;
        Vdlz[1] = xcostheta;
        Vdlz[2] = xcostheta1;
        Vdlz[3] = xcostheta2;
        Vdlz[4] = xphi1*ID;
        Vdlz[5] = xphi2*ID;
        Vdlz[6] = xcosthetaB;
        Vdlz[7] = xcosthetaZ;
        Vdlz[8] = xcosthetaPsi;
        Vdlz[9] = xphiZ*ID;
        Vdlz[10] = xphiPsi*ID;
        Vdlz[11] = xphiMu*ID;
        Vdlz[12] = xcosthetap;
        Vdlz[13] = xalphaMu*ID;
        Vdlz[14] = xmjpsip;                             
        wdp = sig->evaluate(Vdlz); 
      }
//        std::cout << xmkp << " " << wdp << std::endl;
      }
    }
    newBranch->Fill();
    newBranch1->Fill();
    newBranch2->Fill();
       
  }
   fChain->Write("",TObject::kOverwrite);
   flhcbmc->Close();      
    std::cout << "AmplitudeFit() successfully completed!" << std::endl;
 TDatime FinishTime;
 std::cout << "Time(finish)  " << FinishTime.GetYear() << "." << FinishTime.GetMonth() << "." << FinishTime.GetDay() << "    " << FinishTime.GetHour() << ":" << FinishTime.GetMinute() << ":" << FinishTime.GetSecond() << std::endl;
 return 0;            
}

