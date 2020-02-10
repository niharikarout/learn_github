//#ifndef __CINT__
#include "RooGlobalFunc.h"
//#endif
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooBifurGauss.h"
#include "RooAddModel.h"
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "RooBinning.h"
#include "TH1.h"
#include "TH2.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooLandau.h"
#include "TChain.h"
#include<cmath>
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooSuperCategory.h"
#include "RooSimultaneous.h"
#include "RooNLLVar.h"
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace RooFit ;
using namespace std;
//int main(){

void mbcfit(){


/*******************Fit Variables***********************************/

RooRealVar mbc("mbc","M_{bc} (GeV/c^{2})",5.27,5.29);
//RooRealVar md0("md0","md0", 1.70, 2.1);
 
/*******************Input root file**********************************/

TChain* chain=new TChain("bp3");
chain->Add("signal_bcs.root");

Double_t  o_de, o_md0, o_mbc, o_sig;
Int_t nevt=(int)chain->GetEntries();

//chain->SetBranchAddress("B_deltae",&o_de);
chain->SetBranchAddress("Mbc",&o_mbc);
//chain->SetBranchAddress("B_D0_M",&o_md0);
chain->SetBranchAddress("isSignal",&o_sig);
  
RooDataSet* data=new RooDataSet("data","data",RooArgSet(mbc));

//Loading data 

for(int i=0;i<nevt;i++) {
  chain->GetEntry(i);
  mbc.setVal(o_mbc);
  

  if(o_sig ==1 && o_mbc < 5.29 && o_mbc > 5.27)
    data->add(RooArgSet(mbc));
	

}


/*****************************Fit***********************/

// --- Build Gaussian signal PDF ---
RooRealVar sigmean("sigmean","B^{#pm} mass",5.279,5.27,5.29) ;
RooRealVar sigwidth("sigwidth","B^{#pm} width",0.0025,0.001,0.01) ;
RooGaussian gauss("gauss","gaussian PDF",mbc,sigmean,sigwidth) ;
gauss.fitTo(*data) ;

RooPlot* mbcframe = mbc.frame(40) ;
data->plotOn(mbcframe) ;
gauss.plotOn(mbcframe) ;
//gauss.plotOn(mbcframe,Components(argus),LineStyle(kDashed)) ;
gauss.paramOn(mbcframe,data);
Double_t chisq = mbcframe->chiSquare();
RooHist* hpull = mbcframe->pullHist() ;
RooPlot* frame3 = mbc.frame(Title("Pull Distribution")) ;
  frame3->addPlotable(hpull,"P") ;

TCanvas* c1 = new TCanvas() ;
 //TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
//pad1->Draw();             // Draw the upper pad: pad1
 //  pad1->cd();  
  mbcframe->Draw() ;
  

  //c1->cd();          // Go back to the main canvas before defining pad2
  // TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
//pad2->Draw();
//   pad2->cd(); 
// frame3->Draw() ;

cout << "chi2 mbc=" << chisq << endl;



//c1->SaveAs("mbc.root");
//c2->SaveAs("md0_kpi_DR2.root");

//return 0;
}

