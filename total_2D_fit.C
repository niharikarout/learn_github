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

void total_2D_fit(){

/*******************Fit Variables***********************************/

RooRealVar deltae("deltae","#DeltaE (GeV)", -0.13, 0.18);
RooRealVar fbdt("fbdt","C'", 0., 1.);
 
/*******************Input root file**********************************/

TFile *f1 = new TFile("/home/niharika/BELLEII/analysis/B2BII/NewRelease/Dpi/BCS/MCFiles/Transformed/FBDT_generic_stream1.root");
TNtuple *tree1 = (TNtuple*)f1->Get("bp3");
//tree1->AddFriend("variables","/home/niharika/BELLEII/analysis/B2BII/NewRelease/Dpi/expert/expert_generic.root");
   
Double_t de;
tree1->SetBranchAddress("deltaE", &de);
Double_t bdt;
tree1->SetBranchAddress("FBDT_transformed", &bdt);

Int_t nevt1=(int)tree1->GetEntries();
  
RooDataSet* data=new RooDataSet("data","data",RooArgSet(deltae,fbdt));
TH1F *h1 = new TH1F("h1","#DeltaE", 100,0,1.);

//Loading data 

for(int i=0;i<nevt1;i++) {
  tree1->GetEntry(i);
  deltae.setVal(de);
  fbdt.setVal(bdt);

    data->add(RooArgSet(deltae,fbdt));
    //if(sig==1 || (mcErrors & 4) != 0)h1->Fill(sig);
    }



/*****************************Fit signal*************************************************/

RooRealVar meande("Mean_{#DeltaE}","mean gaussian",0.0, -0.005, 0.005);
RooRealVar f_de("f_de", "f", 1., 0., 10.);

RooRealVar frac("Frac1", "frac", 0.759);//, 0., 1.);

RooRealVar sigma1("#sigma1", "sigma", 0.00949);//, 0., 0.1);
RooFormulaVar sigma11("sigma11", "@0*@1", RooArgList(sigma1,f_de));
RooGaussian gauss1("gauss1", "gauss1", deltae, meande, sigma11);  //Gauss1  signal for DeltaE

RooRealVar sigma2("#sigma2", "sigma2", 0.0164);//, 0., 0.1);
RooGaussian gauss2("gauss2", "gauss2", deltae, meande, sigma2); // Gauss2 signal for DeltaE

RooAddPdf sig_de("sig_de","Double gaussian", RooArgList(gauss1,gauss2),frac);	

RooRealVar frac2("Frac2","fraction of the main gaussian",0.9590);//,0.0,1.0);
  
RooRealVar sdm2("#sigma_{L}","width 1_2", 0.0870);//, 0.01, 0.1);
RooRealVar sdm3("#sigma_{R}","width 1_3", 0.0470);//, 0.01, 0.05);
RooBifurGauss gauss3("gauss3","bifurcated gaussian PDF", deltae, meande, sdm2, sdm3); // Gauss3 signal for DeltaE

RooAddPdf signal("signal","signal",RooArgList(sig_de,gauss3),frac2);

/*************************************FBDT signal*******************************************/

RooRealVar a0("a0", "a0", 0.0011, -10., 10.);

RooChebychev sig_fbdt("sig_fbdt","signal_fbdt",fbdt,RooArgSet(a0)) ;

/*******************************************************************************************/

RooProdPdf pdf_signal("pdf_signal","pdf_signal",RooArgList(signal, sig_fbdt));  //2D signal PDF

/****************************************DK component****************************************/

RooRealVar diff("#diff", "diff", 0.05);//, 0., 0.05);
RooFormulaVar mean_de("mean_de", "@0-@1", RooArgList(meande,diff));


RooGaussian gauss1_DK("gauss1_DK", "gauss1_DK", deltae, mean_de, sigma1);  //Gauss1  signal for DeltaE DK component

RooGaussian gauss2_DK("gauss2_DK", "gauss2_DK", deltae, mean_de, sigma2); // Gauss2 signal for DeltaE DK component

RooAddPdf sig_de_DK("sig_de_DK","Double gaussian", RooArgList(gauss1_DK,gauss2_DK),frac);	
  
RooBifurGauss gauss3_DK("gauss3_DK","bifurcated gaussian PDF", deltae, mean_de, sdm2, sdm3); // Gauss3 signal for DeltaE DK component

RooAddPdf signal_DK("signal_DK","signal",RooArgList(sig_de_DK,gauss3_DK),frac2);

/***************************FBDT DK component *********************************************/

RooRealVar c0("c0", "c0", -0.0911);//, -10., 10.);
RooChebychev bkg1("bkg1","Background1",fbdt,RooArgSet(c0)) ;

/********************************************************************************************/

//RooProdPdf pdf_DK("pdf_DK","pdf_DK",RooArgList(signal_DK, DK_fbdt));  //2D DK component PDF
RooProdPdf pdf_DK("pdf_DK","pdf_DK",RooArgList(signal_DK, bkg1));  //2D DK component PDF
  
/**************************************Fit qqbar***********************************************/	

RooRealVar b0("b0", "b0", -0.3063);//, -10., 10.);

RooChebychev bkg("bkg","Background",deltae,RooArgSet(b0)) ;

/********************************** Fit qqbar FBDT********************************************/

RooRealVar d0("d0", "d0", -1.55687);//, -10., 10.);
RooRealVar d1("d1", "d1", 0.686);//, -10., 10.);
RooRealVar d2("d2", "d2", -0.2019);//, -10., 10.);
RooChebychev bkg2("bkg2","Background2",fbdt,RooArgSet(d0,d1,d2)) ;

RooRealVar lambda2("lambda2", "slope2", -19.35);//, -100., 100.);
RooExponential expo2("expo2", "exponential PDF2", fbdt, lambda2);

RooRealVar frac_bkg2("frac_bkg2", "frac_bkg2", 0.654);//, 0., 1.);
RooAddPdf qq_fbdt("qq_fbdt","qq_fbdt",RooArgList(bkg2,expo2),RooArgList(frac_bkg2));


/*********************************************************************************************/

RooProdPdf pdf_qq("pdf_qq","pdf_qq",RooArgList(bkg, qq_fbdt));  //2D qqbar component PDF

/***************************Fit BBbar bkg*******************************/

RooRealVar lambda("lambda", "slope", -11.091);//, -15., 15.);
RooExponential expo("expo", "exponential PDF", deltae, lambda);

/*******************************BBbar bkg FBDT********************************/

RooRealVar e0("e0", "e0", -0.3335);//, -10., 10.);
RooChebychev bkg3("bkg3","Background3",fbdt,RooArgSet(e0)) ;

RooRealVar lambda3("lambda3", "slope3", -18.65);//, -100., 100.);
RooExponential expo3("expo3", "exponential PDF3", fbdt, lambda3);

RooRealVar frac_bkg3("frac_bkg3", "frac_bkg3", 0.975);//, 0., 1.);
RooAddPdf bb_fbdt("bb_fbdt","bb_fbdt",RooArgList(bkg3,expo3),RooArgList(frac_bkg3));

/**************************************************************************************/

RooProdPdf pdf_bb("pdf_bb","pdf_bb",RooArgList(expo, bb_fbdt));  //2D bbbar bkg component PDF

/***************************************************************************************/	 	 
		 
RooRealVar nsig("nsig", "nsig", 5000, -100., 100000.0);
RooRealVar nbb("nbb", "nbb", 2000, -100., 10000.0);
RooRealVar nqq("nqq", "nqq", 5000, -100., 100000.0);
RooRealVar ndk("nsig_DK", "nsig_DK", 39.);//, -5., 100.0);

RooAddPdf sum("sum","sum",RooArgList(pdf_signal,pdf_bb,pdf_qq,pdf_DK),RooArgList(nsig, nbb,nqq,ndk));

sum.fitTo(*data);

//signal region
deltae.setRange("signal1",-0.05,0.05);
fbdt.setRange("signal2",0.65,1.);

RooPlot* deframe = deltae.frame(50) ;
data->plotOn(deframe,CutRange("signal2")) ;
sum.plotOn(deframe,ProjectionRange("signal2")) ;
sum.paramOn(deframe,data);
deframe->getAttText()->SetTextSize(0.06) ;
Double_t chisq = deframe->chiSquare();

TPaveText *box= new TPaveText(0.4, 0.85, 0.1, 0.9,"BRNDC");
box->SetFillColor(10);
box->SetBorderSize(0);
box->SetTextAlign(12);
box->SetTextSize(0.06);
box->SetFillStyle(1002);
TText *text = 0;
Char_t buf[30];
sprintf( buf,  "#chi^{2}/ndf = %f", chisq );
text = box->AddText( buf );
deframe->addObject(box) ;

//frame styles

  deframe->GetXaxis()->SetTitleSize(0.07);
  //deframe->GetXaxis()->SetTitleOffset(1.2);
  deframe->GetXaxis()->SetLabelSize(0.07);
  deframe->GetYaxis()->SetTitleSize(0.07);
  deframe->GetYaxis()->SetLabelSize(0.07);
  
  

RooHist* hpull = deframe->pullHist() ;
RooPlot* frame3 = deltae.frame(Title("Pull Distribution")) ;
frame3->GetYaxis()->SetTitle("Pull") ; 
frame3->GetYaxis()->SetTitleSize(0.2) ;
frame3->GetYaxis()->SetNdivisions(504) ;
frame3->GetYaxis()->SetLabelSize(0.2) ;
frame3->GetXaxis()->SetTitleSize(0.0) ;
frame3->GetXaxis()->SetLabelSize(0.0) ;
frame3->GetYaxis()->SetTitleOffset(0.4) ;
frame3->addPlotable(hpull,"x0 P E1") ; 
frame3->SetMaximum(5.);  
frame3->SetMinimum(-5.);
  
//sig2.plotOn(deframe,Components(signal_de),LineColor(kMagenta),LineStyle(kDashed)) ;
sum.plotOn(deframe,Components(pdf_signal),LineColor(kMagenta),LineStyle(kDashed),ProjectionRange("signal2")) ;
sum.plotOn(deframe,Components(pdf_bb),LineColor(kBlue),LineStyle(kDashed),ProjectionRange("signal2")) ;
sum.plotOn(deframe,Components(pdf_qq),LineColor(kCyan),LineStyle(kDashed),ProjectionRange("signal2")) ;
sum.plotOn(deframe,Components(pdf_DK),LineColor(kGreen),LineStyle(kDashed),ProjectionRange("signal2")) ;

TCanvas* c1 = new TCanvas() ;
deframe->GetYaxis()->SetTitleOffset(1.2) ;
TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
pad1->Draw();             // Draw the upper pad: pad1
pad1->cd();  
deframe->Draw() ;
  
c1->cd();          // Go back to the main canvas before defining pad2
TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
pad2->Draw();
pad2->cd(); 
frame3->Draw() ;
c1->Update();
TLine *l=new TLine(-0.13,0.0,0.18,0.0);
l->SetLineColor(kBlue);
l->SetLineWidth(3);
l->Draw();

//cout << "chi2 de=" << chisq << endl;
//cout << h1->GetEntries() << endl;

/***************************************************************************/

RooPlot* fbdtframe = fbdt.frame(50) ;
data->plotOn(fbdtframe, CutRange("signal1"),RooFit::Name("data")) ;
sum.plotOn(fbdtframe,ProjectionRange("signal1")) ;
//sum.paramOn(fbdtframe,data);
Double_t chisq_fbdt = fbdtframe->chiSquare();


TPaveText *box= new TPaveText(0.4, 0.85, 0.1, 0.9,"BRNDC");
box->SetFillColor(10);
box->SetBorderSize(0);
box->SetTextAlign(12);
box->SetTextSize(0.05F);
box->SetFillStyle(1002);
TText *text = 0;
Char_t buf[30];
sprintf( buf,  "#chi^{2}/ndf = %f", chisq_fbdt );
text = box->AddText( buf );
fbdtframe->addObject(box) ;
fbdtframe->getAttText()->SetTextSize(0.06) ;

fbdtframe->GetXaxis()->SetTitleSize(0.07);
  //deframe->GetXaxis()->SetTitleOffset(1.2);
  fbdtframe->GetXaxis()->SetLabelSize(0.07);
  fbdtframe->GetYaxis()->SetTitleSize(0.07);
  fbdtframe->GetYaxis()->SetLabelSize(0.07);

RooHist* hpull1 = fbdtframe->pullHist() ;
RooPlot* frame4 = fbdt.frame(Title("Pull Distribution")) ;
frame4->GetYaxis()->SetTitle("Pull") ; 
frame4->GetYaxis()->SetTitleSize(0.2) ;
frame4->GetYaxis()->SetNdivisions(504) ;
frame4->GetYaxis()->SetLabelSize(0.2) ;
frame4->GetXaxis()->SetTitleSize(0.0) ;
frame4->GetXaxis()->SetLabelSize(0.0) ;
frame4->GetYaxis()->SetTitleOffset(0.4) ;
frame4->addPlotable(hpull1,"x0 P E1") ; 
frame4->SetMaximum(5.);  
frame4->SetMinimum(-5.);
  
//sig2.plotOn(deframe,Components(signal_de),LineColor(kMagenta),LineStyle(kDashed)) ;
sum.plotOn(fbdtframe,Components(pdf_signal),LineColor(kMagenta),LineStyle(kDashed),ProjectionRange("signal1"),RooFit::Name("signal")) ;
sum.plotOn(fbdtframe,Components(pdf_bb),LineColor(kBlue),LineStyle(kDashed),ProjectionRange("signal1"),RooFit::Name("bbbar_bkg")) ;
sum.plotOn(fbdtframe,Components(pdf_qq),LineColor(kCyan),LineStyle(kDashed),ProjectionRange("signal1"),RooFit::Name("qqbar")) ;
sum.plotOn(fbdtframe,Components(pdf_DK),LineColor(kGreen),LineStyle(kDashed),ProjectionRange("signal1"),RooFit::Name("DK_component")) ;

TCanvas* c2 = new TCanvas() ;
fbdtframe->GetYaxis()->SetTitleOffset(1.2) ;
TPad *pad3 = new TPad("pad3", "pad3", 0, 0.3, 1, 1.0);
pad3->Draw();             // Draw the upper pad: pad1
pad3->cd();  
fbdtframe->Draw() ;

TLegend *leg = new TLegend(0.47,0.7,0.8,0.43,"");
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->SetTextSize(0.065);
leg->SetTextFont(62);
leg->AddEntry(fbdtframe->findObject("data"),"Data","p");
leg->AddEntry(fbdtframe->findObject("signal"),"Signal","l");
leg->AddEntry(fbdtframe->findObject("bbbar_bkg"),"B#bar{B} bkg","l");
leg->AddEntry(fbdtframe->findObject("qqbar"),"q#bar{q} bkg","l");
leg->AddEntry(fbdtframe->findObject("DK_component"),"DK component","l");
//fbdtframe->getAttText()->SetTextSize(0.06) ;
leg->Draw();
  
c2->cd();          // Go back to the main canvas before defining pad2
TPad *pad4 = new TPad("pad4", "pad4", 0, 0.05, 1, 0.3);
pad4->Draw();
pad4->cd(); 
frame4->Draw() ;
c2->Update();
TLine *l1=new TLine(0.0,0.0,1.0,0.0);
l1->SetLineColor(kBlue);
l1->SetLineWidth(3);
l1->Draw();

cout << "chi2 de=" << chisq << endl;
cout << "chi2 fbdt=" << chisq_fbdt << endl;
//cout << h1->GetEntries() << endl;


}

