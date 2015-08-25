#include "TTbarSFbFitTools.h"

#include "TMath.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"

#include "RooHist.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"


using namespace std;

//
TTbarFracFitter::TTbarFracFitter()
{
  using namespace RooFit;

  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().getStream(0).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(ObjectHandling);
  RooMsgService::instance().getStream(1).removeTopic(DataHandling);
  RooMsgService::instance().getStream(1).removeTopic(Fitting);
  RooMsgService::instance().getStream(1).removeTopic(Plotting);
  RooMsgService::instance().getStream(0).removeTopic(InputArguments);
  RooMsgService::instance().getStream(1).removeTopic(InputArguments);
  RooMsgService::instance().getStream(0).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
}

//
TTbarFracFitterResult_t TTbarFracFitter::fit(TObjArray &expTemplates, TH1F *dataH,Int_t idxOfInterest,TString saveResultIn)
{
  using namespace RooFit;

  TTbarFracFitterResult_t result; 
  
  //variable to fit
  RooRealVar x("x",dataH->GetXaxis()->GetTitle(),dataH->GetXaxis()->GetXmin(),dataH->GetXaxis()->GetXmax());
  x.setBins(dataH->GetXaxis()->GetNbins());

  //data
  RooDataHist *data=new RooDataHist("data","data", RooArgSet(x), Import(*dataH) );

  //total events expected
  Float_t totalExp(0);
  for(int i=0; i<expTemplates.GetEntriesFast(); i++)
    totalExp += ((TH1 *)expTemplates.At(i))->Integral();

  //build the pdf components
  RooArgSet expPDFs,expFracs;
  TString poiName("v"); poiName+=idxOfInterest;
  for(int i=0; i<expTemplates.GetEntriesFast(); i++)
    {
      TH1 *h=(TH1 *)expTemplates.At(i);
      TString title(h->GetTitle());
      TString name("v"); name+= i;

      //expected yields for parameter of interest
      Double_t nExpUnc(0);
      float nExp=h->IntegralAndError(1,h->GetXaxis()->GetNbins(),nExpUnc);
      if(name==poiName)
	{
	  result.nExp = nExp;
	  result.nExpUnc = nExpUnc;
	}
      
      if(i<(expTemplates.GetEntriesFast()-1))
	{
	  float frac(nExp/totalExp);
	  RooRealVar *fracVar = new RooRealVar(name,name,frac,0.5*frac,TMath::Min((float)1.0,(float)2*frac));
	  fracVar->SetTitle(title);
	  expFracs.add(*fracVar);
	}
      
      //PDF
      RooDataHist *itempl = new RooDataHist(name+"_hist",name+"_hist", RooArgList(x), Import(*h));
      RooHistPdf *ipdf    = new RooHistPdf (name+"_pdf", name+"_pdf",  RooArgSet(x), *itempl);
      ipdf->SetTitle(title);
      expPDFs.add(*ipdf);
    }

  //create the final pdf
  RooAddPdf *pdf=new RooAddPdf("pdf","pdf", expPDFs, expFracs);

  //fit (using minos has the same behavior of profiling all variables except the parameter of interest)
  RooAbsReal *ll = pdf->createNLL(*data);
  RooMinuit minuit(*ll);
  minuit.setErrorLevel(0.5);
  minuit.migrad();
  minuit.hesse();
  RooArgSet poi(*expFracs.find(poiName));
  result.minuitStatus = minuit.minos(poi);

  //save result
  RooRealVar *fracVar=(RooRealVar *)expFracs.find(poiName);
  result.nObs    = fracVar->getVal()*totalExp;
  result.nObsUnc = fracVar->getError()*totalExp;     
  result.sf      = result.nExp<=0 ? -1 : result.nObs/result.nExp;
  result.sfUnc   = result.nExp<=0 ? -1 : TMath::Sqrt(TMath::Power(result.nObs*result.nExpUnc,2)+TMath::Power(result.nObsUnc*result.nExp,2))/TMath::Power(result.nExp,2);
  
  //
  if(saveResultIn!="")
    {
      TCanvas *canvas=new TCanvas("c","c",500,500);
      TPad *p1 = new TPad("p1","p1",0.0,0.85,1.0,0.0);
      p1->SetRightMargin(0.05);
      p1->SetLeftMargin(0.12);
      p1->SetTopMargin(0.01);
      p1->SetBottomMargin(0.1);
      p1->Draw();
      p1->cd();
      RooPlot *frame=x.frame();
      data->plotOn(frame,RooFit::Name("data"));
      pdf->plotOn(frame,
		  RooFit::Name("totalexp"),
		  RooFit::ProjWData(*data),
		  RooFit::LineColor(kBlue),
		  RooFit::LineWidth(2),
		  RooFit::MoveToBack());
      pdf->plotOn(frame,
		  RooFit::Name("poifit"),
		  RooFit::ProjWData(*data),
		  RooFit::Components(poiName+"_*"),
		  RooFit::FillColor(kOrange+2),
		  RooFit::LineColor(kOrange+2),
		  RooFit::FillStyle(1001),
		  RooFit::DrawOption("f"),
		  RooFit::MoveToBack());
      frame->Draw();
      frame->GetYaxis()->SetTitleOffset(1.5);
      frame->GetXaxis()->SetTitle(x.GetTitle());

      TLatex *label= new TLatex();
      label->SetNDC();
      label->SetTextFont(42);
      label->SetTextSize(0.04);
      label->DrawLatex(0.18,0.94,"#bf{CMS} #it{simulation}");

      TLegend *leg=new TLegend(0.65,0.35,0.95,0.52);
      leg->AddEntry("data",       "Pseudo data",      "p");
      leg->AddEntry("totalexp",   "Total",     "l");
      leg->AddEntry("poifit",     expFracs.find(poiName)->GetTitle(),  "f");
      leg->SetFillStyle(0);
      leg->SetTextFont(43);
      leg->SetTextSize(14);
      leg->SetBorderSize(0);
      leg->Draw();

      canvas->cd();
      TPad *p2 =new TPad("p2","p2",0.0,0.86,1.0,1.0);
      p2->SetBottomMargin(0.05);
      p2->SetRightMargin(0.05);
      p2->SetLeftMargin(0.12);
      p2->SetTopMargin(0.05);
      p2->Draw();
      p2->cd();
      
      RooHist *hpull = frame->pullHist();
      RooPlot *pullFrame = x.frame();
      //hpull->plotOn(pullFrame);
      pullFrame->addPlotable((RooPlotable *)hpull,"P") ;      
      pullFrame->Draw();
      pullFrame->GetYaxis()->SetTitle("Pull");
      pullFrame->GetYaxis()->SetTitleSize(0.2);
      pullFrame->GetYaxis()->SetLabelSize(0.2);
      pullFrame->GetXaxis()->SetTitleSize(0);
      pullFrame->GetXaxis()->SetLabelSize(0);
      pullFrame->GetYaxis()->SetTitleOffset(0.15);
      pullFrame->GetYaxis()->SetNdivisions(4);
      pullFrame->GetYaxis()->SetRangeUser(-3.1,3.1);
      pullFrame->GetXaxis()->SetTitleOffset(0.8);

      canvas->cd();
      TPad *p3 = new TPad("p3","p3",0.7,0.57,0.95,0.82);
      p3->SetRightMargin(0.05);
      p3->SetLeftMargin(0.12);
      p3->SetTopMargin(0.008);
      p3->SetBottomMargin(0.2);
      p3->Draw();
      p3->cd();
      RooPlot *frame2=((RooRealVar *)expFracs.find(poiName))->frame();
      ll->plotOn(frame2,RooFit::ShiftToZero());
      frame2->Draw();
      frame2->GetYaxis()->SetRangeUser(0,12);
      frame2->GetXaxis()->SetRangeUser(165,180);
      frame2->GetYaxis()->SetNdivisions(3);
      frame2->GetXaxis()->SetNdivisions(3);
      frame2->GetXaxis()->SetTitle(expFracs.find(poiName)->GetTitle() + TString(" fraction"));
      frame2->GetYaxis()->SetTitle("LL");
      frame2->GetYaxis()->SetTitleOffset(1.5);
      frame2->GetXaxis()->SetTitleSize(0.08);
      frame2->GetXaxis()->SetLabelSize(0.08);
      frame2->GetYaxis()->SetTitleSize(0.08);
      frame2->GetYaxis()->SetLabelSize(0.08);
      
      canvas->Modified();
      canvas->Update();
      canvas->SaveAs(saveResultIn+".png");
      canvas->SaveAs(saveResultIn+".pdf");
      canvas->Delete();
    }


  //free memory
  delete ll;

  //all done
  return result;
}

//
TTbarFracFitter::~TTbarFracFitter()
{
  
}
