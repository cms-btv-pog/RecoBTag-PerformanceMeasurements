#include "TTbarSFbFitTools.h"

#include "TMath.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"

#include "RooHist.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
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
  std::vector<Float_t> nExp,nExpUnc;
  for(int i=0; i<expTemplates.GetEntriesFast(); i++)
    {
      TH1 *h=((TH1 *)expTemplates.At(i));
      Double_t unc(0);
      nExp[i]=h->IntegralAndError(1,h->GetNbinsX(),unc);
      nExpUnc[i]=unc;
      totalExp += nExp[i];
    }
  //build the pdf components
  RooArgSet expPDFs,expFracs,nonPOIPDFs;
  TString poiName("v"); poiName+=idxOfInterest;
  for(int i=0; i<expTemplates.GetEntriesFast(); i++)
    {
      TH1 *h=(TH1 *)expTemplates.At(i);
      TString title(h->GetTitle());
      TString name("v"); name+= i;
     
      if(i<(expTemplates.GetEntriesFast()-1))
	{
	  float frac(nExp[i]/totalExp);
	  RooRealVar *fracVar = new RooRealVar(name,name,frac,0.5*frac,TMath::Min((float)1.0,(float)2*frac));
	  fracVar->SetTitle(title);
	  expFracs.add(*fracVar);
	}
      
      //PDF
      RooDataHist *itempl = new RooDataHist(name+"_hist",name+"_hist", RooArgList(x), Import(*h));
      RooHistPdf *ipdf    = new RooHistPdf (name+"_pdf", name+"_pdf",  RooArgSet(x), *itempl);
      ipdf->SetTitle(title);
      expPDFs.add(*ipdf);
      if(i!=idxOfInterest) nonPOIPDFs.add(*ipdf);
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

  result.effExp=nExp[idxOfInterest]/totalExp;
  result.effExpUnc=nExpUnc[idxOfInterest]/totalExp;
  Float_t nObs    = fracVar->getVal()*totalExp;
  Float_t nObsUnc = fracVar->getError()*totalExp;  
  result.sf      = nObs/nExp[idxOfInterest];
  result.sfUnc   = TMath::Sqrt(TMath::Power(nObs*nExpUnc[idxOfInterest],2)+TMath::Power(nObsUnc*nExp[idxOfInterest],2))/TMath::Power(nExp[idxOfInterest],2);
  result.eff=result.sf*result.effExp;
  result.effUnc=result.sfUnc*result.effExp;

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
		  RooFit::LineColor(1),
		  RooFit::LineWidth(1),
		  RooFit::MoveToBack());
      

      pdf->plotOn(frame,
		  RooFit::Name("nonpoifit"),
		  RooFit::ProjWData(*data),
		  RooFit::Components(nonPOIPDFs),
		  RooFit::FillColor(kGray),
		  RooFit::LineColor(kGray),
		  //RooFit::Components(poiName+"_*"),
		  //RooFit::FillColor(kOrange+2),
		  //RooFit::LineColor(kOrange+2),
		  RooFit::FillStyle(1001),
		  RooFit::DrawOption("f"),
		  RooFit::MoveToBack());
      frame->Draw();
      frame->GetYaxis()->SetTitleOffset(1.5);
      frame->GetXaxis()->SetTitle(x.GetTitle());
      frame->GetYaxis()->SetRangeUser(0,dataH->GetMaximum()*2.0);

      TLatex *label= new TLatex();
      label->SetNDC();
      label->SetTextFont(42);
      label->SetTextSize(0.04);
      label->DrawLatex(0.18,0.94,"#bf{CMS} #it{simulation}");

      TLegend *leg=new TLegend(0.18,0.7,0.35,0.9);
      leg->AddEntry("data",      "data",                             "p");
      leg->AddEntry("totalexp",  expFracs.find(poiName)->GetTitle(), "l");
      leg->AddEntry("nonpoifit", "others",                           "f");
      leg->SetFillStyle(0);
      leg->SetTextFont(42);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.035);
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
TTbarFracFitterResult_t TTbarFracFitter::fit(TObjArray &passTemplates, TH1F *passDataH,
					     TObjArray &failTemplates, TH1F *failDataH,
					     Int_t idxOfInterest,TString saveResultIn,
					     Float_t lumi)
{
  using namespace RooFit;

  TTbarFracFitterResult_t result; 
  
  //variable to fit
  RooRealVar x("x",passDataH->GetXaxis()->GetTitle(),passDataH->GetXaxis()->GetXmin(),passDataH->GetXaxis()->GetXmax());
  x.setBins(passDataH->GetXaxis()->GetNbins());

  //categories
  RooCategory status("status","status") ;
  status.defineType("pass") ;
  status.defineType("fail");

  //Fix mass values and create a mapped data hist
  std::map<std::string,TH1 *> dataH;
  dataH["pass"]=passDataH;
  dataH["fail"]=failDataH;
  RooDataHist* data = new RooDataHist("data","data",x,status,dataH);

  //compute the expected efficiency
  std::vector<float> expEff(passTemplates.GetEntriesFast(),0.);
  std::vector<float> expEffUnc(passTemplates.GetEntriesFast(),0.);
  for(int i=0; i<passTemplates.GetEntriesFast(); i++)
    {
      TH1 *passH=(TH1 *)passTemplates.At(i);
      TH1 *failH=(TH1 *)failTemplates.At(i);
      if(passH==0 || failH==0) continue;

      Double_t nPassUnc(0);
      Float_t nPass(passH->IntegralAndError(1,passH->GetNbinsX(),nPassUnc)), nFail(failH->Integral());
      Float_t total(nPass+nFail);
      if(total==0) continue;

      expEff[i]=nPass/total;
      expEffUnc[i]=nPassUnc/total;
    }
 
  RooArgSet passPDFs,passCounts,passNonPOIPDFs,failPDFs,failCounts,sfVars,failNonPOIPDFs;
  for(int i=0; i<passTemplates.GetEntriesFast(); i++)
    {
      //build the pass PDF
      TH1 *h=(TH1 *)passTemplates.At(i);
      TString title(h->GetTitle());
      TString name("passv"); name+= i;
      Float_t nExp=h->Integral();
   
      //b-tag efficiency scale factor
      RooRealVar *sfVar= new RooRealVar(Form("sf%d",i),title,1,0,2);      
      sfVars.add(*sfVar);
      if(i!=idxOfInterest) sfVar->setConstant(true);
      
      //npass expected
      RooRealVar *npassVar=new RooRealVar(name,name,nExp,nExp,nExp);
      npassVar->setConstant(true);
      name.ReplaceAll("v","fv");


      //re-scaling of the number of passing jets
      RooFormulaVar *npassFormulaVar = new RooFormulaVar(name,
							 "@0*@1",
							 RooArgSet(*npassVar,*sfVar)
							 );
      passCounts.add( *npassFormulaVar );
      RooDataHist *itempl = new RooDataHist(name+"_hist",name+"_hist", RooArgList(x), Import(*h));
      RooHistPdf *ipdf    = new RooHistPdf (name+"_pdf", name+"_pdf",  RooArgSet(x), *itempl);
      ipdf->SetTitle(title);
      passPDFs.add(*ipdf);
      if(i!=idxOfInterest) passNonPOIPDFs.add(*ipdf);

      //build the fail PDF
      h=(TH1 *)failTemplates.At(i);
      title=h->GetTitle();
      name="failfv"; name+= i;
      nExp=h->Integral();      
      RooFormulaVar *nFailFormulaVar = new RooFormulaVar(name,
							 Form("@1*%f >=0 ? @0*(1-@1*%f)/(@1*%f) : 0",expEff[i],expEff[i],expEff[i]),
							 RooArgSet(*npassFormulaVar,*sfVar)
							 );
      failCounts.add( *nFailFormulaVar );
      itempl = new RooDataHist(name+"_hist",name+"_hist", RooArgList(x), Import(*h));
      ipdf    = new RooHistPdf (name+"_pdf", name+"_pdf",  RooArgSet(x), *itempl);
      ipdf->SetTitle(title);
      failPDFs.add(*ipdf);
      if(i!=idxOfInterest) failNonPOIPDFs.add(*ipdf);
    }

  //create the final pdfs
  RooSimultaneous *pdf = new RooSimultaneous("pdf","pdf",status);
  RooAddPdf *finalPassPDF = new RooAddPdf("passpdf","passpdf", passPDFs, passCounts);
  pdf->addPdf(*finalPassPDF,"pass");
  RooAddPdf *finalFailPDF = new RooAddPdf("failpdf","failpdf", failPDFs, failCounts);
  pdf->addPdf(*finalFailPDF,"fail");
  
  //fit (using minos has the same behavior of profiling all variables except the parameter of interest)
  RooRealVar *sfVar=(RooRealVar *)sfVars.find(Form("sf%d",idxOfInterest));
  RooAbsReal *ll = pdf->createNLL(*data,Extended(true));
  RooMinuit minuit(*ll);
  minuit.setErrorLevel(0.5);
  minuit.migrad();
  minuit.hesse();
  RooArgSet poi(*sfVar);
  result.minuitStatus = minuit.minos(poi);
  
  //save result
  result.effExp    = expEff[idxOfInterest];
  result.effExpUnc = expEffUnc[idxOfInterest];
  result.sf        = sfVar->getVal();
  result.sfUnc     = sfVar->getError();
  result.eff       = result.sf*result.effExp;
  result.effUnc    = TMath::Sqrt( TMath::Power(result.sfUnc*result.effExp,2)+TMath::Power(result.sf*result.effExpUnc,2) );
  
  //
  if(saveResultIn!="")
    {
      TCanvas *canvas=new TCanvas("c","c",1000,500);
      canvas->SetTopMargin(0);
      canvas->SetBottomMargin(0);
      canvas->SetRightMargin(0);
      canvas->SetLeftMargin(0);
      canvas->Divide(2,1);
    
      for(int i=0; i<2; i++)
	{
	  std::string statusName(i==0 ? "pass" : "fail");
	  
	  TPad *c=(TPad *)canvas->cd(i+1);
	  c->SetTopMargin(0);
	  c->SetBottomMargin(0);
	  c->SetRightMargin(0);
	  c->SetLeftMargin(0);
    
	  TPad *p1 = new TPad("p1","p1",0.0,0.85,1.0,0.0);
	  p1->SetRightMargin(0.05);
	  p1->SetLeftMargin(0.12);
	  p1->SetTopMargin(0.01);
	  p1->SetBottomMargin(0.1);
	  p1->Draw();
	  p1->cd();
	  RooPlot *frame=x.frame();
	  data->plotOn(frame,RooFit::Name("data"),RooFit::Cut(("status==status::"+statusName).c_str()));
	  pdf->plotOn(frame,
		      RooFit::Slice(status,statusName.c_str()),
		      RooFit::Name("totalexp"),
		      RooFit::ProjWData(*data),
		      RooFit::LineColor(1),
		      RooFit::LineWidth(1),
		      RooFit::FillStyle(1001),
		      RooFit::FillColor(0),
		      RooFit::MoveToBack());
	  pdf->plotOn(frame,
		      RooFit::Slice(status,statusName.c_str()),
		      RooFit::Name("nonpoifit"),
		      RooFit::ProjWData(*data),
		      RooFit::Components(i==0 ? passNonPOIPDFs : failNonPOIPDFs),
		      RooFit::FillColor(kGray),
		      RooFit::LineColor(kGray),
		      RooFit::FillStyle(1001),
		      RooFit::DrawOption("f"),
		      RooFit::MoveToBack());
	  frame->Draw();
	  frame->GetYaxis()->SetTitleOffset(1.5);
	  frame->GetXaxis()->SetTitle(x.GetTitle());
	  Float_t ymax(TMath::Max(passDataH->GetMaximum(),failDataH->GetMaximum())*1.5);
	  Float_t xmin(frame->GetXaxis()->GetXmin()),xmax(frame->GetXaxis()->GetXmax());
	  frame->GetYaxis()->SetRangeUser(0,ymax);

	  TLatex *label= new TLatex();
	  label->SetNDC();
	  label->SetTextFont(42);
	  label->SetTextSize(0.04);
	  if(i==0)
	    {
	      label->DrawLatex(0.18,0.94,"#bf{CMS} #it{preliminary}");
	      label->DrawLatex(0.70,0.94,Form("13 TeV (%3.1f fb^{-1})",lumi));
	      label->DrawLatex(0.18,0.88,Form("SF=%3.2f#pm%3.2f",sfVar->getVal(),sfVar->getError()));
	    }
	  label->DrawLatex(0.70,0.88,i==0 ? "#it{pass category}" : "#it{fail category}");

	  TLine *line=new TLine();
	  line->SetLineStyle(9);
	  for(Int_t k=1; k<=3; k++)
	    {
	      if(k<3)
		line->DrawLine(xmin+k*(xmax-xmin)/3.,0,xmin+k*(xmax-xmin)/3.,ymax/1.5);
	      label->DrawLatex((2*(k-1)+1.)/6.,0.7,Form("#scale[0.6]{#it{=%d jets}}",k+1));
	    }

	  if(i==0)
	    {
	      TLegend *leg=new TLegend(0.18,0.76,0.6,0.88);
	      leg->SetNColumns(3);
	      leg->AddEntry("data",      "data",            "p");
	      leg->AddEntry("totalexp",  sfVar->GetTitle(), "f");
	      leg->AddEntry("nonpoifit", "others",          "f");
	      leg->SetFillStyle(0);
	      leg->SetTextFont(42);
	      leg->SetBorderSize(0);
	      leg->SetTextSize(0.035);
	      leg->Draw();
	    }
    
	  c->cd();
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
	}

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
