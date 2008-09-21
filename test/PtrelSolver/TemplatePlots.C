void TemplatePlots(const char * filename, Int_t nbins)
{
    gROOT->SetStyle("Plain");

    TFile * file = new TFile(filename);

    MakePlots(file, "n", nbins);

    MakePlots(file, "ntag", nbins, "TCL");
    MakePlots(file, "ntag", nbins, "TCM");
    MakePlots(file, "ntag", nbins, "TCT");

    MakePlots(file, "p", nbins);

    MakePlots(file, "ptag", nbins, "JPL");
    MakePlots(file, "ptag", nbins, "JPM");
    MakePlots(file, "ptag", nbins, "JPT");
}


void MakePlots(TFile * file, const char * sample, Int_t nbins, char * tag = 0, Int_t columns = 3)
{
    char hName[256], clName[256], bName[256];

    sprintf(bName, "%s_%s", sample, tag);
    sprintf(hName, "TemplatePlots %s-sample %s-tag", sample, tag);   

    TCanvas * c1 = new TCanvas(bName, hName);

    c1->Divide(columns, ceil((Double_t)(2*nbins)/columns));

    for (Int_t i=0; i<nbins; i++)
    {
        c1->cd(i+1);
        SetName(hName, sample, "pT", tag, i+1);
        SetName(clName, sample, "pT_cl", tag, i+1);
        SetName(bName, sample, "pT_b", tag, i+1);
        MakePlot(file, hName, clName, bName);
    }

    for (Int_t i=0; i<nbins; i++)
    {
        c1->cd(i+nbins+1);
        SetName(hName, sample, "eta", tag, i+1);
        SetName(clName, sample, "eta_cl", tag, i+1);
        SetName(bName, sample, "eta_b", tag, i+1);
        MakePlot(file, hName, clName, bName);
    }
}


void SetName(char * name, const char * sample, const char * dependency, const char * tag, Int_t i)
{
    if (tag)
    {
        sprintf(name, "%s_%s_%s_%d", sample, dependency, tag, i);
    }
    else
    {
        sprintf(name, "%s_%s_%d", sample, dependency, i);
    }
    std::cout << "Processing : " << name << std::endl;
}


void MakePlot(TFile * file, const char * hName, const char * clName, const char * bName)
{
    char name[256];

    file->cd("templates");

    sprintf(name, "template_%s", clName);
    TH1D * clTemplate = (TH1D*) gDirectory->Get(name);
    sprintf(name, "template_%s", bName);
    TH1D * bTemplate = (TH1D*) gDirectory->Get(name);

    Double_t clNorm = clTemplate->Integral();
    Double_t bNorm = bTemplate->Integral();

    clTemplate->GetListOfFunctions()->Delete();
    clTemplate->SetTitle(hName);
    clTemplate->SetStats(false);
    clTemplate->DrawNormalized();
    bTemplate->GetListOfFunctions()->Delete();
    bTemplate->DrawNormalized("same");

    file->cd("functions");

    sprintf(name, "function_%s", clName);
    TF1 * function = (TF1*) gDirectory->Get(name);
    TH1 * clFunction = function->GetHistogram();
    clFunction->SetLineColor(kRed);
    clFunction->Scale(1./clNorm);
    clFunction->Draw("csame");

    sprintf(name, "function_%s", bName);
    TF1 * function = (TF1*) gDirectory->Get(name);
    TH1 * bFunction = function->GetHistogram();
    bFunction->SetLineColor(kBlue);
    bFunction->Scale(1./bNorm);
    bFunction->Draw("csame");

}

