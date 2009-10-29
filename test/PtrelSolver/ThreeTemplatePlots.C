void ThreeTemplatePlots(const char * filename, Int_t nbins, Int_t columns)
{
    gROOT->SetStyle("Plain");

    TFile * file = new TFile(filename);

    MakePlots(file, "n", nbins, columns);
    MakePlots(file, "ntag", nbins, columns);
    MakePlots(file, "nnoTag", nbins, columns);

    /*MakePlots(file, "ntag", nbins, columns, "TCL");
    MakePlots(file, "ntag", nbins, columns, "TCM");
    MakePlots(file, "ntag", nbins, columns, "TCT");

    MakePlots(file, "ntag", nbins, columns, "JPL");
    MakePlots(file, "ntag", nbins, columns, "JPM");
    MakePlots(file, "ntag", nbins, columns, "JPT");

    MakePlots(file, "ntag", nbins, columns, "JBPL");
    MakePlots(file, "ntag", nbins, columns, "JBPM");
    MakePlots(file, "ntag", nbins, columns, "JBPT");

    MakePlots(file, "ntag", nbins, columns, "SVM");
    MakePlots(file, "ntag", nbins, columns, "SVT");

    MakePlots(file, "ntag", nbins, columns, "CSVL");
    MakePlots(file, "ntag", nbins, columns, "CSVM");
    MakePlots(file, "ntag", nbins, columns, "CSVT");*/

    // MakePlots(file, "ntag", nbins, columns, "SMT");

    MakePlots(file, "p", nbins, columns);
    MakePlots(file, "ptag", nbins, columns);
    MakePlots(file, "pnoTag", nbins, columns);

    /*MakePlots(file, "ptag", nbins, columns, "TCL");
    MakePlots(file, "ptag", nbins, columns, "TCM");
    MakePlots(file, "ptag", nbins, columns, "TCT");

    MakePlots(file, "ptag", nbins, columns, "JPL");
    MakePlots(file, "ptag", nbins, columns, "JPM");
    MakePlots(file, "ptag", nbins, columns, "JPT");

    MakePlots(file, "ptag", nbins, columns, "JBPL");
    MakePlots(file, "ptag", nbins, columns, "JBPM");
    MakePlots(file, "ptag", nbins, columns, "JBPT");

    MakePlots(file, "ptag", nbins, columns, "SVM");
    MakePlots(file, "ptag", nbins, columns, "SVT");

    MakePlots(file, "ptag", nbins, columns, "CSVL");
    MakePlots(file, "ptag", nbins, columns, "CSVM");
    MakePlots(file, "ptag", nbins, columns, "CSVT");*/

    // MakePlots(file, "ptag", nbins, columns, "SMT");

}


void MakePlots(TFile * file, const char * sample, Int_t nbins, Int_t columns = 3, char * tag = 0)
{
    char hName[256], lName[256], cName[256], bName[256];

    if (tag)
    {
        sprintf(bName, "%s_%s", sample, tag);
        sprintf(hName, "TemplatePlots_%s_%s", sample, tag);
    }
    else
    {
        sprintf(bName, "%s", sample);
        sprintf(hName, "TemplatePlots_%s", sample);
    }

    char pName[256];
    sprintf(pName, "%s_InclusiveMu5Pt50.png", hName);

    TCanvas * c1 = new TCanvas(bName, hName);

    c1->Divide(columns, ceil((Double_t)(2*nbins)/columns));

    for (Int_t i=0; i<nbins; i++)
    {
        c1->cd(i+1);
        SetName(hName, sample, "pT", tag, i+1);
        SetName(lName, sample, "pT_l", tag, i+1);
        SetName(cName, sample, "pT_c", tag, i+1);
        SetName(bName, sample, "pT_b", tag, i+1);
        MakePlot(file, hName, lName, cName, bName);
    }

    for (Int_t i=0; i<nbins; i++)
    {
        c1->cd(i+nbins+1);
        SetName(hName, sample, "eta", tag, i+1);
        SetName(lName, sample, "eta_l", tag, i+1);
        SetName(cName, sample, "eta_c", tag, i+1);
        SetName(bName, sample, "eta_b", tag, i+1);
        MakePlot(file, hName, lName, cName, bName);
    }

    TImage *img = TImage::Create();
    img->FromPad(c1);
    img->WriteImage(pName);
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


void MakePlot(TFile * file, const char * hName, const char * lName, const char * cName, const char * bName)
{
    char name[256];

    file->cd("templates");

    sprintf(name, "template_%s", lName);
    TH1D * lTemplate = (TH1D*) gDirectory->Get(name);
    sprintf(name, "template_%s", cName);
    TH1D * cTemplate = (TH1D*) gDirectory->Get(name);
    sprintf(name, "template_%s", bName);
    TH1D * bTemplate = (TH1D*) gDirectory->Get(name);

    if (!lTemplate || !cTemplate || !bTemplate)
    {
        printf("Warning template %s or %s or %s do not exit\n", lName, cName, bName);
        return;
    }

    Double_t lNorm = lTemplate->Integral();
    Double_t cNorm = cTemplate->Integral();
    Double_t bNorm = bTemplate->Integral();

    if ( lNorm > 0. )
    {
        lTemplate->GetListOfFunctions()->Delete();
        lTemplate->SetTitle(hName);
        lTemplate->SetStats(false);
        lTemplate->DrawNormalized();
    }
    else
        printf("Warning template light norm <= 0\n");

    if ( cNorm > 0. )
    {
        cTemplate->GetListOfFunctions()->Delete();
        cTemplate->SetTitle(hName);
        cTemplate->SetStats(false);
        cTemplate->DrawNormalized("same");
    }
    else
        printf("Warning template c norm norm <= 0\n");

    if ( bNorm > 0 )
    {
        bTemplate->GetListOfFunctions()->Delete();
        bTemplate->SetTitle(hName);
        bTemplate->SetStats(false);
        bTemplate->DrawNormalized("same");
    }
    else
        printf("Warning template b norm norm <= 0\n");

    file->cd("functions");

    if ( lNorm > 0. )
    {
        sprintf(name, "function_%s", lName);
        TF1 * function = (TF1*) gDirectory->Get(name);
        TH1 * lFunction = function->GetHistogram();
        lFunction->SetLineColor(kGreen);
        lFunction->Scale(1./lNorm);
        lFunction->Draw("csame");
    }

    if ( cNorm > 0. )
    {
        sprintf(name, "function_%s", cName);
        TF1 * function = (TF1*) gDirectory->Get(name);
        TH1 * cFunction = function->GetHistogram();
        cFunction->SetLineColor(kRed);
        cFunction->Scale(1./cNorm);
        cFunction->Draw("csame");
    }

    if ( bNorm > 0. )
    {
        sprintf(name, "function_%s", bName);
        TF1 * function = (TF1*) gDirectory->Get(name);
        TH1 * bFunction = function->GetHistogram();
        bFunction->SetLineColor(kBlue);
        bFunction->Scale(1./bNorm);
        bFunction->Draw("csame");
    }
}

