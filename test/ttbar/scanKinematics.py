import ROOT

def scanEvents(url=''):
    inF=ROOT.TFile.Open(url)
    tree=inF.Get("btagana/ttree")
    nentries=tree.GetEntriesFast()

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    lpGr=ROOT.TGraph()
    lpGr.SetMarkerStyle(20)
    lpGr.SetFillStyle(0)
    lpGr.SetTitle('l^{+}')
    bmGr=ROOT.TGraph()
    bmGr.SetMarkerStyle(21)
    bmGr.SetFillStyle(0)
    bmGr.SetTitle('b')

    lmGr=ROOT.TGraph()
    lmGr.SetMarkerStyle(24)
    lmGr.SetFillStyle(0)
    lmGr.SetTitle('l^{-}')
    bpGr=ROOT.TGraph()
    bpGr.SetMarkerStyle(25)
    bpGr.SetFillStyle(0)
    bpGr.SetTitle('#bar{b}')
    
    dilGr=ROOT.TGraph()
    dilGr.SetMarkerStyle(29)
    dilGr.SetMarkerColor(ROOT.kGreen)
    dilGr.SetLineColor(ROOT.kGreen)
    dilGr.SetFillStyle(0)
    dilGr.SetTitle('ll')

    otherJetsGr=ROOT.TGraph()
    otherJetsGr.SetMarkerStyle(5)
    otherJetsGr.SetFillStyle(0)
    otherJetsGr.SetTitle('other jets')


    c=ROOT.TCanvas('c','c',500,500)
    frame=ROOT.TH2F('event',';Pseudo-rapidity;Azimuthal angle [rad];',1,-2.5,2.5,1,-3.1415,3.1415)
    
    for i in xrange(0,nentries): 
        tree.GetEntry(i)
        if tree.ttbar_nl<2 : continue


        lp4={}
        for il in xrange(0,2):
            charge=tree.ttbar_lch[il]
            lp4[charge]=ROOT.TLorentzVector()
            lp4[charge].SetPtEtaPhiM(tree.ttbar_lpt[il],tree.ttbar_leta[il],tree.ttbar_lphi[il],0.)
        dilepton=lp4[+1]+lp4[-1]
        if ROOT.TMath.Abs(dilepton.M()-91)<15 : continue         

        bmGr.Set(0)
        bpGr.Set(0)
        otherJetsGr.Set(0)
        selJets=0
        for ij in xrange(0,tree.nJet):

            jp4=ROOT.TLorentzVector()
            jp4.SetPtEtaPhiM(tree.Jet_pt[ij],tree.Jet_eta[ij],tree.Jet_phi[ij],0.)
            
            minDRlj=9999.
            for ilcharge in lp4:
                minDRlj = ROOT.TMath.Min( minDRlj, lp4[ilcharge].DeltaR(jp4) )
            if minDRlj<0.4 : continue

            selJets +=1
            if tree.Jet_flavour[ij]==5 and bmGr.GetN()==0:
                bmGr.SetPoint(0,jp4.Eta(),jp4.Phi())
            elif tree.Jet_flavour[ij]==-5 and bpGr.GetN()==0:
                bpGr.SetPoint(0,jp4.Eta(),jp4.Phi())
            else:
                otherJetsGr.SetPoint(otherJetsGr.GetN(),jp4.Eta(),jp4.Phi())


        if selJets<2 : continue

        c.Clear()
        frame.Draw()
        lpGr.SetPoint(0,lp4[+1].Eta(),lp4[+1].Phi())
        lpGr.Draw('p')
        lmGr.SetPoint(0,lp4[-1].Eta(),lp4[-1].Phi())
        lmGr.Draw('p')
        dilGr.SetPoint(0,dilepton.Eta(),dilepton.Phi())
        dilGr.Draw('p')
        bmGr.Draw('p')
        bpGr.Draw('p')
        otherJetsGr.Draw('p')
        leg=c.BuildLegend(0.1,0.91,0.9,0.99)
        leg.SetFillStyle(0)
        leg.SetBorderSize(1)
        leg.SetNColumns(4)
        c.Modified()
        c.Update()
        raw_input()

scanEvents('root://eoscms//eos/cms/store/group/phys_btag/performance/TTbar/7b810a5/MC13TeV_TTJets/MergedJetTree_0.root')
