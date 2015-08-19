import ROOT

"""
Sources
 CMS AN 022/2015 v15
 https://indico.cern.ch/event/434078/#preview:1614815
"""
def getTriggerEfficiency(pt1,eta1,pt2,eta2,ch):
    if ch==-11*13:
        return 0.91,0.05
    if ch==-11*11:
        return 0.95,0.05
    if ch==-13*13:
        return 0.92,0.05
    raise ValueError('Unknown trigger efficiency for ch=%d' % ch)

"""
Sources
 CMS AN 022/2015 v15
"""
def getLeptonSelectionEfficiencyScaleFactor(pt,eta,pdgId):

    #electrons
    if ROOT.TMath.Abs(pdgId)==11:
        if ROOT.TMath.Abs(eta)<0.8:
            if pt<30 : return 0.927,0.073
            elif pt<40 : return 0.975,0.018
            elif pt<50 : return 0.962,0.036
            else : return 0.955,0.022
        elif ROOT.TMath.Abs(eta)<1.5:
            if pt<30 : return 0.891,0.074
            elif pt<40 : return 0.965,0.020
            elif pt<50 : return 0.968,0.018
            else : return 0.955,0.018
        else:
            if pt<30 : return 0.956,0.059
            elif pt<40 : return 0.995,0.018
            elif pt<50 : return 0.993,0.019
            else : return 0.985,0.023
        raise ValueError('Unknown electron selection efficiency for pt=%3.1f eta=%3.1f' % (pt,eta))

    #muons
    if  ROOT.TMath.Abs(pdgId)==13:
        if ROOT.TMath.Abs(eta)<0.9:
            if pt<30 : return 1.003,0.019
            elif pt<40 : return 1.014,0.015
            elif pt<50 : return 1.001,0.014
            else : return 0.983,0.014
        elif ROOT.TMath.Abs(eta)<1.2:
            if pt<30 : return 0.993,0.019
            elif pt<40 : return 0.994,0.015
            elif pt<50 : return 0.980,0.014
            else : return 0.987,0.015
        else:
            if pt<30 : return 1.023,0.028
            elif pt<40 : return 0.994,0.014
            elif pt<50 : return 0.996,0.014
            else : return 0.979,0.014
        raise ValueError('Unknown muon selection efficiency for pt=%3.1f eta=%3.1f' % (pt,eta))

    #unknown
    raise ValueError('Unknown lepton selection efficiency for pdgId=%d' % pdgId)

"""
Sources
  assuming nominal JES
"""
def getJetEnergyScales(pt,eta,rawsf,area,rho):
    dnSF  = 0.97
    nomSF = 1.0
    upSF  = 1.03
    return dnSF, nomSF, upSF

"""
Sources
  Assuming nominal JER but uncertainties from Run I
  https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
"""
def getJetResolutionScales(pt,eta,genpt,varSign=0):
    eta=ROOT.TMath.Abs(eta)
    ptSF, ptSF_err = 1.000, 0.060
    if eta<0.5 :
        ptSF, ptSF_err = 1.000, ROOT.TMath.Sqrt((0.012**2)+(0.5*(0.062+0.061))**2)
    elif eta<1.1:
        ptSF, ptSF_err = 1.000, ROOT.TMath.Sqrt((0.012**2)+(0.5*(0.056+0.055))**2)
    elif eta<1.7 :
        ptSF, ptSF_err = 1.000, ROOT.TMath.Sqrt((0.017**2)+(0.5*(0.063+0.062))**2)
    elif eta<2.3 :
        ptSF, ptSF_err = 1.000, ROOT.TMath.Sqrt((0.035**2)+(0.5*(0.087+0.085))**2)
    else:
        ptSF, ptSF_err = 1.000, ROOT.TMath.Sqrt((0.127**2)+(0.5*(0.155+0.153))**2)

    dnSF  = ROOT.TMath.Max(0.,(genpt+(ptSF-ptSF_err)*(pt-genpt)))/pt
    nomSF = ROOT.TMath.Max(0.,(genpt+(ptSF)*(pt-genpt)))/pt
    upSF  = ROOT.TMath.Max(0.,(genpt+(ptSF+ptSF_err)*(pt-genpt)))/pt

    return dnSF, nomSF, upSF
