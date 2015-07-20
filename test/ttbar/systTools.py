import ROOT

"""
"""
def smearJetEnergyResolution(pt,eta,genpt,varSign=0):
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

    ptSF_err *= varSign

    return ROOT.TMath.Max(0.,(genpt+(ptSF+ptSF_err)*(pt-genpt)))/pt
