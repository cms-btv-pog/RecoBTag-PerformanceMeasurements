import optparse
import os,sys
import json
import commands
import ROOT
import random
import pickle

from Templated_btagEffFitter import SLICEBINS,SLICEVAR,SYSTVARS

"""
"""
def showFitResult(model,data,frame,pullFrame,nfreePars,labels=[]):

    c = ROOT.TCanvas('c','c',500,500)
    p1 = ROOT.TPad('p1','p1',0.0,0.85,1.0,0.0)
    p1.Draw()
    p1.SetTopMargin(0.01)
    p1.SetRightMargin(0.05)
    c.cd()
    p2 = ROOT.TPad('p2','p2',0.0,0.85,1.0,1.0)
    p2.SetBottomMargin(0.01)
    p2.SetRightMargin(0.05)
    p2.Draw()

    #compute chi2 from comparison of data and model
    p1.cd()
    p1.Clear()
    data.plotOn(frame)
    model.plotOn(frame,ROOT.RooFit.ProjWData(data))
    frame.Draw()
    chi2=frame.chiSquare(nfreePars)
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextFont(42)
    label.SetTextSize(0.035)
    label.DrawLatex(0.2,0.95,'#bf{CMS} #it{simulation}')
    label.DrawLatex(0.2,0.90,' %d variables' % nfreePars)
    label.DrawLatex(0.2,0.85,'#chi^{2}/ndof=%3.2f' % (chi2) )
    for i in xrange(0,len(labels)):
        label.DrawLatex(0.2,0.80-i*5,labels[i])

    #show pull
    p2.cd()
    p2.Clear()
    hpull = frame.pullHist()
    pullFrame.addPlotable(hpull,"P") ;
    pullFrame.Draw()
    pullFrame.GetYaxis().SetTitle("Pull")
    pullFrame.GetYaxis().SetTitleSize(0.2)
    pullFrame.GetYaxis().SetLabelSize(0.2)
    pullFrame.GetXaxis().SetTitleSize(0)
    pullFrame.GetXaxis().SetLabelSize(0)
    pullFrame.GetYaxis().SetTitleOffset(0.15)
    pullFrame.GetYaxis().SetNdivisions(4)
    pullFrame.GetYaxis().SetRangeUser(-3.1,3.1)
    pullFrame.GetXaxis().SetTitleOffset(0.8)
    
    c.Modified()
    c.Update()
    return c,chi2

"""
get distributions from file
"""
def getDiscriminators(opt,islice,syst=''):
    fIn=ROOT.TFile.Open(opt.input)
    othersH=fIn.Get('other_pass0_slice%d%s' % (islice,syst)).Clone('others')
    othersH.Add(fIn.Get('c_pass0_slice%d%s' % (islice,syst)))
    othersH.SetDirectory(0)
    bH=fIn.Get('b_pass0_slice%d%s' % (islice,syst)).Clone('b')
    bH.SetDirectory(0)
    fIn.Close()
    return othersH,bH

"""
"""
def runSystematics(opt,tag,ws):

        print 'Running systematics for',tag

        biasResults={'sig':{},'bkg':{}}
        
        #fix parameters
        ws.loadSnapshot('b_pars')
        ws.loadSnapshot('nonb_pars')
        allVars=ws.allVars()
        varIter = allVars.createIterator()
        var=varIter.Next()
        while var :
            varName=var.GetName()
            if not varName in ['nsig','nbkg'] : 
                ws.var(varName).setConstant()
            var=varIter.Next()
            
        #import data
        variables=ROOT.RooArgSet()
        variables.add(ws.var('KIN'))
        variables.add(ws.factory('bias[0,0,1]'))
        variables.add(ws.factory("wgt[0,0,9999999.]") )

        
        lumi=1000.
        for syst in SYSTLIST:
            othersH,bH=getDiscriminators(opt,sliceCond=tag,syst=syst)
            
            #create pseudo_data
            pseudoData=ROOT.RooDataSet("pseudo_data_%s"%syst,"pseudo_data_%s"%syst,variables,"wgt")
            
            nexp=int(lumi*othersH.Integral())        
            for i in xrange(0,nexp):
                ws.var('KIN').setVal( othersH.GetRandom() )
                ws.var('bias').setVal( random.uniform(0,1) )
                ws.var('wgt').setVal(1.0)
                argset = ROOT.RooArgSet()
                for var in ['KIN','bias','wgt'] : argset.add( ws.var(var) )
                pseudoData.add(argset, ws.var("wgt").getVal())
        
            nexp=int(lumi*bH.Integral())
            for i in xrange(0,nexp):
                ws.var('KIN').setVal( bH.GetRandom())
                ws.var('bias').setVal( random.uniform(0,1) )
                ws.var('wgt').setVal(1.0)
                argset = ROOT.RooArgSet()
                for var in ['KIN','bias','wgt']: argset.add( ws.var(var) )
                pseudoData.add(argset, ws.var("wgt").getVal())

            #unfold distribution with the sPlot technique
            sData = ROOT.RooStats.SPlot("s"+pseudoData.GetName(),"An SPlot from KIN",
                                        pseudoData,ws.pdf('model'),
                                        ROOT.RooArgList(ws.var('nsig'),ws.var('nbkg'))
                                        )
            
            #the weighted data for signal and background species
            for ipop in ['sig','bkg']: 
                ipopData = ROOT.RooDataSet(ipop+'bias',ipop+'bias',pseudoData,pseudoData.get(),'','n%s_sw'%ipop)
                biasH=ROOT.TH1F('biasH',';bias;',25,0,1)
                ipopData.fillHistogram(biasH,ROOT.RooArgList(ws.var('bias')))
                maxProb,maxProbFunc=0.,'pol0'
                for iord in xrange(0,10):
                    func='pol%d'%iord
                    biasH.Fit(func,'0MRQ+')
                    iprob=biasH.GetFunction(func).GetProb()
                    if iprob<maxProb : continue
                    maxProb=iprob
                    maxProbFunc=func
                biasResults[ipop][syst]=biasH.GetFunction(maxProbFunc).Clone('bias_%s_%s'%(ipop,syst))
                biasH.Delete()

            sData.Delete()
            pseudoData.Delete()
                
        fOut=ROOT.TFile.Open('bias_%s.root'%tag,'RECREATE')
        for ipop in biasResults:
                for syst in biasResults[ipop]:
                    biasResults[ipop][syst].Write()
        fOut.Close()

"""
"""        
def unfoldDiscriminatorShape(opt,tag,ws):

        #fix parameters
        ws.loadSnapshot('b_pars')
        ws.loadSnapshot('nonb_pars')
        allVars=ws.allVars()
        varIter = allVars.createIterator()
        var=varIter.Next()
        while var :
            varName=var.GetName()
            if not varName in ['nsig','nbkg'] : 
                ws.var(varName).setConstant()
            var=varIter.Next()
                        
        #import data from trees
        files=os.listdir(opt.data)
        chain=ROOT.TChain('kin')
        for f in files:
                if 'Data13' in f: chain.Add(os.path.join(opt.data,f))
        variables=ROOT.RooArgSet()
        variables.add(ws.var('KIN'))
        variables.add(ws.factory('csv[0,0,1]'))
        variables.add(ws.factory("wgt[0,0,9999999.]") )
        data=ROOT.RooDataSet("data","data",variables,'wgt')

        for i in xrange(0,chain.GetEntries()):
                chain.GetEntry(i)
                if not chain.ttbar_chan in [-11*13]: continue #-13*13, -11*13 -11*11] : continue
                if chain.csv<0 : continue
                ws.var('KIN').setVal( chain.kindisc[0] )
                ws.var('csv').setVal( chain.csv )
                ws.var('wgt').setVal(1.0)
                argset = ROOT.RooArgSet()
                for var in ['KIN','csv','wgt'] : argset.add( ws.var(var) )
                data.add(argset, ws.var("wgt").getVal())

        frame=ws.var('KIN').frame()
        pullFrame=ws.var('KIN').frame()
        c,_ = showFitResult(ws.pdf('model'),data,frame,pullFrame,2)
        c.SaveAs('datafit_%s.png'%tag)
        c.SaveAs('datafit_%s.pdf'%tag)
        
        #unfold distribution with the sPlot technique
        sData = ROOT.RooStats.SPlot("s"+data.GetName(),"An SPlot from KIN",
                                    data,ws.pdf('model'),
                                    ROOT.RooArgList(ws.var('nsig'),ws.var('nbkg'))
                                    )
            
        #the weighted data for signal and background species
        biasF=ROOT.TFile('bias_%s.root'%tag)
        for ipop in ['sig','bkg']:
                ipopdata = ROOT.RooDataSet('%scsv'%ipop,'%scsv'%ipop,data,data.get(),'','n%s_sw'%ipop)
                csvH = ROOT.TH1F('%s_csv'%ipop,';CSV;PDF',15,0,1)        
                ipopdata.fillHistogram(csvH,ROOT.RooArgList(ws.var('csv')))

                histos={}
                for syst in SYSTLIST:
                    corrF=biasF.Get('bias_%s_%s'%(ipop,syst))
                    histos[syst]=csvH.Clone('%s_csv%s'%(ipop,syst))
                    histos[syst].Reset('ICE')
                    for xbin in xrange(1,csvH.GetXaxis().GetNbins()+1):
                        corr=corrF.Eval(csvH.GetXaxis().GetBinCenter(xbin))
                        val=csvH.GetBinContent(xbin)*corr
                        err=csvH.GetBinError(xbin)*corr
                        histos[syst].SetBinContent(xbin,val)
                        histos[syst].SetBinError(xbin,err)
                    histos[syst].Scale(1./histos[syst].Integral())

                errUp=[0]*csvH.GetXaxis().GetNbins()
                errDn=[0]*csvH.GetXaxis().GetNbins()
                for xbin in xrange(1,csvH.GetXaxis().GetNbins()+1):
                    cen=histos[''].GetBinContent(xbin)
                    for syst in SYSTLIST:
                        if syst=='':
                            statUnc=histos[''].GetBinError(xbin)
                            errUp[xbin-1]=errUp[xbin-1]+statUnc**2
                            errDn[xbin-1]=errDn[xbin-1]+statUnc**2
                        else:
                            diff=histos[syst].GetBinContent(xbin)-cen
                            if diff>0 : errUp[xbin-1] = errUp[xbin-1]+diff**2
                            else :      errDn[xbin-1] = errDn[xbin-1]+diff**2
                csvStatGr=ROOT.TGraphAsymmErrors(histos[''])
                csvStatGr.SetName('csvstatunc')
                csvTotalGr=csvStatGr.Clone('csvtotalunc')
                for np in xrange(0,csvTotalGr.GetN()):
                        ex=0.5*csvH.GetXaxis().GetBinWidth(np+1)
                        csvTotalGr.SetPointError(np,ex,ex,ROOT.TMath.Sqrt(errDn[np]),ROOT.TMath.Sqrt(errUp[np]))
                            
                c.Clear()
                csvTotalGr.SetFillStyle(1001)
                csvTotalGr.SetFillColor(ROOT.kGray)
                csvTotalGr.Draw('a2')
                csvTotalGr.GetXaxis().SetTitle('CSV discriminator')
                csvTotalGr.GetYaxis().SetTitle('Events')
                csvTotalGr.GetXaxis().SetRangeUser(0,1)
                csvTotalGr.GetYaxis().SetRangeUser(0,histos[''].GetMaximum()*1.3)
                csvStatGr.SetMarkerStyle(20)
                csvStatGr.SetFillStyle(1001)
                csvStatGr.SetFillColor(ROOT.kOrange+1)
                csvStatGr.Draw('2')
                histos[''].Draw('histsame')              
                label = ROOT.TLatex()
                label.SetNDC()
                label.SetTextFont(42)
                label.SetTextSize(0.035)
                label.DrawLatex(0.12,0.96,'#bf{CMS} #it{preliminary}')
                leg=ROOT.TLegend(0.12,0.91,0.5,0.94)
                leg.SetTextFont(42)
                leg.SetNColumns(2)
                leg.SetBorderSize(0)
                leg.SetTextSize(0.035)
                leg.SetFillStyle(0)
                leg.AddEntry(csvStatGr,'stat','lf')
                leg.AddEntry(csvTotalGr,'total','f')
                leg.Draw()
                c.Modified()
                c.Update()
                raw_input()

"""
"""
def prepareTemplates(opt,islice):

    ws  = ROOT.RooWorkspace('w_slice%d'%islice)
    KIN = ws.factory('KIN[0,-1,1]')
    ws.factory("RooFormulaVar::KIN_shift('@0+1.0',{KIN})")

    othersH,bH = getDiscriminators(opt,islice)
    data={
        'b':ROOT.RooDataHist('b',   'b',   ROOT.RooArgList(KIN),bH),
        'nonb':ROOT.RooDataHist('nonb','nonb',ROOT.RooArgList(KIN),othersH)
        }

    #prepare options for the fit
    rllist = ROOT.RooLinkedList()
    rllist.Add(ROOT.RooFit.Save())
    rllist.Add(ROOT.RooFit.Strategy(2))

    #prepare base models
    ws.factory("SUM::nonb_model("
                "nonb_frac1[0,0,1]*RooArgusBG::nonb_argus(KIN_shift,nonb_endpoint[2.0,1.90,2.1],nonb_shape[-0.5,-10,-0.]),"
                "nonb_frac2[0,0,1]*RooExponential::nonb_expo(KIN_shift,nonb_lambda[-0.5,-10,-0.]),"
                "Chebychev::nonb_pol(KIN,{nonb_a[0,-10,10]})"
                ")")

    ws.factory("SUM::b_model("
                "b_frac1[0,0,1]*RooArgusBG::b_argus(KIN_shift,b_endpoint[2.0,1.90,2.1],b_shape[-0.5,-50,-0.]),"
                "b_frac2[0,0,0.1]*RooGaussian::b_gauss(KIN_shift,b_mean[0.0],b_sigma[0.05,0.01,0.1]),"                
                "Chebychev::b_pol(KIN,{b_a[0,-10,10]})"
                ")")

    #loop to find smallest bias model
    for dist in ['nonb','b']:
        polArgs=''
        ws.factory('%s_polorder[1]'%dist)
        minChi2=999999.
        for i in xrange(1,8):
            if len(polArgs)>0 : polArgs+=','

            #specialize model
            modelName='%s_%d' % (dist,i)
            polArgs+='%s_a%d[0,-10,10]' % (dist,i)
            ws.factory("EDIT::%s(%s_model,%s_pol=Chebychev::%s_pol_%d(KIN,{%s}))" % (modelName,dist,dist,dist,i,polArgs))    

            #fit
            fitRes=ws.pdf(modelName).chi2FitTo(data[dist],rllist)
            floatParsFinal=fitRes.floatParsFinal()
            nPars=len(floatParsFinal)
            c,chi2 = showFitResult(model=ws.pdf(modelName),
                                   data=data[dist],
                                   frame=KIN.frame(),
                                   pullFrame=KIN.frame(),
                                   nfreePars=nPars,
                                   labels=['%s distribution' % dist])

            #save if it minimizes the chi2               
            if chi2>minChi2 : continue
            minChi2=chi2
            ws.var('%s_polorder'%dist).setVal(i)
            parsToFreeze=ROOT.RooArgSet(floatParsFinal)
            ws.saveSnapshot('%s_pars' % dist,parsToFreeze)
            c.SaveAs('%s_%s_model.pdf' % (dist,sliceCond))

    
    #build model to fit data
    ws.factory('nsig[10,0,9999999999999.]')
    ws.factory('nbkg[10,0,9999999999999.]')
    model=ROOT.RooAddPdf("model","b+nonb model",
                         ROOT.RooArgList( ws.pdf("b_%d" % ws.var('b_polorder').getVal()), ws.pdf("nonb_%d" %  ws.var('nonb_polorder').getVal()) ),
                         ROOT.RooArgList( ws.var("nsig"), ws.var("nbkg") )
                         )
    getattr(ws,'import')(model)

    #save workspace    
    ws.writeToFile(os.path.join(opt.outDir,'KINworkspace.root'),False)


"""
Wrapper to be used when run in parallel
"""
def runSystematicsPacked(args):
    opt, tag, ws = args
    return runSystematics(opt=opt,tag=tag,ws=ws)
    
    
"""
steer the script
"""
def main():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.Minimization)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Minimization)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.ObjectHandling)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.DataHandling)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Fitting)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Plotting)
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.Plotting)
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.InputArguments)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.InputArguments)
    ROOT.RooMsgService.instance().getStream(0).removeTopic(ROOT.RooFit.Eval)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Eval)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Integration)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration)
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration)
  
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',       dest='input',       help='input file',                  default=None,        type='string')
    parser.add_option('-o', '--output',      dest='outDir',      help='output directory',            default=None,        type='string')
    parser.add_option(      '--ws',          dest='ws',          help='input file with workspace',   default=None,        type='string')
    parser.add_option(      '--data',        dest='data',        help='input directory with data',   default=None,        type='string')
    parser.add_option(      '--doBias',      dest='doBias',      help='compute bias',                action='store_true')
    parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel',   default=0,           type='int')    
    (opt, args) = parser.parse_args()
  
    #recreate workspace
    if opt.ws is None:
        fOut=ROOT.TFile('KINworkspace.root','RECREATE')
        fOut.Close()
        for islice in xrange(0,len(SLICEBINS)):
            prepareTemplates(opt=opt,islice=islice)
        opt.ws='KINworkspace.root'

    #read workspace from file
    fIn=ROOT.TFile(opt.ws)

    #compute bias, if required
#    if opt.doBias:
#        task_list=[]
#        for sliceCond in SLICES:
#            task_list.append( (opt,sliceCond,fIn.Get('w_%s'%sliceCond) ) )

        #run the analysis jobs
#        if opt.njobs == 0:
#            for opt, tag, ws in task_list:
#                runSystematics(opt,tag,ws)
#        else:
#            from multiprocessing import Pool
#            pool = Pool(opt.njobs)
#            pool.map(runSystematicsPacked, task_list)

    #unfold
#    for sliceCond in SLICES:
#        unfoldDiscriminatorShape(opt,sliceCond,fIn.Get('w_%s'%sliceCond))
         


"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
