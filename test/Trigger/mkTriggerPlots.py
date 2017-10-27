#!/usr/bin/env python

import json
import sys
import ROOT
import optparse
import logging
import os.path
import glob

##Main body of the analysis
if __name__ == '__main__':

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)

    parser.add_option('--Verbose'      , dest='Verbose'      , help='Verbose'     , default=0,              type='int'   )
    parser.add_option('--Dataset'      , dest='Dataset'      , help='Dataset'     , default='BTagMu',       type='string')
    parser.add_option('--Year        ' , dest='Year'         , help='Year'        , default='2017',         type='string')
    parser.add_option('--Period'       , dest='Period'       , help='Period'      , default='',             type='string')
    parser.add_option('--Trigger'      , dest='Trigger'      , help='Trigger'     , default='',             type='string')
    parser.add_option('--inputDir'     , dest='inputDir'     , help='inputDir'    , default='./Histograms', type='string')
    parser.add_option('--outputDir'    , dest='outputDir'    , help='outputDir'   , default='./Plots',      type='string')
    parser.add_option('--plotOption'   , dest='plotOption'   , help='plotOption'  , default='',             type='string')

    (opt, args) = parser.parse_args()

    print 'Make trigger plots for PD', opt.Dataset, 'for years', opt.Year, 'and periods', opt.Period

    if not os.path.isdir(opt.outputDir) :
          os.mkdir(opt.outputDir)

    triggerName  = []

    triggers = {}
    if os.path.exists('triggers.py') :
        handle = open('triggers.py','r')
        exec(handle)
        handle.close()
                    
    samples = {}
    if os.path.exists('samples.py') :
        handle = open('samples.py','r')
        exec(handle)
        handle.close()

    selectedYears    = opt.Year.split(',')
    selectedTriggers = opt.Trigger.split(',')
    selectedPeriods  = opt.Period.split(',')

    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    triggerCanvas = ROOT.TCanvas()

    triggerCanvas.Divide(1, 1)
    triggerCanvas.Range(0,0,1,1)
    triggerCanvas.SetFillColor(10)
    triggerCanvas.SetBorderMode(0)
    triggerCanvas.SetBorderSize(2)
    triggerCanvas.SetTickx(1)
    triggerCanvas.SetTicky(1)
    triggerCanvas.SetLeftMargin(0.14)
    triggerCanvas.SetRightMargin(0.04)
    triggerCanvas.SetTopMargin(0.05)
    triggerCanvas.SetBottomMargin(0.13)
    triggerCanvas.SetFrameFillColor(0)
    triggerCanvas.SetFrameFillStyle(0)
    triggerCanvas.SetFrameBorderMode(0)

    triggerPad = triggerCanvas.GetPad(1)

    triggerPad.SetFillColor(0)
    triggerPad.SetBorderMode(0)
    triggerPad.SetBorderSize(2)
    triggerPad.SetTickx(1)
    triggerPad.SetTicky(1)
    triggerPad.SetLeftMargin(0.14)
    triggerPad.SetRightMargin(0.04)
    #triggerPad.SetTopMargin(0.05)
    #triggerPad.SetBottomMargin(0.31)
    triggerPad.SetTopMargin(0.065)
    triggerPad.SetBottomMargin(0.13)
    triggerPad.SetFrameFillStyle(0)
    triggerPad.SetFrameBorderMode(0)
    triggerPad.SetFrameFillStyle(0)
    triggerPad.SetFrameBorderMode(0)
    if 'LogX'  in opt.plotOption : triggerPad.SetLogx()
    if 'LogY'  in opt.plotOption : triggerPad.SetLogy()
    if 'GridX' in opt.plotOption : triggerPad.SetGridx()
    if 'GridY' in opt.plotOption : triggerPad.SetGridy()
    triggerPad.cd()
    
    titleAxisX = 'Jet p_{T} [GeV]'
    titleAxisY = 'Jets / 10 [GeV]'

    if 'Shape' in opt.plotOption :
        titleAxisY = titleAxisY.replace('Jet', 'Fraction of Jet')

    if opt.Dataset == 'BTagMu' :
        titleAxisX = titleAxisX.replace('Jet', 'muon-Jet')
        titleAxisY = titleAxisY.replace('Jet', 'muon-Jet')

    drawOption = 'histo'

    lineStyle = 1

    histoTrigger = []
    nHistos = 0

    scaleY = 1.1
    
    x1_legy = 0.50; y1_legy = 0.80; x2_legy = 0.80; y2_legy = 0.90;
    x1_legt = 0.50; y1_legt = 0.40; x2_legt = 0.80; y2_legt = 0.70;

    if 'LogX'  in opt.plotOption : 
        x1_legy = 0.60; y1_legy = 0.80; x2_legy = 0.90; y2_legy = 0.90;
        x1_legt = 0.60; y1_legt = 0.40; x2_legt = 0.90; y2_legt = 0.70;

    if 'LogY'  in opt.plotOption : 
        scaleY = 2.0
        x1_legy = 0.75; y1_legy = 0.82; x2_legy = 0.90; y2_legy = 0.90;
        x1_legt = 0.17; y1_legt = 0.17; x2_legt = 0.27; y2_legt = 0.47;
 
    legendYear = ROOT.TLegend(x1_legy, y1_legy, x2_legy, y2_legy)
    legendYear.SetFillColor(0) 
    legendYear.SetBorderSize(0)
    legendYear.SetTextColor(1) 
    legendYear.SetTextSize(0.035)
    legendYear.SetTextFont(42) 
    #legendYear.SetHeader('Triggers') 
    #legendYear.SetMargin(0.2) 

    legendTrigger = ROOT.TLegend(x1_legt, y1_legt, x2_legt, y2_legt)
    legendTrigger.SetFillColor(0) 
    legendTrigger.SetBorderSize(0)
    legendTrigger.SetTextColor(1) 
    legendTrigger.SetTextSize(0.035)
    legendTrigger.SetTextFont(42) 
    #legendTrigger.SetHeader('Triggers') 
    #legendTrigger.SetMargin(0.2)

    for selectedYear in selectedYears :

        lineColor = 1

        for trgName, trigger in triggers.iteritems() :
            if trigger['dataset']==opt.Dataset :

                for selectedTrigger in selectedTriggers :
                    if selectedTrigger in trigger['name'] :

                        for triggerYear in trigger['index'] :
                            if triggerYear == selectedYear :


                                histo = ROOT.TH1F(selectedYear + '_' + trigger['name'], "", 80, 0., 800.)

                                for sampleName, sample in samples.iteritems() :
                                    for selectedPeriod in selectedPeriods :

                                        dataName = opt.Dataset + "_Run" + selectedYear + selectedPeriod
                                        if dataName in sample['name'] :

                                            inputFileName = opt.inputDir + "/" + sample['name'] + '.root'

                                            if os.path.exists(inputFileName) :
                                                inputFile = ROOT.TFile(inputFileName , 'read')
                                                histoPeriod = inputFile.Get( trigger['name'] )
                                                histo.Add( histoPeriod )
                                                inputFile.Close()
                                            else:
                                                print '    Warning:', inputFileName, '.root file not found'

                                histo.SetLineColor( lineColor )
                                histo.SetLineStyle( lineStyle )
                                histo.SetLineWidth( 2 )

                                histo.SetXTitle( titleAxisX )
                                histo.SetYTitle( titleAxisY )

                                histo.GetXaxis().SetLabelFont(42)
                                histo.GetYaxis().SetLabelFont(42)
                                histo.GetXaxis().SetTitleFont(42)
                                histo.GetYaxis().SetTitleFont(42)

                                histo.GetYaxis().SetTitleSize(0.06)
                                histo.GetYaxis().SetLabelSize(0.05)
                                histo.GetXaxis().SetTitleSize(0.06)
                                histo.GetXaxis().SetLabelSize(0.05)

                                histo.GetXaxis().SetTitleOffset(0.95)
                                histo.GetYaxis().SetTitleOffset(1.15)
                                
                                if lineColor == 1 : legendYear.AddEntry(histo, selectedYear, "l")
                                if lineStyle == 1 : legendTrigger.AddEntry(histo, trigger['name'], "l")
                                
                                if 'Shape' in opt.plotOption :
                                    histoIntegral = histo.Integral(-1, -1)
                                    histo.Scale( 1./histoIntegral )

                                if (nHistos>0 and histo.GetMaximum()>histoTrigger[0].GetMaximum()) :
                                    histoTrigger[0].SetMaximum(scaleY*histo.GetMaximum())

                                histoTrigger.append( histo )

                                histoTrigger[nHistos].Draw( drawOption )

                                if 'same' not in drawOption:
                                    drawOption += 'same'

                                nHistos += 1

                        lineColor += 1
                        if lineColor == 10 : lineColor += 1
                        if lineColor == 12 : lineColor = 30

        lineStyle += 1

    legendYear.Draw()
    legendTrigger.Draw()

    plotTitle = opt.outputDir + '/Triggers_' + opt.Dataset + '_' + opt.Year.replace(',','-')
    if opt.Trigger != '' :
        plotTitle += '_' + opt.Trigger.replace(',','-')
    if opt.Period != '' :
        plotTitle += '_' + opt.Period.replace(',','-')
    if opt.plotOption != '' :
        plotTitle += '_' + opt.plotOption

    triggerCanvas.Print(plotTitle + '.png') 
