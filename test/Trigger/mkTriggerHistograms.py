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
    parser.add_option('--outputDir'    , dest='outputDir'    , help='outputDir'   , default='./Histograms', type='string')

    (opt, args) = parser.parse_args()

    print 'Fill trigger histograms for PD', opt.Dataset, 'for periods Run',opt.Year,opt.Period

    if not os.path.isdir(opt.outputDir) :
          os.mkdir(opt.outputDir)

    triggerName  = []
    triggerIndex = []
    trgOffJetPt  = []

    triggers = {}
    if os.path.exists('triggers.py') :
        handle = open('triggers.py','r')
        exec(handle)
        handle.close()

    for trgName, trigger in triggers.iteritems() :
        if trigger['dataset'] == opt.Dataset :
            offJetPt = '0'
            if 'offjetpt' in trigger.keys() :
                offJetPt = trigger['offjetpt']
            for triggerYear in trigger['index'] :
                if triggerYear == opt.Year :
                    triggerName.append( trigger['name'] )
                    triggerIndex.append( int(trigger['index'][opt.Year]) )
                    trgOffJetPt.append( offJetPt )
                    
    samples = {}
    if os.path.exists('samples.py') :
        handle = open('samples.py','r')
        exec(handle)
        handle.close()

    jetPt = 'Jet_pt'
    kinCut = 'Jet_pt[0]>XXX && fabs(Jet_eta)<2.4'
    if opt.Dataset == 'BTagMu' :
        jetPt = 'Jet_pt[PFMuon_IdxJet]'
        kinCut = 'PFMuon_pt>=5 && fabs(PFMuon_eta)<2.4 && PFMuon_GoodQuality>=2'

    Periods = opt.Period.split(',')

    for period in Periods :
        
        dataName = opt.Dataset + "_Run" + opt.Year + period
        dataFound = False

        for sampleName, sample in samples.iteritems() :
            if dataName in sampleName :

                dataFound = True

                if opt.Verbose >= 1 :
                    print '  Processing', sampleName, 'Dataset'

                outputFile = ROOT.TFile.Open( opt.outputDir + "/" + sampleName + ".root", 'recreate')

                chain = ROOT.TChain(sample['treename'])
                
                for inputFile in glob.glob(sample['location'] + sample['filename'] + '*.root') :
                    if opt.Verbose >= 4 :
                        print '      merging file', inputFile 
                    chain.Add('root://eoscms.cern.ch/' + inputFile)

                if opt.Verbose >= 3 :
                    print '      Events in chain = ', chain.GetEntries()

                outputFile.cd()

                for trg in range (0, len(triggerName)) :

                    bitIdx = int(triggerIndex[trg]/32);
                    triggerCut = '(BitTrigger['+str(bitIdx)+'] & ( 1 << ('+str(triggerIndex[trg])+' - '+str(bitIdx)+'*32) ) )'
                    globalCut = triggerCut + ' && ' + kinCut.replace('XXX', trgOffJetPt[trg])

                    if opt.Verbose >= 2 :
                        print '    Filling histogram for ', triggerName[trg], 'cut', globalCut

                    h = ROOT.TH1F(triggerName[trg], "", 80, 0., 800.)
                    chain.Project(h.GetName(), jetPt, globalCut, '') 

                    if opt.Verbose >= 3 :
                        print '      Entries for', triggerName[trg], '=', h.GetEntries()

                    h.Write()

                outputFile.Close()

        if dataFound == False :
            print dataName, "not found"
