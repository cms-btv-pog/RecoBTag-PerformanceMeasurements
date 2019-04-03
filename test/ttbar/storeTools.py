import ROOT
import pickle

"""
Takes a directory on eos (starting from /store/...) and returns a list of all files with 'prepend' prepended
"""
def getEOSlslist(directory, mask='', prepend='root://eoscms//eos/cms'):
    from subprocess import Popen, PIPE
    print 'looking into: '+directory+'...'

    fullPaths=['/eos/cms/'+directory]

    eos_cmd = '/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'
    foundROOTOrEnd=False

    while not foundROOTOrEnd:
        allNewPaths=[]
        for fullPath in fullPaths:
            data = Popen([eos_cmd, 'ls', fullPath],stdout=PIPE)
            out,err = data.communicate()
            items=out.split('\n')
            if items[0].find(".root")!=-1:
                foundROOTOrEnd=True
                print "found ROOT files"
            elif len(items)==0:
                foundROOTOrEnd=True
                print "empty ls"
            else:
                for item in items:
                    if len(item)>0:
                        allNewPaths.append(fullPath+'/'+item)
        if not foundROOTOrEnd:
            fullPaths=allNewPaths


    full_list = []
    for fullPath in fullPaths:
        data = Popen([eos_cmd, 'ls', fullPath],stdout=PIPE)
        out,err = data.communicate()
        items=out.split('\n')

        ### if input file was single root file:
        #if directory.endswith('.root'):
        #    if len(out.split('\n')[0]) > 0:
        #        return [prepend + directory]

        ## instead of only the file name append the string to open the file in ROOT
        for item in items:
            if len(item.split()) == 0: continue
            full_list.append(prepend + fullPath.split("/eos/cms")[-1] + '/' + item)

    print "found n files",len(full_list)

    ## strip the list of files if required
    if mask != '':
        stripped_list = [x for x in full_list if mask in x]
        print "found n files in stripped list",len(stripped_list)
        return stripped_list

    print "found n files",len(full_list)
    returnList=full_list
    return returnList
#    return full_list

"""
Loops over a list of samples and produces a cache file to normalize MC
"""
def produceNormalizationCache(samplesList,inDir,cache,xsecWgts,integLumi):

    #loop over samples
    for tag,sample in samplesList: 

        if sample[1]==1 : 
            xsecWgts[tag]=None
            continue

        if tag in xsecWgts:
            print '[Warning] won\'t override current definition for',tag,'. Use --resetCache option to override'
            continue
        

        input_list=getEOSlslist(directory=inDir+'/'+tag)            
        xsec=sample[0]            
        norigEvents=None
        for f in input_list:
            fIn=ROOT.TFile.Open(f)
            if norigEvents is None:
                norigEvents=fIn.Get('ttbarselectionproducer/wgtcounter').Clone('xsecwgts')
                norigEvents.SetDirectory(0)
                norigEvents.Reset('ICE')
            norigEvents.Add(fIn.Get('ttbarselectionproducer/wgtcounter'))
            fIn.Close()
        try:
            for xbin in xrange(1,norigEvents.GetNbinsX()+1):
                norigEvents.SetBinContent(xbin,xsec/norigEvents.GetBinContent(xbin))
                norigEvents.SetBinError(xbin,0.)
        except:
            print 'No normalization histogram for ',tag
        xsecWgts[tag]  = norigEvents
        integLumi[tag] = 1./norigEvents.GetBinContent(1) if norigEvents else 0.

        if norigEvents:
            print '... %s cross section=%f pb weights sum(initial events)=%3.0f lumi=%3.2f/fb' % (tag,xsec,xsec/norigEvents.GetBinContent(1),integLumi[tag]/1000.)
            
        
    #dump to file    
    cachefile=open(cache,'w')
    pickle.dump(xsecWgts, cachefile, pickle.HIGHEST_PROTOCOL)
    pickle.dump(integLumi, cachefile, pickle.HIGHEST_PROTOCOL)
    cachefile.close()
    print 'Produced normalization cache and pileup weights @ %s'%cache

    return xsecWgts,integLumi
