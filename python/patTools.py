import FWCore.ParameterSet.Config as cms

from copy import deepcopy


def adaptPVs(process, pvCollection=cms.InputTag('offlinePrimaryVertices'), postfix='', sequence='patPF2PATSequence'):

    print "Switching PV collection for " + sequence + postfix + ":", pvCollection
    print "***********************************"

    # PV sources to be exchanged:
    pvExchange = ['Vertices','vertices','pvSrc','primaryVertices','primaryVertex','srcPVs','primaryVertex']

    # exchange the primary vertex source of all relevant modules
    for m in getattr(process,sequence+postfix).moduleNames():
        # only if the module has a source with a relevant name
        for namePvSrc in pvExchange:
            if hasattr(getattr(process,m),namePvSrc):
                setattr(getattr(process,m),namePvSrc,deepcopy(pvCollection))
