import FWCore.ParameterSet.Config as cms

from copy import deepcopy


def adaptPVs(process, pvCollection=cms.InputTag('offlinePrimaryVertices')):

    print "Switching PV collection to ", pvCollection
    print "***********************************"

    # PV sources to be exchanged:
    pvExchange = ['Vertices','vertices','pvSrc','primaryVertices','srcPVs','primaryVertex']

    # exchange the primary vertex source of all relevant modules
    for m in (process.producerNames().split(' ') + process.filterNames().split(' ')):
        # only if the module has a source with a relevant name
        for namePvSrc in pvExchange:
            if hasattr(getattr(process,m),namePvSrc):
                #print m
                setattr(getattr(process,m),namePvSrc,deepcopy(pvCollection))
