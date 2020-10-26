#!/usr/bin/env python
"""
This is a small script that submits CRAB jobs over many datasets
"""
import os
import glob
import argparse

def getOptions() :

  parser = argparse.ArgumentParser(description='CLI for CRAB config file settings')
  parser.add_argument('cmsRun_cfg', type=str, default="runBTagAnalyzer_cfg.py",
      help='The crab script you want to submit')
  parser.add_argument('-i', '--inputFiles', type=str, default='CRAB/input.txt',
      help='Input files that need to be shipped with the job')
  parser.add_argument('-p', '--pyCfgParams', nargs='+', type=str, default='',
      help='Input parameters for config file')
  parser.add_argument('-t', '--maxJobRuntimeMin', type=int, default=2750,
      help='The maximum runtime (in minutes) per job')
  parser.add_argument('-m', '--maxMemoryMB', type=int, default=4000,
      help='Maximum amount of memory (in MB) a job is allowed to use')
  parser.add_argument('-l', '--lumiMask', type=str, default='',
      help='The JSON file containing good lumi list')
  parser.add_argument('-f', '--files', type=str, default='CRAB/tosubmit.txt',
      required=True,
      help='File listing datasets to run over')
  parser.add_argument('-v', '--version', type=str, default='',
      help='BTagAnalyer version')
  parser.add_argument('-s', '--storageSite', type=str, default='T2_CH_CERN',
      help='Storage site')
  parser.add_argument('-o', '--outLFNDirBase', type=str, default='',
      required=True,
      help='EOS path for storage')
  args = parser.parse_args()

  if args.files == None:
    parser.error('-f (sample list) is required')
  if args.outLFNDirBase == None:
    parser.error('-o (EOS LFN path) is required')

  return args

def main():

    args = getOptions()

    from CRABClient.UserUtilities import config
    config = config()

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    config.General.workArea = args.version
    config.General.transferLogs = False
    config.General.transferOutputs = True

    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = args.cmsRun_cfg
    if args.inputFiles != None:
      inFiles = glob.glob( args.inputFiles )
      config.JobType.inputFiles = inFiles
    config.JobType.pyCfgParams = args.pyCfgParams
    #~ config.JobType.sendExternalFolder = True
    config.JobType.maxJobRuntimeMin = args.maxJobRuntimeMin
    config.JobType.maxMemoryMB = args.maxMemoryMB

    config.Data.inputDataset = None
    config.Data.splitting = ''
    config.Data.unitsPerJob = 1
    config.Data.ignoreLocality = True
    config.Data.publication = False
    #~ config.Data.publishDBS = 'phys03'
    config.Data.totalUnits = -1
    config.Site.storageSite = args.storageSite
    if config.Data.ignoreLocality:
       config.Site.whitelist = ['T2_CH_CERN', 'T2_DE_*','T1_US_FNAL*']

    print('Using config  {}'.format(args.cmsRun_cfg))
    print('Writing to versionectory {}'.format(args.version))

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print('Cannot execute command')
            print(hte.headers)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    datasetsFile = open( args.files )
    jobsLines = datasetsFile.readlines()
    jobs = []
    for ijob in jobsLines :
        s = ijob.rstrip()
        if (len(s)==0 or s[0][0]=='#'): continue
        s = ijob.rstrip()
        jobs.append( s )
        print '  --> added ' + s

    for ijob, job in enumerate(jobs) :

        pdname = job.split('/')[1]
        cond = job.split('/')[2]
        datatier = job.split('/')[3]
        requestname = '_'.join(pdname.split('_')[:3]) + '_' + cond
        if len(requestname) > 100: requestname = '_'.join(pdname.split('_')[:2]) + '_' + cond
        config.General.requestName = requestname
        config.Data.inputDataset = job
        if datatier == 'MINIAODSIM':
          config.Data.splitting = 'EventAwareLumiBased'
          config.Data.unitsPerJob = 100
        if datatier == 'FEVT':
          config.Data.splitting = 'EventAwareLumiBased'
          config.Data.unitsPerJob = 100
        elif datatier == 'AODSIM':
          config.Data.splitting = 'FileBased'
        elif datatier == 'MINIAOD':
          config.Data.splitting = 'LumiBased'
          config.Data.unitsPerJob = 40
          config.Data.lumiMask = args.lumiMask
        elif datatier == 'AOD':
          config.Data.splitting = 'LumiBased'
          config.Data.unitsPerJob = 100
          config.Data.lumiMask = args.lumiMask
        if args.outLFNDirBase and not args.outLFNDirBase.isspace():
          config.Data.outLFNDirBase = os.path.join(args.outLFNDirBase,args.version)
        config.Data.outputDatasetTag = cond
        print('Submitting {config.General.requestName} dataset =  {job}'.format(**locals()))
        print('Configuration :')
        print(config)

        try :
          from multiprocessing import Process
          p = Process(target=submit, args=(config,))
          p.start()
          p.join()
          submit(config)
          print('Submitted')
        except :
          print('Not submitted.')

if __name__ == '__main__':
  main()
