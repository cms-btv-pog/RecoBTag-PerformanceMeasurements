#!/usr/bin/env python
import os, argparse, json, fnmatch

from common import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='-- submit crab3 tasks --')

    parser.add_argument('-p', '--production-json', dest='production_json', action='store', default='', required=True,
                        help='path to input production .json file')

    parser.add_argument('--only', dest='only', nargs='+', default=[],
                        help='consider only production-json keys matching one of these strings (via fnmatch)')

    parser.add_argument('--no-submit', dest='no_submit', action='store_true', default=False,
                        help='do not submit crab3 tasks (only crab3 configuration files are produced)')

    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,
                        help='show verbose printouts during execution')

    parser.add_argument('-d', '--dry-run', dest='dry_run', action='store_true', default=False,
                        help='enable dry-run mode')

    opts, opts_unknown = parser.parse_known_args()
    ###

    log_prx = os.path.basename(__file__)+' -- '

    if len(opts_unknown) > 0:
       KILL(log_prx+'unknown command-line arguments: '+str(opts_unknown))

    if not opts.no_submit: which('crab')

    if not os.path.isfile(opts.production_json):
       KILL(log_prx+'invalid path to input production .json [-p]: '+opts.production_json)

    production_dict = json.load(open(opts.production_json))

    for i_dset_key in sorted(production_dict.keys()):

        if len(opts.only) > 0:

           skip_dset = True

           for _tmp in opts.only:

               if fnmatch.fnmatch(i_dset_key, _tmp):

                  skip_dset = False
                  break

           if skip_dset: continue

        if 'crab3' not in production_dict[i_dset_key]:
           KILL(log_prx+'input dictionary for "'+i_dset_key+'" does not contain a key named "crab3" (where are the crab3 configuration parameters?)')

        crab3_conf = production_dict[i_dset_key]['crab3']

        if 'Data.inputDataset' not in crab3_conf:
           KILL(log_prx+'input dictionary for "'+i_dset_key+'" does not contain crab3 configuration parameter "Data.inputDataset"')

        elif not (crab3_conf['Data.inputDataset'].startswith('/') and (len(crab3_conf['Data.inputDataset'].split('/')) == 4)):
           KILL(log_prx+'invalid dataset name for entry "'+i_dset_key+'" in input .json file: '+crab3_conf['Data.inputDataset'])

        for i_key in crab3_conf:
            if isinstance(crab3_conf[i_key], unicode):
               crab3_conf[i_key] = crab3_conf[i_key].encode('UTF8')

            elif isinstance(crab3_conf[i_key], list):
               crab3_conf[i_key] = [(_tmp.encode('UTF8') if isinstance(_tmp, unicode) else _tmp) for _tmp in crab3_conf[i_key]]

        crab3_tag = str(i_dset_key)

        if os.path.exists('crab_'+crab3_tag):
           WARNING(log_prx+'target directory for crab3 task already exists (will be skipped): '+'crab_'+crab3_tag)
           continue

        crab3_cfgfile = 'crab3cfg_'+i_dset_key+'.py'

        if os.path.exists(crab3_cfgfile):
           KILL(log_prx+'target path to crab3 configuration file already exists: '+crab3_cfgfile)

        if not opts.dry_run:

           with open(crab3_cfgfile, 'w') as cfgf:

               cfgf.write("from WMCore.Configuration import Configuration\n")

               cfgf.write("\n")

               cfgf.write("config = Configuration()\n")

               cfgf.write("\n")

               cfgf.write("config.section_('General')\n")
               cfgf.write("config.General.requestName     = \'"+crab3_tag+"\'\n")
               cfgf.write("config.General.transferOutputs = True\n")
               cfgf.write("config.General.transferLogs    = False\n")

               cfgf.write("\n")

               cfgf.write("config.section_('JobType')\n")
               cfgf.write("config.JobType.pluginName  = 'Analysis'\n")
               cfgf.write("config.JobType.allowUndistributedCMSSW = True \n")
#               cfgf.write("config.JobType.outputFiles = ['"+output_filename+"']\n")

               for i_key in crab3_conf:
                   if i_key.startswith('JobType.'):

                      i_out = str(crab3_conf[i_key])
                      if isinstance(crab3_conf[i_key], str): i_out = '\''+i_out+'\''

                      cfgf.write('config.'+i_key+' = '+i_out+'\n')

               cfgf.write("\n")

               cfgf.write("config.section_('User')\n")

               cfgf.write("\n")
               cfgf.write("config.section_('Site')\n")
               cfgf.write("config.Site.storageSite = 'T2_DE_DESY'\n")
#               cfgf.write("config.Site.blacklist   = []\n")

               for i_key in crab3_conf:
                   if i_key.startswith('Site.'):

                      i_out = str(crab3_conf[i_key])
                      if isinstance(crab3_conf[i_key], str): i_out = '\''+i_out+'\''

                      cfgf.write('config.'+i_key+' = '+i_out+'\n')

               cfgf.write("\n")

               cfgf.write("config.section_('Data')\n")
               cfgf.write("config.Data.publication    = False\n")

               if 'Data.ignoreLocality' not in crab3_conf:
                  cfgf.write("config.Data.ignoreLocality = False\n")

               for i_key in crab3_conf:

                   if i_key.startswith('Data.'):

                      if i_key == 'Data.ignoreLocality':
                         i_out = 'True' if int(crab3_conf[i_key]) else 'False'

                      else:
                         i_out = crab3_conf[i_key]
                         if isinstance(i_out, str): i_out = '\''+i_out+'\''

                      cfgf.write('config.'+i_key+' = '+str(i_out)+'\n')

        print colored_text('['+crab3_cfgfile+']', ['1', '93'])

        if not opts.no_submit:

           EXE('crab submit -c '+crab3_cfgfile, verbose=opts.verbose, dry_run=opts.dry_run)

           EXE('mv '+crab3_cfgfile+' crab_'+crab3_tag, verbose=opts.verbose, dry_run=opts.dry_run)
           EXE('rm '+crab3_cfgfile.replace('.py', '.pyc'), verbose=opts.verbose, dry_run=opts.dry_run)
