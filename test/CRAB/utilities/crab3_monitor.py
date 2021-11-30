#!/usr/bin/env python
import argparse, os, json, glob, time

from common import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='-- monitor task(s), resubmit failed jobs, merge outputs --')

    parser.add_argument('-t', '--tasks', dest='crab_tasks', nargs='+', default=[], required=True,
                        help='list of paths to crab3-task directories')

    parser.add_argument('-r', '--resubmit', dest='resubmit', action='store_true', default=False,
                        help='resubmit crab task (if there are failed jobs)')

    parser.add_argument('--hadd', dest='hadd', action='store_true', default=False,
                        help='run ROOT\'s hadd to merge crab3 .root outputs')

    parser.add_argument('--force-hadd', dest='force_hadd', action='store_true', default=False,
                        help='force merging of all available outputs of a crab3 task (even if not all jobs are finished)')

    parser.add_argument('-p', '--production-json', dest='production_json', action='store', default='', required=False,
                        help='path to input production .json file (necessary only when running with --hadd)')

    parser.add_argument('--Tier2-prepath', dest='Tier2_prepath', action='store', default='/pnfs/desy.de/cms/tier2',
                        help='path to /store/ directory on local Tier-2')

    parser.add_argument('--repeat', dest='repeat', action='store_true', default=False,
                        help='repeats resubmitting until all jobs are finished')

    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False,
                        help='show verbose printouts during execution')

    opts, opts_unknown = parser.parse_known_args()

    log_prx = os.path.basename(__file__)+' -- '

    ### ----------

    which('crab')

    # unknown command-line arguments:
    #  - if resubmit=True, additional options for "crab resubmit"
    if opts.resubmit:
       RESUBMIT_OPTS = list(set(opts_unknown))

    elif len(opts_unknown) > 0:
       KILL(log_prx+'unknown command-line arguments: '+str(opts_unknown))

    # hadd
    if opts.hadd:

       if not os.path.isdir(opts.Tier2_prepath):
          KILL(log_prx+'invalid path to Tier2 base directory: '+opts.Tier2_prepath)

       if not os.path.isfile(opts.production_json):
          KILL(log_prx+'invalid path to input production .json [-p] (required by --hadd option): '+opts.production_json)

       production_dict = json.load(open(opts.production_json))

       which('hadd_wrapper.py')

    ### task inspection
    crab3_task_dir_format = '{: <'+str(2+max([len(os.path.basename(os.path.abspath(_tmp))) for _tmp in opts.crab_tasks]))+'}'

    CRAB3_TASK_DIRS_ABSPATH = [os.path.abspath(_tmp) for _tmp in opts.crab_tasks]
    CRAB3_TASK_DIRS_ABSPATH = sorted(list(set(CRAB3_TASK_DIRS_ABSPATH))) # remove duplicates and order

    iTotal = 0
    iFinished = 0
    loopEnd = False

    while (not(loopEnd)):

        if not opts.repeat: loopEnd = True

        for i_crab3_task_dir in CRAB3_TASK_DIRS_ABSPATH:

            if not os.path.isdir(i_crab3_task_dir):
               WARNING(log_prx+'invalid path to crab3 task directory: '+i_crab3_task_dir)
               continue

            iTotal += 1

            crab_status_dict = {}

            status_lines = command_output_lines('crab status --json '+i_crab3_task_dir, stdout=True, stderr=False, permissive=True)

            for i_line in status_lines:
                if i_line.startswith('{') and i_line.endswith('}'): crab_status_dict = json.loads(i_line)

            task_status_dict = {}

            for i_task in sorted(crab_status_dict.keys()):

                i_state = str(crab_status_dict[i_task]['State'])

                if i_state not in task_status_dict: task_status_dict[i_state] = 0

                task_status_dict[i_state] += 1

            jobs_total = sum(task_status_dict.values())
            jobs_finished = (task_status_dict['finished'] if 'finished' in task_status_dict else 0)

            i_crab3_task_dir_basename = os.path.basename(os.path.abspath(i_crab3_task_dir))

            if jobs_total == 0:
               WARNING(log_prx+'no jobs found for '+colored_text(crab3_task_dir_format.format('['+i_crab3_task_dir_basename+']'), ['1']))
               continue

            printout_lines = []

            printout_lines += [colored_text(crab3_task_dir_format.format('['+i_crab3_task_dir_basename+']'), '1')]
            printout_lines += ['| finished = '+colored_text('{:6.2f}%'.format(100. * float(jobs_finished) / float(jobs_total)), ['1', ('92' if jobs_finished == jobs_total else '93')])]

            other_states_desc = []

            for i_state in sorted(task_status_dict.keys()):

                if i_state == 'finished': continue

                other_states_desc += [(i_state+' = {:.2f}%').format(100.* (float(task_status_dict[i_state]) / float(jobs_total)))]

            if len(other_states_desc) != 0: printout_lines += ['('+', '.join(other_states_desc)+')']

            if (jobs_finished == jobs_total): iFinished += 1

            if opts.hadd and (opts.force_hadd or (jobs_finished == jobs_total)):

               if opts.force_hadd and (jobs_finished != jobs_total):
                  WARNING(log_prx+'forcing merging of crab outputs (jobs finished: '+str(jobs_finished)+' / '+str(jobs_total)+')')

               if not i_crab3_task_dir_basename.startswith('crab_'):
                  WARNING(log_prx+'invalid name of crab3 task directory (does not start with "crab_"): '+i_crab3_task_dir_basename)
                  continue

               dset_key = i_crab3_task_dir_basename[len('crab_'):]

               if dset_key not in production_dict:
                  WARNING(log_prx+'dataset key "'+dset_key+'" not found in production JSON (task: '+i_crab3_task_dir+')')
                  continue

               dset_conf = production_dict[dset_key]

               if 'OutputPrePath' not in dset_conf:
                  WARNING(log_prx+'key "OutputPrePath" (str) not found in production JSON entry for dataset "'+dset_key+'"')
                  continue

               if 'OutputSplitFiles' not in dset_conf:
                  WARNING(log_prx+'key "OutputSplitFiles" (int) not found in production JSON entry for dataset "'+dset_key+'"')
                  continue

               if 'crab3' not in dset_conf:
                  WARNING(log_prx+'key "crab3" (dict) not found in production JSON entry for dataset "'+dset_key+'"')
                  continue

               if 'Data.outLFNDirBase' not in dset_conf['crab3']:
                  WARNING(log_prx+'key "crab3":"Data.outLFNDirBase" not found in production JSON entry for dataset "'+dset_key+'"')
                  continue

               OutputPrePath_str = str(dset_conf['OutputPrePath'])

               hadd_output_file   = os.path.abspath(OutputPrePath_str+'.root')
               hadd_output_splitN = str(dset_conf['OutputSplitFiles'])

               if int(hadd_output_splitN) <= 0:
                  WARNING(log_prx+'invalid (non-positive) value for configuration parameter \"OutputSplitFiles\" for dataset "'+dset_key+'"')
                  continue

               tier2_maindir = opts.Tier2_prepath+'/'+dset_conf['crab3']['Data.outLFNDirBase']+'/*/'+i_crab3_task_dir_basename

               crab3_output_files_wildcard = tier2_maindir+'/*/*/*.root'

               crab3_output_files = glob.glob(crab3_output_files_wildcard)
               crab3_output_files = [os.path.abspath(_tmp) for _tmp in crab3_output_files]

               if not opts.force_hadd:

                  if len(crab3_output_files) != jobs_total:
                     WARNING(log_prx+'number of .root files in storage area ('+str(len(crab3_output_files))+') differs from total number of jobs in the crab3 task ('+str(jobs_total)+'): '+tier2_maindir)
                     continue

               else:
                  WARNING(log_prx+'forcing merging of '+str(len(crab3_output_files))+' crab3 outputs (.root files)')

               # expected output files
               exp_output_files = []

               if int(hadd_output_splitN) == 1:
                   exp_output_files = [hadd_output_file]
               else:
                   exp_output_files = [hadd_output_file[:-len('.root')]+'_'+str(tmp_i)+'.root' for tmp_i in range(int(hadd_output_splitN))]

               exp_output_files_exist = 0

               for i_exp_output_files in exp_output_files:
                   if os.path.exists(i_exp_output_files): exp_output_files_exist += 1

               if exp_output_files_exist == len(exp_output_files):

                  printout_lines += [colored_text('[outputs already merged]', ['1', '95'])]
                  print ' '.join(printout_lines)

                  continue

               elif exp_output_files_exist > 0:

                  print ' '.join(printout_lines)

                  for i_exp_output_files in exp_output_files:
                      if os.path.exists(i_exp_output_files):
                         WARNING(log_prx+'target output file already exists (skipping hadd call): '+i_exp_output_files)

                  continue

               printout_lines += [colored_text('[merging outputs...]', ['1', '93'])]

               print ' '.join(printout_lines)

               # output to Tier-2
               if OutputPrePath_str.startswith(opts.Tier2_prepath):

                  if not os.path.isdir(os.path.dirname(hadd_output_file)):
                     WARNING(log_prx+'target Tier-2 directory for output file does not exist: '+os.path.dirname(hadd_output_file))
                     continue

                  # temporary output files in current local directory (later copied to Tier-2 and removed)
                  hadd_output_file_tmp = os.path.abspath(os.path.basename(hadd_output_file))

                  exp_output_files_tmp = []

                  if int(hadd_output_splitN) == 1:
                      exp_output_files_tmp = [hadd_output_file_tmp]
                  else:
                      exp_output_files_tmp = [hadd_output_file_tmp[:-len('.root')]+'_'+str(tmp_i)+'.root' for tmp_i in range(int(hadd_output_splitN))]

                  for i_exp_output_files_tmp in exp_output_files_tmp:
                      if os.path.exists(i_exp_output_files_tmp):
                         WARNING(log_prx+'target path to local (temporary) hadd output already exists: '+i_exp_output_files_tmp)
                         continue

#                   # check available local disk space
#                   statvfs = os.statvfs(os.path.dirname(hadd_output_file_tmp))
#                   local_avail_bytes = int(statvfs.f_frsize * statvfs.f_bavail)
#
#                   exp_output_bytes = sum([os.stat(_tmp_f).st_size for _tmp_f in crab3_output_files])
#
#                   if exp_output_bytes > local_avail_bytes:
#                      WARNING(log_prx+'available local disk space ('+str(local_avail_bytes)+') too small for temporary storage of output files ('+exp_output_bytes+') for dataset "'+dset_key+'"')
#                      continue

                  # create temporary output file(s)
                  EXE('hadd_wrapper.py '+hadd_output_file_tmp+' "'+crab3_output_files_wildcard+'" -s '+hadd_output_splitN, verbose=opts.verbose)

                  # copy temporary output file(s) to Tier-2 and remove them locally
                  for i_tmp in range(len(exp_output_files_tmp)):

                      if not os.path.isfile(exp_output_files_tmp[i_tmp]):
                         WARNING(log_prx+'temporary local output file not found: '+exp_output_files_tmp[i_tmp])
                         continue

                      TIER2_CP_CMD = 'gfal-copy -t 3600'
                      # https://github.com/cms-sw/cmssw/issues/26462
                      # https://github.com/grid-control/grid-control/issues/73
                      EXE('(eval `scram unsetenv -sh`; '+TIER2_CP_CMD+' '+exp_output_files_tmp[i_tmp]+' "srm://dcache-se-cms.desy.de:8443/srm/managerv2?SFN='+exp_output_files[i_tmp]+'"'+')', verbose=opts.verbose)
                      EXE('rm'        +' '+exp_output_files_tmp[i_tmp], verbose=opts.verbose)

               else:

                  if not os.path.isdir(os.path.dirname(hadd_output_file)):
                       EXE('mkdir -p '+os.path.dirname(hadd_output_file), verbose=opts.verbose)

                  EXE('hadd_wrapper.py '+hadd_output_file+' "'+crab3_output_files_wildcard+'" -s '+hadd_output_splitN, verbose=opts.verbose)

            else:

               print ' '.join(printout_lines)

               if opts.resubmit and ('failed' in task_status_dict):

                  cmd_resubmit = 'crab resubmit'
                  if len(RESUBMIT_OPTS) > 0: cmd_resubmit += ' '+(' '.join(set(RESUBMIT_OPTS)))
                  cmd_resubmit += ' '+i_crab3_task_dir

                  EXE(cmd_resubmit, verbose=opts.verbose)

            if opts.repeat:
               loopEnd = (iTotal == iFinished)
               print time.strftime("%H-%M-%S:"), "Sleeping for one hour."
               time.sleep(3600)
    ### ---------------
