#!/usr/bin/env python
import argparse, os, glob, math

from common import *

def hadd_sequence(output_file, input_files, options, remove_inputs=False):

    if os.path.exists(output_file):
       KILL(log_prx+'hadd_sequence: path to target output file already exists: '+output_file)

    OUTPUT_FILE = os.path.abspath(output_file)

    if len(input_files) < 2:
       KILL(log_prx+'hadd_sequence: list of input files contains less than 2 objects: '+str(input_files))

    MAX_INPUT_LENGTH = int(256)

    INPUT_BLOCKS = [input_files[_tmp_idx:_tmp_idx + MAX_INPUT_LENGTH] for _tmp_idx in xrange(0, len(input_files), MAX_INPUT_LENGTH)]

    hadd_cmds = []

    if len(INPUT_BLOCKS) < 2:

       hadd_cmd = 'hadd'
       if options: hadd_cmd += ' \\\n '+options

       hadd_cmd += ' \\\n '+OUTPUT_FILE
       hadd_cmd += ' \\\n '+' \\\n '.join(input_files)

       hadd_cmds += [hadd_cmd]

       if remove_inputs:
          rminp_cmd  = 'rm -f'
          rminp_cmd += ' \\\n '+' \\\n '.join(input_files)

          hadd_cmds += [rminp_cmd]

    else:

        tmp_outputs = []

        for _tmp_block_idx in range(len(INPUT_BLOCKS)):

            i_block = INPUT_BLOCKS[_tmp_block_idx]

            i_out_prepath = os.path.splitext(OUTPUT_FILE)[0]
            i_out_ext     = os.path.splitext(OUTPUT_FILE)[1]

            i_tmp_output = i_out_prepath+'_tmp'+str(_tmp_block_idx)+i_out_ext

            hadd_cmds += hadd_sequence(i_tmp_output, i_block, options, remove_inputs=False)

            tmp_outputs += [i_tmp_output]

        hadd_cmds += hadd_sequence(OUTPUT_FILE, tmp_outputs, options, remove_inputs=True)

    return hadd_cmds
# --

#### main
if __name__ == '__main__':
    ### args
    parser = argparse.ArgumentParser()
    parser.add_argument('--hadd-opts', dest='hadd_opts',
                        action='store', default='',
                        help='string containing command-line options for hadd')
    parser.add_argument('-N', '--N-input', dest='N_input',
                        action='store', type=int, default=-1,
                        help='expected number of input files (will abort if this does not match the actual number of input files)')
    parser.add_argument('-s', '--splitlevel', dest='splitlevel',
                        action='store', type=int, default=1,
                        help='split output in N separate files')
    parser.add_argument('--sge-submit', dest='sge_submit',
                        action='store_true', default=False,
                        help='execute as job on SGE batch system')
    parser.add_argument('--sge-output', dest='sge_output',
                        action='store', default='',
                        help='path to output directory containing SGE submission scripts')
    parser.add_argument('--sge-name', dest='sge_name',
                        action='store', default='',
                        help='basename of SGE submission script')
    parser.add_argument('-v', '--verbose', dest='verbose',
                        action='store_true', default=False,
                        help='enable verbose mode')
    parser.add_argument('-d', '--dry-run', dest='dry_run',
                        action='store_true', default=False,
                        help='enable dry-run mode')
    opts, opts_unknown = parser.parse_known_args()

    log_prx = os.path.basename(__file__)+' -- '
    ###

    which('hadd')

    if len(opts_unknown) < 2:
       KILL(log_prx+'invalid list of input/output files (size='+str(len(opts_unknown))+'): '+str(opts_unknown))

    OUTPUT = opts_unknown[0]

    IFILE_LS = []
    for tmp_if in opts_unknown[1:]:
        IFILE_LS += (glob.glob(tmp_if) if '*' in tmp_if else [tmp_if])

    for tmp_if in IFILE_LS:
        if not os.path.isfile(tmp_if):
           KILL(log_prx+'target input file not found: '+tmp_if)

    if opts.splitlevel <= 0:
       KILL(log_prx+'invalid value for "splitlevel" parameter [-s]: '+str(opts.splitlevel))

    if opts.N_input != -1:
       if len(IFILE_LS) != opts.N_input:
          KILL(log_prx+'unexpected number of input files: (expected={:d}, provided={:d})'.format(opts.N_input, len(IFILE_LS)))
    ###

    CMD_LS = []

    EXT = 'root'

    if not OUTPUT.endswith('.'+EXT):
       KILL(log_prx+'invalid extension for target output file: '+OUTPUT)

    N_BLOCKS = min(len(IFILE_LS), opts.splitlevel)

    output_postfix_format = '{:d}' # '{:0'+str(1+int(math.log10(N_BLOCKS - 1)))+'d}'

    ifile_blocks = [IFILE_LS[tmp_i::N_BLOCKS] for tmp_i in xrange(N_BLOCKS)]

    for i_block in range(len(ifile_blocks)):

        output_postfix = ('' if (N_BLOCKS == 1) else ('_'+output_postfix_format).format(i_block))

        i_OUTPUT = OUTPUT[:-len('.'+EXT)]+output_postfix+'.'+EXT

        CMD_LS += [hadd_sequence(i_OUTPUT, ifile_blocks[i_block], opts.hadd_opts)]
    # ---

    if opts.sge_submit:

       if not opts.sge_output: KILL(log_prx+'unspecified path to output directory for SGE submission scripts [--sge-output]')
       if not opts.sge_name  : KILL(log_prx+'unspecified basename of SGE submission script [--sge-name]')

       if not os.path.isdir(opts.sge_output):
          EXE('mkdir -p '+opts.sge_output, verbose=opts.verbose, dry_run=opts.dry_run)

       SGE_OUTDIR = os.path.abspath(opts.sge_output)

       o_exe = SGE_OUTDIR+'/'+opts.sge_name+'.sh'

       if os.path.exists(o_exe):
          KILL(log_prx+'target path to SGE submission script already exists: '+o_exe)

       o_file = open(o_exe, 'w')

       o_exe_opts = """#!/bin/sh
#$ -l os=sld6
#$ -l site=hh
#$ -cwd
#$ -V
##$ -pe local 8-12
#$ -l h_vmem=8G
#$ -l h_fsize=8G
"""
       o_file.write(o_exe_opts)

       ODIR_LOG = SGE_OUTDIR+'/sge'
       if not os.path.isdir(ODIR_LOG):
          EXE('mkdir -p '+ODIR_LOG, verbose=opts.verbose, dry_run=opts.dry_run)

       SGE_OPTS = [
         '#$ -N '+opts.sge_name,
         '#$ -o '+ODIR_LOG,
         '#$ -e '+ODIR_LOG,
       ]
       for _opt in SGE_OPTS: o_file.write(_opt+'\n')

       for i_cmds in CMD_LS:
           for icmd in i_cmds:
               o_file.write('\n'+icmd+'\n')

       o_file.close()

       EXE('chmod 700 '+o_exe, verbose=opts.verbose, dry_run=opts.dry_run)
       EXE('qsub '     +o_exe, verbose=opts.verbose, dry_run=opts.dry_run)

    else:

       for i_cmds in CMD_LS:
           for icmd in i_cmds:
               EXE(icmd, verbose=opts.verbose, dry_run=opts.dry_run)
#### ----
