#!/usr/bin/env python
import os, subprocess

def KILL(log):
    raise SystemExit('\n '+'\033[1m'+'@@@ '+'\033[91m'+'FATAL'  +'\033[0m'+' -- '+log+'\n')
# --

def WARNING(log):
    print '\n '+'\033[1m'+'@@@ '+'\033[93m'+'WARNING'+'\033[0m'+' -- '+log+'\n'
# --

def EXE(cmd, suspend=True, verbose=False, dry_run=False):
    if verbose: print '\033[1m'+'>'+'\033[0m'+' '+cmd
    if dry_run: return

    _exitcode = os.system(cmd)

    if _exitcode and suspend: raise SystemExit(_exitcode)

    return _exitcode
# --

def get_output(cmd, permissive=False):
    prc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    out, err = prc.communicate()

    if (not permissive) and prc.returncode:
       KILL('get_output -- shell command failed (execute command to reproduce the error):\n'+' '*14+'> '+cmd)

    return (out, err)
# --

def command_output_lines(cmd, stdout=True, stderr=False, permissive=False):

    _tmp_out_ls = []

    if not (stdout or stderr):
       WARNING('command_output_lines -- options "stdout" and "stderr" both set to FALSE, returning empty list')
       return _tmp_out_ls

    _tmp_out = get_output(cmd, permissive=permissive)

    if stdout: _tmp_out_ls += _tmp_out[0].split('\n')
    if stderr: _tmp_out_ls += _tmp_out[1].split('\n')

    return _tmp_out_ls
# --

def rreplace(str__, old__, new__, occurrence__):
    li_ = str__.rsplit(old__, occurrence__)
    return new__.join(li_)
# --

def which(program, permissive=False, verbose=False):

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)

    _exe_ls = []

    if fpath:
        if is_exe(program): _exe_ls += [program]

    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)

            if is_exe(exe_file): _exe_ls += [exe_file]

    _exe_ls = list(set(_exe_ls))

    if len(_exe_ls) == 0:
       log_msg = 'which -- executable not found: '+program

       if permissive:
          if verbose: WARNING(log_msg)
          return None

       else:
          KILL(log_msg)

    if len(_exe_ls) >  1:
       if verbose: WARNING('which -- executable "'+program+'" has multiple matches: \n'+str(_exe_ls))

    return _exe_ls[0]
# --

def is_int(value):

    try: int(value)
    except ValueError: return False

    return True
# --

def is_float(value):

    try: float(value)
    except ValueError: return False

    return True
# --

def colored_text(txt, keys=[]):

    _tmp_out = ''

    for _i_tmp in keys:
        _tmp_out += '\033['+_i_tmp+'m'

    _tmp_out += txt

    if len(keys) > 0: _tmp_out += '\033[0m'

    return _tmp_out
# --
