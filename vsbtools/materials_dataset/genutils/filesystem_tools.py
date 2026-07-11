import os
from subprocess import Popen, PIPE, TimeoutExpired


def sh_execute(cmdstring, add_errors=True):
    proc = Popen(cmdstring.__str__(), stdout=PIPE, stderr=PIPE, shell=True)
    try:
        o, e = proc.communicate(timeout=15)
    except TimeoutExpired:
        proc.kill()
        o, e = proc.communicate()
    if add_errors:
        return (o.decode('utf-8') + e.decode('utf-8')).strip()
    else:
        return o.decode('utf-8').strip()


def add_index(fname, respect_file_extension=True, zerobased=False):

    if respect_file_extension and '.' in fname.split('/')[-1]:
        step = fname.split('.')
        dot_file_extension = '.' + step[-1]
        fname = '.'.join(step[:-1])
    else:
        dot_file_extension = ''
    if zerobased:
        fname_next = fname + '_0'
    else:
        fname_next = fname
    k = 1
    while os.path.exists(fname_next + dot_file_extension):
        fname_next = fname + '_' + str(k)
        k += 1
    return fname_next + dot_file_extension


def mk_new_dir(dirname, zerobased=False):
    newname = add_index(dirname, zerobased=zerobased)
    os.makedirs(newname, exist_ok=True)
    return newname
