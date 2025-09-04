#!python3
import sys
import argparse
import socket
import os
from os.path import join as pj
from datetime import datetime
import pickle
import time
from src.common_tools import cjson_load, add_index, mk_new_dir
from src.stats_database import UspexCalcDB

ap = argparse.ArgumentParser()
ap.add_argument("-f", "--recalc_folder", required=False, help='Calc folder')
ap.add_argument("-t", "--uspex_template", required=True, help="USPEX Template", default='uspex_templates/USPEX_py_001')
ap.add_argument("-n", "--total_calcs", required=False, help='Total calcs', default='2')
ap.add_argument("-p", "--max_parallel_calcs", required=False, help='Maximum parallel calcs', default='1')
ap.add_argument("-o", "--outfile_pattern", required=False, help='Output file pattern', default='log')
ap.add_argument("-s", "--sleep", required=False, help='Sleep time, in minutes', default='5')
ap.add_argument("--maxiter", required=False, help='Maximum number of iteration', default='5')



cwd = os.getcwd()
template = 'USPEX_MATLAB_001'
template_path = pj(cwd, 'uspex_templates', template)
results_path = mk_new_dir(pj(cwd, 'results_' + template))
total_calcs = 4
max_parallel = 2
sleep = 2# minutes
db_file = pj(results_path, 'database.pkl')

log_path = os.path.join(results_path, 'log')
os.makedirs(results_path, exist_ok=True)
sys.stdout = open(log_path, 'a')

calcsDB = UspexCalcDB(results_folder=results_path, template_folder=template_path, total_calcs=total_calcs)

while not os.path.isfile(os.path.join(results_path, 'DONE')) and \
        not os.path.isfile(os.path.join(results_path, 'STOP')):
    calcsDB.update()
    print('*' * 10 + '{:%Y-%m-%d %H:%M}'.format(datetime.now()) + '*' * 10 + '\n')
    calcsDB.submit_jobs(n_par_calcs=max_parallel)
    calcsDB.get_stats(verb=True)
    calcsDB.dump(filename_pkl=db_file, filename_en='energies.txt')
    sys.stdout.flush()
    sleep_sec = int(sleep * 60)
    time.sleep(sleep_sec)




#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# #
# # DEFAULT_RES_FOLDER = 'results'
# #
# ap = argparse.ArgumentParser()
# # ap.add_argument("-c", "--machine", required=False, help='Which cluster')11
# ap.add_argument("-f", "--results_folder", required=False, help='Calc path')
# ap.add_argument("-t", "--template_folder", required=True, help="Source USPEX calc path", default='USPEX_PY_LJ36')
# ap.add_argument("-m", "--maxcalcs", required=False, help='Maximum parallel calcs', default='1')
# # ap.add_argument("-o", "--outfile_pattern", required=False, help='Output file pattern', default='log')
# ap.add_argument("-s", "--sleep", required=False, help='Sleep time, in minutes', default='5')
# ap.add_argument("-n", "--total_calcs", required=False, help='Total uspex calcs to perform', default='5')
# # ap.add_argument("--machines", required=False, help='Machines cjson file', default='machines.cjson')
# # ap.add_argument("--scenarios", required=False, help='Scenarios cjson file', default='scenarios.cjson')
# # ap.add_argument("--min_mult", required=False, help='Use minimal multiplicity', action='store_true')
# input_kwargs = vars(ap.parse_args())
# #
# # if input_kwargs['recalc_folder'] is None:
# #     input_kwargs['recalc_folder'] = add_index(DEFAULT_RES_FOLDER)
# #
# log_path = os.path.join(input_kwargs['results_folder'], 'log')
# os.makedirs(input_kwargs['results_folder'], exist_ok=True)
# sys.stdout = open(log_path, 'a')
# #
# machines_dct = cjson_load(input_kwargs['machines'])
# # Trying to automatically determine the supercomputer if not set
# if not input_kwargs['machine']:
#     hostname = socket.gethostname()
#     input_kwargs['machine'] = 'local'
#     for mach_k, mach_v in machines_dct.items():
#         if mach_v['hostname'] in hostname:
#             input_kwargs['machine'] = mach_k
#     print('Machine is automatically determined as: ' + input_kwargs['machine'])
#     print('In case of error, specify the machine explicitly using -c key and edit the machines.cjson accordingly')
# #
# # Assigning input arguments to variables, minding the defaults
# db_file = os.path.join(input_kwargs['results_folder'], 'database.pkl')
# #
# if os.path.isfile(db_file):
#     with open(db_file, 'rb') as db_fid:
#         database = pickle.load(db_fid)
#     db_file = db_file.split('/')[-1]
# else:
#     database = UspexCalcDB(machines_json=machines_dct,
#                          **{k: input_kwargs[k] for k in ('results_folder', 'template_folder', 'total_calcs')})
#
#
# while not os.path.isfile(os.path.join(input_kwargs['results_folder'], 'DONE')) and \
#         not os.path.isfile(os.path.join(input_kwargs['results_folder'], 'STOP')):
#     database.update()
#     print('*' * 10 + '{:%Y-%m-%d %H:%M}'.format(datetime.now()) + '*' * 10 + '\n')
#     database.submit_jobs(n_par_calcs=int(input_kwargs['maxcalcs']))
#     database.get_stats(verb=True)
#     database.dump(filename_pkl=db_file, filename_en='energies.txt')
#     sys.stdout.flush()
#     sleep_sec = int(float(input_kwargs['sleep']) * 60)
#     time.sleep(sleep_sec)
#
# # sys.stdout.close()