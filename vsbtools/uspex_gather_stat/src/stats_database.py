import re
import pickle
import os
from cclib import parser
from datetime import datetime
from os.path import isfile, isdir, join
# from pathlib import Path
# from .ext_software_io import read_geoms, parse_uspex_calc, dict2gauformat
# from .uspex_calc import UspexSetup
from .common_tools import add_index, sh_execute, checkflag, cjson_load, mk_new_dir

HARTREE = parser.utils.convertor(1, 'hartree', 'eV')  # cclib-compatible

# termination_dict = {'init_gjf':  'initial_config.gjf',
#                     'normal':   {'flags' : 'Normal termination', 'tail': 1}, # To be modified for orb_coeff etc
#                     'restart':  {'corrector': {'route': {'scf': {'yqc': None}}}},
#                     'errors':   {'low_phys_mem': {'flags': 'could not allocate memory', 'tail': 20,
#                                                   'msg': 'Could not allocate memory',
#                                                   'corrector': {'command':{'mem': '4000-100<1000'}}},
#                                  'low_alloc_mem': {'flags': 'Out-of-memory error', 'tail': False,
#                                                    'msg': 'Insufficient memory allocated',
#                                                    'corrector': {'command':{'mem': '1000+1000<10000'}}}},
#                     'fails':    {'bad_input': {'flags': 'Error termination via.*l101', 'tail': 50,
#                                                'msg': 'Bad input file',
#                                                'corrector': {'route': {'scf': {'yqc': None}}}}}}


class UspexCalc:
    # machines = {'rurik': {'jobtemplate': 'job_templates/rurik.sh',
    #                       'runcmd': r'sbatch ', 'jobid': '(?<=job )\d+',
    #                       'queue': "squeue -u `whoami` | grep -Po '(^\s*)\d+' |grep -Po '\d+'"},
    #             'oleg': {'jobtemplate': 'job_templates/oleg.sh',
    #                      'runcmd': r'sbatch ', 'jobid': '(?<=job )\d+',
    #                      'queue': "squeue -u `whoami` | grep -Po '(^\s*)\d+' |grep -Po '\d+'"},
    #             'arkuda': {'jobtemplate': 'job_templates/arkuda.sh',
    #                        'runcmd': r'qsub ', 'jobid': '\d+(?=\.)', 'queue': "qstat | grep -Po '^\s*\d+'"},
    #             'local': {'jobtemplate': 'job_templates/local.sh',
    #                       'runcmd': r'', 'jobid': '\d+', 'queue': "ps -x | grep -Po '^\s*\d+'"}}

    def __init__(self, template_path=None, path=None, status='P', job_id=None, name=None,
                 jobfile=None, machine='uspex', out_fname='results1/Individuals', machine_dict=None):

        if not machine_dict:
            machine_dict = {'jobtemplate': 'job_templates/uspex.sh',
                                      'runcmd': r'', 'jobid': '\d+', 'queue': "ps -x | grep -Po '^\s*\d+'"}
        self.template_path = template_path
        self.job_id = job_id
        self.status = status
        self.path = path  # Full path
        self.name = name  # E.g. USPEX_MATLAB_001_1
        self.jobfile = jobfile
        self.machine = machine
        self.machine_dict = machine_dict
        self.out = out_fname

    def __repr__(self):
        return '< Gaussian task ' + self.name + ' at ' + self.path + '>'
    
    @property
    def gjffilename(self):
        return self.name + '.gjf'

    @property
    def full_gjf_path(self):
        return self.path + '/' + self.gjffilename

    def isrunning(self):
        if self.job_id:
            cur_queue = sh_execute(self.machine_dict['queue'])
            return self.job_id in cur_queue
        else:
            return False

    def copyfiles(self):
        if not os.path.isdir(self.path):
            os.makedirs(self.path)
        sh_execute('cp -r ' + self.template_path + '/* ' + self.path)
        sh_execute('cp ' + self.machine_dict['jobtemplate'] + ' ' + self.jobfile)
        sh_execute("chmod u+x " + self.path + '/' + self.jobfile)

    def submit_job(self):
        print('Submitting ' + self.name)
        submit_cmd = 'cd ' + self.path + ' && ' + self.machine_dict['runcmd'] + ' ' + self.jobfile
        submit_out = sh_execute(submit_cmd, add_errors=True)
        job_id = re.findall(self.machine_dict['jobid'], submit_out)[0]
        self.job_id = job_id
        print(job_id)


class UspexCalcDB(list):

    def __init__(self, machines_json='machines.cjson', template_folder=None, results_folder=None,
                 outfile_pattern='Individuals', machine='uspex', total_calcs=100):

        list.__init__(self)

        self.results_root = results_folder

        template_folder = template_folder

        mach_dcts = cjson_load(machines_json)
        machine_dict = mach_dcts[machine]

        if template_folder and isdir(template_folder):
            for i_calc in range(total_calcs):
                # Creating task entry
                template_name = template_folder.split('/')[-1]
                curr_dest_folder = mk_new_dir(results_folder + '/' + template_name, zerobased=True)
                curr_label = curr_dest_folder.split('/')[-1]
                task = UspexCalc(template_path=template_folder, path=curr_dest_folder, name=curr_label,
                 jobfile=curr_dest_folder + '/jobfile.sh', machine='uspex', out_fname='results1/Individuals')
                self.append(task)
        else:
            raise Exception('Error: destination_folder and (geoms_file or outfile_pattern) must be provided! ')

    def __str__(self):
        return '<UspexCalcDB: ' + str(len(self)) + ' UspexCalc objects in ' + self.results_root + '>'


    def update(self):
        for task in self:
            logfile = task.path + '/' + task.out
            if task.status in ('R', 'F', 'L'):
                # Skip running tasks
                if task.isrunning():
                    continue
                # Update energy and geometry info if at least one SCF is done
                if isfile(join(task.path, 'USPEX_IS_DONE')):
                    print(task.name + ': DONE')
                    task.status = 'D'
                    sh_execute('cd ' + task.path + ' && ./clean')
                    continue

    def submit_jobs(self, n_par_calcs=3):
        submitted = new = 0
        for i, task in enumerate(self):
            if task.status == 'R':
                submitted += 1
            if task.status == 'P':
                task.copyfiles()
                try:
                    task.submit_job()
                except IndexError:
                    print('Failed to submit: ' + task.name)
                else:
                    self[i].status = 'R'
                    submitted += 1
                    new += 1
            if submitted >= n_par_calcs:
                break
        print('Submitted = ' + str(new) + '; Total running = ' + str(submitted))
        return new, submitted

    def get_stats(self, stats_file='stats.txt', verb=False):


        stats_file = join(self.results_root, stats_file)
        stats = {'P': 0, 'R': 0, 'D': 0, 'F': 0, 'L':0}  # Pending, Running, Done, Failed, Loaded
        fullstory = []

        for task in self:
            stats[task.status] += 1
            task_story = '{:16s}'.format(task.name) + '\t' + task.status
            fullstory.append(task_story)

        if verb:
            print('\n'.join(fullstory))

        if stats_file:
            with open(stats_file, 'a+') as stats_fid:
                for tsk in fullstory:
                    stats_fid.write(tsk + '\n')

        if stats['P'] == 0 and stats['R'] == 0 and stats['L'] == 0:
            print('No running or pending tasks!\n')
            if stats['F'] > 0:
                print('We\'ve done all we could with this input\n')
                print('To recover {} failed tasks, try to restart with altered scenarios'.format(stats['F']))
                sh_execute('touch ' + self.results_root + '/STOP')
            else:
                print('Congratulations! All done!')
                sh_execute('touch ' + self.results_root + '/DONE')

        print("Pending: {}, Running: {}, Done: {}, Failed: {}, Loaded: {}".format(
            stats['P'], stats['R'], stats['D'], stats['F'], stats['L']))
        return stats

    def dump(self, filename_pkl='database.pkl', filename_en='energies.txt'):
        filename_en = join(self.results_root, filename_en)
        with open(filename_en, 'w') as en_fid:
            en_fid.write('\n' + '*' * 10 + ' {:%Y-%m-%d %H:%M} '.format(datetime.now()) + '*' * 10 + '\n')
            en_fid.write('task_name\tlast_generation\tbest_energy:\n')
            for task in self:
                if isfile(join(task.path, 'results1', 'BESTIndividuals')):
                    best_en = sh_execute('cd ' + task.path + ' && ./bestenergy.sh')
                    str = '%-10s' % task.name + ': ' + best_en
                    en_fid.write(str + '\n')
                elif task.status == 'D':
                    print('!!! This is weird, task "' + task.name + '" is done, but has no ccdata')
                    # last_ccdata, atoms = parse_gout(task.path + '/' + task.out)
                    # task.ccdata = last_ccdata
                    # en_fid.write('%-10s' % task.name + ': ' + '%-.10f' % (task.ccdata.scfenergies[-1]) + '\n')
        with open(filename_pkl, 'wb') as pkl_fid:
            pickle.dump(self, pkl_fid)


# if __name__ == '__main__':
#     init_gjf = Gjf('../testfolder/complicated.gjf')
#
#     # jobtemplate = '../testfolder/job_arkuda.sh'
#     template_folder = '../testfolder/job_oleg.sh'
#     dest_fold = '../DESTINATION'
#     geoms_file = '../testfolder/POSCARS'
#     database = GauCalcDB(scenarios='../scenarios.cjson', geoms_file=geoms_file, recalc_folder=dest_fold,
#                          machines_json='../machines.cjson', machine='rurik')
#     # database.submitjobs()
#     database.get_stats(verb=True)
#     database.dump()
#
#     del database
#
#     with open('database.pickle', 'rb') as db_fid:
#         db = pickle.load(db_fid)
#
#     # db
#
#     print(db)
#     # with open('pkl.pkl', 'rb') as f:
#     #     x = pickle.load(f)
#     # print(x[0].gjf)
#     # pass
