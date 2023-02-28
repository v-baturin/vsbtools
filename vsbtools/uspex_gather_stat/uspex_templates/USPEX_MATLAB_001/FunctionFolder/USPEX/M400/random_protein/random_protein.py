#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Feb 03, 2014

@author: mrakitin
"""

import os
import sys
# from random import randint
from functions import pseudo_random_angles, unlimited_random_angles, custom_random_angles

# -------------------------------------------------------------------------------
# Arguments processing:
# print len(sys.argv)

prefix = ''
if len(sys.argv) < 5:
    print 'Required parameters are not specified! Exit.'
    exit(1)
else:
    prefix = sys.argv[1] + '.'
    try:
        attempts = int(sys.argv[2])
    except:
        print 'The second argument must be integer! Exit.'
        exit(2)

    random_type = 'pseudo'
    if sys.argv[3] == 'unlimited':
        random_type = 'unlimited'

    if sys.argv[4] == 'memory':
        enable_memory = True
    elif sys.argv[4] == 'nomemory':
        enable_memory = False
    else:
        print 'The fourth argument must be either "memory" or "nomemory"! Exit.'
        exit(3)

if not os.path.exists(prefix + 'make0'):
    print 'File ' + prefix + 'make0 does not exist! Exit.'
    exit(4)
# -------------------------------------------------------------------------------



# -------------------------------------------------------------------------------
# Variables and parameters:
# prefix     = 'ala20.'
make0_file = prefix + 'make0'
make_file = prefix + 'make'
key_file = prefix + 'key'
int_file = prefix + 'int'
int2_file = prefix + 'int_2'
int3_file = prefix + 'int_3'
xyz2_file = prefix + 'xyz_2'
xyz3_file = prefix + 'xyz_3'
pdb_file = prefix + 'pdb'
tmp_file = prefix + 'tmp'
log_file = prefix + 'log'

precision_optirot = 0.001
precision_minimize = 0.00001

marks_list = []
residues_list = []
short_marks = []
carbon_list = []
angles = []
stable_angles = []
stable_energy = 10000000
# attempts      = 10

nums_width = 8
sep = '-------------------------------------------'
# -------------------------------------------------------------------------------


# -------------------------------------------------------------------------------
# Find markers and store it in a list:
f = open(make0_file, 'rb')
content_make0 = f.readlines()
f.close()

for i, row in enumerate(content_make0):
    if row.find(' 0.001 ') > 0:
        start_marks = i
        break

j = 1.0
for i in xrange(start_marks, len(content_make0)):
    j_frac = float(j / 1000.0)
    if content_make0[i].find(str(j_frac)) > 0:
        # short_marks.append(j_frac)    # THIS VERSION IS BUGGY: For 0.10, 0.20, 0.30, etc. it returns 0.1, 0.2, 0.3, etc.! Never use this! Fix is on the next line:
        residues_list.append(content_make0[i].split()[0])
        short_marks.append('%0.3f' % j_frac)
        marks_list.append('%0.4f' % j_frac)
    else:
        break
    j += 1

# print 'Markers list :', marks_list
# print 'Short markers:', short_marks
# print 'residues_list:', residues_list

# Create *.int file with markers and get the corresponding C atom numbers:
print '> Cleaning the directory before execution...'
os.system('rm -f ' + prefix + 'seq* ' + prefix + 'int* ' + prefix + 'xyz* > /dev/null 2>&1')

print '> protein.x test execution to find alpha-carbon numbers...'
os.system('protein.x < ' + make0_file + ' > /dev/null 2>&1')

f = open(int_file, 'rb')
content_int = f.readlines()
f.close()

for i, row in enumerate(content_int):
    for mark in marks_list:
        if row.find(' ' + mark + ' ') > 0:
            carbon_list.append(row.split()[0])
            break

print '\tAlpha-carbon numbers found:', ' '.join(carbon_list)


# -------------------------------------------------------------------------------


# -------------------------------------------------------------------------------
# Function to generate a structure and optimize it.
# Input :
#         (1) passed_angles - if empty, pseudo-random or random angles are used.
#                             Otherwise, customized angles are used.
#         (2) pickup - if set to True, start with existing case_stable.make file.
# Output: energy and a list of angle pairs.

def generate_structure(passed_angles, pickup):
    # print '> Removing files obtained during test execution...'
    os.system('rm -f ' + prefix + 'seq* ' + prefix + 'int* ' + prefix + 'xyz* > /dev/null 2>&1')

    if pickup == False:
        # print '> Generate angles randomly from 7 possible pairs...'
        make_random = ''
        for i in xrange(start_marks):
            make_random += content_make0[i]

        j = 0
        for i in xrange(start_marks, start_marks + len(short_marks)):
            if passed_angles == []:
                if random_type == 'pseudo':
                    rand_angles = pseudo_random_angles()
                elif random_type == 'unlimited':
                    rand_angles = unlimited_random_angles()
            else:
                rand_angles = custom_random_angles(passed_angles)

            make_random += content_make0[i].replace(' ' + str(short_marks[j]) + ' ',
                                                    ' ' + str(rand_angles['phi']) + ' ').replace(
                ' -' + str(short_marks[j]) + ' ', ' ' + str(rand_angles['psi']) + ' ')
            j += 1

        for i in xrange(start_marks + len(short_marks), len(content_make0)):
            make_random += content_make0[i]

        f = open(make_file, 'wb')
        f.writelines(make_random)
        f.close()

    elif pickup == True:
        os.system('cp -f ' + prefix[:-1] + '_stable.make' + ' ' + make_file)

    os.system('protein.x < ' + make_file + ' > /dev/null 2>&1')


    # Optimize angles(!) and get energy and resulted angles for marked atoms:
    # print '> optirot.x execution with random angles...'
    os.system('optirot.x ' + int_file + ' 0 ' + str(precision_optirot) + ' > ' + tmp_file)

    f = open(tmp_file, 'rb')
    content_tmp = f.readlines()
    f.close()

    # Example:
    # Final Function Value :           63.2690

    for i, row in enumerate(content_tmp):
        if row.find('Final Function Value :') >= 0:
            energy = float(row.split(':')[1].strip())
            break

    # print '\tEnergy =', energy, 'kcal/mol'

    # Convert .int_2 to .xyz_2:
    os.system('intxyz.x  ' + int2_file + ' > /dev/null 2>&1')


    # Minimize and get energy and resulted angles for marked atoms:
    os.system('minimize.x ' + xyz2_file + ' ' + str(precision_minimize) + ' > ' + tmp_file)

    f = open(tmp_file, 'rb')
    content_tmp = f.readlines()
    f.close()

    for i, row in enumerate(content_tmp):
        if row.find('Final Function Value :') >= 0:
            energy_min = float(row.split(':')[1].strip())
            break
    if energy_min < energy:
        energy = energy_min

    os.system('xyzint.x ' + xyz3_file + ' A > /dev/null 2>&1')
    os.system('rm -f ' + pdb_file + ' > /dev/null 2>&1')
    os.system('xyzpdb.x ' + xyz3_file + ' ' + key_file + ' > /dev/null 2>&1')

    f = open(int3_file, 'rb')
    content_int3 = f.readlines()
    f.close()

    angles = []
    counter = 0
    for i, row in enumerate(content_int3):
        for atom_num in carbon_list:
            if row.find(' ' + atom_num + ' ') >= 0 and row.find(' ' + atom_num + ' ') < nums_width:
                # print row[1:nums_width]

                phi = float(row.split()[8])  # Column #8 in case.int_2 file
                psi = float(content_int3[i + 1].split()[8])  # Column #8 in case.int_2 file, next line

                if psi > 90.0:
                    psi = psi - 180.0
                elif psi < 0.0:
                    psi = psi + 180.0

                angles_dict = {
                    'residue': residues_list[counter],
                    'phi': phi,
                    'psi': psi,
                }
                counter += 1

                angles.append(angles_dict)
                break

    '''
    print '\tOptimized angles:'
    for i, angle in enumerate(angles):
        print '\t\tAlpha-carbon #%i: phi = %7.2f, psi = %7.2f' % (i+1, angle['phi'], angle['psi'])
    '''
    return_dict = {
        'energy': energy,
        'angles': angles,
    }
    return return_dict


# -------------------------------------------------------------------------------



# -------------------------------------------------------------------------------
# Generate lots of semi-random structures, optimize them and select with the lowest energy.

print '> Generate and optimize ' + str(attempts) + ' different structures using %s-random scheme...' % (random_type)
# result = generate_structure()
# exit(999)

# f = open(log_file, 'wb')

passed_angles = []

pickup = False
if os.path.exists(prefix[:-1] + '_stable.make'):
    pickup = True

for i in xrange(attempts):
    result = generate_structure(passed_angles, pickup)
    energy = result['energy']
    angles = result['angles']
    if energy < stable_energy:
        # os.system('intxyz.x  ' + int3_file + ' > /dev/null 2>&1')

        os.system('cp -f ' + make_file + ' ' + prefix[:-1] + '_stable.make')
        os.system('cp -f ' + int3_file + ' ' + prefix[:-1] + '_stable.int')
        os.system('cp -f ' + xyz3_file + ' ' + prefix[:-1] + '_stable.xyz')
        os.system('cp -f ' + pdb_file + ' ' + prefix[:-1] + '_stable.pdb')

        stable_energy0 = stable_energy
        stable_energy = energy
        stable_angles = angles

        if enable_memory == True:
            passed_angles = stable_angles

        for_write = '\tEnergy lowered at step #%i: %.4f ---> %.4f kcal/mol' % (i + 1, stable_energy0, stable_energy)
        print for_write
        # f.write(for_write + '\n')

        for j, angle in enumerate(stable_angles):
            print '\t\tAlpha-carbon #%3i: phi = %7.2f, psi = %7.2f' % (j + 1, angle['phi'], angle['psi'])

    # Pickup stable structure once on the first step, then ignore it:
    pickup = False


    # Report current stage:
    f = open(prefix + 'current_stage', 'wb')
    current_content = 'Step #%i of %i\n\n' % (i + 1, attempts)
    current_content += 'Current energy: %f kcal/mol\n' % (energy)
    current_content += 'Angles at current energy:\n'
    for i, angle in enumerate(angles):
        current_content += '\tAlpha-carbon #%3i: phi = %7.2f, psi = %7.2f\n' % (i + 1, angle['phi'], angle['psi'])
    f.writelines(current_content)
    f.close()

print '\n' + sep.replace('-', '=')
print 'The lowest energy:', stable_energy, 'kcal/mol'

print 'Angles at the lowest energy:'
for i, angle in enumerate(stable_angles):
    print '\tAlpha-carbon #%3i: phi = %7.2f, psi = %7.2f' % (i + 1, angle['phi'], angle['psi'])
print ''

f = open('input_stable.angles', 'wb')
f.write('%-8s %7s %7s\n' % ('residue', 'phi', 'psi'))
for i, angle in enumerate(stable_angles):
    f.write('%-8s %7.2f %7.2f\n' % (angle['residue'], angle['phi'], angle['psi']))
f.close()

'''
#-------------------------------------------------------------------------------
# Convert .int file to .xyz for MOLDEN processing:

print ''
print '> intxyz.x execution to convert '+ int2_file +' file to '+ xyz2_file +' for MOLDEN...'
os.system('intxyz.x  ' + int2_file + ' > /dev/null 2>&1')
print ''
#-------------------------------------------------------------------------------
'''

exit(0)
