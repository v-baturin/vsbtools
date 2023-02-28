#!/usr/bin/env python
# -*- coding: utf-8 -*-

import random, math, string,sys

print "==================================================="
print "==================================================="
print "=== Transition Path Sampling Shooter for LAMMPS ==="
print "==================================================="
print "==================================================="
import random, math, string,sys

aut = 2.4188843265e-17   # second
aul = 5.2917720859e-11   # meter
auT = 3.1577464e+5       # kelvin
aue = 4.35974417e-18     # J
kb1 = 1.3806488e-23      # J.K-1
kb2 = 1.987204118e-3     # Kcal.mol-1.K-1
avg = 6.0221415e23       # Avogadro number
aum = 9.10938291e-31     # Kg

class atom:
	def set_xyz(self, value):
		id,tag,type,charge,x,y,z,bin1,bin2,bin3 = value.split()
		#id,type,x,y,z,bin1,bin2,bin3 = value.split()
		self.id = string.atoi(id)
		self.tag = 1 #string.atoi(tag)
		self.type = string.atoi(type)
		self.charge = 1 #string.atof(charge)
		self.x,self.y,self.z = string.atof(x), string.atof(y), string.atof(z)
		self.bin1 = bin1
		self.bin2 = bin2
		self.bin3 = bin3
	def set_vel(self, value):
		id,vx,vy,vz = value.split()
		self.id = string.atoi(id)
		self.vx,self.vy,self.vz = string.atof(vx), string.atof(vy), string.atof(vz)

def VectorModule(v):
	return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def RandomGaussianNumber(deviation, mean):	# Box-Muller method
	ransq = 1.1
	while ransq>1.0:
		ran1 = 2 * random.uniform(0,1) - 1.0
		ran2 = 2 * random.uniform(0,1) - 1.0
		ransq = ran1*ran1 + ran2*ran2
	return ran1 * math.sqrt(-2.0 * math.log(ransq)/ransq) * deviation + mean # vel * s + m // m: mean && s: standard deviation, here m=0 && s=1

def GaussianVelocity(mass, temperature, m_factor): # mass in g/mol && temp. in K ---- J = kg.m**2.s**-2
	return RandomGaussianNumber(m_factor, 0.0) * math.sqrt((kb1 * temperature / mass) * 1.0e3 * avg) * (1.0e+10/1.0e+15)
#       A.sec**-1                                J.K**-1   K           g       kg   Avogadro   = A.sec**-1 

def RotateTwoVectors(q1, q2, v1, v2, mass1, mass2, dvx, dvy, dvz):
	dx = q2[0] - q1[0]
	dy = q2[1] - q1[1]
	dz = q2[2] - q1[2]

	distance = math.sqrt(dx*dx + dy*dy + dz*dz)	# distance between atoms

	if distance == 0.0:
		print "!!!!!!!!!!! same atom ... something wrong"
		sys.exit()

	v1_new = [0.0,0.0,0.0]
	v2_new = [0.0,0.0,0.0]

	v1_new[0] = v1[0] + (dvx * (dx/distance))	# act +dv on 1st atom
	v1_new[1] = v1[1] + (dvy * (dy/distance))
	v1_new[2] = v1[2] + (dvz * (dz/distance))

	mass_scale = mass1/mass2

	v2_new[0] = v2[0] - (dvx * mass_scale * (dx/distance))	# act -dv on 2nd atom
	v2_new[1] = v2[1] - (dvy * mass_scale * (dy/distance))
	v2_new[2] = v2[2] - (dvz * mass_scale * (dz/distance))

	return v1_new, v2_new	# the velocities returned here need to be rescaled by sqrt[T_desired/T_calculated]

def ScaleVelocity(To,Tn,velocity):
	return math.sqrt(To/Tn) * velocity

######################### input_shoot FILE #########################
print "     ----------------------"
print "---> Open file: input_shoot"
print "     ----------------------"

input_shoot = open('input_shoot','r')	# file contains shooting parameters
types = int(input_shoot.readline())
print types, "atomic type(s)"
ntypes = []
masses = []
labels = []
natoms = 0
for i in range(types):
	a,b,c = input_shoot.readline().split()
	ntypes.append(string.atoi(a))
	masses.append(string.atof(b))
	labels.append(c)
	print "\t", ntypes[i], labels[i], "atoms / mass =", masses[i]
	natoms = natoms + ntypes[i]
print "\t", natoms, "atoms in total\n"

print "     ------------------"
print "---> Open file: lastrun"
print "     ------------------"
print ".... check shooting direction"
lastrun = open('lastrun','r')
sucfail, direction = lastrun.readline().split()	# (1:suc) (-1:fail) (A:B->A) (B:A->B)
sucfail = string.atoi(sucfail)
lastrun.close()
print "     -------------------"
print "---> Close file: lastrun"
print "     -------------------"
print " "
amplitudeAB, magnitudeAB = input_shoot.readline().split()	# for A-->B
amplitudeBA, magnitudeBA = input_shoot.readline().split()	# for B-->A
amplitudeAB, magnitudeAB = string.atof(amplitudeAB),string.atof(magnitudeAB)
amplitudeBA, magnitudeBA = string.atof(amplitudeBA),string.atof(magnitudeBA)
if direction == "A":
	print ".... last path was B--->A"
	if sucfail == 1:
		print ".... last path was successful\n                   ----------"
		amplitude, magnitude = amplitudeAB, magnitudeAB
		print ".... new path is A--->B"
	elif sucfail == -1:
		print ".... last path was failed\n                   ------"
		amplitude, magnitude = amplitudeBA, magnitudeBA
		print ".... new path is B--->A"
	else:
		print "... something wrong !!!! check lastrun file ..."
		sys.exit()
elif direction == "B":
	print ".... last path was A--->B"
	if sucfail == 1:
                print ".... last path was successful\n                   ----------"
                amplitude, magnitude = amplitudeBA, magnitudeBA
		print ".... new path is B--->A"
        elif sucfail == -1:
                print ".... last path was failed\n                   ------"
                amplitude, magnitude = amplitudeAB, magnitudeAB
		print ".... new path is A--->B"
        else:
                print "... something wrong !!!! check lastrun file ..."
                sys.exit()
else:
	"... something wrong !!!! check lastrun file ..."
	sys.exit()
print " "
print "shooting parameters: amplitude = ", amplitude
print "                     magnitude = ", magnitude
print " "
myseed = string.atoi(input_shoot.readline().split()[0])
print ".... Random Number Generator Seed:", myseed
print ".... seeding RNG"
random.seed(myseed)     # RNG seed
input_shoot.close()
print "     -----------------------"
print "---> Close file: input_shoot"
print "     -----------------------"
print " "
######################### ascii.data FILE #########################
print "     ---------------------"
print "---> Open file: ascii.data (ASCII restart LAMMPS file)"
print "     ---------------------"
ascii_data = open('ascii.data','r')

list_atoms=[]
for i in range(natoms):
	item = atom()
	list_atoms.append(item)

######################### GET X Y Z ... #########################
print ".... collecting data: x y z"
while 1:
	line = ascii_data.readline()
	if "Atoms" in line:
		ascii_data.readline()	# skip blank line
		for i in range(natoms):
			item = atom()
			line = ascii_data.readline()	# read line that contains x y z
			item.set_xyz(line)		# set x y z
			list_atoms[item.id - 1] = item	# id starts at 1, python array index at 0
		break					# x y z read ... stop
######################### GET VX VY VZ #########################
print ".... collecting data: vx vy vz"
while 1:
	line = ascii_data.readline()
	if not line: break
	if "Velocities" in line:
		ascii_data.readline()   # skip blank line
		for i in range(natoms):
			line = ascii_data.readline()	# read line that contains vx vy vz
			a,b,c,d = line.split()
			list_atoms[string.atoi(a) - 1].set_vel(line)	# set vx vy vz
ascii_data.close()
print "     ----------------------"
print "---> Close file: ascii.data"
print "     ----------------------"
print " "
######################### KINETIC ENERGY  #########################
print ".... calculate kinetic energy"
ktot = 0.0
for item in list_atoms:
	velocity2 = item.vx**2 + item.vy**2 + item.vz**2
	ktot = ktot + 0.5 * masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * velocity2 * (1.0e-10/1.0e-15)**2 * (aut/aul)**2
#       a.u.                g/mol                   Kg       a.u mass    Avogadro    A/fs        (m/sec)**2             a.u. length/time
print "Kinetic energy =", ktot * aue, "[Joule]"
print "               =", ktot * aue * 6.24150974e18, "[eV]"
######################### TEMPERATURE #########################
print ".... calculate temperature"
temp = 2.0 * ktot * aue / (kb1 * (3*natoms - 3))    # T = 2 * E / (kB * (3N - 3)) --- system temperature
print "                              ==========="
print "Temperature (Instantaneous) =", temp, "[K] --- before shoot move"
print "                              ==========="
########################## TOTAL LINEAR MOMENTUM #########################
print ".... check for total linear momentum conservation (before shooting)"
P = 0.0	# total linear momentum
Px = 0.0
Py = 0.0
Pz = 0.0
CMx = 0.0
CMy = 0.0
CMz = 0.0
CM = 0.0
Mtot = 0.0
average_v = 0.0
average_p = 0.0
for item in list_atoms:
	average_v = average_v + VectorModule([item.vx,item.vy,item.vz])	# Ang.fs-1

	Px = Px + masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * item.vx * (1.0e-10/1.0e-15) * (aut/aul)
	Py = Py + masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * item.vy * (1.0e-10/1.0e-15) * (aut/aul)
	Pz = Pz + masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * item.vz * (1.0e-10/1.0e-15) * (aut/aul)

	average_p = average_p + VectorModule([Px*Px, Py*Py, Pz*Pz]) * (aum*aul/aut)

	CMx = CMx + masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * item.vx * (1.0e-10/1.0e-15) * (aut/aul)
	CMy = CMy + masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * item.vy * (1.0e-10/1.0e-15) * (aut/aul)
	CMz = CMz + masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * item.vz * (1.0e-10/1.0e-15) * (aut/aul)

	Mtot = Mtot + masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg)

P = Px + Py + Pz

print "Total Linear Momentum =", P, "[a.u.] --- before shoot move"
print "                      =", P * aum * aul / aut, "[Kg.m.sec-1]"
print "Px, Py, Pz            =", Px, Py, Pz, "[a.u.]"
print "Px, Py, Pz            =", Px * aum * aul / aut, Py * aum * aul / aut, Pz * aum * aul / aut, "[Kg.m.sec-1]"
CM = (CMx + CMy + CMz) / Mtot
print "Velocity of Center-of-Mass Vcm =", CM * aut/aul,"[m.s-1]"
print "                               =", CM * aut/aul * 1.0e+10/1.0e15,"[Ang.fs-1]"

print "<v> =", average_v/natoms, "[Ang.fs-1]"
print "<p> =", average_p/natoms, "[Kg.m.sec-1]"

print "\treference atom (vx,vy,vz):",list_atoms[0].vx,list_atoms[0].vy,list_atoms[0].vz
print "                                   -----------------------------------------------------"

########################## SHOOT MOVE ##########################
print "------------"
print "Shooting ..."
print "------------"
print " "

pair = 0

print ".... apply +dp on atom i and -dp on atom j\n"
i = 0				# atom i
j = len(list_atoms) - 1		# atom j
while pair < int(natoms/2):			# loop over all atoms to form i-j pairs
	# print "Treating Pair: ", pair+1, " ", i, "---",j
	if i == j: break # only 1 atom or the same atom was picked

	dvx = GaussianVelocity(masses[list_atoms[i].type - 1], temp, amplitude)
	dvy = GaussianVelocity(masses[list_atoms[i].type - 1], temp, amplitude)
	dvz = GaussianVelocity(masses[list_atoms[i].type - 1], temp, amplitude)

	qo1 = [list_atoms[i].x,list_atoms[i].y,list_atoms[i].z]
	qo2 = [list_atoms[j].x,list_atoms[j].y,list_atoms[j].z]

	vo1 = [list_atoms[i].vx,list_atoms[i].vy,list_atoms[i].vz]
	vo2 = [list_atoms[j].vx,list_atoms[j].vy,list_atoms[j].vz]

	vn1 = [0.0,0.0,0.0]
	vn2 = [0.0,0.0,0.0]

	vn1,vn2 = RotateTwoVectors(qo1, qo2, vo1, vo2, masses[list_atoms[i].type - 1], masses[list_atoms[j].type - 1], dvx, dvy, dvz)

	###### replace old velocities in original list ######
	list_atoms[i].vx = vn1[0]	# atom i
	list_atoms[i].vy = vn1[1]
	list_atoms[i].vz = vn1[2]
	list_atoms[j].vx = vn2[0]	# atom j
	list_atoms[j].vy = vn2[1]
	list_atoms[j].vz = vn2[2]

	i = i + 1
	j = j - 1

	pair = pair + 1
print "--------"
print "Done ..."
print "--------"
print " "
######################### KINETIC ENERGY  (AFTER SHOOT) #########################
print ".... calculate kinetic energy"
ktotn = 0.0
for item in list_atoms:
        velocity2 = item.vx**2 + item.vy**2 + item.vz**2
        ktotn = ktotn + 0.5 * masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * velocity2 * (1.0e-10/1.0e-15)**2 * (aut/aul)**2
#       a.u.                g/mol                   Kg       a.u mass    Avogadro    A/fs        (m/sec)**2             a.u. length/time
print "Kinetic energy =", ktotn * aue, "[Joule]"
print "               =", ktotn * aue * 6.24150974e18, "[eV]"
######################### TEMPERATURE (AFTER SHOOT) #########################
print ".... calculate temperature"
tempn = 2.0 * ktotn * aue / (kb1 * (3*natoms - 3))    # T = 2 * E / (kB * (3N - 3)) --- system temperature
print "                            ============="
print "Temperature (Perturbated) =", tempn, "[K] --- after shoot move"
print "                            ============="

print "\treference atom (vx,vy,vz):",list_atoms[0].vx,list_atoms[0].vy,list_atoms[0].vz

######################### RESCALE VELOCITIES  #########################
print " "
print "------------------------"
print "Rescaling Velocities ..."
print "------------------------"
for item in list_atoms:
	item.vx = ScaleVelocity(temp,tempn,item.vx)
	item.vy = ScaleVelocity(temp,tempn,item.vy)
	item.vz = ScaleVelocity(temp,tempn,item.vz)
######################### KINETIC ENERGY (AFTER SHOOT + RESCALE) #########################
print ".... calculate kinetic energy"
ktotn = 0.0
for item in list_atoms:
        velocity2 = item.vx**2 + item.vy**2 + item.vz**2
        ktotn = ktotn + 0.5 * masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * velocity2 * (1.0e-10/1.0e-15)**2 * (aut/aul)**2
#       a.u.                g/mol                   Kg       a.u mass    Avogadro    A/fs        (m/sec)**2             a.u. length/time
print "Kinetic energy =", ktotn * aue, "[Joule]"
print "               =", ktotn * aue * 6.24150974e18, "[eV]"
######################### TEMPERATURE (AFTER SHOOT + RESCALE) #########################
print ".... calculate temperature"
tempn = 2.0 * ktotn * aue / (kb1 * (3*natoms - 3))    # T = 2 * E / (kB * (3N - 3)) --- system temperature
print "                         ==========="
print "Temperature (Rescaled) =", tempn, "[K] --- after shoot move and velocity scaling"
print "                         ==========="
########################## TOTAL LINEAR MOMENTUM #########################
print ".... check for total linear momentum conservation (after shooting)"
P = 0.0 # total linear momentum
Px = 0.0
Py = 0.0
Pz = 0.0
CMx = 0.0
CMy = 0.0
CMz = 0.0
CM = 0.0
Mtot = 0.0
average_v = 0.0
average_p = 0.0
for item in list_atoms:
	average_v = average_v + VectorModule([item.vx,item.vy,item.vz]) # Ang.fs-1

        Px = Px + masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * item.vx * (1.0e-10/1.0e-15) * (aut/aul)
        Py = Py + masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * item.vy * (1.0e-10/1.0e-15) * (aut/aul)
        Pz = Pz + masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * item.vz * (1.0e-10/1.0e-15) * (aut/aul)

	average_p = average_p + VectorModule([Px*Px, Py*Py, Pz*Pz]) * (aum*aul/aut)

        CMx = CMx + masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * item.vx * (1.0e-10/1.0e-15) * (aut/aul)
        CMy = CMy + masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * item.vy * (1.0e-10/1.0e-15) * (aut/aul)
        CMz = CMz + masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg) * item.vz * (1.0e-10/1.0e-15) * (aut/aul)

        Mtot = Mtot + masses[item.type - 1] * 1.0e-3 * (1.0/aum) * (1.0/avg)

P = Px + Py + Pz

print "<v> =", average_v/natoms, "[Ang.fs-1]"
print "<p> =", average_p/natoms, "[Kg.m.sec-1]"

print "Total Linear Momentum =", P, "[a.u.] --- after shoot move"
print "                      =", P * aum * aul / aut, "[Kg.m.sec-1]"
print "Px, Py, Pz            =", Px, Py, Pz, "[a.u.]"
print "Px, Py, Pz            =", Px * aum * aul / aut, Py * aum * aul / aut, Pz * aum * aul / aut, "[Kg.m.sec-1]"
CM = (CMx + CMy + CMz) / Mtot
print "Velocity of Center-of-Mass Vcm =", CM * aut/aul,"[m.s-1]"
print "                               =", CM * aut/aul * 1.0e+10/1.0e15,"[Ang.fs-1]"

print "\treference atom (vx,vy,vz):",list_atoms[0].vx,list_atoms[0].vy,list_atoms[0].vz
print "                                   -----------------------------------------------------"
########################## NEW LAMMPS RESTART FILE #########################
print "     ----------------------"
print "---> Create file: ascii.new (ASCII restart LAMMPS file - after shoot move)"
print "     ----------------------"
ascii_new = open('ascii.new','w')
ascii_data = open('ascii.data','r')
while True:
	line = ascii_data.readline()
	if not line: break
	if "Velocities" in line:
		ascii_new.write(line)
		line = ascii_data.readline()
		print ".... write new velocities to restart file"
		for item in list_atoms:
			ascii_data.readline()
			ascii_new.write("\n%d %2.16e %2.16e %2.16e" % (item.id,item.vx,item.vy,item.vz))
	ascii_new.write(line)
ascii_data.close()
ascii_new.close()
print "     ---------------------"
print "---> Close file: ascii.new"
print "     ---------------------"
print " "

