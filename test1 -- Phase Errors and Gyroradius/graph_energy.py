import matplotlib.pyplot as plt
import csv
import numpy as np
import ast
from statistics import mean
import configparser


config = configparser.RawConfigParser()
config.read('config.txt')
config_dict = dict(config.items('variables'))

pos_data = 'pos_'+config_dict['file']+'.csv'
vel_data = 'vel_'+config_dict['file']+'.csv'

positions = []
velocities = []
time = []
with open(pos_data) as csvfile1:
	read_position = csv.reader(csvfile1, delimiter = ',')

	for row in read_position:

		time.append( float(row[-1]) )

		pos = []

		for i in range((len(row)-1) // 3):

			pos.append( [ float(row[i]), float(row[i+1]), float(row[i+2]) ] )

		positions.append(pos)

with open(vel_data) as csvfile1:
	read_velocity = csv.reader(csvfile1, delimiter = ',')

	for row in read_velocity:

		vel = []

		for i in range((len(row)-1) // 3):

	  		vel.append([float(row[i]),float(row[i+1]),float(row[i+2])])

		velocities.append(vel)


def total_kinetic_energy(velocities):

	total_KE = 0

	for v in velocities:

		v_sq = v[0]**2 + v[1]**2 + v[2]**2

		total_KE += (0.5)*v_sq

	return total_KE


def total_potential_energy(positions):

	total_PE = 0

	for pos in positions:

		r_sq = pos[0]**2 + pos[1]**2 + pos[2]**2

		total_PE += 0.5 * r_sq

	return total_PE


def total_energy(positions,velocities):

	position = positions[0]
	velocity = velocities[0]

	vel_sq = velocity[0]**2 + velocity[1]**2 + velocity[2]**2

	return 0.5*(vel_sq)

E_0 = 0
delH = []

for i in range(len(time)):
	#total_e = total_kinetic_energy(velocities[i]) + total_potential_energy(positions[i])
	total_e = total_energy(positions[i],velocities[i])

	if i == 0:
		E_0 = total_e

	delH.append( (total_e - E_0) )

print(max([abs(h) for h in delH]))
plt.title("Relative Change in Energy over 100,000 orbits using dt = T/8")
plt.xlabel("orbits")
plt.ylabel("Relative Change in Energy")
time = [(t/(2*np.pi)) for t in time]
plt.plot(time,delH)
plt.show()



