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


actual_phase = []

for i in range(len(time)):

	angle = np.arctan2(positions[i][0][1],positions[i][0][0])
	if (angle < 0):
		angle = 2*np.pi + angle

	actual_phase.append(angle)

predicted_phase = []
x = np.cos(time)
y = -np.sin(time)

for i in range(len(time)):
	angle = np.arctan2(y[i],x[i])

	if (angle < 0):
		angle = 2*np.pi + angle

	predicted_phase.append(angle)

plt.xlabel("orbits")
plt.ylabel("Error in phase")
plt.title("Phase Error for constant B-field and no E-field\nfor dt = T/8 over " + str(int(np.ceil(max(time)/(2*np.pi)))) + " orbits")
time = [(t/(2*np.pi)) for t in time]

plt.plot(time,[predicted_phase[i]-actual_phase[i] for i in range(len(actual_phase))])
plt.show()