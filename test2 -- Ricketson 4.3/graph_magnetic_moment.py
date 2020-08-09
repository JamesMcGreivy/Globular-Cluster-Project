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

mag = lambda x: np.sqrt(np.dot(x,x))


def KE_perp(velocity, B):

	V_parr = [(np.dot(velocity, B) / (mag(B)**2))*b for b in B]
	V_perp = [velocity[i]-V_parr[i] for i in range(3)]
	return (0.5)*(mag(V_perp)**2)


def B(position):

	return [0,0,np.sqrt(1 + 4*((position[0]-12)**2) )]

mag_moment = [KE_perp(velocities[i][0],B(positions[i][0]))/mag(B(positions[i][0])) for i in range(len(time))]

plt.xlabel("t")
plt.ylabel("μ")
plt.title('Magnetic Moment (μ) Over Time for\nTest Case 4.3 with ΩcΔt = 0.78')
plt.plot(time,mag_moment)
plt.show()


