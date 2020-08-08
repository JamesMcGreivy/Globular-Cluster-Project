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
time = []
with open(pos_data) as csvfile1:
	read_position = csv.reader(csvfile1, delimiter = ',')

	for row in read_position:

		time.append( float(row[-1]) )

		pos = []

		for i in range((len(row)-1) // 3):

			pos.append( [ float(row[i]), float(row[i+1]), float(row[i+2]) ] )

		positions.append(pos)

radius = [positions[i][0][0]**2 + positions[i][0][1]**2 for i in range(len(positions))]
print(radius[:10])
del_radius = [1-r for r in radius]

plt.xlabel("orbits")
plt.ylabel("Error in gyroradius")
plt.title("Gyroradius Error for constant B-field and no E-field\nfor dt = T/8 over " + str(int(np.ceil(max(time)/(2*np.pi)))) + " orbits")
time = [(t/(2*np.pi)) for t in time]

plt.plot([time[i] for i in range(len(time)) if i % 6 == 0],[del_radius[i] for i in range(len(time)) if i % 6 == 0])
plt.show()