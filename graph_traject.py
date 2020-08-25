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

positions = []
time = []
with open('output/'+pos_data) as csvfile1:
	read_position = csv.reader(csvfile1, delimiter = ',')

	for row in read_position:

		time.append( float(row[-1]) )

		pos = []

		for i in range((len(row)-1) // 3):

			pos.append( [ float(row[i]), float(row[i+1]), float(row[i+2]) ] )

		positions.append(pos)

x =[]
y = []
for i in range(len(positions)):
	x.append(positions[i][0][0])
	y.append(positions[i][0][1])
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.set_xlabel("x position")
ax1.set_ylabel("y position")
plt.title("Trajectory for\nTest Case 4.3 with ΩcΔt = 0.78")
ax1.set_xlim(-1,1)
ax1.set_ylim(-1,1)
ax1.plot(x,y)
plt.show()