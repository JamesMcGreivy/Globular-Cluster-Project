import matplotlib.pyplot as plt
from matplotlib import cm
import csv
import numpy as np

plt.gcf().subplots_adjust(left=0.15)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 2
# Generate 2 colors from the 'tab10' colormap
colors = cm.get_cmap('tab10', 2)
ax = plt.subplot(111)
ax.xaxis.set_tick_params(which='major', size=7, width=2, direction='in')
ax.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in')
ax.yaxis.set_tick_params(which='major', size=7, width=2, direction='in')
ax.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in')


pos_data_everhart = 'Everhart 0.63/pos_data.csv'
vel_data_everhart = 'Everhart 0.63/vel_data.csv'

pos_everhart = []
vel_everhart = []
time_everhart = []
with open(pos_data_everhart) as csvfile1:
	read_position = csv.reader(csvfile1, delimiter = ',')

	for row in read_position:

		time_everhart.append( float(row[-1]) )

		pos = []

		for i in range((len(row)-1) // 3):

			pos.append( [ float(row[i]), float(row[i+1]), float(row[i+2]) ] )

		pos_everhart.append(pos)

with open(vel_data_everhart) as csvfile1:
	read_velocity = csv.reader(csvfile1, delimiter = ',')

	for row in read_velocity:

		vel = []

		for i in range((len(row)-1) // 3):

	  		vel.append([float(row[i]),float(row[i+1]),float(row[i+2])])

		vel_everhart.append(vel)

pos_data_boris = 'Boris 0.001/pos_data.csv'
vel_data_boris = 'Boris 0.001/vel_data.csv'

pos_boris = []
vel_boris = []
time_boris = []
with open(pos_data_boris) as csvfile1:
	read_position = csv.reader(csvfile1, delimiter = ',')

	for row in read_position:

		time_boris.append( float(row[-1]) )

		pos = []

		for i in range((len(row)-1) // 3):

			pos.append( [ float(row[i]), float(row[i+1]), float(row[i+2]) ] )

		pos_boris.append(pos)

with open(vel_data_boris) as csvfile1:
	read_velocity = csv.reader(csvfile1, delimiter = ',')

	for row in read_velocity:

		vel = []

		for i in range((len(row)-1) // 3):

	  		vel.append([float(row[i]),float(row[i+1]),float(row[i+2])])

		vel_boris.append(vel)


radius_everhart = [pos_everhart[i][0][0]**2 + pos_everhart[i][0][1]**2 for i in range(len(pos_everhart))]
del_radius_everhart = [abs(1-r) for r in radius_everhart]
time_everhart = [(t/(2*np.pi)) for t in time_everhart]

radius_boris = [pos_boris[i][0][0]**2 + pos_boris[i][0][1]**2 for i in range(len(pos_boris))]
del_radius_boris = [abs(1-r) for r in radius_boris]
time_boris = [(t/(2*np.pi)) for t in time_boris]



plt.xlabel("# of Orbits")
plt.ylabel("Gyroradius Error")

ax.plot([time_boris[i] for i in range(len(time_boris)) if i % 41 == 0],[del_radius_boris[i] for i in range(len(time_boris)) if i % 41 == 0],label = "Boris, Ωc∆t = 0.001",linewidth=0.85)
ax.plot([time_everhart[i] for i in range(len(time_everhart))],[del_radius_everhart[i] for i in range(len(time_everhart))],label = "Everhart, Ωc∆t = 0.63",linewidth=1)

ax.legend(loc = "upper left", frameon=False, fontsize=12)
plt.savefig('radius_errors.eps', format='eps')


