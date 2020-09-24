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

pos_data_boris_001 = 'Boris 0.001/pos_data.csv'
vel_data_boris_001 = 'Boris 0.001/vel_data.csv'

pos_boris_001 = []
vel_boris_001 = []
time_boris_001 = []
with open(pos_data_boris_001) as csvfile1:
	read_position = csv.reader(csvfile1, delimiter = ',')

	for row in read_position:

		time_boris_001.append( float(row[-1]) )

		pos = []

		for i in range((len(row)-1) // 3):

			pos.append( [ float(row[i]), float(row[i+1]), float(row[i+2]) ] )

		pos_boris_001.append(pos)

with open(vel_data_boris_001) as csvfile1:
	read_velocity = csv.reader(csvfile1, delimiter = ',')

	for row in read_velocity:

		vel = []

		for i in range((len(row)-1) // 3):

	  		vel.append([float(row[i]),float(row[i+1]),float(row[i+2])])

		vel_boris_001.append(vel)

pos_data_boris_1 = 'Boris 0.1/pos_data.csv'
vel_data_boris_1 = 'Boris 0.1/vel_data.csv'

pos_boris_1 = []
vel_boris_1 = []
time_boris_1 = []
with open(pos_data_boris_1) as csvfile1:
	read_position = csv.reader(csvfile1, delimiter = ',')

	for row in read_position:

		time_boris_1.append( float(row[-1]) )

		pos = []

		for i in range((len(row)-1) // 3):

			pos.append( [ float(row[i]), float(row[i+1]), float(row[i+2]) ] )

		pos_boris_1.append(pos)

with open(vel_data_boris_1) as csvfile1:
	read_velocity = csv.reader(csvfile1, delimiter = ',')

	for row in read_velocity:

		vel = []

		for i in range((len(row)-1) // 3):

	  		vel.append([float(row[i]),float(row[i+1]),float(row[i+2])])

		vel_boris_1.append(vel)



radius_everhart = [pos_everhart[i][0][0]**2 + pos_everhart[i][0][1]**2 for i in range(len(pos_everhart))]
del_radius_everhart = [abs(1-r) for r in radius_everhart]
time_everhart = [(t/(2*np.pi)) for t in time_everhart]

radius_boris_001 = [pos_boris_001[i][0][0]**2 + pos_boris_001[i][0][1]**2 for i in range(len(pos_boris_001))]
del_radius_boris_001 = [abs(1-r) for r in radius_boris_001]
time_boris_001 = [(t/(2*np.pi)) for t in time_boris_001]


radius_boris_1 = [pos_boris_1[i][0][0]**2 + pos_boris_1[i][0][1]**2 for i in range(len(pos_boris_1))]
del_radius_boris_1 = [abs(1-r) for r in radius_boris_1]
time_boris_1 = [(t/(2*np.pi)) for t in time_boris_1]

time_boris_001 = [time_boris_001[i] for i in range(len(time_boris_001)) if i%5000 == 0]
del_radius_boris_001 = [del_radius_boris_001[i] for i in range(len(del_radius_boris_001)) if i%5000 == 0]

time_boris_1 = [time_boris_1[i] for i in range(len(time_boris_1)) if i%50 == 0]
del_radius_boris_1 = [del_radius_boris_1[i] for i in range(len(del_radius_boris_1)) if i%50 == 0]

time_everhart = [time_everhart[i] for i in range(len(time_everhart)) if i%50 == 0]
del_radius_everhart = [del_radius_everhart[i] for i in range(len(del_radius_everhart)) if i%50 == 0]


plt.xlabel("# of Orbits")
plt.ylabel("Gyroradius Error")

ax.plot(time_everhart,del_radius_everhart,label = "Everhart, Ωc∆t = 0.63",linewidth=0.70)
ax.plot(time_boris_001,del_radius_boris_001,label = "Boris, Ωc∆t = 0.001",linewidth=0.70)
ax.plot(time_boris_1,del_radius_boris_1,label = "Boris, Ωc∆t = 0.1",linewidth=0.70)


ax.legend(loc = "upper left", frameon=False, fontsize=12)
plt.savefig('radius_errors.eps', format='eps')


