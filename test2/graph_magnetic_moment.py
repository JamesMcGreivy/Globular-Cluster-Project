import matplotlib.pyplot as plt
from matplotlib import cm
import csv
import numpy as np

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 2
# Generate 2 colors from the 'tab10' colormap
colors = cm.get_cmap('tab10', 2)
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = fig.add_subplot(221)
ax1.xaxis.set_tick_params(which='major', size=7, width=2, direction='in')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in')
ax1.yaxis.set_tick_params(which='major', size=7, width=2, direction='in')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.yaxis.set_ticks([])
ax2.xaxis.set_ticks([])

pos_data_boris= 'Boris 0.001/pos_data.csv'
vel_data_boris = 'Boris 0.001/vel_data.csv'


pos_data_ours = 'Everhart 0.63/pos_data.csv'
vel_data_ours = 'Everhart 0.63/vel_data.csv'

pos_boris = []
vel_boris = []
time_boris = []

with open(pos_data_boris) as csvfile:
	read_position = csv.reader(csvfile, delimiter = ',')

	for row in read_position:

		time_boris.append( float(row[3]) )
		pos = [ float(row[0]), float(row[1]), float(row[2]) ]
		pos_boris.append(pos)

with open(vel_data_boris) as csvfile:
	read_velocity = csv.reader(csvfile, delimiter = ',')

	for row in read_velocity:
		vel = [ float(row[0]),float(row[1]),float(row[2]) ]
		vel_boris.append(vel)


pos_ours = []
vel_ours = []
time_ours = []
with open(pos_data_ours) as csvfile:
	read_position = csv.reader(csvfile, delimiter = ',')

	for row in read_position:

		time_ours.append( float(row[3]) )
		pos = [ float(row[0]), float(row[1]), float(row[2]) ]
		pos_ours.append(pos)

with open(vel_data_ours) as csvfile:
	read_velocity = csv.reader(csvfile, delimiter = ',')

	for row in read_velocity:
		vel = [ float(row[0]),float(row[1]),float(row[2]) ]
		vel_ours.append(vel)


mag = lambda x: np.sqrt(np.dot(x,x))


def KE_perp(velocity, B):

	V_parr = [(np.dot(velocity, B) / (mag(B)**2))*b for b in B]
	V_perp = [velocity[i]-V_parr[i] for i in range(3)]
	return (0.5)*(mag(V_perp)**2)


def B(position):
	return [0,0,np.sqrt(1 + 4*((position[0]-12)**2) )]

mag_moment_ours = [KE_perp(vel_ours[i],B(pos_ours[i]))/mag(B(pos_ours[i])) for i in range(len(time_ours))]
mag_moment_boris = [KE_perp(vel_boris[i],B(pos_boris[i]))/mag(B(pos_boris[i])) for i in range(len(time_boris))]

plt.setp(ax1, xlabel='t')
plt.setp(ax1, ylabel='μ')
#mag_moment_boris = np.interp(time_ours,time_boris,mag_moment_boris)

ax1.plot(time_ours,mag_moment_ours, label = "Everhart, Ωc∆t = 0.63",linewidth=0.8)
ax1.plot(time_boris,mag_moment_boris, label = "Boris, Ωc∆t = 0.01",linewidth=0.8)

dist1 = len(time_ours)
dist2 = len(time_boris)
ax2.plot(time_ours[dist1//2 -2: 1003*dist1//2000 -1],mag_moment_ours[dist1//2 - 2: 1003*dist1//2000 - 1],linewidth=0.85)
ax2.plot(time_boris[dist2//2 : 1003*dist2//2000],mag_moment_boris[dist2//2 : 1003*dist2//2000],linewidth=0.75)

ax1.legend(loc = "upper right", frameon=False, fontsize=12)
plt.savefig('mag_moment.eps', format='eps')

