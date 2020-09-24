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
ax1.xaxis.set_tick_params(which='major', size=7, width=2, direction='in')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in')
ax1.yaxis.set_tick_params(which='major', size=7, width=2, direction='in')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

pos_data_boris_001= 'Boris 0.001/pos_data.csv'
vel_data_boris_001 = 'Boris 0.001/vel_data.csv'


pos_data_ours = 'Everhart 0.63/pos_data.csv'
vel_data_ours = 'Everhart 0.63/vel_data.csv'

pos_boris_001 = []
vel_boris_001 = []
time_boris_001 = []

with open(pos_data_boris_001) as csvfile:
	read_position = csv.reader(csvfile, delimiter = ',')

	for row in read_position:

		time_boris_001.append( float(row[3]) )
		pos = [ float(row[0]), float(row[1]), float(row[2]) ]
		pos_boris_001.append(pos)

with open(vel_data_boris_001) as csvfile:
	read_velocity = csv.reader(csvfile, delimiter = ',')

	for row in read_velocity:
		vel = [ float(row[0]),float(row[1]),float(row[2]) ]
		vel_boris_001.append(vel)


pos_data_boris_1= 'Boris 0.1/pos_data.csv'
vel_data_boris_1 = 'Boris 0.1/vel_data.csv'

pos_boris_1 = []
vel_boris_1 = []
time_boris_1 = []

with open(pos_data_boris_1) as csvfile:
	read_position = csv.reader(csvfile, delimiter = ',')

	for row in read_position:

		time_boris_1.append( float(row[3]) )
		pos = [ float(row[0]), float(row[1]), float(row[2]) ]
		pos_boris_1.append(pos)

with open(vel_data_boris_1) as csvfile:
	read_velocity = csv.reader(csvfile, delimiter = ',')

	for row in read_velocity:
		vel = [ float(row[0]),float(row[1]),float(row[2]) ]
		vel_boris_1.append(vel)



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
mag_moment_boris_001 = [KE_perp(vel_boris_001[i],B(pos_boris_001[i]))/mag(B(pos_boris_001[i])) for i in range(len(time_boris_001))]
mag_moment_boris_1 = [KE_perp(vel_boris_1[i],B(pos_boris_1[i]))/mag(B(pos_boris_1[i])) for i in range(len(time_boris_1))]

plt.setp(ax1, xlabel='t')
plt.setp(ax1, ylabel='μ')
#mag_moment_boris = np.interp(time_ours,time_boris,mag_moment_boris)

ax1.plot(time_ours,mag_moment_ours, label = "Everhart, Ωc∆t = 0.63",linewidth=0.8,color="mediumblue")
ax1.plot(time_boris_001,mag_moment_boris_001, label = "Boris, Ωc∆t = 0.001",linewidth=0.8,linestyle="dashed",color="indianred")
ax1.plot(time_boris_1,mag_moment_boris_1, label = "Boris, Ωc∆t = 0.1",linewidth=0.8,color="orange")


from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


axins = inset_axes(ax1, 2, 2, loc=2)
x1, x2 = 155, 160
y1, y2 = .125, 0.23
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

axins.plot(time_ours,mag_moment_ours, label = "Everhart, Ωc∆t = 0.63",linewidth=0.8,color="mediumblue")
axins.plot(time_boris_001,mag_moment_boris_001, label = "Boris, Ωc∆t = 0.001",linewidth=0.8,linestyle="dashed",color="indianred")

plt.yticks(visible=False)
plt.xticks(visible=False)

mark_inset(ax1, axins, loc1=1, loc2=4)


ax1.legend(loc = "upper right", frameon=False, fontsize=12)
plt.savefig('mag_moment.eps', format='eps')

