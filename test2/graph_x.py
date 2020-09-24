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

#Get Everhart Data
data_everhart = 'Everhart 0.63/pos_data.csv'

pos_everhart = []
time_everhart = []
with open(data_everhart) as csvfile:
	read_position = csv.reader(csvfile, delimiter = ',')

	for row in read_position:

		time_everhart.append( float(row[3]) )
		pos = [ float(row[0]), float(row[1]), float(row[2]) ]
		pos_everhart.append(pos)

#Get Boris_001 Data
data_boris_001 = 'Boris 0.001/pos_data.csv'

pos_boris_001 = []
time_boris_001 = []
with open(data_boris_001) as csvfile:
	read_position = csv.reader(csvfile, delimiter = ',')

	for row in read_position:

		time_boris_001.append( float(row[3]) )
		pos = [ float(row[0]), float(row[1]), float(row[2]) ]
		pos_boris_001.append(pos)

#Get Boris_01 Data
data_boris_1 = 'Boris 0.1/pos_data.csv'

pos_boris_1 = []
time_boris_1 = []
with open(data_boris_1) as csvfile:
	read_position = csv.reader(csvfile, delimiter = ',')

	for row in read_position:

		time_boris_1.append( float(row[3]) )
		pos = [ float(row[0]), float(row[1]), float(row[2]) ]
		pos_boris_1.append(pos)




x_everhart =[]
y_everhart = []
for i in range(len(pos_everhart)):
	x_everhart.append(pos_everhart[i][0])
	y_everhart.append(pos_everhart[i][1])

x_boris_1 =[]
y_boris_1 = []
for i in range(len(pos_boris_1)):
	x_boris_1.append(pos_boris_1[i][0])
	y_boris_1.append(pos_boris_1[i][1])

x_boris_001 =[]
y_boris_001 = []
for i in range(len(pos_boris_001)):
	x_boris_001.append(pos_boris_001[i][0])
	y_boris_001.append(pos_boris_001[i][1])




ax1.plot(time_everhart,x_everhart, label = "Everhart, Ωc∆t = 0.63",linewidth=0.75,color="mediumblue")
ax1.plot(time_boris_1,x_boris_1, label = "Boris, Ωc∆t = 0.1",linewidth=0.75,color="orange")
ax1.plot(time_boris_001,x_boris_001, label = "Boris, Ωc∆t = 0.001",linewidth=0.75,color="indianred",linestyle="dashed")


ax1.set_xlim(450,600)
ax1.set_ylim(-1.5,1.5)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

axins = inset_axes(ax1, 1.5, 1.5, loc=2)
x1, x2 = 587, 600
y1, y2 = -0.9, 0.3
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

axins.plot(time_everhart,x_everhart, label = "Everhart, Ωc∆t = 0.63",linewidth=0.75,color="mediumblue")
axins.plot(time_boris_1,x_boris_1, label = "Boris, Ωc∆t = 0.1",linewidth=0.75,color="orange")
axins.plot(time_boris_001,x_boris_001, label = "Boris, Ωc∆t = 0.001",linewidth=0.75,color="indianred",linestyle="dashed")

plt.yticks(visible=False)
plt.xticks(visible=False)

mark_inset(ax1, axins, loc1=1, loc2=4)


plt.setp(ax1, xlabel='t')
plt.setp(ax1, ylabel='x')

ax1.legend(loc = "upper right", frameon=True, fontsize=12)
plt.savefig('x_trajectory.eps', format='eps')

plt.show()