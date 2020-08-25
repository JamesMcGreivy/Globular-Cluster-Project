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

#Get Boris Data
data_boris = 'Boris 0.001/pos_data.csv'

pos_boris = []
time_boris = []
with open(data_boris) as csvfile:
	read_position = csv.reader(csvfile, delimiter = ',')

	for row in read_position:

		time_boris.append( float(row[3]) )
		pos = [ float(row[0]), float(row[1]), float(row[2]) ]
		pos_boris.append(pos)


x_everhart =[]
y_everhart = []
for i in range(len(pos_everhart)):
	x_everhart.append(pos_everhart[i][0])
	y_everhart.append(pos_everhart[i][1])

x_boris =[]
y_boris = []
for i in range(len(pos_boris)):
	x_boris.append(pos_boris[i][0])
	y_boris.append(pos_boris[i][1])



ax1.plot(x_boris,y_boris, label = "Boris, Ωc∆t = 0.001",linewidth=0.75)
ax1.plot(x_everhart,y_everhart, label = "Everhart, Ωc∆t = 0.63",linewidth=0.75)

ax2.plot(x_boris[len(time_boris)//2:1008*len(time_boris)//2000],y_boris[len(time_boris)//2:1008*len(time_boris)//2000],linewidth=0.75)
ax2.plot(x_everhart[len(time_everhart)//2 -2:1008*len(time_everhart)//2000 -1],y_everhart[len(time_everhart)//2 - 2:1008*len(time_everhart)//2000 -1],linewidth=0.75)

plt.setp(ax1, xlabel='x')
plt.setp(ax1, ylabel='y')

ax1.legend(loc = "upper right", frameon=True, fontsize=12)
plt.savefig('trajectory.eps', format='eps')