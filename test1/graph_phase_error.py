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


actual_phase_everhart = []
for i in range(len(time_everhart)):

	angle = np.arctan2(pos_everhart[i][0][1],pos_everhart[i][0][0])
	if (angle < 0):
		angle = 2*np.pi + angle

	actual_phase_everhart.append(angle)

predicted_phase_everhart = []
x = np.cos(time_everhart)
y = -np.sin(time_everhart)

for i in range(len(time_everhart)):
	angle = np.arctan2(y[i],x[i])

	if (angle < 0):
		angle = 2*np.pi + angle

	predicted_phase_everhart.append(angle)



actual_phase_boris = []
for i in range(len(time_boris)):

	angle = np.arctan2(pos_boris[i][0][1],pos_boris[i][0][0])
	if (angle < 0):
		angle = 2*np.pi + angle

	actual_phase_boris.append(angle)

predicted_phase_boris = []
x = np.cos(time_boris)
y = -np.sin(time_boris)

for i in range(len(time_boris)):
	angle = np.arctan2(y[i],x[i])

	if (angle < 0):
		angle = 2*np.pi + angle

	predicted_phase_boris.append(angle)


plt.xlabel("# of Orbits")
plt.ylim([10e-15,1])
plt.ylabel("Phase Error")
time_everhart = [(t/(2*np.pi)) for t in time_everhart]
time_boris = [(t/(2*np.pi)) for t in time_boris]

everhart_phase_error = [abs(predicted_phase_everhart[i]-actual_phase_everhart[i]) for i in range(len(actual_phase_everhart))]
boris_phase_error = [abs(predicted_phase_boris[i]-actual_phase_boris[i]) for i in range(len(actual_phase_boris))]

boris_phase_error = [p if p < (2*np.pi-0.01) else abs(p-(2*np.pi)) for p in boris_phase_error]
ax.plot(time_everhart,everhart_phase_error,label = "Everhart, Ωc∆t = 0.63",linewidth=1)
ax.plot(time_boris,boris_phase_error,label = "Boris, Ωc∆t = 0.01",linewidth=1)
plt.yscale("log")

ax.legend(loc = "upper left", frameon=False, fontsize=12)


plt.savefig('phase_errors.eps', format='eps')
