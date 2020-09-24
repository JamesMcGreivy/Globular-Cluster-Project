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


pos_data_boris_2_1 = 'Boris 2 0.1/pos_data.csv'
vel_data_boris_2_1 = 'Boris 2 0.1/vel_data.csv'

pos_boris_2_1 = []
vel_boris_2_1 = []
time_boris_2_1 = []
with open(pos_data_boris_2_1) as csvfile1:
	read_position = csv.reader(csvfile1, delimiter = ',')

	for row in read_position:

		time_boris_2_1.append( float(row[-1]) )

		pos = []

		for i in range((len(row)-1) // 3):

			pos.append( [ float(row[i]), float(row[i+1]), float(row[i+2]) ] )

		pos_boris_2_1.append(pos)

with open(vel_data_boris_2_1) as csvfile1:
	read_velocity = csv.reader(csvfile1, delimiter = ',')

	for row in read_velocity:

		vel = []

		for i in range((len(row)-1) // 3):

	  		vel.append([float(row[i]),float(row[i+1]),float(row[i+2])])

		vel_boris_2_1.append(vel)


pos_data_boris_2_001 = 'Boris 2 0.001/pos_data.csv'
vel_data_boris_2_001 = 'Boris 2 0.001/vel_data.csv'

pos_boris_2_001 = []
vel_boris_2_001 = []
time_boris_2_001 = []
with open(pos_data_boris_2_001) as csvfile1:
	read_position = csv.reader(csvfile1, delimiter = ',')

	for row in read_position:

		time_boris_2_001.append( float(row[-1]) )

		pos = []

		for i in range((len(row)-1) // 3):

			pos.append( [ float(row[i]), float(row[i+1]), float(row[i+2]) ] )

		pos_boris_2_001.append(pos)

with open(vel_data_boris_2_001) as csvfile1:
	read_velocity = csv.reader(csvfile1, delimiter = ',')

	for row in read_velocity:

		vel = []

		for i in range((len(row)-1) // 3):

	  		vel.append([float(row[i]),float(row[i+1]),float(row[i+2])])

		vel_boris_2_001.append(vel)


def arctan2(y,x):
	angle = np.arctan2(y,x)
	if (angle < 0):
		angle += 2*np.pi
	return angle


everhart_phase_error = []
actual_phase, actual_angle = 0,0
predicted_phase, predicted_angle = 0,0 
x = np.cos(time_everhart)
y = -np.sin(time_everhart)

for i in range(len(time_everhart)):

	actual_angle = arctan2(pos_everhart[i][0][1],pos_everhart[i][0][0])
	
	predicted_angle = arctan2(y[i],x[i])

	everhart_phase_error.append(abs(actual_angle-predicted_angle))


boris_phase_error_001 = []
actual_phase, actual_angle = 0,0
predicted_phase, predicted_angle = 0,0 
x = np.cos(time_boris_001)
y = -np.sin(time_boris_001)

for i in range(len(time_boris_001)):

	if arctan2(pos_boris_001[i][0][1],pos_boris_001[i][0][0]) < (actual_angle - (actual_phase * 2*np.pi)):
		actual_phase +=1
	actual_angle = (actual_phase * 2*np.pi) + arctan2(pos_boris_001[i][0][1],pos_boris_001[i][0][0])
	
	if arctan2(y[i],x[i]) < (predicted_angle - (predicted_phase*2*np.pi)):
		predicted_phase +=1
	predicted_angle = (predicted_phase * 2*np.pi) + arctan2(y[i],x[i])

	boris_phase_error_001.append(abs(actual_angle-predicted_angle))


boris_phase_error_1 = []
actual_phase, actual_angle = 0,0
predicted_phase, predicted_angle = 0,0 
x = np.cos(time_boris_1)
y = -np.sin(time_boris_1)

for i in range(len(time_boris_1)):

	if arctan2(pos_boris_1[i][0][1],pos_boris_1[i][0][0]) < (actual_angle - (actual_phase * 2*np.pi)):
		actual_phase +=1
	actual_angle = (actual_phase * 2*np.pi) + arctan2(pos_boris_1[i][0][1],pos_boris_1[i][0][0])
	
	if arctan2(y[i],x[i]) < (predicted_angle - (predicted_phase*2*np.pi)):
		predicted_phase +=1
	predicted_angle = (predicted_phase * 2*np.pi) + arctan2(y[i],x[i])

	boris_phase_error_1.append(abs(actual_angle-predicted_angle))

boris_phase_error_2_1 = []
actual_phase, actual_angle = 0,0
predicted_phase, predicted_angle = 0,0 
x = np.cos(time_boris_2_1)
y = -np.sin(time_boris_2_1)

for i in range(len(time_boris_2_1)):

	if arctan2(pos_boris_2_1[i][0][1],pos_boris_2_1[i][0][0]) < (actual_angle - (actual_phase * 2*np.pi)):
		actual_phase +=1
	actual_angle = (actual_phase * 2*np.pi) + arctan2(pos_boris_2_1[i][0][1],pos_boris_2_1[i][0][0])
	
	if arctan2(y[i],x[i]) < (predicted_angle - (predicted_phase*2*np.pi)):
		predicted_phase +=1
	predicted_angle = (predicted_phase * 2*np.pi) + arctan2(y[i],x[i])

	boris_phase_error_2_1.append(abs(actual_angle-predicted_angle))


boris_phase_error_2_001 = []
actual_phase, actual_angle = 0,0
predicted_phase, predicted_angle = 0,0 
x = np.cos(time_boris_2_001)
y = -np.sin(time_boris_2_001)

for i in range(len(time_boris_2_001)):

	if arctan2(pos_boris_2_001[i][0][1],pos_boris_2_001[i][0][0]) < (actual_angle - (actual_phase * 2*np.pi)):
		actual_phase +=1
	actual_angle = (actual_phase * 2*np.pi) + arctan2(pos_boris_2_001[i][0][1],pos_boris_2_001[i][0][0])
	
	if arctan2(y[i],x[i]) < (predicted_angle - (predicted_phase*2*np.pi)):
		predicted_phase +=1
	predicted_angle = (predicted_phase * 2*np.pi) + arctan2(y[i],x[i])

	boris_phase_error_2_001.append(abs(actual_angle-predicted_angle))




plt.xlabel("# of Orbits")
plt.ylabel("Phase Error")
time_everhart = [(t/(2*np.pi)) for t in time_everhart]
time_boris_1 = [(t/(2*np.pi)) for t in time_boris_1]
time_boris_001 = [(t/(2*np.pi)) for t in time_boris_001]


ax.plot(time_everhart,everhart_phase_error,label = "Everhart, Ωc∆t = 0.63",linewidth=1)
ax.plot(time_boris_001,boris_phase_error_001,label = "Boris, Ωc∆t = 0.001",linewidth=1)
ax.plot(time_boris_1,boris_phase_error_1,label = "Boris, Ωc∆t = 0.1",linewidth=1)
ax.plot(time_boris_2_001,boris_phase_error_2_001,label = "Boris (tan), Ωc∆t = 0.001",linewidth=1)
ax.plot(time_boris_2_1,boris_phase_error_2_1,label = "Boris (tan), Ωc∆t = 0.1",linewidth=1)



plt.yscale("log")
plt.ylim([1e-15,100000000])

print(max(boris_phase_error_1))

ax.legend(loc = "upper left", frameon=False, fontsize=12)


plt.savefig('phase_errors_test.eps', format='eps')
