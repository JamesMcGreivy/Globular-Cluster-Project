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

pos_boris_001 = []
vel_boris_001 = []
time_boris_001 = []
with open(pos_data_boris) as csvfile1:
	read_position = csv.reader(csvfile1, delimiter = ',')

	for row in read_position:

		time_boris_001.append( float(row[-1]) )

		pos = []

		for i in range((len(row)-1) // 3):

			pos.append( [ float(row[i]), float(row[i+1]), float(row[i+2]) ] )

		pos_boris_001.append(pos)

with open(vel_data_boris) as csvfile1:
	read_velocity = csv.reader(csvfile1, delimiter = ',')

	for row in read_velocity:

		vel = []

		for i in range((len(row)-1) // 3):

	  		vel.append([float(row[i]),float(row[i+1]),float(row[i+2])])

		vel_boris_001.append(vel)


pos_data_boris = 'Boris 0.1/pos_data.csv'
vel_data_boris = 'Boris 0.1/vel_data.csv'

pos_boris_1 = []
vel_boris_1 = []
time_boris_1 = []
with open(pos_data_boris) as csvfile1:
	read_position = csv.reader(csvfile1, delimiter = ',')

	for row in read_position:

		time_boris_1.append( float(row[-1]) )

		pos = []

		for i in range((len(row)-1) // 3):

			pos.append( [ float(row[i]), float(row[i+1]), float(row[i+2]) ] )

		pos_boris_1.append(pos)

with open(vel_data_boris) as csvfile1:
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



def total_kinetic_energy(velocities):

	total_KE = 0

	for v in velocities:

		v_sq = v[0]**2 + v[1]**2 + v[2]**2

		total_KE += (0.5)*v_sq

	return total_KE


def total_potential_energy(positions):

	total_PE = 0

	for pos in positions:

		r_sq = pos[0]**2 + pos[1]**2 + pos[2]**2

		total_PE += 0.5 * r_sq

	return total_PE


def total_energy(positions,velocities):

	position = positions[0]
	velocity = velocities[0]

	vel_sq = velocity[0]**2 + velocity[1]**2 + velocity[2]**2

	return 0.5*(vel_sq)

E_0 = 0
delH_everhart = []

for i in range(len(time_everhart)):
	#total_e = total_kinetic_energy(velocities[i]) + total_potential_energy(positions[i])
	total_e = total_energy(pos_everhart[i],vel_everhart[i])

	if i == 0:
		E_0 = total_e

	delH_everhart.append( abs((total_e - E_0)/E_0) )
time_everhart = [(t/(2*np.pi)) for t in time_everhart]

E_0 = 0
delH_boris_001 = []

for i in range(len(time_boris_001)):
	#total_e = total_kinetic_energy(velocities[i]) + total_potential_energy(positions[i])
	total_e = total_energy(pos_boris_001[i],vel_boris_001[i])

	if i == 0:
		E_0 = total_e

	delH_boris_001.append( abs((total_e - E_0)/E_0) )
time_boris_001 = [(t/(2*np.pi)) for t in time_boris_001]

E_0 = 0
delH_boris_1 = []

for i in range(len(time_boris_1)):
	#total_e = total_kinetic_energy(velocities[i]) + total_potential_energy(positions[i])
	total_e = total_energy(pos_boris_1[i],vel_boris_1[i])

	if i == 0:
		E_0 = total_e

	delH_boris_1.append( abs((total_e - E_0)/E_0) )
time_boris_1 = [(t/(2*np.pi)) for t in time_boris_1]


E_0 = 0
delH_boris_2_1 = []

for i in range(len(time_boris_2_1)):
	#total_e = total_kinetic_energy(velocities[i]) + total_potential_energy(positions[i])
	total_e = total_energy(pos_boris_2_1[i],vel_boris_2_1[i])

	if i == 0:
		E_0 = total_e

	delH_boris_2_1.append( abs((total_e - E_0)/E_0) )
time_boris_2_1 = [(t/(2*np.pi)) for t in time_boris_2_1]


E_0 = 0
delH_boris_2_001 = []

for i in range(len(time_boris_2_001)):
	#total_e = total_kinetic_energy(velocities[i]) + total_potential_energy(positions[i])
	total_e = total_energy(pos_boris_2_001[i],vel_boris_2_001[i])

	if i == 0:
		E_0 = total_e

	delH_boris_2_001.append( abs((total_e - E_0)/E_0) )
time_boris_2_001 = [(t/(2*np.pi)) for t in time_boris_2_001]



plt.plot(time_everhart,delH_everhart,label = "Everhart, Ωc∆t = 0.63",linewidth=0.75)
plt.plot(time_boris_001,delH_boris_001,label = "Boris, Ωc∆t = 0.001",linewidth=0.75)
plt.plot(time_boris_1,delH_boris_1,label = "Boris, Ωc∆t = 0.1",linewidth=0.75)
plt.plot(time_boris_2_001,delH_boris_2_001,label = "Boris (tan), Ωc∆t = 0.001",linewidth=0.75)
plt.plot(time_boris_2_1,delH_boris_2_1,label = "Boris (tan), Ωc∆t = 0.1",linewidth=0.75)


plt.ylim([0,8e-14])
plt.xlabel("orbits")
plt.ylabel("Relative Change in Energy")

ax.legend(loc = "upper left", frameon=False, fontsize=12)
plt.savefig('energy_errors.eps', format='eps')





