import matplotlib.pyplot as plt
from matplotlib import cm
import csv
import numpy as np

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



class Data:

	def __init__(self,positions,velocities,times):
		
		assert(len(positions) == len(velocities) == len(times))
		
		self.times = times
		self.positions  = positions
		self.velocities = velocities


#Get all of the data for boris 0.001
pos_boris = 'boris 0.001/pos_data.csv'
vel_boris = 'boris 0.001/vel_data.csv'

positions = []
velocities = []
times = []

with open(pos_boris) as csvfile:
	
	read_position = csv.reader(csvfile, delimiter = ',')
	for row in read_position:

		times.append( float(row[3]) )
		pos = [ float(row[0]), float(row[1]), float(row[2]) ]
		positions.append(pos)

with open(vel_boris) as csvfile:
	
	read_velocity = csv.reader(csvfile, delimiter = ',')
	for row in read_velocity:
		vel = [ float(row[0]),float(row[1]),float(row[2]) ]
		velocities.append(vel)

boris_data_001 = Data(positions,velocities,times)

#Get all of the data for boris 0.001
pos_boris = 'boris 0.1/pos_data.csv'
vel_boris = 'boris 0.1/vel_data.csv'

positions = []
velocities = []
times = []

with open(pos_boris) as csvfile:
	
	read_position = csv.reader(csvfile, delimiter = ',')
	for row in read_position:

		times.append( float(row[3]) )
		pos = [ float(row[0]), float(row[1]), float(row[2]) ]
		positions.append(pos)

with open(vel_boris) as csvfile:
	
	read_velocity = csv.reader(csvfile, delimiter = ',')
	for row in read_velocity:
		vel = [ float(row[0]),float(row[1]),float(row[2]) ]
		velocities.append(vel)

boris_data_1 = Data(positions,velocities,times)


#Get all of the data for Everhart 0.63
pos_everhart = 'Everhart 0.63/pos_data.csv'
vel_everhart = 'Everhart 0.63/vel_data.csv'

positions = []
velocities = []
times = []

with open(pos_everhart) as csvfile:
	
	read_position = csv.reader(csvfile, delimiter = ',')
	for row in read_position:

		times.append( float(row[3]) )
		pos = [ float(row[0]), float(row[1]), float(row[2]) ]
		positions.append(pos)

with open(vel_everhart) as csvfile:
	
	read_velocity = csv.reader(csvfile, delimiter = ',')
	for row in read_velocity:
		vel = [ float(row[0]),float(row[1]),float(row[2]) ]
		velocities.append(vel)

everhart_data = Data(positions,velocities,times)

# Functions for determining change in energy

def total_energy(position,velocity):

	vel_sq = velocity[0]**2 + velocity[1]**2 + velocity[2]**2
	return (0.5*(vel_sq)) - position[1]

def get_energy_error(data):

	data.E0 = total_energy(data.positions[0],data.velocities[0])
	delEnergy = []

	for i in range( len( data.times) ):

		energy = total_energy(data.positions[i],data.velocities[i])
		delEnergy.append(abs(energy - data.E0))

	return delEnergy

def get_normalized_energy_error(data):
	
	delEnergy = get_energy_error(data)
	delEnergy = [E/2.5 for E in delEnergy]
	return delEnergy



plt.xlabel("t")
plt.ylabel("Change in Energy")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.plot(boris_data_001.times, get_energy_error(boris_data_001), label = "Boris, Ωc∆t = 0.001",color="indianred",linewidth=0.75)
ax.plot(boris_data_1.times, get_energy_error(boris_data_1), label = "Boris, Ωc∆t = 0.1",color="orange",linewidth=0.75)
ax.plot(everhart_data.times, get_energy_error(everhart_data), label = "Everhart, Ωc∆t = 0.63",color="mediumblue",linewidth=0.75)



ax.legend(loc = "upper left", frameon=False, fontsize=12)

plt.savefig('DeltaE.eps', format='eps')




