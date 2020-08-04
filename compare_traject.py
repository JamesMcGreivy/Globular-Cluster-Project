import matplotlib.pyplot as plt
import csv
import numpy as np

time = []
positions = []
with open('NewFiles/Output/positions.csv') as csvfile1:
	read_position = csv.reader(csvfile1, delimiter = ',')

	for row in read_position:

		time.append( float(row[-1]) )

		pos = []

		for i in range((len(row)-1) // 3):

			pos.append( [ float(row[i]), float(row[i+1]), float(row[i+2]) ] )

		positions.append(pos)

velocities = []
with open('NewFiles/Output/velocities.csv') as csvfile1:
	read_velocity = csv.reader(csvfile1, delimiter = ',')

	for row in read_velocity:

		vel = []

		for i in range((len(row)-1) // 3):

  			vel.append([float(row[i]),float(row[i+1]),float(row[i+2])])

		velocities.append(vel)


#c1 = 1,0,0
#c2 = 0,0,1
pos_0 = [1,1,1]
vel_0 = [1,1,1]
g = [1,0,0]
omega = [0,1,0]

def x_perp(t):
	#y direction
	return vel_0[1]*t + pos_0[1]

def x_parallel(t):

	return [np.sin(t) + pos_0[0], -np.cos(t) + t + (1 + pos_0[2])]

pos = []
for t in time:
	pos.append([x_parallel(t)[0],x_perp(t),x_parallel(t)[1]])




plt.title("Change in Position over \n100,000 orbits using dt = T/4")
plt.xlabel("Number of Orbits")
plt.ylabel("Difference in Position")
time = [(t/(2*np.pi)) for t in time]
#print(max([p[0][0] - p[1][0][0] for p in zip(pos,positions)]))
#plt.plot(time,[p[0][0] - p[1][0][0] for p in zip(pos,positions)])
plt.plot(time,[p[0][0] for p in positions])
plt.show()

