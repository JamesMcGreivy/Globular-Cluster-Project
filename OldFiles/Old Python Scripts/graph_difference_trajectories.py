
import csv
import numpy as np
x1 = []
x2 = []

x = []

with open('Data/Output/Data1.csv') as csvfile1:
  read_position = csv.reader(csvfile1, delimiter = ',')

  for row in read_position:
    x1.append(row)

with open('Data/Output/Data2.csv') as csvfile2:
  read_position = csv.reader(csvfile2, delimiter = ',')

  for row in read_position:
    x2.append(row)

for k in range((len(x1[0]) - 2)//4):
	
	y = []
	
	for i in range(len(x1)):
		distance = np.sqrt((float(x1[i][4*k]) - float(x2[i][4*k]))**2 + (float(x1[i][(4*k)+1]) - float(x2[i][(4*k)+1]))**2 + (float(x1[i][(4*k)+2]) - float(x2[i][(4*k)+2]))**2)
		y.append([x1[i][-2],distance])

	x.append(y)

from matplotlib import pyplot as plt

fig = plt.figure()
ax1 = fig.add_subplot(111)


for i in range(len(x)):
	dist = []
	time = []

	for data in x[i]:
		dist.append(data[1])
		time.append(float(data[0]))

	line, = ax1.plot(time,dist)
	line.set_label("Particle " + str(i+1))

ax1.set_xlabel("Time")
ax1.set_ylabel("Distance ")
ax1.set_title("Distance Between Particles Simulated With and Without Adaptive Time-Stepping\nUsing Leap Frogging")
ax1.legend()

plt.show()
