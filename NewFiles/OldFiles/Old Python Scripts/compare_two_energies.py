import csv
from matplotlib import pyplot as plt

energy_1, time_1, energy_2, time_2 = [],[], [], []

with open('Data/Output/Energy1.csv') as csvfile:
  read_energy_1 = csv.reader(csvfile, delimiter = ",")
  
  for row in read_energy_1:
    energy_1.append(float(row[7])/float(row[5]))
    time_1.append(float(row[9]))

with open('Data/Output/Energy2.csv') as csvfile:
  read_energy_2 = csv.reader(csvfile, delimiter = ",")
  
  for row in read_energy_2:
    energy_2.append(float(row[7])/float(row[5]))
    time_2.append(float(row[9]))

fig = plt.figure()
ax1 = fig.add_subplot(111)

line1, = ax1.plot(time_1,energy_1,color="blue")
line2, = ax1.plot(time_2,energy_2,color="red")
	
line1.set_label("Method 1")
line2.set_label("Method 2")

ax1.legend()

ax1.set_title("Algorithm Energy Conservation")
ax1.set_xlabel("Time")
ax1.set_ylabel("Relative Change in Total Energy of System")
plt.show()


