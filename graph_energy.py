import csv
from matplotlib import pyplot as plt

energy, time = [],[]


with open('NewFiles/Output/energy.csv') as csvfile:
  read_energy = csv.reader(csvfile, delimiter = ",")
  
  for row in read_energy:
    energy.append(float(row[7]))
    time.append(float(row[9]))

fig = plt.figure()
ax1 = fig.add_subplot(111)

line1, = ax1.plot(time,energy,color="blue")

ax1.set_title("Energy Conservation")
ax1.set_xlabel("Time")
ax1.set_ylabel("Relative Change in Total Energy of System")
plt.show()


