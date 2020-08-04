import csv;

x = []

bound = 150

with open('Data/Output/timestep_eval_data.csv') as csvfile:
  readCSV = csv.reader(csvfile, delimiter = ',')

  for row in readCSV:
    x.append(row)
   # print(row)


from mpl_toolkits.mplot3d import axes3d
import numpy as np 
from matplotlib import animation, cm, pyplot as plt



fig = plt.figure()
ax = plt.subplot(111)
print(x[1][4])
ax.hist([float(y[4]) for y in x])
ax.set_xlim(0,1)
ax.set_title("Size of Timesteps Taken During Simulation")
ax.set_xlabel("Size of Timestep")
ax.set_ylabel("Number Taken")
ax.set_xticks([0.0625,0.125,0.25,0.5,1])
plt.show()