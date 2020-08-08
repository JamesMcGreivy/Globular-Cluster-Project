import csv;

position, velocity = [], []


with open('NewFiles/Output/positions.csv') as csvfile1:
  read_position = csv.reader(csvfile1, delimiter = ',')

  for row in read_position:
    position.append(row)

with open('NewFiles/Output/velocities.csv') as csvfile1:
  read_velocity = csv.reader(csvfile1, delimiter = ',')

  for row in read_velocity:
    velocity.append(row)


from mpl_toolkits.mplot3d import axes3d
import numpy as np 
from matplotlib import animation, cm, pyplot as plt



fig = plt.figure()
ax1 = fig.add_subplot(111)

data1, data2 = [], []
prev_x = 0
prev_y = 0
prev_y_vel = 0

for zipped_row in zip(position,velocity):
  x = float(zipped_row[0][0])
  y = float(zipped_row[0][1])
  y_vel = float(zipped_row[1][1])
  
  if abs(x) == 0:
    p_y = y_vel * np.sqrt( 1 / (1 - (y_vel**2)))  
    
    data1.append(y)
    data2.append(p_y)

  if prev_x < 0 and x > 0:
    y = y - (x*(y - prev_y)/(x-prev_x))
    y_vel = y_vel - (x*(y_vel - prev_y_vel)/(x-prev_x))
    p_y = y_vel * np.sqrt( 1 / (1 - (y_vel**2)))
  
    data1.append(y)
    data2.append(p_y)


  prev_x = float(zipped_row[0][0])
  prev_y = float(zipped_row[0][1])
  prev_y_vel = float(zipped_row[1][1])




plt.scatter(data1, data2, s = 0.25)
plt.xlabel("y")
plt.ylabel("p_y")

plt.show()