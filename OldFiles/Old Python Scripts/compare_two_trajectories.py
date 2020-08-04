
import csv;

x1 = []
x2 = []

with open('Data/Output/Data1.csv') as csvfile1:
  read_position = csv.reader(csvfile1, delimiter = ',')

  for row in read_position:
    x1.append(row)

with open('Data/Output/Data2.csv') as csvfile2:
  read_position = csv.reader(csvfile2, delimiter = ',')

  for row in read_position:
    x2.append(row)


from mpl_toolkits.mplot3d import axes3d
import numpy as np 
from matplotlib import animation, cm, pyplot as plt



fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')

bound = 500
i=0
m1 = 20
m2 = 20
def animate(i):
  global bound
  ax1.clear()

  i1 = m1*i
  i2 = m2*i

  for k in range((len(x1[i]) - 2)//4):
    points1 = ax1.scatter(float(x1[i1][4*k]),float(x1[i1][(4*k)+1]),float(x1[i1][(4*k)+2]),color="orange", s = ((float(x1[0][(4*k)+3]))/4))
    points2 = ax1.scatter(float(x2[i2][4*k]),float(x2[i2][(4*k)+1]),float(x2[i2][(4*k)+2]),color="red", s = ((float(x2[0][(4*k)+3]))/4))

  
  ax1.set_xlim(-bound,bound)
  ax1.set_ylim(-bound,bound)
  ax1.set_zlim(-(bound/2),(bound/2))
  ax1.set_xlabel("Time : " + str(x1[i][-2]))

  print("Rendered: " , i1)

  return points1

frames = len(x1) // (m1)
print(len(x1),m1)
print(len(x2),m2)
ani = animation.FuncAnimation(fig, animate, interval=100, blit=False, frames = frames)
ani.save('Data/Animations/orbit.gif', fps=25)

