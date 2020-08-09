import sys
import csv;

pos_data = 'pos_'+sys.argv[1]+'.csv'

x = []


with open('Output/'+pos_data) as csvfile1:
  read_position = csv.reader(csvfile1, delimiter = ',')

  for row in read_position:
    x.append(row)


from mpl_toolkits.mplot3d import axes3d
import numpy as np 
from matplotlib import animation, cm, pyplot as plt



fig = plt.figure()
ax1 = fig.add_subplot(111)#, projection='3d')

bound = 10
i=0
m = 25
def animate(i):
  global bound
  ax1.clear()

  i = m*i

  for k in range((len(x[i]) - 1)//3):
    points1 = ax1.scatter(float(x[i][3*k]),float(x[i][(3*k)+1]))#,float(x[i][(3*k)+2]),color="blue",s=2)
  
  ax1.set_xlim(0,25)
  ax1.set_ylim(-2,4)
  ax1.set_ylabel("y-axis")
  #ax1.set_zlim(-bound,bound)
 #ax1.set_zlabel("z-axis")
  ax1.set_xlabel("Time : " + str(x[i][-1]))

  print("Rendered: " , i)

  return points1

print(len(x))
frames = len(x) // (m)
ani = animation.FuncAnimation(fig, animate, interval=100, blit=False, frames = frames)
ani.save('Output/orbit.gif', fps=25)

