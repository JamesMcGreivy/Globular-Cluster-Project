import csv;

x = []

bound = 150

with open('histogram.csv') as csvfile:
  readCSV = csv.reader(csvfile, delimiter = ',')

  for row in readCSV:
    x.append(row)
    print(row)


from mpl_toolkits.mplot3d import axes3d
import numpy as np 
from matplotlib import animation, cm, pyplot as plt



fig = plt.figure()
ax = plt.subplot(111)
ax.hist([float(y[4]) for y in x])



plt.show()