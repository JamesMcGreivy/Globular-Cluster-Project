import matplotlib.pyplot as plt
import csv
import numpy as np
from statistics import mean


convergence_with_b = []
convergence_without_b =[]
variation_vs_gyroradius_with_b = []
variation_vs_gyroradius_without_b = []

with open('NewFiles/Output/convergence.csv') as csvfile1:
	read_data = csv.reader(csvfile1, delimiter = ',')

	i = 0
	for row in read_data:
		if (float(row[3]) == 0):
			convergence_without_b.append( float(row[0]) )
			variation_vs_gyroradius_without_b.append( float(row[1]))
		else:
			convergence_with_b.append( float(row[0]) )
			variation_vs_gyroradius_with_b.append( float(row[1]))

v = [100*(v[0]-v[1])/(v[0]) for v in zip(convergence_without_b,convergence_with_b)]


plt.plot(variation_vs_gyroradius_without_b,v)

plt.xlabel("Scale difference between variation in E/B and gyroradius")
plt.ylabel("% less convergences required")

plt.title("Number of Convergences Required for non-constant E and B field\nUsing Predicted b-values")
plt.xscale('log')
plt.show()