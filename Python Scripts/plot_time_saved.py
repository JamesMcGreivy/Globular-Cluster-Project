import matplotlib.pyplot as plt
import csv
import numpy as np
from statistics import mean


E = [0.1,0.01,0.001,0.0001,0.00001,0.000001]
speed = [10,17,20,24,30,40]

diff_in_e_to_radius = [e/2 for e in E]

plt.plot(diff_in_e_to_radius,speed)

plt.title("% Increase in Computational Speed using Predicted B-Values\nOver 200 Timesteps (dt = T/10)")
plt.ylabel("% increase in computation speed")
plt.xscale('log')
plt.xlabel("scale length difference between change in E and gyro-radius")
#time = [(t/(2*np.pi)) for t in time]
plt.show()



