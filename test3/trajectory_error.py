import numpy as np 
import csv
from matplotlib import cm
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 20
plt.rcParams['axes.linewidth'] = 2


# Generate 2 colors from the 'tab10' colormap
colors = cm.get_cmap('tab10', 2)
ax = plt.subplot(111)
ax.xaxis.set_tick_params(which='major', size=7, width=2, direction='in')
ax.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in')
ax.yaxis.set_tick_params(which='major', size=7, width=2, direction='in')
ax.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in')


plt.xlabel("t")
plt.ylabel("Estimated Error in Trajectory")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)


#Get all of the data for boris 10
boris_data_10 = pd.read_csv('test3/data/boris10.csv').values
pos_boris_10 = boris_data_10[:,:3]
vel_boris_10 = boris_data_10[:,3:6]
time_boris_10 = boris_data_10[:,6]

#Get all of the data for boris 100
boris_data_100 = pd.read_csv('test3/data/boris100.csv').values
pos_boris_100 = boris_data_100[:,:3]
vel_boris_100 = boris_data_100[:,3:6]
time_boris_100 = boris_data_100[:,6]

#Get all of the data for boris 1000
boris_data_1000 = pd.read_csv('test3/data/boris1000.csv').values
pos_boris_1000 = boris_data_1000[:,:3]
vel_boris_1000 = boris_data_1000[:,3:6]
time_boris_1000 = boris_data_1000[:,6]

#Get all of the data for collocation 2
c_data_2 = pd.read_csv('test3/data/collocation2.csv').values
pos_c_2 = c_data_2[:,:3]
vel_c_2 = c_data_2[:,3:6]
time_c_2 = c_data_2[:,6]

#Get all of the data for collocation 5
c_data_5 = pd.read_csv('test3/data/collocation5.csv').values
pos_c_5 = c_data_5[:,:3]
vel_c_5 = c_data_5[:,3:6]
time_c_5 = c_data_5[:,6]

#Get all of the data for collocation 10
c_data_10 = pd.read_csv('test3/data/collocation10.csv').values
pos_c_10 = c_data_10[:,:3]
vel_c_10 = c_data_10[:,3:6]
time_c_10 = c_data_10[:,6]

def interp(time, pos):
	diffX = pos[:,0] - np.interp(time, time_boris_1000, pos_boris_1000[:,0])
	diffY = pos[:,1] - np.interp(time, time_boris_1000, pos_boris_1000[:,1])
	return np.sqrt(diffX**2 + diffY**2)

ax.plot(time_boris_10, interp(time_boris_10, pos_boris_10), label="Boris, T/10", linewidth=1.2,color="orange")
ax.plot(time_boris_100, interp(time_boris_100, pos_boris_100), label="Boris, T/100", linewidth=1.2)

ax.plot(time_c_10, interp(time_c_10, pos_c_10), label="Collocation, T/10", linewidth=1.2,color="black")
ax.plot(time_c_5, interp(time_c_5, pos_c_5), label="Collocation, T/5", linewidth=1.2,color="lightblue")
ax.plot(time_c_2, interp(time_c_2, pos_c_2), label="Collocation, T/2", linewidth=1.2,color="indianred")

ax.set_ylim(1e-11, 1e1)

ax.legend(loc = "best", frameon=False, fontsize=12)
ax.set_yscale("log")

plt.tight_layout()
plt.savefig('test3/trajectory_error.eps', format='eps')


