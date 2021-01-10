import matplotlib.pyplot as plt
from matplotlib import cm
import csv
import numpy as np
import pandas as pd

plt.gcf().subplots_adjust(left=0.15)
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

#Get all of the data for boris 0.01
boris_data_01 = pd.read_csv('test1/data/boris100.csv').values[::10]
pos_boris_01 = boris_data_01[:,:3]
vel_boris_01 = boris_data_01[:,3:6]
time_boris_01 = boris_data_01[:,6]

#Get all of the data for boris 0.001
boris_data_001 = pd.read_csv('test1/data/boris1000.csv').values[::100]
pos_boris_001 = boris_data_001[:,:3]
vel_boris_001 = boris_data_001[:,3:6]
time_boris_001 = boris_data_001[:,6]

#Get all of the data for collocation 0.5
c_data_2 = pd.read_csv('test1/data/collocation2.csv').values
pos_c_2 = c_data_2[:,:3]
vel_c_2 = c_data_2[:,3:6]
time_c_2 = c_data_2[:,6]

def get_phase_error(pos, vel, time): 
	
	def arctan2(y,x):
		angle = np.arctan2(y,x)
		if (angle < 0):
			angle += 2*np.pi
		return angle

	phase_error = []
	actual_phase, actual_angle = 0,0
	predicted_phase, predicted_angle = 0,0 
	x = np.cos(time)
	y = -np.sin(time)

	for i in range(len(time)):

		if arctan2(pos[i][1],pos[i][0]) < (actual_angle - (actual_phase * 2*np.pi)):
			actual_phase +=1
		actual_angle = (actual_phase * 2*np.pi) + arctan2(pos[i][1],pos[i][0])
	
		if arctan2(y[i],x[i]) < (predicted_angle - (predicted_phase*2*np.pi)):
			predicted_phase +=1
		predicted_angle = (predicted_phase * 2*np.pi) + arctan2(y[i],x[i])

		phase_error.append(abs(actual_angle-predicted_angle))
	
	return phase_error

plt.xlabel("# of Orbits")
plt.ylabel("Phase Error (radians)")
orbits_c = [(t/(2*np.pi)) for t in time_c_2]
orbits_boris_01 = [(t/(2*np.pi)) for t in time_boris_01]
orbits_boris_001 = [(t/(2*np.pi)) for t in time_boris_001]


ax.plot(orbits_c, get_phase_error(pos_c_2, vel_c_2, time_c_2), label = "Mixed-Basis Collocation, ∆t = T/2",linewidth=1)
ax.plot(orbits_boris_001, get_phase_error(pos_boris_001, vel_boris_001, time_boris_001), label = "Boris, ∆t = T/1000",linewidth=1)
ax.plot(orbits_boris_01, get_phase_error(pos_boris_01, vel_boris_01, time_boris_01), label = "Boris, ∆t = T/100",linewidth=1)

plt.yscale("log")
plt.ylim([1e-15,100000000])

ax.legend(loc = "upper left", frameon=False, fontsize=12)
plt.tight_layout()
plt.savefig('test1/phase_errors.eps', format='eps')
