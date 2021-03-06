import matplotlib.pyplot as plt
from matplotlib import cm
import csv
import numpy as np
import pandas as pd

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 14
plt.rcParams['axes.linewidth'] = 2
# Generate 2 colors from the 'tab10' colormap
colors = cm.get_cmap('tab10', 2)
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.xaxis.set_tick_params(which='major', size=7, width=2, direction='in')
ax1.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in')
ax1.yaxis.set_tick_params(which='major', size=7, width=2, direction='in')
ax1.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)


#Get all of the data for boris 50
boris_data_10 = pd.read_csv('test3/data/boris10.csv').values
pos_boris_10 = boris_data_10[:,:3]
vel_boris_10 = boris_data_10[:,3:6]
time_boris_10 = boris_data_10[:,6]

#Get all of the data for boris 500
boris_data_100 = pd.read_csv('test3/data/boris100.csv').values
pos_boris_100 = boris_data_100[:,:3]
vel_boris_100 = boris_data_100[:,3:6]
time_boris_100 = boris_data_100[:,6]

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


mag = lambda x: x[0]**2 + x[1]**2 + x[2]**2

def KE_perp(position, velocity):
	magMoment = []
	for z in zip(position,velocity):
		pos = z[0]
		vel = z[1]
		magMoment.append((mag(vel) / (2*np.linalg.norm(B(pos)))))
	return np.array(magMoment)

def B(pos):
	return np.array([0,0,np.sqrt(1 + 4*((pos[0] - 12))**2)])


mag_moment_c_2 = KE_perp(pos_c_2, vel_c_2)
mag_moment_c_5 = KE_perp(pos_c_5, vel_c_5)
mag_moment_boris_10 = KE_perp(pos_boris_10, vel_boris_10)
mag_moment_boris_100 = KE_perp(pos_boris_100, vel_boris_100)

plt.setp(ax1, xlabel='t')
plt.setp(ax1, ylabel='μ')
ax1.set_xlim(0,300)

ax1.plot(time_boris_10,mag_moment_boris_10, label = "Boris, ∆t = T/10",linewidth=0.9,linestyle='dashed')
ax1.plot(time_boris_100,mag_moment_boris_100, label = "Boris, ∆t = T/100",linewidth=0.9)
ax1.plot(time_c_5,mag_moment_c_5, label = "Collocation, ∆t = T/5",linewidth=1.0)
ax1.plot(time_c_2,mag_moment_c_2, label = "Collocation, ∆t = T/2",linewidth=1.0)

ax1.legend(loc = "upper right", frameon=False, fontsize=12)
plt.savefig('test3/mag_moment.eps', format='eps')
