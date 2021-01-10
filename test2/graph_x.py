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

#Get all of the data for boris 0.01
boris_data_01 = pd.read_csv('test2/data/boris100.csv').values
pos_boris_01 = boris_data_01[:,:3]
vel_boris_01 = boris_data_01[:,3:6]
time_boris_01 = boris_data_01[:,6]

#Get all of the data for boris 0.001
boris_data_001 = pd.read_csv('test2/data/boris1000.csv').values
pos_boris_001 = boris_data_001[:,:3]
vel_boris_001 = boris_data_001[:,3:6]
time_boris_001 = boris_data_001[:,6]


#Get all of the data for collocation 0.2
c_data_5 = pd.read_csv('test2/data/collocation5.csv').values
pos_c_5 = c_data_5[:,:3]
vel_c_5 = c_data_5[:,3:6]
time_c_5 = c_data_5[:,6]

#Get all of the data for collocation 0.1
c_data_10 = pd.read_csv('test2/data/collocation10.csv').values
pos_c_10 = c_data_10[:,:3]
vel_c_10 = c_data_10[:,3:6]
time_c_10 = c_data_10[:,6]


ax1.plot(time_boris_01, pos_boris_01[:,0], label = "Boris, ∆t = T/100",linewidth=1.0,color="orange")
ax1.plot(time_boris_001, pos_boris_001[:,0], label = "Boris, ∆t = T/1000",linewidth=1.0,color="indianred")
ax1.plot(time_c_10, pos_c_10[:,0], label = "Mixed Basis Collocation, ∆t = T/10",linewidth=1.0,color="mediumblue")
ax1.plot(time_c_5, pos_c_5[:,0], label = "Mixed Basis Collocation, ∆t = T/5",linewidth=1.0)

ax1.set_xlim(5850,6000)
ax1.set_ylim(-1.0,2.5)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

axins = inset_axes(ax1, 1.30, 1.6, loc=2)
x1, x2 = 5958, 5962
y1, y2 = -0.10, 0.9
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

axins.plot(time_boris_01, pos_boris_01[:,0], label = "Boris, ∆t = T/100",linewidth=1.0,color="orange")
axins.plot(time_boris_001, pos_boris_001[:,0], label = "Boris, ∆t = T/1000",linewidth=1.0,color="indianred")
axins.plot(time_c_10, pos_c_10[:,0], label = "Mixed Basis Collocation, ∆t = T/10",linewidth=1.0,color="mediumblue")
axins.plot(time_c_5, pos_c_5[:,0], label = "Mixed Basis Collocation, ∆t = T/5",linewidth=1.0)


plt.yticks(visible=False)
plt.xticks(visible=False)

mark_inset(ax1, axins, loc1=1, loc2=4)


plt.setp(ax1, xlabel='t')
plt.setp(ax1, ylabel='x')

ax1.legend(loc = "upper right", frameon=True, fontsize=12)
plt.savefig('test2/x_trajectory.eps', format='eps')
