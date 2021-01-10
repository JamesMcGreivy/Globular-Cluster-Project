import matplotlib.pyplot as plt
from matplotlib import cm
import csv
import numpy as np
import pandas as pd

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 20
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

ax1.plot(pos_boris_100[:,0],pos_boris_100[:,1], label = "Boris, ∆t = T/100",linewidth=1.0,color="darkgreen",linestyle='dashed')
ax1.plot(pos_boris_10[:,0],pos_boris_10[:,1], label = "Boris, ∆t = T/10",linewidth=1.0,color="orange")
ax1.plot(pos_c_5[:,0],pos_c_5[:,1], label = "Collocation, ∆t = T/5",linewidth=1.0,color="lightblue")
ax1.plot(pos_c_2[:,0],pos_c_2[:,1], label = "Collocation, ∆t = T/2",linewidth=1.0,color="indianred")


ax1.set_ylim(-1, 4)
ax1.set_xlim(0, 25)

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

axins = inset_axes(ax1, 1.5, 1.5, loc=2)
x1, x2 = 15.2, 15.8
y1, y2 = 0.4, 1.1
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

axins.plot(pos_boris_100[:,0],pos_boris_100[:,1], label = "Boris, ∆t = T/100",linewidth=1.0,color="darkgreen",linestyle='dashed')
axins.plot(pos_boris_10[:,0],pos_boris_10[:,1], label = "Boris, ∆t = T/10",linewidth=1.0,color="orange")
axins.plot(pos_c_5[:,0],pos_c_5[:,1], label = "Collocation, ∆t = T/5",linewidth=1.0,color="lightblue")
axins.plot(pos_c_2[:,0],pos_c_2[:,1], label = "Collocation, ∆t = T/2",linewidth=1.0,color="indianred")

plt.yticks(visible=False)
plt.xticks(visible=False)

mark_inset(ax1, axins, loc1=1, loc2=4)


plt.setp(ax1, xlabel='x')
plt.setp(ax1, ylabel='y')

ax1.legend(loc = "upper right", frameon=True, fontsize=12)
plt.tight_layout()
plt.savefig('test3/trajectory.eps', format='eps')