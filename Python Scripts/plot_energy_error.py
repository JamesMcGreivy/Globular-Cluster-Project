import matplotlib.pyplot as plt
import csv
import numpy as np
from statistics import mean


E = [0.1,0.01,0.001,0.0001,0.00001,0.000001]
error_predicted_b = [1.32e-14,1.44e-14,1.495e-14,1.26e-14,1.518e-14,1.52e-14]
error_normally =[1.288e-14,1.49e-14,1.490e-14,1.31e-14,1.524e-14,1.5e-14]

diff_in_e_to_radius = [e/2 for e in E]

predict, =plt.plot(diff_in_e_to_radius,error_predicted_b,color="red",label="Using predicted b-values")
normal, =plt.plot(diff_in_e_to_radius,error_normally,color="blue",label="Using Previous b-values")


plt.legend()
plt.title("Relative Error in Energy Using Predicted B-Values\nOver 200 Timesteps (dt = T/10)")
plt.ylabel("Maximum Relative Energy Error")
plt.xscale('log')
plt.xlabel("scale length difference between change in E and gyro-radius")
#time = [(t/(2*np.pi)) for t in time]
plt.show()



