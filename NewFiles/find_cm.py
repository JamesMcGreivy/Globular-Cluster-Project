import numpy as np


x = np.zeros((8,8))

h = [ 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626];

def find_coeff(j,k):

	if (j < 0 or k < 0 or j > 7 or k > 7):
		print("Error")
		return
	
	elif (j < k):
		x[k][j] = 0

	elif (j == k):
		x[k][j] = 1

	elif (k == 1 and j > 1):
		x[k][j] = -h[j] * x[1,j-1]

	elif (k < j):
		x[k][j] = x[k-1,j-1] - h[j] * x[k, j-1]

	else:
		print("Error")
		return

for k in range(0,8):
	for j in range(0,8):
		find_coeff(j,k)



for row in x:
	string = ""
	for elem in row:
		string = string + (str(elem)+"00000000000000000")[:15] + "  "
	print(string)

