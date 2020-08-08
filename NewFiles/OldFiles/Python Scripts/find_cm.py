import numpy as np


x = np.zeros((7,7))

h = [ 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626];

def find_coeff(j,k):

	if (j < 0 or k < 0 or j > 6 or k > 6):
		print("Error")
		return

	elif (j == k):
		x[k][j] = 1

	elif (k == 0 and j > 0):
		x[k][j] = -h[j] * x[0,j-1]

	elif (k < j):
		x[k][j] = x[k-1,j-1] - h[j] * x[k, j-1]

	else:
		print("Error")
		return

for k in range(0,7):
	for j in range(k,7):
		find_coeff(j,k)



for row in x:
	string = ""
	for elem in row:
		string = string + (str(elem)+"000000000000000000000000000000000000000000000000000")[:22] + "  "
	print(string)

for k in range(0,7):
	for j in range(k,7):
		if k == j:
			assert(x[k][j] == 1)

		elif k == 0 and j > 0:
			assert(x[k][j] == -h[j]*x[k][j-1])

		elif k < j:
			assert(x[k][j] == x[k-1][j-1] - h[j]*x[k][j-1] )

