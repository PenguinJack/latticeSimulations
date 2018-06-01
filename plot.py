#!/usr/bin/env python
#
# Generates a contourplot by reading a csv file
# with columns x, y, f(x,y)
#
# Author: Jack Jenkins

def main():

 	import numpy as np
	import matplotlib.pyplot as plt
	import matplotlib.tri as tri

	# Load data from CSV
	data = np.genfromtxt("scan_data.csv", delimiter=',',skip_header=1)

	T = data[:, 0]
	h = data[:, 1]
	m = data[:, 2]

	triang = tri.Triangulation(T, h)
	plt.tricontourf(T,h,m, 20, linewidth=0.5, extend='both')
	plt.axis("tight")
	plt.colorbar()
	plt.xlabel('Temperature (T)')
	plt.ylabel('External Field (h)')
	plt.title('Magnetization')

	plt.savefig('test.pdf')
	plt.show()
	return 0

if ( __name__ == '__main__'):
	main()