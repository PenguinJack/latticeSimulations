#!/usr/bin/env python
#
# The main part of the program
# Scans over temperature and external field parameters in the Ising Model,
# and writes the magnetization to a csv file
#
# Author: Jack Jenkins



# Declare global variables
GRID_SIZE = 10	
STEPS = 10000			# Number of monte carlo steps.
NEIGHBOR_COUPLING = -1.	# Parameter J in the Ising Hamiltonian H=-JS_i S_j-h S_i

def main():
	import numpy as np
	import matplotlib.pyplot as plt
	import csv
	# inform user that the program is running
	print "Running monte carlo."

	# specify the temperatures and external fields
	# at which the magnetization and energy density will be calculated
	temperatures = np.arange(0.1, 5.0, 0.1)
	externalFields = np.arange(-1.025, 1.025, .05)

	# scan over temperatures and external fields
	data = []
	for T in temperatures:
		for h in externalFields:
			print "On temperature %f and external field %f" % (T, h)
			result = runMetropolis(GRID_SIZE, STEPS, T, NEIGHBOR_COUPLING, h)
			# store temperature, external field, magnetization, its uncertainty,
			# energy density, and its uncertainty
			data.append([T,h,result])
			#data.append([T,h,result[0],result[1],result[2],result[3]])

	# write to file
	with open("scan_data.csv", "w") as csvfile:
		writer = csv.writer(csvfile)
		[writer.writerow(["Temperature","External Field","Magnetization"])]
		[writer.writerow(r) for r in data]

	print "Data has been written to 'scan_data.csv'"

	return 0


def flip(x):
	return -x

def monteCarloStep(grid, gridSize, T, J, h):
# one step of the metropolis algorithm
# T is temperature, J is nearest-neighbor coupling, h is external magnetic field
	import math
	import random
	# pick a point on the lattice
	x = random.randint(0,gridSize-1)
	y = random.randint(0,gridSize-1)
	# calculate the energy cost for flipping that spin
	energy_difference = \
			2*J * grid[x][y] * (grid[x][(y+1) % gridSize] + grid[x][(y-1) % gridSize] \
			+ grid[(x+1) % gridSize][y] + grid[(x-1) % gridSize][y]) \
			+2*h * grid[x][y]
	# define probability for flipping a spin
	if energy_difference<0:
		prob_flip = 1
	else:
		prob_flip = math.exp(-energy_difference / T)
	
	if prob_flip > random.random():
		grid[x][y] = flip(grid[x][y]) 

	return grid

def runMetropolis(gridSize, steps, T, J, h):
# run the metropolis algorithm
	import numpy as np
	grid = initializeSpins(gridSize)
	magnetizationTimeseries = []
	energyDensityTimeseries = []
	for step in range(steps):
		grid=monteCarloStep(grid,gridSize, T, J, h)
		# control the sampling rate
		if step % (steps / 100) == 0 and step>steps/2:
			magnetizationTimeseries.append(magnetization(grid, gridSize))
	return np.average(magnetizationTimeseries)

# Initialize spins
def initializeSpins(gridSize):
	import numpy as np
	return np.random.choice((-1,1),size=(gridSize, gridSize))

# Calculate the magnetization
def magnetization(grid, gridSize):
	import numpy as np
	return np.sum(grid)/(gridSize**2.0)

# Calculate the correlation length (to do)

# Calculate the energy density
def energyDensity(grid, gridSize, J, h):
	result=0.
	for x in range(gridSize):
		for y in range(gridSize):
			result += -J * grid[x][y] * (grid[x][(y+1) % gridSize] + grid[(x+1)%gridSize][y]) \
					 - h*grid[x][y]
	return result/(gridSize**2.0)

if ( __name__ == '__main__'):
	main()