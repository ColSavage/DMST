#!/usr/bin/python3
import queue

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

from TSPClasses import *
import time
import numpy as np
from queue import PriorityQueue

import heapq
import itertools


class TSPSolver:
	def __init__(self, gui_view):
		self._scenario = None

	def setupWithScenario(self, scenario):
		self._scenario = scenario

	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution, 
		time spent to find solution, number of permutations tried during search, the 
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

	def defaultRandomTour(self, time_allowance=60.0):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time() - start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation(ncities)
			route = []
			# Now build the route using the random permutation
			for i in range(ncities):
				route.append(cities[perm[i]])
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results

	''' <summary>
		This is the entry point for the greedy solver, which you must implement for 
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

	def greedy(self, time_allowance=60.0):

		pass

	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints: 
		max queue size, total number of states created, and number of pruned states.</returns> 
	'''

	def branchAndBound(self, time_allowance=60.0):

		results = {}
		cities = self._scenario.getCities()
		nCities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()

		# Get initial bssf value using the randomTour function provided
		randomTourResults = self.defaultRandomTour(60.0)
		bssf = randomTourResults['cost']

		# TODO: Make a PQ, use bssf to prune states based on their lowerbound produced with the function below

		# This creates a square matrix array that is n x n cities big, all values initialized to zero
		costMatrix = np.zeros((nCities, nCities))

		# This for loop will add the cost to visit each city using the function costTo() found in the City class in TSPClasses.py
		for i in range(nCities):  # This would be the rows
			for j in range(nCities):  # This would be the columns
				cost = City.costTo(cities[i], cities[j])
				if cost != np.inf:
					costMatrix[i][j] = cost
				else:
					costMatrix[i][j] = np.inf

		# reduceCostMatrix returns a dictionary that contains, the reduced matrix and the bound found
		matrixData = self.reduceCostMatrix(nCities, costMatrix)
		currentMatrix = matrixData['reduced']
		# Create the parentState(Based off of the proj. specs), Consists of a matrix, bound, and partial tour
		parentState = {}
		parentState['bound'] = matrixData['bound']
		parentState['matrix'] = matrixData['reduced']
		parentState['tour'] = [0, 0]

		startCity = 0
		heapq = PriorityQueue()
		heapq.put(parentState)

		while heapq.empty() is not True:
			parentState = heapq.get()
			# Expand
			for i in range(1, nCities):
				if currentMatrix[startCity, i] != np.inf:
					childBound = parentState['bound'] + currentMatrix[startCity, i]

					# Since (startCity, i) was visited, we need to mark the matrix to say that we can no longer visit these cities, Make both row and column inf, using a copy of the parent matrix
					childMatrix = np.copy(parentState['matrix'])
					childMatrix[startCity:] = np.inf
					childMatrix[:i] = np.inf
					childMatrix[i:startCity] = np.inf

					# According to proj specs, after reducing matrix again, we create a new child state(Containing a matrix, bound, and partial tour), the tour is concatenated with the parentstate
					newData = self.reduceCostMatrix(nCities, childMatrix)
					childState = {}
					childState['bound'] = childBound			# TODO: May need to set this to the bound computed in the reduceCostMatrix function
					childState['matrix'] = newData['reduced']
					childState['tour'] = np.concatenate((parentState['tour'], [startCity, i]))

					# Compare bound with the bssf to see if we should add the state to the queue
					if childState['bound'] < bssf:
						heapq.put(childState)

			# Dequeue the most promising child state and follow the same process until a solution is found











		# heapq.put(stateOne)
		# while heapq.empty() is not True:
		# 	state = heapq.get()
		# 	# TODO: make an expand function that will get the points to visit in the row
		# 	if state['bound'] < bssf:
		# 		rowOfCities = state['reduced'][0:]
		# 		for city in rowOfCities:
		# 			# TODO: make a test function that will calculate the bound value



	# while not foundTour and time.time() - start_time < time_allowance:

	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''

	def fancy(self, time_allowance=60.0):
		pass

	def reduceCostMatrix(self, numCities, matrixToReduce):
		lowestCost = np.inf
		bound = 0
		data = {}
		rowValue = 0
		# This for loop will reduce the matrix by rows, finding the smallest value and then subtracting that value from each element in the row
		# TODO: Modify to use np.min() function
		# for k in range(numCities):
		# 	listOfLowest = matrixToReduce.min(axis=1)
		# 	for x in listOfLowest:
		# 		rowValue = np.linspace(x, x, numCities)
		# 		matrixToReduce[k] = np.subtract(matrixToReduce[k], rowValue)
		# 		lowestCost = np.inf
		# 		bound +=
		for i in range(numCities):
			for j in range(numCities):
				if matrixToReduce[i][j] < lowestCost:
					lowestCost = matrixToReduce[i][j]
					bound += lowestCost
			rowValue = np.linspace(lowestCost, lowestCost, numCities)
			matrixToReduce[i] = np.subtract(matrixToReduce[i], rowValue)
			lowestCost = np.inf

		# This for loop will reduce the matrix column wise if there are no zeros in then column
		for j in range(numCities):
			zeroFound = False
			lowestCost = np.inf
			for i in range(numCities):
				if matrixToReduce[i][j] == 0:
					zeroFound = True
					break
				elif matrixToReduce[i][j] != np.inf:
					if matrixToReduce[i][j] < lowestCost:
						lowestCost = matrixToReduce[i][j]
						bound += lowestCost
			if not zeroFound:
				columnValue = np.linspace(lowestCost, lowestCost, numCities)

				matrixColum = matrixToReduce[:, j]
				matrixToReduce[:, j] = np.subtract(matrixColum, columnValue)

		data["reduced"] = matrixToReduce
		data["bound"] = bound
		return data
