#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time
import numpy as np
from TSPClasses import *
import heapq
import itertools



class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
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
	
	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
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

	def greedy( self,time_allowance=60.0 ):
		pass
	
	
	
	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints: 
		max queue size, total number of states created, and number of pruned states.</returns> 
	'''
		
	def branchAndBound( self, time_allowance=60.0 ):
		start_time = time.time()
		cities = self._scenario.getCities()
		chosenRes = None
		currMin = float('inf')
		numResults = 0
		stateTotal = 0
		for city in cities:
			ncities = len(cities)
			x = [[cities[i].costTo(cities[j]) for j in range(ncities)] for i in range(ncities)]
			# we've built the matrix, now reduce it

			response = self.reduce2DArray(x)
			reducedX = response[0]
			reducedCost = response[1]

			# This response should have an array on 0 should have array of city indicis path and at 1 it should have the cost
			response = self.pathFinder(ncities, reducedX, set(), 0, reducedCost, [0], 1)
			stateTotal = stateTotal + response[2]
			if currMin == response[1]:
				numResults = numResults + 1
			if currMin > response[1]:
				numResults = 1
				currMin = response[1]
				chosenRes = response
			cities.append(cities.pop(0))
		# turn indices into array of cities
		routeArray = chosenRes[0]
		bssfCount = chosenRes[2]
		route = []
		for i in range(len(routeArray)):
			route.append(cities[routeArray[i]])
		# 	this result is like the bssf in the example default
		result = TSPSolution(route)
		result.cost = chosenRes[1]
		end_time = time.time()


		results = {}
		results['cost'] = result.cost
		results['time'] = end_time - start_time
		results['count'] = numResults
		results['soln'] = result
		results['max'] = bssfCount
		results['total'] = stateTotal
		results['pruned'] = stateTotal - bssfCount



		return results



	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''
		
	def fancy( self,time_allowance=60.0 ):
		pass
		

	# Returns a tuple of [reducedX, reducedCost]
	def reduce2DArray(self, x):
		# loop through i and reduce rows
		reducedCost = 0
		for i in range(len(x)):
			myMin = min(x[i])
			reducedCost = reducedCost + myMin
			for j in range(len(x)):
				x[i][j] = x[i][j] - myMin

		# now loop through the j using a transpose of x
		xTranspose = [[x[j][i] for j in range(len(x))] for i in range(len(x))]
		for i in range(len(x)):
			myMin = min(xTranspose[i])
			reducedCost = reducedCost + myMin
			for j in range(len(x)):
				xTranspose[i][j] = xTranspose[i][j] - myMin

		reducedX = [[xTranspose[j][i] for j in range(len(x))] for i in range(len(x))]
		return [reducedX, reducedCost]

	# routeArray should start with one index already at 0 and startIndex starts at 0
	def pathFinder(self, ncities, reducedX, pathSet, startIndex, currCost, routeArray, bssf):
		minDict = {}
		for n in range(ncities):
			if (n == startIndex):
				continue
			if (n in pathSet):
				continue
			minDict[n] = self.reducedCostTo(startIndex, n, reducedX, pathSet, currCost)
			bssf = bssf + 1
		winnerIndex = min(minDict, key=minDict.get)
		winnerCost = minDict[winnerIndex]
		pathSet.add(startIndex)
		pathSet.add(winnerIndex)
		routeArray.append(winnerIndex)
		if len(pathSet) != ncities:
			return self.pathFinder(ncities, reducedX, pathSet, winnerIndex, winnerCost, routeArray, bssf)
		else:
			return [routeArray, winnerCost, bssf]



	def reducedCostTo(self, startIndex, destIndex, reducedX, pathSet, prevCost):
		newX = [[0 for j in range(len(reducedX))] for i in range(len(reducedX))]
		for i in range(len(reducedX)):
			for j in range(len(reducedX)):
				if (startIndex == i):
					newX[i][j] = float('inf')
				elif(destIndex == j):
					newX[i][j] = float('inf')
				elif((startIndex == j) and (destIndex == i)):
					newX[i][j] = float('inf')
				elif((j in pathSet) and (destIndex == i)):
					newX[i][j] = float('inf')
				else:
					newX[i][j] = reducedX[i][j]
		response = self.reduce2DArray(reducedX)
		totalCost = reducedX[startIndex][destIndex] + prevCost + response[1]
		return totalCost
