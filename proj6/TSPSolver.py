#!/usr/bin/python3
import copy
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
import copy
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
        results = {}

        cities = self._scenario.getCities()
        ncities = len(cities)
        count = 1
        bssf = None
        start_time = time.time()

        for i in range(ncities):
            if time.time() - start_time >= time_allowance:
                break
            currentCity = cities[i]
            visited = set()
            route = []

            visited.add(currentCity)
            route.append(currentCity)
            breakEarly = False
            for k in range(ncities - 1):
                shortestPathSoFar = np.inf
                shortestPathSoFarIndex = -1
                for j in range(ncities):
                    # if city at j is not in visited, and if a path exists between currentCity and the cities[j], then
                    # compare it and update if needed.
                    if not cities[j] in visited:
                        if currentCity.costTo(cities[j]) < shortestPathSoFar:
                            shortestPathSoFar = currentCity.costTo(cities[j])
                            shortestPathSoFarIndex = j
                if shortestPathSoFarIndex == -1:  # no path starting at city i
                    breakEarly = True
                    break
                else:
                    currentCityIndex = shortestPathSoFarIndex
                    currentCity = cities[currentCityIndex]
                    visited.add(currentCity)
                    route.append(currentCity)
            if (breakEarly):
                continue

            # check to see if it connects back! If it does, you can exit. If not, try again
            if (currentCity.costTo(cities[i]) < np.inf):
                tempSol = TSPSolution(route)
                if (bssf == None or tempSol.cost < bssf.cost):
                    bssf = tempSol

            i += 1

        end_time = time.time()
        results['cost'] = bssf.cost
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    ''' <summary>
        This is the entry point for the branch-and-bound algorithm that you will implement
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution, 
        time spent to find best solution, total number solutions found during search (does
        not include the initial BSSF), the best solution found, and three more ints: 
        max queue size, total number of states created, and number of pruned states.</returns> 
    '''

    def branchAndBound(self, time_allowance=60.0):
        # Branch and bound TSP algorithm, branch to cover the search
        # space, but prune branches that will be less than optimal. Priority
        # of what to expand is determined using depth first, then bound.

        results = {}

        count = 0  # total number solutions found during search
        max = 1  # max queue size
        total = 1  # total number of states created
        pruned = 0  # number of pruned states
        updates = 0  # number of bssf updates

        # Initial bssf using the greedy algorithm
        greedyResults = self.greedy()  # Best case O(n^3) time, O(n^2) space, worst case O(n!)
        bssf = greedyResults['soln']

        start_time = time.time()

        # Prepare initial state - simplify initial matrix,
        # get inital bound, add first city to path
        matrix = self.initialMatrix()  # O(n^2) time and space
        bound, matrix = self.simplifyMatrix(0, matrix)  # O(n^2) time and space
        city = list(matrix.keys())[0]
        city_path = [city]
        index_path = [city._index]

        # Create state object, add to list and heapify
        state = State(bound, matrix, city_path, index_path)  # O(1) time, O(n^2) space
        PQ = [state]
        heapq.heapify(PQ)  # O(1) time and space (because we only put one element in)

        # Run the algorithm until we find optimal or run out of time
        # Loop is O(n^3) time and space for each iteration
        # Number of iterations is O(b^n) on average, O((n-1)!) at worst (or until time limit)
        while len(PQ) != 0 and time.time() - start_time < time_allowance:

            # Pop the highest priority state
            state = heapq.heappop(PQ)
            cost = state.bound

            # Prune state if its bound is already worse than the bssf
            if cost >= bssf.cost:
                pruned += 1
                continue

            # Get all the information about the state to expand new states
            matrix = state.matrix
            city_path = state.path
            index_path = state.index_path
            city = city_path[-1]

            # Loop over all the cities to consider all possible new states
            # O(n^3) time and space
            for i, c in enumerate(matrix[city]):  # loop runs O(n) times

                # If we have already visited the city, move on
                if i in index_path:  # O(n)
                    continue

                # If the city is unreachable, don't expand that path (prune)
                if c == math.inf:
                    pruned += 1
                    total += 1
                    continue

                # Here the city is unvisited and reachable, let's explore
                newCity = list(matrix.keys())[i]
                newCost = cost + c
                newPath = city_path + [newCity]  # O(1) time, O(n) space
                newIndex = index_path + [i]  # O(1) time, O(n) space

                # If this is the last city, and it makes a cycle, it's a solution
                # O(1) time and space
                if len(newPath) == len(matrix.keys()) and self.fullCycle(newPath, matrix):

                    # Add the cost of returning to the first city
                    newCost += matrix[newPath[-1]][newPath[0]._index]
                    count += 1

                    # Update the bssf if the new solution is better
                    if newCost < bssf.cost:
                        bssf = TSPSolution(newPath)
                        updates += 1
                    continue

                # If still a partial path, update cost and matrix with city in path
                newMatrix = self.StrikeOut(matrix, city, newCity)  # O(n^2) time and space
                newCost, newMatrix = self.simplifyMatrix(newCost, newMatrix)  # O(n^2) time and space

                # If the cost of the partial path is still less than the bssf
                # create a new state and add it to the queue
                if newCost < bssf.cost:
                    newState = State(newCost, newMatrix, newPath, newIndex)  # O(1) time, O(n^2) space
                    heapq.heappush(PQ, newState)
                    total += 1
                    if len(PQ) > max:
                        max = len(PQ)

                # Otherwise, prune because it won't be optimal
                else:
                    pruned += 1
                    total += 1

        end_time = time.time()

        # Prune the remaining items on the queue
        pruned += len(PQ)

        # Tally up the results and return them
        results['cost'] = bssf.cost
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = max
        results['total'] = total
        results['pruned'] = pruned

        return results

    def initialMatrix(self):  # O(n^2) time and space
        # Create the initial matrix by getting the value of each edge in the
        # graph and putting it in a dictionary with city objects as the keys
        # and a list of distances (in order of the cities) as the values

        cities = self._scenario.getCities()
        matrix = {}
        for city in cities:
            matrix[city] = []
            for city2 in cities:
                matrix[city].append(city.costTo(city2))  # O(n^2) - nested for loops
        return matrix

    def simplifyMatrix(self, bound, matrix):  # O(n^2) time and space
        # Simplify the matrix by subtracting the the minimum of each row from itself
        # then the minumum of each column from itself, and add each minimum to the bound
        # because its going to take at least those minumum values for the whole tour

        newMatrix = {}
        for city, weights in matrix.items():
            minimum = min(weights)
            if minimum == math.inf:
                newWeights = [math.inf] * len(weights)
            else:
                bound += minimum
                newWeights = []
                for weight in weights:
                    newWeights.append(weight - minimum)  # O(n^2) - nested for loops
            newMatrix[city] = newWeights

        for city in newMatrix:
            column = []
            for city2, weights in newMatrix.items():
                column.append(weights[city._index])  # O(n^2) - nested for loops
            minimum = min(column)
            if minimum == math.inf:
                continue
            bound += minimum
            for city2, weights in newMatrix.items():
                newMatrix[city2][city._index] -= minimum

        return bound, newMatrix

    def StrikeOut(self, matrix, city_from, city_to):  # O(n^2) time and space
        # Change the row of the "from" city and the column of the "to"
        # city to infinity to prevent coming from/to those cities again
        # (also the cell that goes right back to the "from" city)

        newMatrix = {}
        for city, weights in matrix.items():
            newWeights = []
            for i, weight in enumerate(weights):
                if city == city_from or i == city_to._index:
                    newWeights.append(math.inf)
                    continue
                if city == city_to and i == city_from._index:
                    newWeights.append(math.inf)
                    continue
                newWeights.append(weight)

            newMatrix[city] = newWeights

        return newMatrix

    def fullCycle(self, path, matrix):  # O(1) time and space
        # Checks if the path makes a full cycle (does the last city
        # have a path back to the first city), returns a boolean
        return matrix[path[-1]][path[0]._index] != math.inf

    ''' <summary>
        This is the entry point for the algorithm you'll write for your group project.
        </summary>
        <returns>results dictionary for GUI that contains three ints: cost of best solution, 
        time spent to find best solution, total number of solutions found during search, the 
        best solution found.  You may use the other three field however you like.
        algorithm</returns> 
    '''

    '''
    1. Let 1 be the starting and ending point for salesman. 
    2. Construct MST from with 1 as root using Prim’s Algorithm.
    3. List vertices visited in preorder walk of the constructed MST and add 1 at the end. 
    
    '''

    def fancy(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities()

        currentBest = self.greedy(time_allowance)  # Get the results from greedy, use this to compare
        bssf = currentBest['soln']
        changed = True

        results = currentBest
        route = copy.deepcopy(bssf.route)
        count = 0
        startTime = time.time()

        while changed is True and time.time() - startTime < time_allowance:
            changed = False
            for i in range(0, len(route) - 1):                                                                          # O(n^2)
                shouldBreak = False
                for j in range(i + 1, len(route)):
                    route = copy.deepcopy(bssf.route)

                    if i != j:
                        temp = route[i]
                        route[i] = route[j]
                        route[j] = temp

                        tempSolution = TSPSolution(route)
                        if tempSolution.cost < bssf.cost:
                            bssf = tempSolution
                            shouldBreak = True
                            changed = True
                            count = count + 1
                            break
                if shouldBreak:
                    break

                    #11 0 8 1

        endTime = time.time()

        results['cost'] = bssf.cost
        results['time'] = endTime - startTime
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None

        return results


class State:
    # A state consists of a bound, matrix, and partial path (cities and indices)
    # The priority of the state is a tuple with the depth (how many more cities need to be
    # added to the path, i.e. lower number is better) and the bound (again, lower is better)
    def __init__(self, bound, matrix, path, index_path):  # O(1) time, O(n^2) space
        self.bound = bound
        self.matrix = matrix
        self.path = path
        self.index_path = index_path
        self.priority = (len(matrix) - len(path), bound)

    # A state is compared based on its priority (described above)
    def __lt__(self, other):  # O(1) time and space
        return self.priority < other.priority
