'''CMAES_evolStrat.py
 Carries out the process to find the quantities associated with
  the fitness and feasibility of solutions for the afib project
 For UWyo Ind. Study in Genetic Algorithms Spring 2018
'''

import numpy as np
from subprocess import run
from collections import Counter
from operator import itemgetter
from CMAES_xmlMake import getRuntimes, xmlwriter_checkQ
from createAblation import createAblationFile, setOnes

from matplotlib import pyplot as plt

import time


def analyzeSolutions(solutions, tissueWidth, ablationFile, paramFiles,
                     logFiles, batchToolPath, tissuePath, outDir):
    ''' For each solution, find:
         the number of ablated cells,
         number of quarantined cells,
         number of cells connected to the boundary,
         and the average run time of each solution over each paramfile
         then append each to a list to be returned'''
    ablations = len(solutions) * [0]  # Number cells ablated per solution
    quarantines = len(solutions) * [0]  # Num cells quarantined per solution
    connecteds = len(solutions) * [0]  # Num cells connected to edge
    runTimes = len(solutions) * [0]  # The mean activations averaged over
    #                                      all paramfiles

    for i, x in enumerate(solutions):
        # Generate the ablation set file
        createAblationFile(x, ablationFile)
        # count the number ablated and add to list
        t1 = time.time()
        cells = np.loadtxt(ablationFile)
        ablations[i] = np.sum(cells)
        t2 = time.time()
        print('Finding Ablated took', t2 - t1)
        quarantines[i] = getQuarantined(cells, tissueWidth)
        t3 = time.time()
        print('Finding quarantined took', t3 - t2)
        connecteds[i] = getContiguous(x)
        t4 = time.time()
        print('Finding connected took', t4 - t3)

        # List to hold the activations of each paramFile
        tempTimes = len(paramFiles) * [0]
        # Each solution needs to generalize over all the tissues presented
        for j, xmlFile in enumerate(paramFiles):
            # For redirecting batchtool output
            t5 = time.time()
            with open(logFiles[j], 'w') as file:
                #  Ensure the output csv file exists
                outputFile = os.path.join(
                    outDir, 'apTime_' + str(i) + '_' + str(j) + '.csv')
                run(['touch', outputFile])

                # List of arguments for subprocess.run
                runcommand = [batchToolPath +
                              '/VisibleEP_BatchTool_1', xmlFile]
                run(runcommand, stdout=file)  # Running tissue in the simulator
            t6 = time.time()
            print('Running in sim took', t6 - t5)
            tempTimes[j] = getRuntimes(logFiles[j])  # simulation steps
            # Relocate the output to make sure in the correct location
            run(['mv', outDir + 'simdata__sim0_input0_ApTimeFlag.csv',
                 outputFile])
        runTimes[i] = np.sum(tempTimes) / len(paramFiles)
        print('Finding runtime took', time.time() - t4)

    return ablations, quarantines, connecteds, runTimes


def is_feasible(numAblated, numQuarantined, tissueSize):
    ''' If any numAblated > 18% of the total cells or
        if any numQuarantined > 20% of the total cells
        return False.
        Otherwise, return True.
    '''
    percentAblated = numAblated / tissueSize
    percentQuarantined = numQuarantined / tissueSize
    if((percentAblated > 0.18) or (percentQuarantined > 0.2)):
        return False
    # Only return true if none of the values violate the conditions
    return True


def getContiguous(solution, tissueWidth):
    ''' Finds the number of cells connected to the boundary of the tissue
        using principles established for graph traversal.
        Input is a afib solution, 1-d list of 24 entries in range [0,1]
    '''
    # Create the point pairs from solutions
    xs = []
    ys = []
    for i in range(0, len(solution)):
        if i % 2 == 0:
            xs.append(solution[i])
        else:
            ys.append(solution[i])
    # Map all points from [0,1] to [0, tissueWidth-1] end inclusive
    for i in range(len(xs)):
        x = np.floor(tissueWidth-1 * xs[i])
        if x > tissueWidth-1:
            x = tissueWidth-1
        elif x < 0:
            x = 0
        else:
            x = int(x)
        xs[i] = x

    for i in range(len(ys)):
        y = np.floor(tissueWidth-1 * ys[i])
        if y > tissueWidth-1:
            y = tissueWidth-1
        elif y < 0:
            y = 0
        else:
            y = int(y)
        ys[i] = y

    # Check for points on the boundary
    xHas0 = xs.count(0)
    xHas79 = xs.count(tissueWidth-1)
    yHas0 = ys.count(0)
    yHas79 = ys.count(tissueWidth-1)

    if(xHas0 or xHas79 or yHas0 or yHas79):
        # Create the x-y point pairs
        points = np.array([(xs[i], ys[i])
                           for i in range(0, len(xs))])

        # Create the cell set
        cells = np.zeros((tissueWidth, tissueWidth))

        # Iterate over point pairs and set overlaid values to 1s
        for i in np.arange(0, len(points), 2):
            setOnes(cells, points[int(i)], points[int(i + 1)])

        # cells need to be transposed so x,y coordinates match with
        #  indices in the matrix
        cells = np.transpose(cells)

        frontier = []  # List of points to explore
        visited = {}  # Dict of points already visited
        numConnected = 0

        # Add every point on the boundary to the frontier
        for i, x in enumerate(xs):
            if((x == 0) or (x == tissueWidth-1)):
                frontier.append(points[i])
        for i, y in enumerate(ys):
            if((y == 0) or (y == tissueWidth-1)):
                frontier.append(points[i])

        while(len(frontier) > 0):
            # Pop the last element from the list into focus
            focus = frontier.pop(-1)

            # Skip over any points already visited
            while((str((focus[0], focus[1])) in visited.keys())):
                if(len(frontier) > 0):
                    focus = frontier.pop(-1)
                else:
                    # If there are no more points in the frontier,
                    #  return the number connected
                    return numConnected

            # Increment the number connected
            numConnected += 1
            # Add the current focus to the visited list
            visited[str((focus[0], focus[1]))] = True

            # Add points around the focus to frontier if the point is in the
            #  square defined by [0,tissueWidth-1], they are 1 in cells, and
            #  not visited
            # Move left
            if(focus[0] - 1 >= 0):
                testX = focus[0] - 1
                testY = focus[1]
                if((cells[testX, testY] == 1) and
                   not(str((testX, testY)) in visited.keys())):
                    frontier.append((testX, testY))
            # Move right
            if(focus[0] + 1 <= tissueWidth-1):
                testX = focus[0] + 1
                testY = focus[1]
                if((cells[testX, testY] == 1) and
                   not(str((testX, testY)) in visited.keys())):
                    frontier.append((testX, testY))
            # Move down
            if(focus[1] - 1 >= 0):
                testX = focus[0]
                testY = focus[1] - 1
                if((cells[testX, testY] == 1) and
                   not(str((testX, testY)) in visited.keys())):
                    frontier.append((testX, testY))
            # Move up
            if(focus[1] + 1 >= 0):
                testX = focus[0]
                testY = focus[1] + 1
                if((cells[testX, testY] == 1) and
                   not(str((testX, testY)) in visited.keys())):
                    frontier.append((testX, testY))

        return numConnected

    else:  # If no points are on the boundary, no cells are connected
        return 0


def getQuarantined(cells,  tissueWidth):
    '''
    Changes a matrix based on the ablated cells to find the separate groups
     and returns the number of cells that aren't ablated and aren't in the
     largest group
    '''
    # Ensure square
    cells = cells.reshape(tissueWidth, tissueWidth)

    unGrouped = []  # The cells we need to assign to a group
    ablated = []  # Cells ablated in ablnFile
    explored = []  # Cells we have alread "focused" on
    frontier = []  # Cells we are considering to group
    for i in range(0, tissueWidth):
        for j in range(0, tissueWidth):
            if(cells[i, j] == 0):
                unGrouped.append((i, j))
            else:
                ablated.append((i, j))

    groupId = 1  # Start at 1, the groupId for ablated cells
    while(len(unGrouped) > 0):
        groupId += 1  # Increment groupId to an unused value
        frontier.append(unGrouped[0])
        while(len(frontier) > 0):
            # Consider a new cell
            focus = frontier.pop()
            # Ensure not in explored already (shouldn't happen)
            while(focus in explored):
                if(len(frontier) > 0):
                    focus = frontier.pop()
                else:
                    break

            explored.append(focus)  # Now seen this cell
            cells[focus[0], focus[1]] = groupId  # Assign group
            if(focus in unGrouped):  # Safety
                unGrouped.remove(focus)  # Now grouped
            # Add points around the focus to frontier if the point is in the
            #  square defined by [0,tissueWidth-1], they are unablated,
            #  and not explored
            # Move left
            if(focus[0] - 1 >= 0):
                testX = focus[0] - 1
                testY = focus[1]
                if(not((testX, testY) in ablated) and
                   not((testX, testY) in explored)):
                    frontier.append((testX, testY))
            # Move right
            if(focus[0] + 1 <= tissueWidth-1):
                testX = focus[0] + 1
                testY = focus[1]
                if(not((testX, testY) in ablated) and
                   not((testX, testY) in explored)):
                    frontier.append((testX, testY))
            # Move down
            if(focus[1] - 1 >= 0):
                testX = focus[0]
                testY = focus[1] - 1
                if(not((testX, testY) in ablated) and
                   not((testX, testY) in explored)):
                    frontier.append((testX, testY))
            # Move up
            if(focus[1] + 1 <= tissueWidth-1):
                testX = focus[0]
                testY = focus[1] + 1
                if(not((testX, testY) in ablated) and
                   not((testX, testY) in explored)):
                    frontier.append((testX, testY))

    # Count the number in each group
    count = Counter(cells.flatten().tolist())
    # Remove the number ablated
    del count[1]
    # Get a list of (element, count) tuples
    count = list(count.items())
    # Ensure the list is sorted by count
    count.sort(key=itemgetter(1))
    count.reverse()
    # Quarantined are the total cells that can't be
    #  reached by the largest group
    quarantined = 0
    for i in range(1, len(count)):
        quarantined += count[i][1]

    return quarantined
