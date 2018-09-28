''' createAblation.py
Creates an ablation file for a 80x80 grid of simulated heart
 cells that are flattened to 1-d given a 24 entry 1-d numpy
 array with entries between [0,1]
 maps to [0, 79] and pairs up points.
 Output is a .txt file with a column of 1's, ablated, or 0's
'''

import numpy as np
import random


def setOnes(List, p1, p2, thickness=0.000001):
    '''
    Given a np.array to modify, 2 (x,y) points, and a thickness
     set the entries in List overlapped by the line created by
     the points to 1 if in the thickness
    '''
    theList = np.transpose(List)
    theList[int(p1[0]), int(p1[1])] = 1
    theList[int(p2[0]), int(p2[1])] = 1
    # find the min x-value
    if(p2[0] == p1[0]):  # vertical slope
        # Set all the yvalues in this range to 1
        for i in range(int(min(p1[1], p2[1])), int(max(p1[1], p2[1]))):
            theList[int(p1[0]), i] = 1
        List = np.transpose(theList)
        return  # Return early to avoid division by 0 in slope
    # calculate slope and y-intercept
    elif(p2[1] == p1[1]):  # horizontal slope
        for i in range(int(min(p1[0], p2[0])), int(max(p1[0], p2[0]))):
            theList[i, int(p2[1])] = 1
        List = np.transpose(theList)
        return
    elif(p1[0] < p2[0]):
        minp = p1
        maxp = p2
    elif(p2[0] < p1[0]):
        minp = p2
        maxp = p1
    slope = float(maxp[1] - minp[1]) / (maxp[0] - minp[0])
    if((abs(slope) < 2.0) & (abs(slope) >= 1)):
        thickness = 0.05
    elif ((abs(slope) >= 0.4) & (abs(slope) < 1)):
        thickness = 0.1
    elif(abs(slope) < 0.4):
        thickness = 0.5
    b = minp[1] - minp[0] * slope  # b = y0 - m*x0
    # iterate over the x-values between the two points with 1000 x-coordinates
    for i in np.arange(int(minp[0]), int(maxp[0]),
                       (int(maxp[0] - minp[0]) / 1000)):
        yval = slope * i + b  # Calculate yvalue for given x
        # Define the margin yvalues with safety ternary operation
        ymin = yval - thickness if (yval - thickness > 0) else 0
        ymax = yval + thickness if (yval + thickness < 79) else 79
        # Iterate over the y-values between the two points
        for j in range(int(min(minp[1], maxp[1])), int(max(minp[1], maxp[1]))):
            # Turn i into an appropriate integer
            x = int(round(i))
            # Set the cells overlapped by the margin to 1
            if(round(yval) == j):
                theList[x, j] = 1
            if(round(ymin) == j):
                theList[x, j] = 1
            if(round(ymax) == j):
                theList[x, j] = 1
    List = np.transpose(theList)


def createAblationFile(solution, filename='ablation.txt'):
    '''
    Given a CMA-ES solution as a 24 entry long bitstring
     write an ablation file based on the points defined
     by the solutoin
    '''
    # Create the cell set
    cells = np.zeros((80, 80))
    # Create the point pairs from solutions
    xs = []
    ys = []
    for i in range(0, len(solution)):
        if i % 2 == 0:
            xs.append(solution[i])
        else:
            ys.append(solution[i])
    # Map all points from [0,1] to [0,79] end inclusive
    for i in range(len(xs)):
        x = np.floor(79 * xs[i])
        if x > 79:
            x = 79
        elif x < 0:
            x = 0
        else:
            x = int(x)
        xs[i] = x

    for i in range(len(ys)):
        y = np.floor(79 * ys[i])
        if y > 79:
            y = 79
        elif y < 0:
            y = 0
        else:
            y = int(y)
        ys[i] = y

    points = np.array([(xs[i], ys[i])
                       for i in range(0, len(xs))])

    # Iterate over point pairs and set overlaid values to 1s
    for i in np.arange(0, len(points), 2):
        xs = [points[int(i)][0], points[int(i) + 1][0]]
        ys = [points[int(i)][1], points[int(i) + 1][1]]
        setOnes(cells, points[int(i)], points[int(i + 1)])

    # Write the resulting file to an output
    file = open(filename, 'w')
    cells = cells.flatten()
    for c in cells:
        file.write(str(int(c)))
        file.write('\n')
    file.close()


def createDisconnAblnFile(solution, file1='ablation_diss.txt',
                          file2='ablation_num.txt'):
    '''
    Creates two kinds of ablation files where any ablations that would be
    within 3 cells from the edge, are instead unablated. The two kinds of
    ablation files made are:
        1) All ablations within 3 cells are unablated
        2) Same as number one, but the number of cells ablated is made the
           same as the number of cells ablated normally by the solution
    '''
    # Create the cell set
    cells = np.zeros((80, 80))
    # Create the point pairs from solutions
    xs = []
    ys = []
    for i in range(0, len(solution)):
        if i % 2 == 0:
            xs.append(solution[i])
        else:
            ys.append(solution[i])
    # Map all points from [0,1] to [0,79] end inclusive
    for i in range(len(xs)):
        x = np.floor(79 * xs[i])
        if x > 79:
            x = 79
        elif x < 0:
            x = 0
        else:
            x = int(x)
        xs[i] = x
    for i in range(len(ys)):
        y = np.floor(79 * ys[i])
        if y > 79:
            y = 79
        elif y < 0:
            y = 0
        else:
            y = int(y)
        ys[i] = y

    # Create an array of points
    points = np.array([(xs[i], ys[i])
                       for i in range(0, len(xs))])

    # Iterate over point pairs and set overlaid values to 1s
    for i in np.arange(0, len(points), 2):
        xs = [points[int(i)][0], points[int(i) + 1][0]]
        ys = [points[int(i)][1], points[int(i) + 1][1]]
        setOnes(cells, points[int(i)], points[int(i + 1)])

    # Count the number of ones in cells
    totalAblated = np.sum(cells)
    disconnCells = np.zeros((80, 80))
    for i in range(3, 77):
        for j in range(3, 77):
            disconnCells[i][j] = cells[i][j]
    disconnAbln = np.sum(disconnCells)
    # Reatach cells ablated in disconnAbln to ablations in the valid area
    diff = totalAblated - disconnAbln
    sameNum = np.copy(disconnCells)
    # This loop searches for neighbors of ablated cells in the valid area
    #  to ablate while there is still a difference between the number of
    #  cells ablated in sameNum and totalAblated
    while(diff > 0):
        i, j = np.random.random_integers(4, 76, 2)
        if disconnCells[i, j] == 1:
            while(disconnCells[i, j] == 1):
                delI = random.choice([-1, 1])
                delJ = random.choice([-1, 1])
                i = i + delI
                j = j + delJ
            sameNum[i, j] = 1
            diff = diff - 1

    # Write the results to output files
    with open(file1, 'w') as f:
        t = disconnCells.flatten()
        for c in t:
            f.write('{}\n'.format(c))
    with open(file2, 'w') as f:
        t = sameNum.flatten()
        for c in t:
            f.write('{}\n'.format(c))
