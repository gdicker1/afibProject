''' CMAES_fileWrite.py
  Creates the resistance and apd files for the afib project
   for both the tissue with a patch (homogeneous) and without (heterogeneous)
  For UWyo Ind. Study in Genetic Algorithms Spring 2018
'''

import numpy as np
from numpy import random as npr
import csv
from subprocess import run
from createAblation import createAblationFile
from CMAES_evolStrat import analyzeSolutions, getContiguous
from CMAES_xmlMake import xmlwriter_checkQ


def editResistanceXML2(tissueFile):
    """
    Get the resistance values in the .vepxml file given by tissueFile.
    """
    R = []
    N = 0
    EDIT = False
    with open(tissueFile, 'r') as file:
        for line in file:
            if("</connections>" in line):
                return R
            elif(EDIT):
                if(N == 1):
                    # First entry is prefaced by <connections>
                    ind = line.index('>')
                    strR = line[ind + 1:]
                    strR = strR.split(',')
                    strR[2] = strR[2].split(';')[0]
                    R.append([float(s) for s in strR])
                    N = 0
                else:
                    strR = line.split(",")  # Resistances are comma separated
                    # Remove the semi-colon on the last value
                    strR[2] = strR[2].split(';')[0]
                    # Append a list of floats to R
                    R.append([float(s) for s in strR])
            if("</cells>" in line):
                EDIT = True
                N = 1
    return R


def createHomogResistanceFile(ResList, Resistance, ResFile):
    """
    Create the Resistance text file for a homogenous set of tissue
    """
    # The resistance to be replaced in the connection is sampled from a uniform
    #  distribution centered around Resistance and going as low as 3
    # Can't subtract less than 3 or ** Stability Error ** occurs
    # Print to resfile
    with open(ResFile, 'w') as file:
        for i in range(len(ResList)):
            R = ResList[i]
            # Replace the connection resistance with 9
            R[2] = Resistance
            if i == len(ResList) - 1:
                # The last line doesn't have a newline at the end
                theString = str(int(R[0])) + "," + str(int(R[1])) + "," + str(int(R[2])) + ";"
            else:
                theString = str(int(R[0])) + "," + str(int(R[1])) + "," + str(int(R[2])) + ";\n"
            file.write(theString)
    return


def createHeteroResistanceFile(ResList, Resistance, tissueHeight,
                               PatchDiam, PatchRes, ResFile):
    """
    Create the Resistance text file for a heterogeneous set of
    tissue with a patch of different Resistance (40% higher)
    in the middle
    """
    # Create a square array to represent the tissue
    tempRC = np.zeros((tissueHeight, tissueHeight))
    # Identify where the patch will start
    patchStart = np.floor(tissueHeight / 2) - 1
    # Set the values in the patch to NaN to make them easy to find
    for i in np.arange(patchStart, patchStart + PatchDiam, 1):
        for j in np.arange(patchStart, patchStart + PatchDiam, 1):
            tempRC[int(i), int(j)] = np.NaN
    # Flatten the array so the indices will match up with those in
    #  ResList (a flat array)
    tempRC = np.reshape(tempRC, (-1, 1))
    # Save the indices that would be in the patch
    patchFlatPos = []
    for i in range(len(tempRC)):
        if(tempRC[i, 0] == np.NaN):
            patchFlatPos.append(i)

    # Print to resfile
    with open(ResFile, 'w') as file:
        for i in range(len(ResList)):
            R = ResList[i]
            # Value to replace the resistance in this connection
            if R[0] in patchFlatPos:
                R[2] = PatchRes
            else:
                R[2] = Resistance
            if i == len(ResList) - 1:
                # The last line doesn't have a newline at the end
                theString = str(int(R[0])) +"," + str(int(R[1]))+ "," + str(int(R[2])) + ";"
            else:
                theString = str(int(R[0])) +"," + str(int(R[1])) + "," + str(int(R[2])) + ";\n"
            file.write(theString)
    return


def createHomogAPD(tissueSize, APDrefernce, APDfile):
    """
    Write a set of tissueSize APD values to the APDfile given a
    base line APD as APDreference
    """
    # Create a tissueSize number of random integers in APD-25,APD+25
    SetAPD = npr.random_integers(APDrefernce - 5, APDrefernce + 5,
                                 tissueSize)
    with open(APDfile, 'w') as file:
        for apd in SetAPD:
            theString = "{!s}\n"
            file.write(theString.format(int(apd)))
    return


def createHeteroAPD(tissueHeight, APDrefernce, PatchDiam,
                    PatchAPDref, APDfile):
    """
    Write a set of tissueHeight^2 APD values to APDfile with a patch of
    different APD in the middle
    """
    # Create a square array to represent the tissue
    setAPD = np.zeros((tissueHeight, tissueHeight))
    # Identify where the patch will start
    patchStart = np.floor(tissueHeight / 2) - 1
    # Set the values in the patch to NaN to make them easy to find
    for i in np.arange(patchStart, patchStart + PatchDiam, 1):
        for j in np.arange(patchStart, patchStart + PatchDiam, 1):
            setAPD[int(i), int(j)] = np.NaN
    # Flatten the array so the indices match up to ResList's (a flat array)
    setAPD = np.reshape(setAPD, (-1, 1))
    for i in range(len(setAPD)):
        if(setAPD[i, 0] == np.NaN):  # in the patch
            setAPD[i] = npr.random_integers(PatchAPDref - 5, PatchAPDref + 5)
        else:
            setAPD[i] = npr.random_integers(APDrefernce - 5, APDrefernce + 5)
    with open(APDfile, 'w') as file:
        for apd in setAPD:
            theString = "{!s}\n"
            file.write(theString.format(int(apd[0])))
    return


def addPatchAPD(homAdpFile, heterApdFile, tissueWidth,
                patchWidth, patchPos, patchAPD):
    ''' Add a patch of different APD to an existing
        homogeneous tissue and save to the heterApdFile
    '''
    # Set up patch start and stopping indices
    #  Note: the patch and tissue are assumed square
    patchXStart = patchPos
    patchXStop = patchPos + patchWidth
    patchYStart = patchPos
    patchYStop = patchPos + patchWidth

    # Load the apds and put in 2D form
    apdArr = np.loadtxt(homAdpFile)
    apdArr = np.reshape(apdArr, (tissueWidth, tissueWidth))

    # Iterate over the values in the patch randomly sampling
    #  from uniform [patchAPD-5, patchAPD+5]
    for i in range(patchXStart, patchXStop):
        for j in range(patchYStart, patchYStop):
            newApd = npr.random_integers(patchAPD - 5, patchAPD + 5)
            apdArr[i, j] = newApd

    # Flatten the apdArr again and write to heterApdFile, overwritting contents
    apdArr = apdArr.flatten()
    with open(heterApdFile, 'w') as f:
        for apd in apdArr:
            f.write(str(apd) + '\n')


def addPatchRes(homResFile, heterResFile, tissueWidth,
                patchWidth, patchPos, patchRes, nonPatchRes):
    ''' Add a patch of different resistance to an existing
        homogeneous tissue and save to heterResFile
    '''
    # Set up patch start and stopping indices
    #  Note: the patch and tissue are assumed square
    patchXStart = patchPos
    patchXStop = patchPos + patchWidth
    patchYStart = patchPos
    patchYStop = patchPos + patchWidth

    resList = np.genfromtxt(homResFile, delimiter=',', )

    # Create a square array to represent the tissue
    tempRC = np.zeros((tissueWidth, tissueWidth))
    for i in range(patchXStart, patchXStop):
        for j in range(patchYStart, patchYStop):
            tempRC[i][j] = np.NaN

    tempRC = tempRC.flatten()
    patchFlatPos = []
    for i in range(len(tempRC)):
        if(tempRC[i] == np.NaN):
            patchFlatPos.append(i)

    # Print to resfile
    with open(heterResFile, 'w') as file:
        for i in range(len(resList)):
            R = resList[i]
            # Value to replace the resistance in this connection
            if R[0] in patchFlatPos:
                R[2] = patchRes
            else:
                R[2] = nonPatchRes
            if i == len(resList) - 1:
                # The last line doesn't have a newline at the end
                theString = str(int(R[0])) +"," + str(int(R[1]))+ "," + str(int(R[2])) + ";"
            else:
                theString = str(int(R[0])) +"," + str(int(R[1])) + "," + str(int(R[2])) + ";\n"
            file.write(theString)
    return


def writeSolutionsHeader(file):
    with open(file, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        theRow = ['iteration', 'numAblated', 'numQuarintened',
                  'numConnected', 'fitness']
        for i in range(24):
            theRow.append('sol_' + str(i))
        writer.writerow(theRow)
    return


def writeSolutions(file, iteration, ablations, quarantines, connecteds,
                   fitnesses, solutions, numToWrite=None):
    '''Append the solutions to the given file as a row
       the first element is the iteration number, then fitness of solution,
       and finally the solution itself.
       If numToWrite is None, all are written. If numToWrite is finite,
        the top numToWrite solutions are written based on fitness.
    '''
    # Setup numToWrite for the default case
    if numToWrite is None:
        numToWrite = len(fitnesses)
    # Open the file
    with open(file, 'a', newline='') as f:
        # Create a csv writer
        writer = csv.writer(f, delimiter=',')
        # Create a list of fitnesses combined with rows sorted by fitness
        indices = list(range(len(fitnesses)))
        # Hacky way to sort indices by fitness
        newIndices = np.array(sorted(zip(fitnesses, indices)), dtype=int).T[1]
        numWritten = 0  # Track how many solutions written for early stop
        for i in newIndices:
            theRow = [iteration]
            theRow.extend([ablations[i], quarantines[i], connecteds[i],
                           fitnesses[i]])
            theRow.extend(list(solutions[i]))
            writer.writerow(theRow)  # Write the row
            numWritten += 1
            if(numWritten >= numToWrite):
                break  # Break if numToWrite has been reached
    return


def writeDataFile(dataFile, solutionsFile,
                  batchToolPath, tissuePath,
                  tissueSize):
    '''Copy iteration and fitness information from the solutionFile,
       use the solutions to generate the percent ablated,
       quarantined, and connected with the boundary,
       then write iteration, fitness, and percentages to dataFile.
    '''
    solutions = np.loadtxt(solutionsFile, delimiter=',')
    iterations = solutions.T[0].T
    fitnesses = solutions.T[1].T
    outDir = 'checkQ/'
    ablnFile = outDir + 'testAbln.txt'
    run(['touch', ablnFile])
    xmlFiles = []
    logFiles = []
    connecteds = []

    sols = [solutions[i][2:] for i in range(solutions.shape[0])]
    for i, s in enumerate(sols):
        createAblationFile(s, ablnFile)
        xmlFiles.append(outDir + 'checkQ_' + str(i) + '.csv')
        run(['touch', xmlFiles[-1]])
        logFiles.append(outDir + 'checkQ_' + str(i) + '.log')
        run(['touch', logFiles[-1]])

        xmlwriter_checkQ(ablnFile, outDir, tissuePath, xmlFiles[-1])

        connecteds.append(getContiguous(s, ablnFile))

    ablations, quarantines, _ = analyzeSolutions(sols, ablnFile, xmlFiles,
                                                 logFiles, batchToolPath,
                                                 outDir)

    with open(dataFile, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        indices = list(range(len(fitnesses)))
        for i in indices:
            theRow = [iterations[i], fitnesses[i],
                      ablations[i] / tissueSize, quarantines[i] / tissueSize,
                      connecteds[i] / tissueSize]
            writer.writerow[theRow]

    return
