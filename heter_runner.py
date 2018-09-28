'''
heter_runner.py
Starts the python process for evolving ablation sets on simulated heart
    tissue with a patch of altered tissue. For UWyo Ind. Study in
    Genetic Algorithms Spring 2018
'''

import numpy as np
import cma
from subprocess import run
import argparse
import time
import os
import CMAES_fileWrite as fW
import CMAES_xmlMake as xmlMake
import CMAES_evolStrat as evS

if __name__ == '__main__':
    # Read in information from command-line and assign them to variables
    parser = argparse.ArgumentParser()
    parser.add_argument('batchtoolpath',
                        help='Absolute path to the VisibleEP_BatchTool_1 file',
                        type=str)
    parser.add_argument('tissuepath',
                        help='Absolute path to the directory containing the reference tissue and homogeneous tissue files',
                        type=str)
    parser.add_argument('start',
                        help='Starting index of files from tissueDir to use for the run. Must be in [0,299]',
                        type=int)
    parser.add_argument('end',
                        help='Ending index of files from tissueDir to use for the run. Must be in [0,299] and greater than start',
                        type=int)
    argDict = vars(parser.parse_args())
    batchtoolpath = argDict['batchtoolpath']
    tissuePath = argDict['tissuepath']
    startInd = argDict['start']
    stopInd = argDict['end']

    # Error checking of the command-line arguments
    if(not(0 <= startInd) or not(300 > startInd)):
        raise argparse.ArgumentTypeError('start not in the correct range [0,299]')
    if(not(0 <= stopInd) or not(300 > stopInd)):
        raise argparse.ArgumentTypeError('end not in the correct range [0,299]')
    if(not(startInd <= stopInd)):
        raise argparse.ArgumentTypeError('start not less than or equal to end')
    # Ensure paths don't end in /
    if(batchtoolpath[-1] == '/'):
        batchtoolpath = batchtoolpath[:-1]
    if(tissuePath[-1] == '/'):
        tissuePath = tissuePath[:-1]

    times = []  # Used when debugging to figure out how long each part takes
    times.append(time.time())
    print('Starting heterogeneous CMA-ES_runnner for the 2018 Afib Project')
    print('Geo. Dylan Dickerson at', times[0])
    print('batchtoolpath=', batchtoolpath)
    print('tissuepath=', tissuePath)
    print('startInd=', startInd)
    print('stopInd=', stopInd)
    # set the tissue and patch parameters
    APD = 130  # APD = 100+/-25
    RES = 9  # For uniform RC in homogeneous tissue
    PAPD = 80  # Patch APD = 50+/-25
    PRES = 13  # Patch RES is 40% higher than RES
    duration = 10000  # Let simulations run for at most 10 seconds
    MWRduration = 2500  # Minimum time for multi-wavelet rentry
    startRecord = 1000  # When to start the cell record
    stopRecord = 1500  # When to stop recording
    tissueHeight = 80  # We are 80 cells tall
    tissueWidth = 80  # 80 cells wide
    tissueSize = tissueHeight * tissueWidth  # Total number of cells in tissue
    # Files that es.logger.add() will create.
    # These need to be moved into Text_Files at the end of a run
    cmaesOutputFiles = ['outcmaesaxlencorr.dat', 'outcmaesaxlen.dat',
                        'outcmaesfit.dat', 'outcmaesstddev.dat',
                        'outcmaesxmean.dat', 'outcmaesxrecentbest.dat']

    # Create the file names and dirs
    times.append(time.time())
    print('Creating filenames at', times[-1])
    # Directory to put new files in, makes sure the current working dir is used
    runDir = os.getcwd() + '/'
    print('runDir=', runDir)
    apDir = os.path.join(runDir, 'APTimes')
    print('apDir=', apDir)
    apInitDir = os.path.join(apDir, 'initial')
    apMidDir = os.path.join(apDir, 'midPoint')
    apFinDir = os.path.join(apDir, 'final')
    # Ensure the checkQ directory exists
    if(not('checkQ' in os.listdir(runDir))):
        run(['mkdir', 'checkQ'])
    # Make directories to save the info to create gifs
    if(not('APTimes' in os.listdir(runDir))):
        run(['mkdir', 'APTimes'])
        run(['mkdir', apInitDir])
        run(['mkdir', apMidDir])
        run(['mkdir', apFinDir])
    tissueFile = os.path.join(tissuePath, '80x80_RC10.vepxml')
    ablnFile = 'abln.txt'
    apdFile = 'apd_'
    resFile = 'res_'
    consoleLogFile = 'consoleLog_'
    xmlFile = 'paramfile_'

    heterabln = runDir + ablnFile
    heterXmlFiles = []
    heterConsoleFiles = []
    hetCmaLoggingFile = 'cmaConsoleLog.txt'
    solutionFile = runDir + 'solutions.csv'
    # times.append(time.time())
    # print('Creating files took: ', # times[-2] - # times[-1])

    # print('Creating tissues at', # times[-1])
    # # Generate the data for heterogeneous tissues
    # Generate a patch position either at 10, 20, 30, or 40 from 0,0
    patchPos = 10 * np.random.random_integers(1, 5)
    print('Patch is at ({}, {})\n'.format(patchPos, patchPos))
    for i in range(startInd, stopInd + 1):
        # Setup the specific filenames for this run
        heterapd = apdFile + str(i) + '.txt'
        heterres = resFile + str(i) + '.txt'
        runtimeFile = 'xmlRT.txt'
        heterconsole = runDir + consoleLogFile + str(i) + '.txt'
        heterxml = runDir + xmlFile + str(i) + '.txt'
        # Get the reference APD and Res files in the homogeneous dir
        homogapd = os.path.join(tissuePath, 'homogeneous',
                                apdFile + str(i) + '.txt')
        homogres = os.path.join(tissuePath, 'homogeneous',
                                resFile + str(i) + '.txt')

        # See if run time goes long enough,  and
        #  if not regen APDfile,  startpos,  and numstims until it does
        # Randomize a new starting position for the stimulation
        with open(heterabln, 'w') as f:
            # Write tissueSize zeros to ablation file
            for j in range(tissueSize):
                f.write("0\n")

        # Create heterres file
        fW.addPatchRes(homogres, heterres, tissueWidth,
                       20, patchPos, PRES, RES)

        # Create heterapd
        fW.addPatchAPD(homogapd, heterapd, tissueWidth,
                       20, patchPos, PAPD)

        # Randomize a new starting position for the stimulation
        startpos = [np.random.random_integers(0, tissueWidth),
                    np.random.random_integers(0, tissueHeight)]
        # Randomize the number of stimulation events at startpos
        numstims = 950 + round(10 * np.random.random()) * 10
        # Generate the xml file to find run# times
        # print('Making heter xml for runtime')
        xmlMake.xmlwriter_runtime(heterabln, heterapd, heterres,
                                  tissueFile, runDir, startpos,
                                  runtimeFile, duration, numstims)
        # Find the runtime for this xml file
        # print('Finding runtime')
        runtime = xmlMake.findRuntime(heterconsole, runtimeFile, batchtoolpath)

        # print('runtime:', runtime)
        # Repeat the finding a patch position, creating res&apd files,
        #  randomizing startpos, and num stims until
        #  we get the tissue to run for long enough
        while(runtime < MWRduration):
            # print('Looping runtime')

            # Create heterres file
            fW.addPatchRes(homogres, heterres, tissueWidth,
                           20, patchPos, PRES, RES)

            # Create heterapd
            fW.addPatchAPD(homogapd, heterapd, tissueWidth,
                           20, patchPos, PAPD)

            # Randomize a new starting position for the stimulation
            startpos = [np.random.random_integers(0, tissueWidth),
                        np.random.random_integers(0, tissueHeight)]
            # Randomize the number of stimulation events at startpos
            numstims = 950 + round(10 * np.random.random()) * 10
            # Generate the xml file to find run# times
            # print('Making heter xml for runtime')
            xmlMake.xmlwriter_runtime(heterabln, heterapd, heterres,
                                      tissueFile, os.getcwd(), startpos,
                                      runtimeFile, duration, numstims)
            # Find the runtime for this xml file
            # print('Finding runtime')
            runtime = xmlMake.findRuntime(heterconsole, runtimeFile,
                                          batchtoolpath)
            # print('runtime:', runtime)
        xmlMake.runtimeToVM(runtimeFile, heterxml, duration, tissueFile,
                            'ApTimeFlag', runDir, startRecord, stopRecord,
                            heterabln, heterapd, heterres)
        # print('Appending xmlfile', i)
        heterXmlFiles.append(heterxml)
        heterConsoleFiles.append(heterconsole)
    # # times.append(time.time())
    # print('Creating tissues took', # times[-2] - # times[-1])
    # print('Starting CMA-ES strategy at', # times[-1])

    # CMA-ES loop
    # Create an evolutionary strategy with
    #  a genome 24 entries long sampled from uniform [0.3,0.7],
    #  a sigma of 0.2,
    #  number of children made per gen = 15,
    #  number of parent slots = 7
    #  Restrict the values to the range [0,1]

    # Ensure the initial sample is within quarantine and ablation limits
    x0 = np.random.uniform(0.3, 0.7, 24)
    ablns, quars, _, _ = evS.analyzeSolutions([x0], ablnFile, heterXmlFiles,
                                              heterConsoleFiles,
                                              batchtoolpath, tissueFile,
                                              runDir)
    while(not evS.is_feasible(ablns[0], quars[0], tissueSize)):
        x0 = np.random.uniform(0.3, 0.7, 24)
        ablns, quars, _, _ = evS.analyzeSolutions([x0], ablnFile,
                                                  heterXmlFiles,
                                                  heterConsoleFiles,
                                                  batchtoolpath, tissueFile,
                                                  runDir)
    for i in range(len(heterXmlFiles)):
        apFile = 'apTime_0_' + str(i) + '.csv'
        run(['mv', apFile, apInitDir])
    # Create an evolutionary strategy with
    #  a genome 24 entries long,  a sigma of 0.1,
    #  number of children made per gen = 15,
    #  number of parent slots = popsize // 2 = 7
    #  Restrict the values to the range [0,1]
    es = cma.CMAEvolutionStrategy(x0, 0.2,
                                  {'popsize': 15,
                                   'bounds': [0, 1],
                                   'maxiter': 10000})
    # times.append(time.time())
    # print('Defining evolutionary strategy took', # times[-2] - # times[-1])

    # This loop contains the main CMA-ES steps.
    #  New solutions (children) are generated,
    #  then their fitnesses are calculated,
    #  and then relevant info is logged both by cmaes and in a
    #  custom file that logs iteration, fitness, solution for
    #  each of the 15 individuals in the combined parent, child population
    iteration = 0
    midPoint = 200
    fW.writeSolutionsHeader(solutionFile)
    while not es.stop():
        loopStartTime = time.time()
        times.append(loopStartTime)
        solutions = es.ask()
        times.append(time.time())
        ablns, quars, conns, runTimes = evS.analyzeSolutions(
                                             solutions, ablnFile,
                                             heterXmlFiles, heterConsoleFiles,
                                             batchtoolpath, tissueFile,
                                             runDir)
        times.append(time.time())
        fitnesses = len(solutions) * [1]

        for i in range(len(fitnesses)):
            if(evS.is_feasible(ablns[i], quars[i], tissueSize)):
                # Only calculate actual fitness if the solution is feasible
                fitnesses[i] = runTimes[i] / duration
            else:
                # If not feasible bias solution away by assigning high fitness
                fitnesses[i] = 1
        es.tell(solutions, fitnesses)
        # Only log info and update iteration if all feasible
        if(not(len(solutions) == fitnesses.count(1))):
            if iteration % 10 == 0:
                # Log the solution with my custom code
                fW.writeSolutions(solutionFile, iteration,
                                  ablns, quars, conns,
                                  fitnesses, solutions)
            # If at midPoint, record info to create gifs
            if iteration == midPoint:
                for i in range(len(solutions)):
                    for j in range(len(heterXmlFiles)):
                        apFile = 'apTime_' + str(i) + '_' + str(j) + '.csv'
                        run(['mv', apFile, apMidDir])
            es.logger.add()  # Log with the built in CMAES logger
            iteration += 1  # Increase the iteration
        times.append(time.time())
        print('CMAES Loop #{} took {} \n'.format(iteration,
                                                 times[-1] - loopStartTime))
    for i in range(len(solutions)):
        for j in range(len(heterXmlFiles)):
            apFile = 'apTime_' + str(i) + '_' + str(j) + '.csv'
            run(['mv', apFile, apFinDir])
    print('Entire Process took:', times[-1] - times[0])
