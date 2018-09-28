''' hom_fileMake.py
Creates apd and resistance files offline for homogeneous tissue
For UWyo Ind. Study in Genetic Algorithms Spring 2018
'''

import os
import numpy as np
import subprocess
import CMAES_fileWrite as fW
import CMAES_xmlMake as xmlMake
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('start', metavar='str', type=int,
                        default=0,
                        help='Starting index for file names')
    parser.add_argument('stop', metavar='stp', type=int,
                        default=300,
                        help='Final index for file names (non-inclusive)')
    parser.add_argument('runDir', type=str)

    argDict = vars(parser.parse_args())
    startInd = int(argDict['start'])
    stopInd = int(argDict['stop'])
    runDir = argDict['runDir']

    APD = 130  # APD = 100+/-25
    RES = 9  # For uniform RC in homogeneous tissue
    duration = 2500  # We want 2.5 seconds of activation
    numTissues = 1  # We want 10 homogeneous tissues
    tissueHeight = 80  # We are 80 cells tall
    tissueWidth = 80  # 80 cells wide
    tissueSize = tissueHeight * tissueWidth  # Total number of cells in tissue

    # File creation & path formatting
    if(not os.path.exists(runDir)):
        subprocess.run(['mkdir', '-p', runDir])
    tissueFile = os.path.join('Tissues', '80x80_RC10.vepxml')
    ablnFile = 'abln.txt'
    apdFile = 'apd_'
    resFile = 'res_'
    consoleLogFile = 'consoleLog_'
    xmlFile = 'paramfile_'
    batchtoolpath = os.path.join('Batchtool', 'VepCore')
    homogXmlFiles = []
    homogabln = os.path.join(runDir, ablnFile)

    # Ablation File
    #  Write a column of 0s to the ablation file
    with open(homogabln, 'w') as f:
        # Write tissueSize zeros to ablation file
        for j in range(tissueSize):
            f.write("0\n")

    for i in range(startInd, stopInd):
        # print('Creating tissueNum:', i)
        # Setup the specific filenames for this run
        hapdF = apdFile + str(i) + '.txt'
        homogapd = os.path.join(runDir, hapdF)
        hresF = resFile + str(i) + '.txt'
        homogres = os.path.join(runDir, hresF)
        homogconsole = os.path.join(runDir, consoleLogFile + str(i) + '.txt')
        homogxml = os.path.join(runDir, xmlFile + str(i) + '.txt')

        # Create the initial files

        # Resistance file
        R = fW.editResistanceXML2(tissueFile)
        fW.createHomogResistanceFile(R, RES, homogres)
        # APD File
        fW.createHomogAPD(tissueSize, APD, homogapd)

        # See if run time goes long enough,  and
        #  if not regen APDfile,  startpos,  and numstims until it does
        # Randomize a new starting position for the stimulation
        startpos = [np.random.random_integers(0, tissueWidth / 2),
                    np.random.random_integers(0, tissueHeight / 2)]
        # Ranndomize the number of stimulation events at startpos
        numstims = 950 + round(10 * np.random.random()) * 10
        # Generate the xml file to find runtimes
        # # print('Calling xmlwriter_runtime')
        xmlMake.xmlwriter_runtime(homogabln, homogapd, homogres,
                                  tissueFile, runDir, startpos,
                                  homogxml, duration, numstims)
        # Find the runtime for this xml file
        # print('Calling findRuntime')
        runtime = xmlMake.findRuntime(homogconsole, homogxml, batchtoolpath)
        # print('runtime:', runtime)
        # Repeat the previous 4 commands until
        #  it runs for long enough
        while(runtime < duration):
            # # print("Looping runtime")
            fW.createHomogAPD(tissueSize, APD, homogapd)
            fW.createHomogResistanceFile(R, RES, homogres)
            startpos = [np.random.random_integers(0, tissueWidth),
                        np.random.random_integers(0, tissueHeight)]
            numstims = 950 + round(10 * np.random.random()) * 10
            # print('xmlmake')
            xmlMake.xmlwriter_runtime(homogabln, homogapd, homogres,
                                      tissueFile, runDir, startpos,
                                      homogxml, duration, numstims)
            # print('runtime')
            runtime = xmlMake.findRuntime(homogconsole, homogxml,
                                          batchtoolpath)
            # print('runtime:', runtime)
