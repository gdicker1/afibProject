''' CMAES_xmlMake.py
Uses lxml to create 3 xml file types to be called by the BatchTool
    creates one for total data collection of activations
            one without saving data output to find how long a tissue runs for
            one to attempt to identify which cells are quarintined and never
             activate.
For UWyo Ind. Study in Genetic Algorithms Spring 2018
'''

from lxml import etree
import numpy as np
from contextlib import redirect_stdout
from subprocess import run


def xmlwriter_VM(albntextfile, apdtextfile, restextfile, tissuefile, outdir,
                 startposition, xmlfile, duration, numstims, output):
    ''' Writes an xml file for the batch tool that modifies the
        ADP and resitance of a reference tissue and save the output
        to a directory
    '''
    CardEP = etree.Element("cardiac-EP-model")
    Simulations = etree.SubElement(CardEP, "simulations",
                                   concurrency="1",)
    Simulation = etree.SubElement(Simulations, "simulation",
                                  iterations=str(duration),
                                  triggeredStop="true")
    # Define the input and output files
    Inputs = etree.SubElement(Simulation, "inputs")
    Input = etree.SubElement(Inputs, "input", format="xml")
    Input.text = tissuefile
    Outputs = etree.SubElement(Simulation, "outputs")
    Output = etree.SubElement(Outputs, "output", type=output,
                              dir=outdir, startIndex="1000",
                              stopIndex=str(duration))
    Events = etree.SubElement(Simulation, "events")
    # Move Electrode Event
    MoveEvent = etree.SubElement(Events, "event",
                                 occurAt=str(0), type="move-electrode")
    MoveIndex = etree.SubElement(MoveEvent, "index")
    MoveIndex.text = str(0)
    Xpos = etree.SubElement(MoveEvent, "x")
    Xpos.text = str(startposition[0])
    Ypos = etree.SubElement(MoveEvent, "y")
    Ypos.text = str(startposition[1])
    Zpos = etree.SubElement(MoveEvent, "z")
    Zpos.text = str(0.01)

    # Ablation file
    AblEvent = etree.SubElement(Events, "event",
                                occurAt=str(0), type="ablation")
    AblFile = etree.SubElement(AblEvent, "file")
    AblFile.text = albntextfile

    # Nrepol (APD) file
    ADPEvent = etree.SubElement(Events, "event",
                                occurAt=str(0), type="nrepol")
    APDFile = etree.SubElement(ADPEvent, "file")
    APDFile.text = apdtextfile

    # Resistance file
    ResEvent = etree.SubElement(Events, "event",
                                occurAt=str(0), type="change-resistance")
    ResFile = etree.SubElement(ResEvent, "file")
    ResFile.text = restextfile

    # Stimulation events
    for i in np.arange(0, numstims, 10):
        StimEvent = etree.SubElement(Events, "event",
                                     occurAt=str(i),
                                     type="stimulate-electrode")
        StimIndex = etree.SubElement(StimEvent, "index")
        StimIndex.text = str(0)
        Magnitude = etree.SubElement(StimEvent, "magnitude")
        Magnitude.text = str(0.2)

    # Write to the xmlfile
    tree = etree.ElementTree(CardEP)
    tree.write(xmlfile, pretty_print=True,
               xml_declaration=True, encoding="utf-8")
    return


def xmlwriter_runtime(albntextfile, apdtextfile, restextfile, tissuefile,
                      outdir, startposition, xmlfile, duration, numstims):
    ''' Writes an xml file for the batch tool that modifies the
        ADP and resitance of a reference tissue and does not save
        the output to a directory
    '''
    CardEP = etree.Element("cardiac-EP-model")
    Simulations = etree.SubElement(CardEP, "simulations",
                                   concurrency="1",)
    Simulation = etree.SubElement(Simulations, "simulation",
                                  iterations=str(duration),
                                  triggeredStop="true")
    # Define the input and output files
    Inputs = etree.SubElement(Simulation, "inputs")
    Input = etree.SubElement(Inputs, "input", format="xml")
    Input.text = tissuefile
    Outputs = etree.SubElement(Simulation, "outputs")
    Events = etree.SubElement(Simulation, "events")
    # Move Electrode Event
    MoveEvent = etree.SubElement(Events, "event",
                                 occurAt=str(0), type="move-electrode")
    MoveIndex = etree.SubElement(MoveEvent, "index")
    MoveIndex.text = str(0)
    Xpos = etree.SubElement(MoveEvent, "x")
    Xpos.text = str(startposition[0])
    Ypos = etree.SubElement(MoveEvent, "y")
    Ypos.text = str(startposition[1])
    Zpos = etree.SubElement(MoveEvent, "z")
    Zpos.text = str(0.01)

    # Ablation file
    AblEvent = etree.SubElement(Events, "event",
                                occurAt=str(0), type="ablation")
    AblFile = etree.SubElement(AblEvent, "file")
    AblFile.text = albntextfile

    # Nrepol (APD) file
    ADPEvent = etree.SubElement(Events, "event",
                                occurAt=str(0), type="nrepol")
    APDFile = etree.SubElement(ADPEvent, "file")
    APDFile.text = apdtextfile

    # Resistance file
    ResEvent = etree.SubElement(Events, "event",
                                occurAt=str(0), type="change-resistance")
    ResFile = etree.SubElement(ResEvent, "file")
    ResFile.text = restextfile

    # Stimulation events
    for i in np.arange(0, numstims, 10):
        StimEvent = etree.SubElement(Events, "event",
                                     occurAt=str(i),
                                     type="stimulate-electrode")
        StimIndex = etree.SubElement(StimEvent, "index")
        StimIndex.text = str(0)
        Magnitude = etree.SubElement(StimEvent, "magnitude")
        Magnitude.text = str(0.2)

    # Write to the xmlfile
    tree = etree.ElementTree(CardEP)
    tree.write(xmlfile, pretty_print=True,
               xml_declaration=True, encoding="utf-8")
    return


def xmlwriter_checkQ(albntextfile, outdir, tissuefile, xmlfile):
    '''Create a xml file for the batch tool that doesn't modify
       ADP or resistance in the tissue file and runs for a shorter
       duration with only 1 stimulation event in the middle of one
       edge
    '''
    CardEP = etree.Element("cardiac-EP-model")
    Simulations = etree.SubElement(CardEP, "simulations",
                                   concurrency="1",)
    Simulation = etree.SubElement(Simulations, "simulation",
                                  iterations=str(1000), triggeredStop="true")
    # Define the input and output files
    Inputs = etree.SubElement(Simulation, "inputs")
    Input = etree.SubElement(Inputs, "input", format="xml")
    Input.text = tissuefile
    Outputs = etree.SubElement(Simulation, "outputs")
    Output = etree.SubElement(Outputs, "output", type='ApTimeFlagSum',
                              dir=outdir, startIndex="0", stopIndex=str(1000))
    Events = etree.SubElement(Simulation, "events")
    # Move Electrode Event
    MoveEvent1 = etree.SubElement(Events, "event",
                                  occurAt=str(0), type="move-electrode")
    MoveIndex1 = etree.SubElement(MoveEvent1, "index")
    MoveIndex1.text = str(0)
    Xpos1 = etree.SubElement(MoveEvent1, "x")
    Xpos1.text = str(0)
    Ypos1 = etree.SubElement(MoveEvent1, "y")
    Ypos1.text = str(39)
    Zpos1 = etree.SubElement(MoveEvent1, "z")
    Zpos1.text = str(0.01)

    # Ablation file
    AblEvent = etree.SubElement(Events, "event",
                                occurAt=str(0), type="ablation")
    AblFile = etree.SubElement(AblEvent, "file")
    AblFile.text = albntextfile

    # Stimulation events
    StimEvent1 = etree.SubElement(Events, "event",
                                  occurAt=str(0),
                                  type="stimulate-electrode")
    StimIndex1 = etree.SubElement(StimEvent1, "index")
    StimIndex1.text = str(0)
    Magnitude1 = etree.SubElement(StimEvent1, "magnitude")
    Magnitude1.text = str(0.2)

    # Write to the xmlfile
    tree = etree.ElementTree(CardEP)
    open(xmlfile, "w")
    tree.write(xmlfile, pretty_print=True,
               xml_declaration=True, encoding="utf-8")
    return


def runtimeToVM(runtimeFile, vmFile,
                duration,
                inputFile,
                outputType, outputDir, start, stop,
                albntextfile, apdtextfile, restextfile):
    '''Create a copy of previously existing runtime xml file,
       add an output directory and modify the ablation, apd change,
       and resistance change events so they point to the correct files.
    '''
    # Read in the runtimeFile as an ElementTree
    with open(runtimeFile, 'r') as f:
        tree = etree.parse(f)
    root = tree.getroot()
    # The simulation element is the "grandchild" of the root
    simulationElement = root[0][0]
    # input is the first "grandchild" of simulation
    inputElement = simulationElement[0][0]
    # outputs is the second child of simulation
    outputsElement = simulationElement[1]
    # events is the third child of simultion
    #  the second child of events is the ablation event
    #  the third child of events is the apd (nrepol) event
    #  the fourth child is the resistance event
    ablnEventElement = simulationElement[2][1]
    apdEventElement = simulationElement[2][2]
    resEventElement = simulationElement[2][3]

    # update the iterations to the desired duration
    simulationElement.set('iterations', str(duration))

    # update the text of the input element to inputFile
    inputElement.text = inputFile

    # Add an output event to outputs
    Output = etree.SubElement(outputsElement, "output", type=outputType,
                              dir=outputDir, startIndex=str(start),
                              stopIndex=str(stop))

    # update the file subelement's text of the *EventElements to the files
    ablnEventElement[0].text = albntextfile
    apdEventElement[0].text = apdtextfile
    resEventElement[0].text = restextfile

    # Write to the vmFile
    tree = etree.ElementTree(root)
    tree.write(vmFile, pretty_print=True,
               xml_declaration=True, encoding="utf-8")
    return


def runtimeToVM_AblnOnly(runtimeFile, vmFile,
                         duration,
                         outputType, outputDir, start, stop,
                         albntextfile):
    # Read in the runtimeFile as an ElementTree
    with open(runtimeFile, 'r') as f:
        tree = etree.parse(f)
    root = tree.getroot()
    # The simulation element is the "grandchild" of the root
    simulationElement = root[0][0]
    # outputs is the second child of simulation
    outputsElement = simulationElement[1]
    # events is the third child of simultion
    #  the second child of events is the ablation event
    #  the third child of events is the apd (nrepol) event
    #  the fourth child is the resistance event
    ablnEventElement = simulationElement[2][1]

    # update the iterations to the desired duration
    simulationElement.set('iterations', str(duration))

    # Add an output event to outputs
    Output = etree.SubElement(outputsElement, "output", type=outputType,
                              dir=outputDir, startIndex=str(start),
                              stopIndex=str(stop))

    # update the file subelement's text of the *EventElements to the files
    ablnEventElement[0].text = albntextfile
    # Write to the vmFile
    tree = etree.ElementTree(root)
    tree.write(vmFile, pretty_print=True,
               xml_declaration=True, encoding="utf-8")
    return


def findRuntime(consoleLogFile, xmlfile, batchtoolpath):
    '''Runs the batch tool on the given xml file and uses
       getRuntimes to examine the consolLogFile to determine
       how long the tissue ran in the simulator
    '''
    with open(consoleLogFile, 'w') as f:
        with redirect_stdout(f):  # Redirect console outuput to consoleLogFile
            runcommand = [batchtoolpath + '/VisibleEP_BatchTool_1', xmlfile]
            runcommandString = runcommand[0] + ' ' + \
                               runcommand[1] + ' > ' + consoleLogFile
            try:
                run(runcommand, stdout=f, timeout=None, check=True)
            except:
                print('Unexpected error running \"' + runcommandString + '\"')
                raise

    return getRuntimes(consoleLogFile)


def getRuntimes(consoleLogFile):
    """
    Searches the given consoleLogFile to find the line for the two lines
    containing __Runtime summary__ and a: b and returns b if found.
    If the pattern isn't matched numpy.NaN is returned
    """
    with open(consoleLogFile, 'r') as file:
        nextContainsRuntime = False
        for line in file:
            if(('__Runtime summary__' in line) &
               (nextContainsRuntime is False)):
                nextContainsRuntime = True
            elif(nextContainsRuntime):
                if(':' in line):
                    splitLine = line.split(" ")
                    return int(splitLine[1][:-1])
                else:
                    nextContainsRuntime = False
    return np.NaN
