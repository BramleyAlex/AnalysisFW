import sys
import os
import ROOT as r

# Define the colors for the output
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
def DEBUG(text):
    return bcolors.OKGREEN+text+bcolors.ENDC
def TITLE(text):
    return bcolors.OKBLUE+text+bcolors.ENDC
def HEADER(text):
    return bcolors.HEADER+bcolors.BOLD+bcolors.OKBLUE+bcolors.UNDERLINE+text+bcolors.ENDC
def ERROR(text):
    return bcolors.FAIL+text+bcolors.ENDC


# Load the databases
from infofile import infos
from dataSets import dataSets, totRealLum, realList, dataCombos, dirs, processes
from DatasetsPaths import findMainPath

# For parallel processing and timing
import multiprocessing
from itertools import product
import time

# For parsing the script arguments
import argparse

# Load the C++ library
sys.path.append(findMainPath()+'/build/python/.')
from AnalysisFW import CLoopConfig, LM, MJ, HM, BDT, BDT_Cut

# This is a hack to allow the C++ class to be pickled, see:
# https://www.boost.org/doc/libs/1_84_0/libs/python/doc/html/reference/topics/pickle_support.html
def CLoopConfig_getinitargs(self):
  return (self.m_saveHistograms,self.m_saveEvents,self.m_reweightMjj,self.m_bdtWeightsPath, self.m_region, self.m_massRegion)
# now inject __getinitargs__ (Python is a dynamic language!)
CLoopConfig.__getinitargs__ = CLoopConfig_getinitargs


def getLuminosity(sampleName):
    if "2018" in sampleName:
        return 58.4501
    elif "2017" in sampleName:
        return 43.5873
    else :
        return 36.2369

def getSampleID(sampleName):
    z_sample=0
    if ("Zee" in sampleName) or ("Zmumu" in sampleName) or ("Ztautau" in sampleName):
        z_sample=1
        if "sherpa" in sampleName or "Sherpa" in sampleName:
            z_sample=2
            if "NLO" in sampleName:
                z_sample=4
        if "VBF" in sampleName:
            z_sample=0
        if "MG" in sampleName:
            z_sample=3
            if "NLO" in sampleName:
                z_sample=5
    return z_sample

def getNormalise(sampleName,lumFactor):
    if sampleName in realList:
        return 1.0
    else:
        weight = (lumFactor * 1000 * infos[sampleName]["xsec"] * infos[sampleName]["fil_eff"] * infos[sampleName]["kfac"])/infos[sampleName]["sumw"]
        return weight

def getAbsoluteFilePath(sampleName):
    correctPath = ""
    # Get the filename
    try:
        filename = dataSets[sampleName]
    except KeyError:
        print("Sample ", sampleName, " not defined. Exiting.")
        sys.exit(1)
    # search through several directories to find where the input file is located
    for path in dirs:
        if filename in os.listdir(path):
            correctPath = path
            break

    if correctPath == "":
        print("File for sample ", sampleName, " not found. Exiting.")
        sys.exit(1)

    return correctPath + filename

def getSamplesToRun(option,allData,allMC):
    if option=="All":
        for key,value in dataCombos.items():
            if "data" in key:
                allData+=value
            else:
                if "truth" not in key and "sys" not in key and "old" not in key and "VV_EW_Semi" not in key:
                    allMC+=value
        return
    elif option=="Data":
        for key,value in dataCombos.items():
            if "data" in key:
                allData+=value
        return
    elif option=="MC":
        for key,value in dataCombos.items():
            if "data" not in key:
                if "truth" not in key and "sys" not in key and "old" not in key and "VV_EW_Semi" not in key:
                    allMC+=value
        return
    
def createOutputDirectory(outputPath,treeName,verbosity):
    if not os.path.exists(outputPath):
        if verbosity=="DEBUG":
            print(DEBUG("Creating output directory at: "), os.path.abspath(outputPath))
        os.makedirs(outputPath)
    if not os.path.exists(outputPath+"/"+treeName):
        if verbosity=="DEBUG":
            print(DEBUG("Creating output for tree at: "), os.path.abspath(outputPath)+"/"+treeName)
        os.makedirs(outputPath+"/"+treeName)

def createConfigObject(jobTypeArgument,verbosity,region,massRegion):
    makeHistograms = 'h' in jobTypeArgument
    if verbosity=="DEBUG" and makeHistograms:
        print(DEBUG("Making histograms!"))
    makeNTuples = 'n' in jobTypeArgument
    if verbosity=="DEBUG" and makeNTuples:
        print(DEBUG("Making ntuples!"))
    makeReweighting = 'r' in jobTypeArgument
    if verbosity=="DEBUG" and makeReweighting:
        print(DEBUG("Making reweighting!"))
    mvaWeightsPath = findMainPath()+"/data/MVA-Weights/10Folds_BDT-0.4.6.weights.xml"
    if verbosity=="DEBUG":
        print(DEBUG("MVA weights path: "), mvaWeightsPath)
    if massRegion=="low":
        massRegion_int = 0
    elif massRegion=="medium":
        massRegion_int = 1
    elif massRegion=="high":
        massRegion_int = 2
    return CLoopConfig(makeHistograms,makeNTuples,makeReweighting,mvaWeightsPath,region,massRegion_int)

def str2bool(str):
    if isinstance(str, bool):
        return str
    if str.lower() in ("true", "yes", "t", "y", 1):
        return True
    elif str.lower() in ("false", "no", "f", "n", 0):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean Expected")

def createArgumentParser():
    # Parse the script arguments
    parser = argparse.ArgumentParser()
    executionMode = parser.add_mutually_exclusive_group()
    executionMode.add_argument("--samples", help="Group of samples to run over.",type=str,choices=["MC","Data","All"],default="All")
    executionMode.add_argument("--singleSample", help="Run over a single sample.",type=str,default="")
    executionMode.add_argument("--inputFile", help="Input txt file with a list of samples to run over.",type=str,default="")
    parser.add_argument("--verbosity", help="Verbosity level.",type=str,default="INFO",choices=["INFO","DEBUG"])
    parser.add_argument("--treeName", help="Name of the tree to run over.",type=str,default="T_s2thh_NOMINAL")
    parser.add_argument("--jobType", help="Type of job to run.",type=str,default="h",choices=["h","n","hn","hr","hnr","nr"])
    parser.add_argument("--outputDir", help="Path of to the directory used to store the processed samples.",type=str,default=findMainPath()+"/Results")
    parser.add_argument("--j", help="Number of cores to use.",type=int,default=1)
    parser.add_argument("--region", help="Which region to run over.",type=str,default="all",choices=["all","SR","CR","CRa","CRb","CRc"]) #Phase this out/overhaul
    parser.add_argument("--massRegion", help="low: <116, medium: 116-160, high: >160",type=str,default="low",choices=["low","medium","high"]) #Phase this out
    parser.add_argument("--analysisType",help="Type of physics analysis to run",type=str,default="LM",choices=["LM","MJ","HM","BDT","BDTCUT"])
    parser.add_argument("--histAdd",help="Add histograms into groups defined in dataSets.py. Only runs if using --samples",type=str2bool,nargs='?',default=False,const=True)
    return parser

def getArgumentTupleForSampleGroup(treeName,sampleGroup,verbosity,outputPath,analysisConfig):
    return product([treeName],sampleGroup,[verbosity],[outputPath],[analysisConfig])

def mergeFiles(treeName,verbosity,outputPath,cores):
    # Output Path
    path = outputPath+"/"+treeName+"/"
    for key, value in processes.items():
        file_output = key + ".root"
        file_input = ""
        if verbosity=="DEBUG":
            print(DEBUG("Creating file "+file_output))
        for i in value:
            file_input += (" "+ path + i + ".root")

        success = os.system(("hadd -f -j {0} "+outputPath+"/"+file_output+file_input).format(cores))
        if success != 0:
            print(ERROR("Error creating "+file_output))
            sys.exit(1)
    success = os.system("rm -r "+path)
    if success !=0:
        print(ERROR("Error deleting {0} folder".format(treeName)))
        sys.exit(1)

def runAnalysis(treeName,sampleName,verbosity,outputPath,analysisConfig):
    # Get the absolute path of the file
    filePath = getAbsoluteFilePath(sampleName)
    if verbosity=="DEBUG":
        print(TITLE("Sample name: "), sampleName)
        print(DEBUG("File path: "), filePath)

    # Open root file and get the tree
    file = r.TFile.Open(filePath, "READ")
    tree = r.TTree()
    tree = file.Get(treeName)

    # Get the sample ID
    sampleID = getSampleID(sampleName)
    if verbosity=="DEBUG":
        print(DEBUG("Sample ID: "), sampleID)
    # Get the luminosity
    lumFactor = getLuminosity(sampleName)
    if verbosity=="DEBUG":
        print(DEBUG("Luminosity: "), lumFactor)
    # Get the normalisation weight
    weight = getNormalise(sampleName,lumFactor)
    if verbosity=="DEBUG":
        print(DEBUG("Normalisation weight: "), weight)
    print("Tree address: {0}".format(tree))
    #analysis = CLoop(r.addressof(tree), sampleName)

    if args.analysisType == "LM": analysis = LM(r.addressof(tree), sampleName)
    elif args.analysisType == "MJ": analysis = MJ(r.addressof(tree), sampleName)
    elif args.analysisType == "HM": analysis = HM(r.addressof(tree), sampleName)
    elif args.analysisType == "BDT": analysis = BDT(r.addressof(tree), sampleName)
    elif args.analysisType == "BDTCUT": analysis = BDT_Cut(r.addressof(tree), sampleName)
    else:
        print(ERROR("Cloop derived class not found."))
        sys.exit(1)
    analysis.Loop(weight, sampleID, sampleName,analysisConfig)
    
    del analysis
    file.Close()
    success = os.system("mv "+sampleName+".root "+outputPath+"/"+treeName+"/"+sampleName+treeName+".root")
    if success != 0:
        print(ERROR("Error moving file to Results directory."))
        sys.exit(1)
    if verbosity=="DEBUG":
        print(DEBUG("==========================================================================================================================================="))

if __name__ == "__main__":
    # Parse the script arguments
    parser = createArgumentParser()
    args = parser.parse_args()

    # Start timer and greet the user
    initTime = time.time()
    print(HEADER("Welcome to the AnalysisFW"))

    # Define debug mode
    verbosity = args.verbosity

    # Define number of cores
    nCPU = args.j if verbosity=="INFO" else 1
    if verbosity=="DEBUG" and args.j>1:
        print(DEBUG("Only one core will be used when in DEBUG mode."))

    # Create the output directory
    createOutputDirectory(args.outputDir,args.treeName,verbosity)

    # Create the config object to pass to CLoop
    config = createConfigObject(args.jobType,verbosity,args.region,args.massRegion)

    # If a single sample is chosen, run just over that
    if args.singleSample != "":
        runAnalysis(args.treeName,args.singleSample,verbosity,args.outputDir,config)
    
    # If an input file is given, run over the samples in the file
    elif args.inputFile != "":
        print(TITLE("Running over samples in input file: "+args.inputFile))
        listOfSamples = []
        with open(args.inputFile) as f:
            for line in f:
                listOfSamples.append(line.strip())
        samplesTuple = getArgumentTupleForSampleGroup(args.treeName,listOfSamples,verbosity,args.outputDir, config)
        print(TITLE("Running over "+str(len(listOfSamples))+" samples\n"))
        with multiprocessing.Pool(processes=nCPU) as pool:
            pool.starmap(runAnalysis, samplesTuple)
    # If a group of samples is chosen, run over all of them
    elif args.samples != "":
        # Select the correct datasets
        allData = []
        allMC = []
        getSamplesToRun(args.samples,allData,allMC)
        dataTuple = getArgumentTupleForSampleGroup(args.treeName,allData,verbosity,args.outputDir, config)
        mcTuple = getArgumentTupleForSampleGroup(args.treeName,allMC,verbosity,args.outputDir, config)
        print(TITLE("Running over "+str(len(allData))+" DATA samples\n"))
        with multiprocessing.Pool(processes=nCPU) as pool:
            pool.starmap(runAnalysis, dataTuple)   
        print(TITLE("Running over "+str(len(allMC))+" MC samples\n"))
        with multiprocessing.Pool(processes=nCPU) as pool:
            pool.starmap(runAnalysis, mcTuple)
        # Merge files if true
        if args.histAdd:
            print(TITLE("Merging files"))
            mergeFiles(args.treeName,verbosity,args.outputDir,nCPU)
    # Say goodbye and print the time taken
    print(HEADER("Analysis done"))
    print(DEBUG("Time taken: "+str(time.time()-initTime)))
