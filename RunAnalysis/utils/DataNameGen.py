from os import listdir
import sys
import os
import ROOT as r
import csv

from MetadataGen import *

def getYear(directoryName):
  if "data15" in directoryName:
    return "_2015"
  elif "data16" in directoryName:
    return "_2016"
  elif "data17" in directoryName:
    return "_2017"
  elif "data18" in directoryName:
    return "_2018"

def getName(folderName):
  return folder_name[18:20] + "_" +folder_name[27]
   

def getDSID(rootFiles):
  DSID = -1
  #print(rootFiles)
  firstFile = r.TFile.Open(os.path.join(samplePath,rootFiles[0]))
  DSID = int(firstFile.Get("sumOfWeights").GetTitle())
  firstFile.Close()
  if DSID == -1:
      print("WARNING.. DSID not found in the root file ", rootFiles[0])
  return DSID
  
def checkInput():
    # First check that exactly two arguments are passed
    if len(sys.argv)!=2:
        print("Please provide one arguments: <pathToData>")
        print("<pathToData> is the root directory where all the MC and Data samples are stored.")
        sys.exit(1)

    pathToDatasets = sys.argv[1]

    # Check that the first is a directory
    if (not os.path.isdir(pathToDatasets)):
        print("The first argument ", pathToDatasets, " is not a directory. Exiting.")
        sys.exit(1)

if __name__ == "__main__":
  # Check that the script inputs are correct
  checkInput()
  pathToDatasets = sys.argv[1]

  samplesPaths = [os.path.join(pathToDatasets,f) for f in listdir(pathToDatasets) if os.path.isdir(os.path.join(pathToDatasets, f))]

  # Create the relevant dictionries to store metadata
  NameROOTFile_dict = {} # {"INDIVIDUAL_SAMPLE_NAME": "FILE_NAME.root"}
  Combo_dict = {} # {"GROUP_SAMPLE_NAME": [INDIVIDUAL_SAMPLE_NAME-1,INDIVIDUAL_SAMPLE_NAME-2,...] }


  # Loop over all the directories.
  for samplePath in samplesPaths:
    folder_name = samplePath.split("/")[-1]
    if "data" in folder_name:
      rootFiles = [f for f in listdir(samplePath) if f.endswith('.root')]
      if len(rootFiles) != 0:
        DSID = getDSID(rootFiles)
        #PERIOD_SUFFIX = getYear(samplePath)
        PERIOD_SUFFIX = getName(folder_name)
        #print("{0}  {1}".format(DSID, PERIOD_SUFFIX))
        NameROOTFile_dict.update(buildSampleNameROOTFileDict(DSID,rootFiles,PERIOD_SUFFIX))
        Combo_dict.update(buildSamplesComboDict(DSID,rootFiles,PERIOD_SUFFIX))
      


    
  outputFile = open('rootFileNamesDATA.txt','a+')
  for key, value in NameROOTFile_dict.items():
    outputFile.write(f'"{key}": "{value}",\n')
  outputFile.close()

  outputFile = open('all_data.txt', "a+")
  outputFile.write("[")
  for key, value in NameROOTFile_dict.items():
    outputFile.write(f'"{key}",')
  outputFile.write("]")
  outputFile.close()

  outputFile = open('samplesComboDATA.txt','a+')
  for key, value in Combo_dict.items():
    outputFile.write(f'"{key}": {value},\n')
  outputFile.close()