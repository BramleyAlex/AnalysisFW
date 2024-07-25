# This script should be used as:
# python3 MetadataGen.py <pathToData> <pathToMetadataFile>
# <pathToData> is the absolute path to the root directory where all the MC and Data samples are stored.
# <pathToMetadataFile> is the absolute path to the metadata file that contains the DSIDs and the cross-sections, etc.

from os import listdir
import sys
import os
import ROOT as r
import csv

# Dictionary to map DSIDs into user-defined sample names
sampleNameDSID_dict = {
700360: 'VBF_Ztautau',
410649: 'Wt_DR_dilepton_atop', # single top
410648: 'Wt_DR_dilepton_top', # single top
410645: 'st_schan_atop', # single top
410644: 'st_schan_top', # single top
410658: 'st_tchan_top', # single top
410472: 'ttbar_dil', # tt bar
410470: 'ttbar_nonallhad', # tt bar
345123: 'tautauh30h20', 
345120: 'tautaul13l7',  #Higgs
345121: 'ggHttlm15hp20', #Higgs
364128: 'Ztautau_sherpa1',
364129: 'Ztautau_sherpa2',
364130: 'Ztautau_sherpa3',
364131: 'Ztautau_sherpa4',
364132: 'Ztautau_sherpa5',
364133: 'Ztautau_sherpa6',
364134: 'Ztautau_sherpa7',
364135: 'Ztautau_sherpa8',
364136: 'Ztautau_sherpa9',
364137: 'Ztautau_sherpa10',
364138: 'Ztautau_sherpa11',
364139: 'Ztautau_sherpa12',
364140: 'Ztautau_sherpa13',
364141: 'Ztautau_sherpa14',
700338: 'Sh_2211_Wenu_maxHTpTV2_BFilter', # W_Sherpa - W+jets?
700339: 'Sh_2211_Wenu_maxHTpTV2_CFilterBVeto', # W_Sherpa
700340: 'Sh_2211_Wenu_maxHTpTV2_CVetoBVeto', # W_Sherpa
700488: 'VV_QCD1',
700490: 'VV_QCD2',
700489: 'VV_QCD3',
700601: 'VV_QCD4',
700602: 'VV_QCD5',
700603: 'VV_QCD6',
700600: 'VV_QCD7', 
700604: 'VV_QCD8',
700492: 'VV_QCD9',
700494: 'VV_QCD10',
700493: 'VV_QCD11',
700363: 'W_EWK_sherpa1',
700362: 'W_EWK_sherpa2',
700364: 'W_EWK_sherpa3',
700341: 'Sh_2211_Wmunu_maxHTpTV2_BFilter', # W_Sherpa
700342: 'Sh_2211_Wmunu_maxHTpTV2_CFilterBVeto', # W_Sherpa
700343: 'Sh_2211_Wmunu_maxHTpTV2_CVetoBVeto', # W_Sherpa
700347: 'Sh_2211_Wtaunu_H_maxHTpTV2_BFilter', # W_Sherpa
700348: 'Sh_2211_Wtaunu_H_maxHTpTV2_CFilterBVeto', # W_Sherpa
700349: 'Sh_2211_Wtaunu_H_maxHTpTV2_CVetoBVeto', # W_Sherpa
700344: 'Sh_2211_Wtaunu_L_maxHTpTV2_BFilter', # W_Sherpa
700345: 'Sh_2211_Wtaunu_L_maxHTpTV2_CFilterBVeto', # W_Sherpa
700346: 'Sh_2211_Wtaunu_L_maxHTpTV2_CVetoBVeto', # W_Sherpa
700358: 'VBF_Zee_sherpa',
700359: 'VBF_Zmumu_sherpa',
700361: 'Sh_2211_Znunu2jets_Min_N_TChannel', # Z_EWK
700592: 'VV_EWK1',
700594: 'VV_EWK2',
700589: 'VV_EWK3',
700587: 'VV_EWK4',
700593: 'VV_EWK5',
700590: 'VV_EWK6',
700591: 'VV_EWK7',
700588: 'VV_EWK8',
410659: 'st_tchan_atop', #single top
345122: 'ggHttlp15hm20', #Higgs
1: 'data',
}

# THIS FUNCTION MIGHT CHANGE FROM ONE NTUPLE FORMAT TO ANOTHER
def getDSID(rootFiles):
    DSID = -1
    firstFile = r.TFile.Open(os.path.join(samplePath,rootFiles[0]))
    DSID = int(firstFile.Get("sumOfWeights").GetTitle())
    firstFile.Close()
    if DSID == -1:
        print("WARNING.. DSID not found in the root file ", rootFiles[0])
    return DSID

def getPeriodSuffix(directoryName):
    if "MC16a" in directoryName:
        return "_2015"
    elif "MC16d" in directoryName:
        return "_2017"
    elif "MC16e" in directoryName:
        return "_2018"

# This function assumes that all the root files in the list have the same DSID
def buildSampleNameROOTFileDict(DSID,rootFilesList,periodSuffix):
    sampleNameROOTFile_dict = {}
    sampleName = sampleNameDSID_dict[DSID]+periodSuffix
    counter = 0
    for rootFile in rootFilesList:
        sampleNameROOTFile_dict[sampleName+'_'+str(counter)] = rootFile
        counter += 1
    return sampleNameROOTFile_dict

def buildSamplesComboDict(DSID,rootFilesList,periodSuffix):
    samplesCombo_dict = {}
    sampleGroupName = sampleNameDSID_dict[DSID]+periodSuffix
    samplesCombo_dict[sampleGroupName] = [sampleGroupName+'_'+str(i) for i in range(len(rootFilesList))]
    return samplesCombo_dict

def getMetadataFromFile(pathToMetadataFile,DSID,sumWeights):
    metadata_dict = {'DSID' : DSID ,'events' : sumWeights,'red_eff' : 1,'sumw' : sumWeights,'xsec' : 0.0,'kfac' : 0.0,'fil_eff' :0.0}
    # Look for the DSID and extract the metadata
    with open(pathToMetadataFile,"r") as csv_file:
        csv_reader=csv.reader(csv_file,delimiter=',')
        for row in csv_reader:
            if str(DSID) in row: # DSID is the first element of the row and is a string
                metadata_dict['xsec']=float(row[2])
                metadata_dict['kfac']=float(row[3])
                metadata_dict['fil_eff']=float(row[4])
                break
    return metadata_dict

def getSumOfWeights(rootFiles,samplePath):
    sumWeights = 0
    for rootFile in rootFiles:
        file = r.TFile.Open(os.path.join(samplePath,rootFile))
        sumWeights += file.Get("sumOfWeights").GetBinContent(4)
        file.Close()
    return sumWeights

def buildSampleNameMetadataDict(DSID,periodSuffix,rootFilesList,metadata_dict):
    nameMetadata_dict = {}
    sampleGroupName = sampleNameDSID_dict[DSID]+periodSuffix
    for i in range(len(rootFilesList)):
        nameMetadata_dict[sampleGroupName+'_'+str(i)] = metadata_dict
    return nameMetadata_dict

def buildFolderNameDict(DSID,samplePath,periodSuffix):
    folderName_dict = {}
    folder = samplePath.split("/")
    sampleGroupName = sampleNameDSID_dict[DSID]+periodSuffix
    folderName_dict[sampleGroupName] = folder[-1]
    return folderName_dict

def checkInputs():
    # First check that exactly two arguments are passed
    if len(sys.argv)!=3:
        print("Please provide two arguments: <pathToData> <pathToMetadataFile>")
        print("<pathToData> is the root directory where all the MC and Data samples are stored.")
        print("<pathToMetadataFile> is the absolute path to the metadata file that contains the DSIDs and the cross-sections, etc.")
        sys.exit(1)

    pathToDatasets = sys.argv[1]
    pathToMetadata = sys.argv[2]

    # Check that the first is a directory and the second is a file
    if (not os.path.isdir(pathToDatasets)):
        print("The first argument ", pathToDatasets, " is not a directory. Exiting.")
        sys.exit(1)
    if (not os.path.isfile(pathToMetadata)):
        print("The second argument ", pathToMetadata, " is not a file. Exiting.")
        sys.exit(1)

if __name__ == "__main__":
    # Check that the script inputs are correct
    checkInputs()
    pathToDatasets = sys.argv[1]
    pathToMetadata = sys.argv[2]

    samplesPaths = [os.path.join(pathToDatasets,f) for f in listdir(pathToDatasets) if os.path.isdir(os.path.join(pathToDatasets, f))]

    # Create the relevant dictionries to store metadata
    sampleNameROOTFile_dict = {} # {"INDIVIDUAL_SAMPLE_NAME": "FILE_NAME.root"}
    samplesCombo_dict = {} # {"GROUP_SAMPLE_NAME": [INDIVIDUAL_SAMPLE_NAME-1,INDIVIDUAL_SAMPLE_NAME-2,...] }
    sampleNameMetadata_dic = {} # {"INDIVIDUAL_SAMPLE_NAME": {'DSID' : val ,'events' : val,'red_eff' : 1,'sumw' : val,'xsec' : val,'kfac' : val,'fil_eff' :val} }
    folderName_dict = {}

    # Loop over all the directories.
    for samplePath in samplesPaths:
        rootFiles = [f for f in listdir(samplePath) if f.endswith('.root')]

        if len(rootFiles) != 0:
            DSID = getDSID(rootFiles)
            PERIOD_SUFFIX = getPeriodSuffix(samplePath)
            SUM_WEIGHTS = getSumOfWeights(rootFiles,samplePath)
            METADATA_DIC = getMetadataFromFile(pathToMetadata,DSID,SUM_WEIGHTS)
            # Build the objects of interest
            sampleNameROOTFile_dict.update(buildSampleNameROOTFileDict(DSID,rootFiles,PERIOD_SUFFIX))
            samplesCombo_dict.update(buildSamplesComboDict(DSID,rootFiles,PERIOD_SUFFIX))
            sampleNameMetadata_dic.update(buildSampleNameMetadataDict(DSID,PERIOD_SUFFIX,rootFiles,METADATA_DIC))
            folderName_dict.update(buildFolderNameDict(DSID,samplePath,PERIOD_SUFFIX))

    outputFile = open('rootFileNames.txt','a+')
    for key, value in sampleNameROOTFile_dict.items():
        outputFile.write(f'"{key}": "{value}",\n')
    outputFile.close()

    outputFile = open('samplesCombo.txt','a+')
    for key, value in samplesCombo_dict.items():
        outputFile.write(f'"{key}": {value},\n')
    outputFile.close()

    outputFile = open('metadata.txt','a+')
    for key, value in sampleNameMetadata_dic.items():
        outputFile.write(f'"{key}": {value},\n')
    outputFile.close()

    outputFile = open('folderName.txt', 'a+')
    for key, value in folderName_dict.items():
        outputFile.write(f'"{key}": "{value}",\n')
    outputFile.close()
    


