import time
from classes.Dataset import *

'''
Function: createInitialDataSets
Argument list: 
Purpose: Create dictionary of datasets
Return types, values:
Dependencies:  h5py, Dataset Class
Creator: James Stauder
Date created: 1/31/18
Last edited: 1/31/18
'''


def createInitialDatasets():
    print "Creating data sets"
    t0 = time.time()

    datasetDict = {}

    dataFile = h5py.File(dataFileName, 'r')
    map['x1'] = len(dataFile['bed'][:][0])
    map['y1'] = len(dataFile['bed'][:])
    map['proj_x1'] = dataFile['x'][:][-1]
    map['proj_y1'] = dataFile['y'][:][-1]
    velocity = Dataset('velocity')
    datasetDict['velocity'] = velocity

    smb = Dataset('smb')
    datasetDict['smb'] = smb

    bed = Dataset('bed')
    datasetDict['bed'] = bed

    surface = Dataset('surface')
    datasetDict['surface'] = surface

    thickness = Dataset('thickness')
    datasetDict['thickness'] = thickness

    t2m = Dataset('t2m')
    datasetDict['t2m'] = t2m

    dataFile.close()

    print "Loaded all data sets in ", time.time() - t0, " seconds"
    return datasetDict
