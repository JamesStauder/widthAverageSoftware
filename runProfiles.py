import sys
import glob, os
import numpy as np
import h5py
import time 
from matplotlib.path import Path
from helperFiles.classes.FlowIntegrator import *
import matplotlib.pyplot as plt
import math

def main(argv):
    allDataFile = h5py.File('AllDataSets0.h5', 'r')
    smallerDataFile = h5py.File('GreenlandInBedCoord_V2.h5', 'r')
    os.chdir(argv[1])
    myFiles = []

    
    for file in glob.glob("*.h5"):
        myFiles.append(file)
        
    
    t0 = time.time()
    print "starting VX"
    datasetDict = {}
    datasetDict['VX'] = Dataset('VX', allDataFile)
    print "VX",  time.time() - t0
    print "starting VY"
    t0 = time.time()
    datasetDict['VY'] = Dataset('VY', allDataFile)
    print "VY",  time.time() - t0
    print "starting thickness"
    t0 = time.time()
    datasetDict['thickness'] = Dataset('thickness', allDataFile)
    print "thickness",  time.time() - t0
    t0 = time.time()
    
    smallerDataFile.close()
    
    integrator = createIntegrator(datasetDict['VX'], datasetDict['VY'])
    
    
    for file in myFiles:

        tempFile = h5py.File(file, 'r')
        shear1X = tempFile['Shear_0_x'][:]
        shear1Y = tempFile['Shear_0_y'][:]
        shear2X = tempFile['Shear_1_x'][:]
        shear2Y = tempFile['Shear_1_y'][:]
        

        myShears = [list(zip(shear1X, shear1Y)), list(zip(shear2X,shear2Y))]
        
        
        tempFile.close()
        numTrials = range(1,3)
        widthList = getWidths(myShears)
        
        method1TotalValues = []
        method2TotalValues = []
        method1StepValues = []
        method2StepValues = []
        
        for x in numTrials:

            method1TotalValue, method1StepValue, method2TotalValue, method2StepValue = 0,0,0,0
            
            try:
                method1TotalValue, method1StepValue = method1(myShears,x,widthList, 1000, datasetDict, integrator)
            except:
                e = sys.exc_info()[0]
                print "error on method1 file: ", file, "\nerror trial: ", x, "\nerror: ", e
            method1TotalValues.append(method1TotalValue)
            method1StepValues.append(method1StepValue)
            
            try:
                method2TotalValue, method2StepValue = method2(myShears,x,widthList, 1000, datasetDict)
            except:
                e = sys.exc_info()[0]
                print "error on method2 file: ", file, "\nerror trial: ", x, "\nerror: ", e
            method2TotalValues.append(method2TotalValue)
            method2StepValues.append(method2StepValue)
            

        try:
            method3Value = method3(myShears, allDataFile)
        except:
            e = sys.exc_info()[0]
            print "error on method3 file: ", file,"\nerror: ", e
        
        
        outFileName = "widthAverage" + file
        outDataFile = h5py.File(outFileName, 'w')
        
        outDataFile.create_dataset("method1TotalVolume", data=method1TotalValues)
        outDataFile.create_dataset("method1StepVolume", data=method1StepValues)
        outDataFile.create_dataset("method2TotalVolume", data=method2TotalValues)
        outDataFile.create_dataset("method2StepVolume", data=method2StepValues)
        
        
        method3ValueArray = []
        for i in range(0, len(method1TotalValues)):
            method3ValueArray.append(method3Value)
        outDataFile.create_dataset("method3Value", data=method3ValueArray)
        
        
        
        outDataFile.close()
        
    
def getWidths(shears):
    shear1 = np.asarray(shears[0])
    shear2 = np.asarray(shears[1])
    width = sqrt((shear1[:,0] - shear2[:,0])**2 + (shear1[:,1] - shear2[:,1])**2)
    return np.asarray(width)

def createIntegrator(vx,vy):
    return FlowIntegrator(vx,vy)
    
def method1(shears, numSamples, widthList, resolution,datasetDict, integrator):
    shear1 = np.asarray(shears[0])
    shear2 = np.asarray(shears[1])
    
    dx = (shear2[0][0] - shear1[0][0])/(numSamples+1)
    dy = (shear2[0][1] - shear1[0][1])/(numSamples+1)
    
    currX = shear1[0][0]
    currY = shear1[0][1]
    
    flowlineSamples = []
    for _ in range(0, numSamples):
        currX = currX + dx
        currY = currY + dy
        
        newLine = []
        
        for i in range(len(shear1)):
            newLine.append(None)
        newLine[0] = [currX, currY]
        a = integrator.integrate(currX, currY, newLine, 0, resolution)

        
        if None not in newLine:
            flowlineSamples.append(newLine)
        else:
            print "integration error on averaging method. Ommitting line"
        
    print "Ommitted ", len(flowlineSamples) - numSamples, " lines. Out of a possible ", numSamples
    
    thicknessList = []
    for step in range(len(flowlineSamples[0])):
        tempThickness = []
        for sample in range(len(flowlineSamples)):
            samplePoint = flowlineSamples[sample][step]
            tempThickness.append(datasetDict['thickness'].getInterpolatedValue(samplePoint[0], samplePoint[1]))
        thicknessList.append(np.mean(np.asarray(tempThickness)))
    thicknessList = np.asarray(thicknessList)
    
    volume = (thicknessList[:] * widthList[:]) * resolution
    
    return sum(volume), volume
        
        
        
    
    

def method2(shears, numSamples, widthList, resolution, datasetDict):

    shear1 = np.asarray(shears[0])
    shear2 = np.asarray(shears[1])
    
    avgThicknessList = []
    
    for x in range(0, len(shear1)):
        totalThickness = 0
        shearPoints = [[shear1[x][0], shear1[x][1]], [shear2[x][0], shear2[x][1]]]
        
        dx = (shearPoints[1][0] - shearPoints[0][0]) / (numSamples + 1)
        dy = (shearPoints[1][1] - shearPoints[0][1]) / (numSamples + 1)
        
        currPoint = shearPoints[0]
        
        for y in range(numSamples):
            currPoint[0] = currPoint[0] + dx
            currPoint[1] = currPoint[1] + dy
            
            totalThickness = totalThickness + datasetDict['thickness'].getInterpolatedValue(currPoint[0], currPoint[1])
            
        avgThicknessList.append(totalThickness/numSamples)
    
    volume = ((avgThicknessList * widthList) * resolution)
    
    return sum(volume), volume
        
def method3(shears, allDataFile):

    shear1 = np.asarray(shears[0])
    shear2 = np.asarray(shears[1])
    
    myShears = []
    myShears.extend(shear1)
    myShears.extend(list(reversed(shear2)))
    myShears.append(shear1[0])
    myShears = np.asarray(myShears)
    
    myPath = Path(myShears, closed=True)
    
    maxY, minY = max(myShears[:,1]), min(myShears[:,1])
    minY = minY - 150  #buffer
    maxY = maxY + 150 #buffer

    maxX, minX = max(myShears[:,0]), min(myShears[:,0])
    maxX = maxX + 150
    minX = minX - 150

    
    xs = allDataFile['x'][:]
    ys = allDataFile['y'][:]

    thickness = allDataFile['thickness'][:]
    
    allDataFile.close()

    xx, yy = np.meshgrid(xs, ys)
    XY = np.dstack((xx,yy))
    XYFlat = XY.reshape(-1, 2)

    thicknessFlat = thickness.reshape(-1)





    dataMaxY = XYFlat[0][1]
    dataMinY = XYFlat[-1][1]


    YToCutMax = abs(dataMaxY - maxY)/150
    YToCutMin = abs(dataMinY - minY)/150

    numX = 10018
    samplePoints = XYFlat[int(math.floor(YToCutMax)*numX):-(int(math.floor(YToCutMin) * numX))]
    thicknessFlat = thicknessFlat[int(math.floor(YToCutMax)*numX):-(int(math.floor(YToCutMin) * numX))]



    newSamplePoints = []
    newThickness = []
    for i in range(numX):
        newSamplePoints.extend(samplePoints[i::numX])
        newThickness.extend(thicknessFlat[i::numX])

    newSamplePoints = np.asarray(newSamplePoints)
    newThickness = np.asarray(newThickness)





    dataMaxX = newSamplePoints[-1][0]
    dataMinX = newSamplePoints[0][0]


    XToCutMax = int(math.floor(abs(dataMaxX - maxX)/150))
    XToCutMin = int(math.floor(abs(dataMinX - minX)/150))

    newMaxY, newMinY = newSamplePoints[0][1], newSamplePoints[-1][1]

    numY = int(abs(newMaxY - newMinY)/150) +1










    finalSamplePoints = newSamplePoints[XToCutMin*numY:-(XToCutMax*numY)]
    finalThicknessPoints = newThickness[XToCutMin*numY:-(XToCutMax*numY)]




    thicknessValuesInShape = []



    pointsInShape = []

    time0 = time.time()
    for i in range(len(finalSamplePoints)):

        point = [int(finalSamplePoints[i][0]), int(finalSamplePoints[i][1])]

        if myPath.contains_point(point):
            thicknessValuesInShape.append(finalThicknessPoints[i])
            pointsInShape.append(point)


    trueVolume = sum(thicknessValuesInShape) * 150**2

    return trueVolume

#pulled from software and modified. 
class Dataset:
    def __init__(self, name, dataFile):
        self.name = name
        self.dataFile = dataFile
        
        
        map = {'x0': 0, 'y0': 0,  'x1': 10018, 'y1': 17946,
       'proj_x0': -638000, 'proj_x1': 864550,
       'proj_y0': -657600, 'proj_y1': -3349350}
       
        map['x1'] = len(dataFile['bed'][:][0])
        map['y1'] = len(dataFile['bed'][:])
        map['proj_x1'] = dataFile['x'][:][-1]
        map['proj_y1'] = dataFile['y'][:][-1]
       

        bed_xarray = linspace(map['proj_x0'], map['proj_x1'], map['x1'], endpoint=True)
        bed_yarray = linspace(map['proj_y1'], map['proj_y0'], map['y1'], endpoint=True)

        if self.name == 'velocity':
            self.data, self.vx, self.vy = self.setData(name)

            self.vxInterp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(self.vx).transpose())
            self.vyInterp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(self.vy).transpose())

        else:
            self.data = self.setData(name)

        self.interp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(self.data).transpose())

    def setData(self, name):
        if name == 'velocity':
            vx = self.dataFile['VX'][:]
            vy = self.dataFile['VY'][:]
            data = sqrt(vx ** 2 + vy ** 2)
            return data, vx, vy
        else:
            data = self.dataFile[name][:]
            return data

    def getInterpolatedValue(self, xPosition, yPosition):
        return self.interp(xPosition, yPosition)[0][0]






if __name__ == '__main__':
    main(sys.argv)