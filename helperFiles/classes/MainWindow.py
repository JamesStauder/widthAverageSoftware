from PyQt4.QtGui import *
import math
from Instructions import *
from FlowIntegrator import *
from Dataset import *
from Marker import *
from ModelGUI import *
from ..caching.cachingFunctions import *
import time
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as mpatches
import matplotlib.patches as patches
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
from decimal import Decimal

'''
Class: MainWindow
Argument list:
Purpose: Create main window for GUI. Has many functions based on what user does
Return types, values:
Dependencies: pyQT, dolfin, math
Creator: James Stauder
Date created: 1/31/18
Last edited: 5/29/18
'''


class MainWindow(QMainWindow):
    def __init__(self, parent=None):

        super(MainWindow, self).__init__(parent)

        self.setWindowTitle("Greenland")
        self.setMinimumHeight(1000)
        self.setMinimumWidth(1200)

        self.centralWidget = QtGui.QWidget()
        self.setCentralWidget(self.centralWidget)
        self.mainLayout = QtGui.QHBoxLayout()
        self.centralWidget.setLayout(self.mainLayout)

        # index of current map
        self.currentMap = 0

        # marker selected variables
        self.isMarkerSelected = False
        self.whichMarkerSelected = None
        self.selectedMarkerPosition = None
        self.whichIndexOfFlowlineSelected = None

        # Flowline information
        self.flowlineDistance = 100000
        self.lengthOfFlowline = 1
        self.flowlines = []
        self.flowlineMarkers = []
        self.integratorPerMarker = 10

        '''
        Side widget with button
        '''
        self.maxWidth = 300

        self.buttonBoxWidget = QtGui.QWidget()
        self.buttonBox = QtGui.QVBoxLayout()
        self.buttonBoxWidget.setLayout(self.buttonBox)

        self.mapList = QtGui.QComboBox()
        self.maps = ['Velocity', 'Bed', 'Surface', 'SMB', 'Thickness', 't2m']
        self.mapList.addItems(self.maps)
        self.mapList.setMaximumWidth(self.maxWidth)
        self.buttonBox.addWidget(self.mapList)

        self.spatialResolutionWidget = QtGui.QWidget()
        self.spatialResolutionLayout = QtGui.QHBoxLayout()
        self.spatialResolutionWidget.setLayout(self.spatialResolutionLayout)
        self.spatialResolutionLabel = QtGui.QLabel('Spatial Resolution(m)')
        self.spatialResolutionLineEdit = QtGui.QLineEdit('1000')
        self.spatialResolutionLayout.addWidget(self.spatialResolutionLabel)
        self.spatialResolutionLayout.addWidget(self.spatialResolutionLineEdit)
        self.buttonBox.addWidget(self.spatialResolutionWidget)

        self.distanceWidget = QtGui.QWidget()
        self.distanceLayout = QtGui.QHBoxLayout()
        self.distanceWidget.setLayout(self.distanceLayout)
        self.distanceLabel = QtGui.QLabel('distance(km)')
        self.distanceLineEdit = QtGui.QLineEdit('100')
        self.spatialResolutionLayout.addWidget(self.distanceLabel)
        self.spatialResolutionLayout.addWidget(self.distanceLineEdit)
        self.buttonBox.addWidget(self.distanceWidget)

        self.averageWidget = QtGui.QWidget()
        self.averageLayout = QtGui.QHBoxLayout()
        self.averageWidget.setLayout(self.averageLayout)
        self.widthAverageButton = QtGui.QCheckBox('Use Width Average')
        self.widthAverageButton.setTristate(False)
        self.averageLayout.addWidget(self.widthAverageButton)
        self.buttonBox.addWidget(self.averageWidget)

        self.profileWidget = QtGui.QWidget()
        self.profileLayout = QtGui.QHBoxLayout()
        self.profileWidget.setLayout(self.profileLayout)
        self.profileLabel = QtGui.QLabel('output file name')
        self.profileLineEdit = QtGui.QLineEdit('myProfile.h5')
        self.profileLayout.addWidget(self.profileLabel)
        self.profileLayout.addWidget(self.profileLineEdit)
        self.buttonBox.addWidget(self.profileWidget)

        self.instructionButton = QtGui.QPushButton('Instructions')
        self.instructionButton.setEnabled(True)
        self.instructionButton.setMaximumWidth(self.maxWidth)
        self.buttonBox.addWidget(self.instructionButton)

        self.plotPathButton = QtGui.QPushButton('Plot Path')
        self.plotPathButton.setEnabled(False)
        self.plotPathButton.setMaximumWidth(self.maxWidth)
        self.buttonBox.addWidget(self.plotPathButton)

        self.runModelButton = QtGui.QPushButton('Run Model')
        self.runModelButton.setEnabled(False)
        self.runModelButton.setMaximumWidth(self.maxWidth)
        self.buttonBox.addWidget(self.runModelButton)

        self.resetButton = QtGui.QPushButton('Reset')
        self.resetButton.setEnabled(True)
        self.resetButton.setMaximumWidth(self.maxWidth)
        self.buttonBox.addWidget(self.resetButton)

        self.velocityWidthButton = QtGui.QPushButton('Create Profile')
        self.velocityWidthButton.setEnabled(False)
        self.velocityWidthButton.setMaximumWidth(self.maxWidth)
        self.buttonBox.addWidget(self.velocityWidthButton)

        self.textOut = QtGui.QTextBrowser()
        self.textOut.setMaximumWidth(self.maxWidth)
        self.buttonBox.addWidget(self.textOut)

        self.leftSideWidget = QtGui.QWidget()
        self.leftSide = QtGui.QVBoxLayout()
        self.leftSideWidget.setLayout(self.leftSide)

        self.imageItemContainer = QtGui.QStackedWidget()

        self.leftSide.addWidget(self.imageItemContainer)

        self.mainLayout.addWidget(self.leftSideWidget)
        self.mainLayout.addWidget(self.buttonBoxWidget)

        self.buttonBoxWidget.setMaximumWidth(self.maxWidth + 12)

        self.connectButtons()

    '''
    Function: addToImageItemContainer
    Argument list: datasetDict
    Purpose: add the different dataset widgets to the imageItemContainer
    Return types, values: 
    Dependencies: 
    Creator: James Stauder
    Date created: 2/25/18
    Last edited: 3/5/18
    '''

    def addToImageItemContainer(self, datasetDict):
        self.imageItemContainer.addWidget(datasetDict['velocity'].plotWidget)
        self.imageItemContainer.setCurrentWidget(datasetDict['velocity'].plotWidget)

        self.imageItemContainer.addWidget(datasetDict['bed'].plotWidget)
        self.imageItemContainer.addWidget(datasetDict['surface'].plotWidget)
        self.imageItemContainer.addWidget(datasetDict['thickness'].plotWidget)
        self.imageItemContainer.addWidget(datasetDict['t2m'].plotWidget)
        self.imageItemContainer.addWidget(datasetDict['smb'].plotWidget)

    '''
    Function: changeMap
    Argument list: 
        index: index of which map to use
    Purpose: Changes the map to a different colormap
    Return types, values: 
    Dependencies: 
    Creator: James Stauder
    Date created: 2/25/18
    Last edited: 3/5/18
    '''

    def changeMap(self, index):
        vr = self.imageItemContainer.currentWidget().getPlotItem().getViewBox().viewRange()
        indexToDatasetDict = {
            0: 'velocity',
            1: 'bed',
            2: 'surface',
            3: 'smb',
            4: 'thickness',
            5: 't2m'}
        if index != self.currentMap:
            oldMap = self.currentMap
            self.currentMap = index

        self.imageItemContainer.setCurrentWidget(self.datasetDict[indexToDatasetDict[self.currentMap]].plotWidget)
        self.datasetDict[indexToDatasetDict[self.currentMap]].imageItem.hoverEvent = self.mouseMove
        self.datasetDict[indexToDatasetDict[self.currentMap]].imageItem.mouseClickEvent = self.mouseClick

        self.datasetDict[indexToDatasetDict[self.currentMap]].plotWidget.getPlotItem().getViewBox().setRange(
            xRange=vr[0],
            yRange=vr[1],
            padding=0.0)
        for line in self.flowlineMarkers:
            for marker in line:
                marker.plotWidget = self.datasetDict[indexToDatasetDict[self.currentMap]]
                self.datasetDict[indexToDatasetDict[oldMap]].plotWidget.removeItem(marker.cross[0])
                self.datasetDict[indexToDatasetDict[oldMap]].plotWidget.removeItem(marker.cross[1])
                self.datasetDict[indexToDatasetDict[self.currentMap]].plotWidget.addItem(marker.cross[0])
                self.datasetDict[indexToDatasetDict[self.currentMap]].plotWidget.addItem(marker.cross[1])

                if marker.lines[0]:
                    self.datasetDict[indexToDatasetDict[self.currentMap]].plotWidget.addItem(marker.lines[0])

    '''
     Function: mouseClick
     Argument list: 
        e: event trigger from mouse being clicked
     Purpose: 
        Create a new flowline or move a previous flowline
     Return types, values: None
     Dependencies: None
     Creator: James Stauder
     Date created: 2/25/18
     Last edited: 5/29/18
     '''

    def mouseClick(self, e):

        # If no marker is selected
        if self.isMarkerSelected is False:

            # Check to see if click selects a marker. If so memoize the marker and the flowline Position
            # This block checks every marker along flowline and reintegrates up. Commented out with new
            # averaging method
            '''
            for i in range(len(self.flowlineMarkers)):
                for j in range(len(self.flowlineMarkers[i])):
                    if self.flowlineMarkers[i][j].checkClicked(e.pos()):
                        self.isMarkerSelected = True
                        self.whichMarkerSelected = self.flowlineMarkers[i][j]
                        self.selectedMarkerPosition = [i, j]

                        self.displayMarkerVariables()
                        tempX, tempY = self.whichMarkerSelected.dx, self.whichMarkerSelected.dy
                        for k in range(len(self.flowlines[i])):
                            if self.flowlines[i][k] == [tempX, tempY]:
                                self.whichIndexOfFlowlineSelected = [i, k]
                        break
            '''

            # Checks to see only if first marker in each flowline is detected.
            for i in range(len(self.flowlineMarkers)):
                if self.flowlineMarkers[i][0].checkClicked(e.pos()):
                    self.isMarkerSelected = True
                    self.whichMarkerSelected = self.flowlineMarkers[i][0]
                    self.selectedMarkerPosition = [i, 0]

                    self.displayMarkerVariables()
                    self.whichIndexOfFlowlineSelected = [i, 0]
                    break

            # If no marker selected previously or currently create new flowline. Also cannot create more
            # then 2 flowlines.
            if (len(self.flowlines) < 2) and self.isMarkerSelected is False:
                self.spatialResolutionLineEdit.setReadOnly(True)
                self.distanceLineEdit.setReadOnly(True)
                self.flowlineDistance = int(self.distanceLineEdit.text()) * 1000
                self.lengthOfFlowline = int(self.flowlineDistance / float(self.spatialResolutionLineEdit.text()))
                self.integratorPerMarker = int(math.ceil(10000 / (float(self.spatialResolutionLineEdit.text()))))
                xClickPosition = e.pos().x()
                yClickPosition = e.pos().y()

                dx, dy = colorToProj(xClickPosition, yClickPosition)

                # Create new flowline
                newFlowline = []
                for x in range(0, self.lengthOfFlowline):
                    newFlowline.append(None)
                newFlowline[0] = [dx, dy]

                newFlowline = self.flowIntegrator.integrate(dx, dy, newFlowline, 0,
                                                            float(self.spatialResolutionLineEdit.text()))

                if None in newFlowline:
                    print "Integration Error. Try Again"
                    return

                # Create a flowline of markers spaced out based on the IntegratorPerMarker
                newFlowlineMarkers = newFlowline[::self.integratorPerMarker]

                for i in range(len(newFlowlineMarkers)):
                    dx = newFlowlineMarkers[i][0]
                    dy = newFlowlineMarkers[i][1]
                    cx, cy = colorCoord(dx, dy)
                    newFlowlineMarkers[i] = Marker(cx, cy, dx, dy, self.imageItemContainer.currentWidget())

                self.displayMarkers(newFlowlineMarkers)

                self.flowlines.append(newFlowline)
                self.flowlineMarkers.append(newFlowlineMarkers)

                if len(self.flowlines) == 2:
                    self.velocityWidthButton.setEnabled(True)

        # Release the marker that was previously held
        else:
            self.isMarkerSelected = False
            self.whichMarkerSelected = None
            self.textOut.clear()

    '''
    Function: mouseMove
    Argument list: 
    Purpose: This function is used to move the marker that is selected and create a new integration path. 
    Return types, values: 
    Dependencies: 
    Creator: James Stauder
    Date created: 2/25/18
    Last edited: 3/9/18
    TODO: 
        This can be a bit confusing to read. The code is kind of wordy. We could possibly redraw flowline with the 
        display Markers function but that would require some changes to the Markers function to take an index.
    '''

    def mouseMove(self, e):

        if self.isMarkerSelected:

            # change the x , y values of the marker at the selected index
            xPositionOfCursor = e.pos().x()
            yPositionOfCursor = e.pos().y()
            self.whichMarkerSelected.cx = xPositionOfCursor
            self.whichMarkerSelected.cy = yPositionOfCursor
            self.whichMarkerSelected.updateCross()

            # change the x, y values of the flowline at the selected index
            whichFlowlineSelected = self.whichIndexOfFlowlineSelected[0]
            indexSelected = self.whichIndexOfFlowlineSelected[1]
            self.flowlines[whichFlowlineSelected][indexSelected] = [self.whichMarkerSelected.dx,
                                                                    self.whichMarkerSelected.dy]

            self.flowlines[whichFlowlineSelected] = self.flowIntegrator.integrate(
                self.whichMarkerSelected.dx, self.whichMarkerSelected.dy,
                self.flowlines[whichFlowlineSelected], indexSelected,
                float(self.spatialResolutionLineEdit.text()))

            # Remove every marker past the one we selected
            for i in range(self.selectedMarkerPosition[1] + 1, len(self.flowlineMarkers[0])):
                self.imageItemContainer.currentWidget().removeItem(
                    self.flowlineMarkers[self.selectedMarkerPosition[0]][i])

                # get the flowline position of the new marker
                newPosition = self.flowlines[whichFlowlineSelected][i * self.integratorPerMarker]
                cx, cy = colorCoord(newPosition[0], newPosition[1])

                # Create new marker with new data
                self.flowlineMarkers[self.selectedMarkerPosition[0]][i] = Marker(
                    cx, cy, newPosition[0], newPosition[1],
                    self.imageItemContainer.currentWidget())
            # This section redraws the new markers
            for i in range(self.selectedMarkerPosition[1] + 1, len(self.flowlineMarkers[0])):
                self.imageItemContainer.currentWidget().addItem(
                    self.flowlineMarkers[self.selectedMarkerPosition[0]][i].getCross()[0])
                self.imageItemContainer.currentWidget().addItem(
                    self.flowlineMarkers[self.selectedMarkerPosition[0]][i].getCross()[1])

                xa = [self.flowlineMarkers[self.selectedMarkerPosition[0]][i - 1].cx,
                      self.flowlineMarkers[self.selectedMarkerPosition[0]][i].cx]
                ya = [self.flowlineMarkers[self.selectedMarkerPosition[0]][i - 1].cy,
                      self.flowlineMarkers[self.selectedMarkerPosition[0]][i].cy]
                self.flowlineMarkers[self.selectedMarkerPosition[0]][i].setLine(
                    pg.PlotDataItem(xa, ya, connect='all', pen=skinnyBlackPlotPen), 0)
                self.flowlineMarkers[self.selectedMarkerPosition[0]][i - 1].setLine(
                    self.flowlineMarkers[self.selectedMarkerPosition[0]][i].lines[0], 1)

                self.imageItemContainer.currentWidget().addItem(
                    self.flowlineMarkers[self.selectedMarkerPosition[0]][i].lines[0])

            self.displayMarkerVariables()

            # Connect lines between marker selected and previous marker
            if self.whichMarkerSelected.lines[0] is not None:
                self.whichMarkerSelected.lines[0].setData(
                    [self.whichMarkerSelected.lines[0].getData()[0][0], self.whichMarkerSelected.cx],
                    [self.whichMarkerSelected.lines[0].getData()[1][0], self.whichMarkerSelected.cy])

    '''
    Function: calcVelocityWidth
    Argument list: 
    Purpose: 
        Calculates velocity width by connecting the ith marker of each shear margin. This displays lines between
        each displayed marker as well. This also does an averaging scheme where we create lines between the two shear
        margins starting at the terminus. The number of lines is determined by numberOfLines. The averaging is done
        by averaging the ith index of each line together.
        TODO: Reference paper when written
    Return types, values: 
    Dependencies: Two selected shear margins. This is susceptible to user errors. TODO: Fix the errors
    Creator: James Stauder
    Date created: 2/25/18
    Last edited: 3/9/18
    '''

    def calcVelocityWidth(self):
        t0 = time.time()

        self.flowlines = self.flowlines[0:2]
        midValues, naiveValues, sampleValues = self.createThreeProfiles()

        oldVolumesStep = []
        newVolumesStep = []
        oldVolumes = []
        newVolumes = []
        trueVolumes = []
        x1, y1 = self.flowlineMarkers[0][0].cx, self.flowlineMarkers[0][0].cy
        x2, y2 = self.flowlineMarkers[1][0].cx, self.flowlineMarkers[1][0].cy

        numbersToTest = range(1, 200, 1)

        timeTrue = time.time()

        print "calculating true volume"
        for x in range(0, len(self.flowlines[0]) - 1):
            print "solving true for step ", x
            trueVolumes.append(self.calcTrueVolume(x))

        print "calculated true volume in time: ", time.time() - timeTrue

        for numberOfLines in numbersToTest:
            self.flowlines = self.flowlines[0:2]
            dx = (x2 - x1) / (numberOfLines + 1)
            dy = (y2 - y1) / (numberOfLines + 1)
            currX = x1
            currY = y1

            # Create center flowlines

            for _ in range(0, numberOfLines):
                currX = currX + dx
                currY = currY + dy

                xProj, yProj = colorToProj(currX, currY)

                newLine = []

                for i in range(self.lengthOfFlowline):
                    newLine.append(None)
                newLine[0] = [xProj, yProj]
                newLine = self.flowIntegrator.integrate(xProj, yProj, newLine, 0,
                                                        float(self.spatialResolutionLineEdit.text()))

                # This checks to see if integration worked
                if None not in newLine:
                    self.flowlines.append(newLine)
                else:
                    print "integration error on averaging method. Ommitting line"

            print "Ommitted ", numberOfLines - len(self.flowlines) + 2, " lines. Out of a possible ", numberOfLines
            tempOldS, tempNewS, tempOld, tempNew, = self.calcVolume()
            oldVolumes.append(tempOld)
            newVolumes.append(tempNew)
            oldVolumesStep.append(tempOldS)
            newVolumesStep.append(tempNewS)

        widthList = self.getWidths()
        trueDX = self.getDX()
        sampleFlowlines = self.flowlines[2:]
        midFlowline = sampleFlowlines[len(sampleFlowlines) / 2]

        thicknessList = []

        for step in midFlowline:
            thickness = self.datasetDict['thickness'].getInterpolatedValue(step[0], step[1])[0][0]
            thicknessList.append(thickness)
        thicknessList = np.asarray(thicknessList)

        midVolumeStep = (
            ((thicknessList[:-1] + thicknessList[1:]) / 2) *
            ((widthList[:-1] + widthList[1:]) / 2) *
            trueDX)
        print "Profile creation took :", time.time() - t0

        trueVolume = sum(trueVolumes)
        midVolume = sum(midVolumeStep)

        horizLineDataTrue = np.array([trueVolume for i in numbersToTest])
        horizLineDataMid = np.array([midVolume for i in numbersToTest])
        plt.plot(numbersToTest, (oldVolumes), label='Naive way', color='red')
        plt.plot(numbersToTest, (newVolumes), label='New way', color='blue')
        plt.plot(numbersToTest, horizLineDataMid, label='Midline', color='green')
        plt.plot(numbersToTest, horizLineDataTrue, label='True Volume?', color='black')
        plt.legend()
        plt.show()

        plt.plot(oldVolumesStep[1], label='Old volume by Step', color='red')
        plt.plot(newVolumesStep[1], label='New volume by Step', color='blue')
        plt.plot(midVolumeStep, label='Mid Volume by Step', color='green')
        plt.plot(trueVolumes, label='True volume by Step', color='black')
        plt.legend()
        plt.show()

        outFile = h5py.File('.data/' + str(self.profileLineEdit.text()), 'w')

        outFile.create_dataset("/true_volume_step", data=trueVolumes)
        outFile.create_dataset("/old_volume_step", data=oldVolumesStep)
        outFile.create_dataset("/new_volume_step", data=newVolumesStep)
        outFile.create_dataset("/mid_volume_step", data=midVolumeStep)
        outFile.create_dataset("/shear1", data=self.flowlines[0])
        outFile.create_dataset("/shear2", data=self.flowlines[1])

        outFile.create_dataset('mid_bed', data = midValues['bed'])
        outFile.create_dataset('mid_surface', data=midValues['surface'])
        outFile.create_dataset('mid_velocity', data=midValues['velocity'])

        outFile.create_dataset('naive_bed', data = naiveValues['bed'])
        outFile.create_dataset('naive_surface', data=naiveValues['surface'])
        outFile.create_dataset('naive_velocity', data=naiveValues['velocity'])

        outFile.create_dataset('sample_bed', data = sampleValues['bed'])
        outFile.create_dataset('sample_surface', data=sampleValues['surface'])
        outFile.create_dataset('sample_velocity', data=sampleValues['velocity'])


        outFile.close()

    def createThreeProfiles(self):

        numSamples = 50

        midFlowlineValues = {}
        naiveAverageValues = {}
        sampleAverageValues = {}

        midFlowlineValues['bed'] = []
        midFlowlineValues['surface'] = []
        midFlowlineValues['velocity'] = []

        naiveAverageValues['bed'] = []
        naiveAverageValues['surface'] = []
        naiveAverageValues['velocity'] = []

        sampleAverageValues['bed'] = []
        sampleAverageValues['surface'] = []
        sampleAverageValues['velocity'] = []

        x1, y1 = self.flowlineMarkers[0][0].cx, self.flowlineMarkers[0][0].cy
        x2, y2 = self.flowlineMarkers[1][0].cx, self.flowlineMarkers[1][0].cy
        midX = (x2 + x1) / 2
        midY = (y2 + y1) / 2
        xProj, yProj = colorToProj(midX, midY)

        midFlowline = []
        for i in range(self.lengthOfFlowline):
            midFlowline.append(None)
        midFlowline[0] = [xProj, yProj]
        midFlowline = self.flowIntegrator.integrate(xProj, yProj, midFlowline, 0,
                                                    float(self.spatialResolutionLineEdit.text()))

        if None in midFlowline:
            print "midFlowline error. Exiting"
            sys.exit()

        for i in midFlowline:
            midFlowlineValues['bed'].append(self.datasetDict['bed'].getInterpolatedValue(i[0], i[1])[0][0])
            midFlowlineValues['surface'].append(self.datasetDict['surface'].getInterpolatedValue(i[0], i[1])[0][0])
            midFlowlineValues['velocity'].append(self.datasetDict['velocity'].getInterpolatedValue(i[0], i[1])[0][0])

        shear1 = self.flowlines[0]
        shear2 = self.flowlines[1]

        for i in range(len(self.flowlines[0])):
            P1 = shear1[i]
            P2 = shear2[i]

            xValues = np.linspace(P1[0], P2[0], num=numSamples, endpoint=True)
            yValues = np.linspace(P1[1], P2[1], num=numSamples, endpoint=True)

            bedValues = []
            surfaceValues = []
            velocityValues = []
            for j in range(len(xValues)):
                bedValues.append(self.datasetDict['bed'].getInterpolatedValue(xValues[j], yValues[j])[0][0])
                surfaceValues.append(self.datasetDict['surface'].getInterpolatedValue(xValues[j], yValues[j])[0][0])
                velocityValues.append(self.datasetDict['velocity'].getInterpolatedValue(xValues[j], yValues[j])[0][0])

            naiveAverageValues['bed'].append(np.average(bedValues))
            naiveAverageValues['surface'].append(np.average(surfaceValues))
            naiveAverageValues['velocity'].append(np.average(velocityValues))

        xValues = np.linspace(x1, x2, num=numSamples, endpoint=True)
        yValues = np.linspace(y1, y2, num=numSamples, endpoint=True)

        # cut off beginning and endPoint
        for i in range(1, len(xValues) - 1):
            xProj, yProj = colorToProj(xValues[i], yValues[i])

            newLine = []

            for j in range(self.lengthOfFlowline):
                newLine.append(None)

            newLine[0] = [xProj, yProj]

            newLine = self.flowIntegrator.integrate(xProj, yProj, newLine, 0,
                                                    float(self.spatialResolutionLineEdit.text()))

            if None in newLine:
                print "error in sampling. Exiting"
                sys.exit()

            self.flowlines.append(newLine)

        for i in range(len(self.flowlines[0])):

            bedValues = []
            surfaceValues = []
            velocityValues = []

            for j in range(len(self.flowlines)):
                bedValues.append(self.datasetDict['bed'].getInterpolatedValue(self.flowlines[j][i][0],
                                                                              self.flowlines[j][i][1])
                                 [0][0])

                surfaceValues.append(self.datasetDict['surface'].getInterpolatedValue(self.flowlines[j][i][0],
                                                                                      self.flowlines[j][i][1])
                                     [0][0])

                velocityValues.append(self.datasetDict['velocity'].getInterpolatedValue(self.flowlines[j][i][0],
                                                                                        self.flowlines[j][i][1])
                                      [0][0])

            sampleAverageValues['bed'].append(np.average(bedValues))
            sampleAverageValues['surface'].append(np.average(surfaceValues))
            sampleAverageValues['velocity'].append(np.average(velocityValues))



        return midFlowlineValues, naiveAverageValues, sampleAverageValues

    def calcTrueVolume(self, step):

        '''
        trueVolume = 0
        resolution = 50

        for i in range(len(self.flowlines[0]) - 1):
            shearPoints1 = [[self.flowlines[0][i][0], self.flowlines[0][i][1]],
                            [self.flowlines[1][i][0], self.flowlines[1][i][1]]]

            totalWidth = sqrt(
                (shearPoints1[0][0] - shearPoints1[1][0]) ** 2 + (shearPoints1[0][1] - shearPoints1[1][1]) ** 2)
            numPoints = math.ceil(totalWidth / resolution)

            dx1 = (shearPoints1[1][0] - shearPoints1[0][0]) / (numPoints + 1)
            dy1 = (shearPoints1[1][1] - shearPoints1[0][1]) / (numPoints + 1)

            shearPoints2 = [[self.flowlines[0][i + 1][0], self.flowlines[0][i + 1][1]],
                            [self.flowlines[1][i + 1][0], self.flowlines[1][i + 1][1]]]

            dx2 = (shearPoints2[1][0] - shearPoints2[0][0]) / (numPoints + 1)
            dy2 = (shearPoints2[1][1] - shearPoints2[0][1]) / (numPoints + 1)

            line1Values = [shearPoints1[0]]
            line2Values = [shearPoints2[0]]

            for _ in range(int(numPoints) + 1):
                line2Values.append([line2Values[-1][0] + dx2, line2Values[-1][1] + dy2])
                line1Values.append([line1Values[-1][0] + dx1, line1Values[-1][1] + dy1])

            for j in range(len(line1Values) - 1):
                # Really confusing way to make a quadrilateral set of points to calculate the volume
                points = []


                points.append([line1Values[j][0], line1Values[j][1],
                               self.datasetDict['surface'].getInterpolatedValue(line1Values[j][0], line1Values[j][1])[0][
                                   0] -
                               self.datasetDict['thickness'].getInterpolatedValue(line1Values[j][0], line1Values[j][1])[
                                   0][0]])
                points.append([line1Values[j][0], line1Values[j][1],
                               self.datasetDict['surface'].getInterpolatedValue(line1Values[j][0], line1Values[j][1])[
                                   0][0]])

                points.append([line1Values[j + 1][0], line1Values[j + 1][1],
                               self.datasetDict['surface'].getInterpolatedValue(line1Values[j + 1][0],
                                                                            line1Values[j + 1][1])[0][0] -
                               self.datasetDict['thickness'].getInterpolatedValue(line1Values[j + 1][0],
                                                                                  line1Values[j + 1][1])[0][0]])
                points.append([line1Values[j + 1][0], line1Values[j + 1][1],
                               self.datasetDict['surface'].getInterpolatedValue(line1Values[j + 1][0],
                                                                                line1Values[j + 1][1])[0][0]])



                points.append([line2Values[j][0], line2Values[j][1],
                               self.datasetDict['surface'].getInterpolatedValue(line2Values[j][0], line2Values[j][1])[0][
                                   0] -
                               self.datasetDict['thickness'].getInterpolatedValue(line2Values[j][0], line2Values[j][1])[
                                   0][0]])
                points.append([line2Values[j][0], line2Values[j][1],
                               self.datasetDict['surface'].getInterpolatedValue(line2Values[j][0], line2Values[j][1])[
                                   0][0]])

                points.append([line2Values[j + 1][0], line2Values[j + 1][1],
                               self.datasetDict['surface'].getInterpolatedValue(line2Values[j + 1][0],
                                                                            line2Values[j + 1][1])[0][0] -
                               self.datasetDict['thickness'].getInterpolatedValue(line2Values[j + 1][0],
                                                                            line2Values[j + 1][1])[0][0]])
                points.append([line2Values[j + 1][0], line2Values[j + 1][1],
                               self.datasetDict['surface'].getInterpolatedValue(line2Values[j + 1][0],
                                                                                line2Values[j + 1][1])[0][0]])



                dt = Delaunay(np.asarray(points), qhull_options='QbB Pp')
                tets = dt.points[dt.simplices]

                trueVolume = trueVolume + np.sum(self.tetrahedron_volume(tets[:, 0], tets[:, 1],
                                                                         tets[:, 2], tets[:, 3]))
        testVolume = trueVolume
        '''

        shear1 = np.asarray(self.flowlines[0][step:step + 2])
        shear2 = np.asarray(self.flowlines[1][step:step + 2])

        myShears = []
        myShears.extend(shear1)
        myShears.extend(list(reversed(shear2)))
        myShears.append(shear1[0])
        myShears = np.asarray(myShears)

        myPath = Path(myShears, closed=True)

        maxY, minY = max(myShears[:, 1]), min(myShears[:, 1])
        minY = minY - 150  # buffer
        maxY = maxY + 150  # buffer

        maxX, minX = max(myShears[:, 0]), min(myShears[:, 0])
        maxX = maxX + 150
        minX = minX - 150

        allDataFile = h5py.File(fullDataFileName, 'r')

        xs = allDataFile['x'][:]
        ys = allDataFile['y'][:]

        thickness = allDataFile['thickness'][:]

        xx, yy = np.meshgrid(xs, ys)
        XY = np.dstack((xx, yy))
        XYFlat = XY.reshape(-1, 2)

        thicknessFlat = thickness.reshape(-1)

        dataMaxY = XYFlat[0][1]
        dataMinY = XYFlat[-1][1]

        YToCutMax = abs(dataMaxY - maxY) / 150
        YToCutMin = abs(dataMinY - minY) / 150

        numX = 10018
        samplePoints = XYFlat[int(math.floor(YToCutMax) * numX):-(int(math.floor(YToCutMin) * numX))]
        thicknessFlat = thicknessFlat[int(math.floor(YToCutMax) * numX):-(int(math.floor(YToCutMin) * numX))]

        newSamplePoints = []
        newThickness = []
        for i in range(numX):
            newSamplePoints.extend(samplePoints[i::numX])
            newThickness.extend(thicknessFlat[i::numX])

        newSamplePoints = np.asarray(newSamplePoints)
        newThickness = np.asarray(newThickness)

        dataMaxX = newSamplePoints[-1][0]
        dataMinX = newSamplePoints[0][0]

        XToCutMax = int(math.floor(abs(dataMaxX - maxX) / 150))
        XToCutMin = int(math.floor(abs(dataMinX - minX) / 150))

        newMaxY, newMinY = newSamplePoints[0][1], newSamplePoints[-1][1]

        numY = int(abs(newMaxY - newMinY) / 150) + 1

        finalSamplePoints = newSamplePoints[XToCutMin * numY:-(XToCutMax * numY)]
        finalThicknessPoints = newThickness[XToCutMin * numY:-(XToCutMax * numY)]

        thicknessValuesInShape = []

        pointsInShape = []

        time0 = time.time()
        for i in range(len(finalSamplePoints)):

            point = [int(finalSamplePoints[i][0]), int(finalSamplePoints[i][1])]

            if myPath.contains_point(point):
                thicknessValuesInShape.append(finalThicknessPoints[i])
                pointsInShape.append(point)

        '''
        pointsInShape = np.asarray(pointsInShape)
        plt.scatter(finalSamplePoints[:,0], finalSamplePoints[:,1],1, color='red')
        plt.scatter(pointsInShape[:,0], pointsInShape[:,1],1, color='green')
        plt.plot(myShears[:,0], myShears[:,1], color='blue')
        plt.show()
        '''

        trueVolume = sum(thicknessValuesInShape) * 150 ** 2

        return trueVolume

    def tetrahedron_volume(self, a, b, c, d):
        return np.abs(np.einsum('ij,ij->i', a - d, np.cross(b - d, c - d))) / 6

    def getWidths(self):

        shear1 = np.asarray(self.flowlines[0])
        shear2 = np.asarray(self.flowlines[1])

        width = sqrt((shear1[:, 0] - shear2[:, 0]) ** 2 + (shear1[:, 1] - shear2[:, 1]) ** 2)
        return np.asarray(width)

    def getDX(self):
        dx = []

        for i in range(0, len(self.flowlines[0]) - 1):
            P1 = self.flowlines[0][i]
            P2 = self.flowlines[1][i]
            P3 = self.flowlines[0][i + 1]
            P4 = self.flowlines[1][i + 1]

            d1 = (abs(
                (P3[1] - P4[1]) * P1[0] +
                (P4[0] - P3[0]) * P1[1] +
                (P3[0] * P4[1] - P4[0] * P3[1])
            )) / (sqrt(
                (P3[1] - P4[1]) ** 2 +
                (P4[0] - P3[0]) ** 2
            ))

            d2 = (abs(
                (P3[1] - P4[1]) * P2[0] +
                (P4[0] - P3[0]) * P2[1] +
                (P3[0] * P4[1] - P4[0] * P3[1])
            )) / (sqrt(
                (P3[1] - P4[1]) ** 2 +
                (P4[0] - P3[0]) ** 2
            ))

            dx.append((d1 + d2) / 2)

        return dx

    def calcVolume(self):

        widthList = self.getWidths()
        trueDX = self.getDX()

        thicknessList = []
        for step in range(0, len(self.flowlines[0])):
            thickness = []

            for sample in range(2, len(self.flowlines)):
                samplePoint = self.flowlines[sample][step]
                thickness.append(
                    self.datasetDict['thickness'].getInterpolatedValue(samplePoint[0], samplePoint[1])[0][0])

            thicknessList.append(np.mean(np.asarray(thickness)))
        thicknessList = np.asarray(thicknessList)

        newMethodVolumeStep = (
            ((thicknessList[:-1] + thicknessList[1:]) / 2) *
            ((widthList[:-1] + widthList[1:]) / 2) *
            trueDX)

        thicknessList = []
        for step in range(0, len(self.flowlines[0])):

            shearPoints = [[self.flowlines[0][step][0], self.flowlines[0][step][1]],
                           [self.flowlines[1][step][0], self.flowlines[1][step][1]]]

            dx = (shearPoints[1][0] - shearPoints[0][0]) / (len(self.flowlines) - 1)
            dy = (shearPoints[1][1] - shearPoints[0][1]) / (len(self.flowlines) - 1)

            currPoint = shearPoints[0]

            thickness = []
            for sample in range(2, len(self.flowlines)):
                currPoint[0] = currPoint[0] + dx
                currPoint[1] = currPoint[1] + dy

                thickness.append(self.datasetDict['thickness'].getInterpolatedValue(
                    currPoint[0], currPoint[1])[0][0])

            thicknessList.append(np.mean(np.asarray(thickness)))

        thicknessList = np.asarray(thicknessList)

        oldMethodVolumeStep = (
            ((thicknessList[:-1] + thicknessList[1:]) / 2) *
            ((widthList[:-1] + widthList[1:]) / 2) *
            trueDX)

        return oldMethodVolumeStep, newMethodVolumeStep, sum(oldMethodVolumeStep), sum(newMethodVolumeStep)

        '''
        Function: displayMarkers
        Argument list: 
            flowline: flowline in which to display
        Purpose: Takes a flowline of markers and displays them on the gui
        Return types, values: None
        Dependencies: None
        Creator: James Stauder
        Date created: 3/2/18
        Last edited: 3/2/18
        '''

    def displayMarkers(self, flowline):

        # Add first marker. This needs to be peeled because the for loop
        # connects the markers backwards
        self.imageItemContainer.currentWidget().addItem(flowline[0].getCross()[0])
        self.imageItemContainer.currentWidget().addItem(flowline[0].getCross()[1])

        for i in range(1, len(flowline)):
            self.imageItemContainer.currentWidget().addItem(flowline[i].getCross()[0])
            self.imageItemContainer.currentWidget().addItem(flowline[i].getCross()[1])

            xValuesOfMarkers = [flowline[i - 1].cx, flowline[i].cx]
            yValuesOfMarkers = [flowline[i - 1].cy, flowline[i].cy]

            # Create lines from each marker
            flowline[i].setLine(
                pg.PlotDataItem(xValuesOfMarkers, yValuesOfMarkers, connect='all', pen=skinnyBlackPlotPen), 0)
            flowline[i - 1].setLine(flowline[i].lines[0], 1)

            self.imageItemContainer.currentWidget().addItem(flowline[i].lines[0])

    '''
    Function: displayMarkerVariables
    Argument list: None
    Purpose: Displays the marker variables of the marker selected
    Return types, values: None
    Dependencies: Marker to be selected
    Creator: James Stauder
    Date created: 2/25/18
    Last edited: 2/25/18
    '''

    def displayMarkerVariables(self):
        self.textOut.clear()
        selectedMarkerX = self.whichMarkerSelected.dx
        selectedMarkerY = self.whichMarkerSelected.dy

        self.textOut.append(str((self.whichMarkerSelected.dx, self.whichMarkerSelected.dy)))

        for x in self.maps:
            stringOut = str(self.datasetDict[x.lower()].getInterpolatedValue(selectedMarkerX, selectedMarkerY))
            self.textOut.append(x + ": " + stringOut[2:-2])

    '''
    Function: createIntegrator
    Argument list: Nones
    Purpose: Create integrator class. This will allow us to integrate up the ice flow
    Return types, values: None
    Dependencies: None
    Creator: James Stauder
    Date created: 2/5/18
    Last edited: 2/5/18
    '''

    # TODO: Does this have to be tied to mw? Can this be changed in some way?
    def createIntegrator(self):
        vx = Dataset('VX')
        vy = Dataset('VY')
        self.flowIntegrator = FlowIntegrator(vx, vy)

    def runModel(self):
        m = ModelGUI(self)

    def reset(self):
        del self.flowlines[:]
        del self.flowlineMarkers[:]
        for x in self.datasetDict:
            self.datasetDict[x].pathData = None
        self.runModelButton.setEnabled(False)
        self.spatialResolutionLineEdit.setReadOnly(False)
        self.distanceLineEdit.setReadOnly(False)

    '''
    Function: connectButtons
    Argument list: None
    Purpose: connect the buttons of the gui to various functions
    Return types, values: None
    Dependencies: None
    Creator: James Stauder
    Date created: 2/5/18
    Last edited: 2/5/18
    '''

    def connectButtons(self):
        self.mapList.currentIndexChanged.connect(self.changeMap)
        self.instructionButton.clicked.connect(self.showInstructions)
        self.velocityWidthButton.clicked.connect(self.calcVelocityWidth)
        self.runModelButton.clicked.connect(self.runModel)
        self.resetButton.clicked.connect(self.reset)

    def showInstructions(self):
        Instructions(self)
