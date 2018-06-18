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

        oldVolumes = []
        newVolumes = []
        x1, y1 = self.flowlineMarkers[0][0].cx, self.flowlineMarkers[0][0].cy
        x2, y2 = self.flowlineMarkers[1][0].cx, self.flowlineMarkers[1][0].cy

        numbersToTest = range(20, 100)


        '''
        trueVolume = self.calcTrueVolume()
        print trueVolume
        '''
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
            tempOld, tempNew = self.calcVolume()
            oldVolumes.append(tempOld)
            newVolumes.append(tempNew)
        print "Profile creation took :", time.time() - t0
        print oldVolumes
        print newVolumes
        #print trueVolume

        plt.plot(oldVolumes, label='Old bad way', color='red')
        plt.plot(newVolumes, label='Jimmys awesome way', color='blue')
        # plt.plot(trueVolume, label='True Volume?', color='black')
        plt.legend()
        plt.show()

    def calcTrueVolume(self):
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
        return trueVolume

    def tetrahedron_volume(self, a, b, c, d):
        return np.abs(np.einsum('ij,ij->i', a - d, np.cross(b - d, c - d))) / 6

    def calcVolume(self):
        newMethodVolume = 0
        oldMethodVolume = 0

        for i in range(len(self.flowlines)):
            for j in range(len(self.flowlines[0])):
                self.flowlines[i][j].append(self.datasetDict['thickness'].getInterpolatedValue(self.flowlines[i][j][0], self.flowlines[i][j][1])[0][0]
                                            )


        myFlowlines = np.asarray(self.flowlines)
        myFlowlines = np.rot90(myFlowlines,k=3)

        






        averageThicknessList = []
        widthList = []

        for i in range(0, len(self.flowlines[0])):

            totalThickness = 0
            shearPoints = [[self.flowlines[0][i][0], self.flowlines[0][i][1]],
                           [self.flowlines[1][i][0], self.flowlines[1][i][1]]]

            dx = (shearPoints[1][0] - shearPoints[0][0]) / (len(self.flowlines) - 1)
            dy = (shearPoints[1][1] - shearPoints[0][1]) / (len(self.flowlines) - 1)

            totalWidth = sqrt(
                (shearPoints[0][0] - shearPoints[1][0]) ** 2 + (shearPoints[0][1] - shearPoints[1][1]) ** 2)
            currPoint = shearPoints[0]

            for j in range(0, len(self.flowlines)):
                totalThickness = totalThickness + self.datasetDict['thickness'].getInterpolatedValue(
                    currPoint[0], currPoint[1])[0][0]

                currPoint[0] = currPoint[0] + dx
                currPoint[1] = currPoint[1] + dy

            widthList.append(totalWidth)
            averageThicknessList.append((totalThickness / len(self.flowlines)))

        for i in range(0, len(averageThicknessList) - 1):
            oldMethodVolume = oldMethodVolume + \
                              (((averageThicknessList[i] * widthList[i]) + (
                                  averageThicknessList[i + 1] * widthList[i + 1])) / 2
                               * float(self.spatialResolutionLineEdit.text()))

        return oldMethodVolume, newMethodVolume

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
