from pyqtgraph.Qt import QtGui
from PyQt4 import QtCore
from BlackBox import *
from ModelPlotter import *

'''
Class: ModelGUI
Argument list: 
    parent -> MainWindow parent
Purpose: Create gui for the model run
Return types, values:
Dependencies: pyqtGraph
Creator: James Stauder
Date created: 3/29/18
Last edited: 5/29/18
'''


class ModelGUI(QtGui.QMainWindow):
    def __init__(self, parent):
        self.hdfName = None
        self.plotter = None
        self.parent = parent

        QtGui.QMainWindow.__init__(self, self.parent)

        self.centerWidget = QtGui.QWidget()
        self.setCentralWidget(self.centerWidget)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.centerWidget.setLayout(self.horizontalLayout)

        self.createRightPanel()
        self.createLeftPanel()

        self.horizontalLayout.addWidget(self.leftPanelWidget)
        self.horizontalLayout.addWidget(self.rightPanelWidget)

        self.runButton.clicked.connect(self.runModelEvent)
        self.pauseButton.clicked.connect(self.pauseModel)

        runAverage = self.parent.widthAverageButton.checkState() == 2

        self.runModel = BlackBox('.data/latestProfile.h5', float(self.timeEndLineEdit.text()),
                                 float(self.timeStepLineEdit.text()), average=runAverage)
        self.plots = ModelPlotter(self.runModel.strs, self.runModel.mesh, self.runModel.Bhat, self.plot1, self.plot2,
                                  self.plot3)

        self.showMaximized()
        self.run = True
        self.show()
        self.closeEvent = self.windowClosed

    # Pause model run when window closed
    def windowClosed(self, e):
        if self.run:
            self.run = False

    def createRightPanel(self):
        self.rightPanelWidget = QtGui.QWidget()
        self.rightPanelWidget.setMaximumWidth(300)
        self.rightPanelLayout = QtGui.QGridLayout()
        self.rightPanelLayout.setAlignment(QtCore.Qt.AlignTop)
        self.rightPanelLayout.setSpacing(4)

        self.rightPanelWidget.setLayout(self.rightPanelLayout)

        self.timeEndLabel = QtGui.QLabel('time end(yr):')
        self.timeEndLineEdit = QtGui.QLineEdit('20000')
        self.timeStepLabel = QtGui.QLabel('time step(yr):')
        self.timeStepLineEdit = QtGui.QLineEdit('10')
        self.timeCurrent = QtGui.QLabel('Current year: ')

        # Buttons
        self.runButton = QtGui.QPushButton('Run Model')

        self.pauseButton = QtGui.QPushButton('Pause Model')
        self.pauseButton.setEnabled(False)

        self.rightPanelLayout.addWidget(self.timeEndLabel, 0, 0)
        self.rightPanelLayout.addWidget(self.timeEndLineEdit, 0, 1)
        self.rightPanelLayout.addWidget(self.timeStepLabel, 1, 0)
        self.rightPanelLayout.addWidget(self.timeStepLineEdit, 1, 1)
        self.rightPanelLayout.addWidget(self.timeCurrent, 3, 0, 1, 2)
        self.rightPanelLayout.addWidget(self.runButton, 4, 0, 1, 2)
        self.rightPanelLayout.addWidget(self.pauseButton, 5, 0, 1, 2)

    def createLeftPanel(self):
        self.leftPanelWidget = QtGui.QWidget()
        self.leftPanelLayout = QtGui.QVBoxLayout()
        self.leftPanelWidget.setLayout(self.leftPanelLayout)

        self.plot1 = pg.PlotWidget()
        self.plot2 = pg.PlotWidget()
        self.plot3 = pg.PlotWidget()

        self.slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.slider.setEnabled(False)
        self.sliderLabel = QtGui.QLabel('Year: ')
        self.leftPanelLayout.addWidget(self.plot1)
        self.leftPanelLayout.addWidget(self.plot2)
        self.leftPanelLayout.addWidget(self.plot3)
        self.leftPanelLayout.addWidget(self.slider)
        self.leftPanelLayout.addWidget(self.sliderLabel)

    #run Model
    def runModelEvent(self):
        self.runButton.setEnabled(False)
        self.pauseButton.setEnabled(True)

        self.run = True
        while self.runModel.t < self.runModel.timeEnd and self.run is True:
            self.plots.refreshPlot(self.runModel)
            pg.QtGui.QApplication.processEvents()

    # pause Model
    def pauseModel(self):
        if self.run:
            self.run = False
            self.pauseButton.setText('Resume')
        else:
            self.run = True
            self.pauseButton.setText('Pause')
            while self.runModel.t < self.runModel.timeEnd and self.run is True:
                self.plots.refreshPlot(self.runModel)
                pg.QtGui.QApplication.processEvents()
