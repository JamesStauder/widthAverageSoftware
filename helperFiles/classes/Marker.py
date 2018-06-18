from ..math_functions import *

'''
Class: Marker
Argument list: 
    cx - Color coordinate x position
    cy - Color coordinate y position
    dx - Projection x position
    dy - Projection y position
plotWidget - Which widget to plot to
plotCross - Boolean on whether or not to plot

Purpose: This is the marker class that is used to create the x markers on the GUI
Dependencies: pyQT
Creator: Patrick Kreitzberg
Date created: Unknown
Last edited: 3/7/18 by James Stauder
'''


class Marker:
    def __init__(self, cx, cy, dx, dy, plotWidget, plotCross=True):
        self.plotWidget = plotWidget
        self.plotCross = plotCross
        self.cx, self.cy = cx, cy  # color coordinates
        self.dx, self.dy = dx, dy  # data coordinates
        self.px, self.py = colorToProj(self.cx, self.cy)
        self.pen = blackPlotPen
        self.pen.setWidth(2)
        self.c = 3
        if plotCross:
            self.cross = [
                pg.PlotDataItem([self.cx - self.c, self.cx + self.c], [self.cy - self.c, self.cy + self.c],
                                connect='all', pen=self.pen)
                , pg.PlotDataItem([self.cx - self.c, self.cx + self.c], [self.cy + self.c, self.cy - self.c],
                                  connect='all', pen=self.pen)
            ]
        self.lines = [None] * 2
        self.intLine = None

    def __del__(self):
        if self.lines[0]:
            self.plotWidget.removeItem(self.lines[0])
        if self.lines[1]:
            self.plotWidget.removeItem(self.lines[1])
        if self.plotCross:
            self.plotWidget.removeItem(self.cross[0])
            self.plotWidget.removeItem(self.cross[1])
        if self.intLine:
            self.plotWidget.removeItem(self.intLine)

    '''
    Function: updateCross
    Argument list: 
    Purpose: 
        updates the dx and dy coordinates of the cross. This is called when the cx and cy coordinates change(when
        the user moves the marker.) This then also moves the bars around
    Return types, values: 
    Dependencies: 
    Creator: Patrick Kreitzberg
    Date created: Unknown
    Last edited: 3/7/18
    '''

    def updateCross(self):
        self.dx, self.dy = colorToProj(self.cx, self.cy)
        self.cross[0].setData([self.cx - self.c, self.cx + self.c], [self.cy - self.c, self.cy + self.c], connect='all',
                              pen=self.pen)
        self.cross[1].setData([self.cx - self.c, self.cx + self.c], [self.cy + self.c, self.cy - self.c], connect='all',
                              pen=self.pen)
        self.cross[0].updateItems()
        self.cross[1].updateItems()

    def setPlotWidget(self, pw):
        self.plotWidget = pw

    '''
    Function: checkClicked
    Argument list: 
        pos - x and y coordinates of the clicked position
    Purpose: checks to see if the pos argument is on the marker
    Return types, values: 
        Boolean
        True - if marker is clicked
        False - if marker is not clicked
    Dependencies: 
    Creator: Patrick Kreitzberg
    Date created: Unknown
    Last edited: 3/7/18
    '''
    def checkClicked(self, pos):
        if self.cx - self.c <= pos.x() <= self.cx + self.c and self.cy - self.c <= pos.y() <= self.cy + self.c:
            return True
        else:
            return False

    '''
    Function: setLine
    Argument list: 
        line - the black line between markers
        index - which index does the line belong to (index 0 is marker previous, index 1 is marker next)
    Purpose: set the lines of the marker to nearby markers in an organized format
    Return types, values: 
    Dependencies: 
    Creator: Patrick Kreitzberg
    Date created: Unknown
    Last edited: 3/7/18
    '''
    def setLine(self, line, i):
        self.lines[i] = line

    def getLine(self):
        return self.lines

    def getCross(self):
        return self.cross[0], self.cross[1]
