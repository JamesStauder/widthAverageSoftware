from pylab import sqrt, linspace
from scipy.interpolate import RectBivariateSpline

from ..createColorMaps import *

import numpy as np

'''
Class: Dataset
Argument list: name of dataset, pen type(used for plotting)
Purpose: This is the class of datasets. This will store velocity, smb, etc. This takes the Velocity in X and Y direction
and makes one dataset of just Velocity. This velocity dataset ONLY stores the magnitude but not direction.

Dependencies: pylabs sqrt and linspace, RectBivariateSplint, numpy
Creator: James Stauder
Date created:2/23/18
Last edited: 3/2/18
'''


class Dataset:
    def __init__(self, name):
        self.name = name

        bed_xarray = linspace(map['proj_x0'], map['proj_x1'], map['x1'], endpoint=True)
        bed_yarray = linspace(map['proj_y1'], map['proj_y0'], map['y1'], endpoint=True)

        if self.name == 'velocity':
            self.data, self.vx, self.vy = self.setData(name)

            self.vxInterp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(self.vx).transpose())
            self.vyInterp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(self.vy).transpose())

        else:
            self.data = self.setData(name)

        self.interp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(self.data).transpose())

        # Only create color map if we wish to render possible in the future where we wish to create
        # Dataset for backend web
        if name != 'VX' and name != 'VY':
            createColorMap(self)

    def setData(self, name):
        dataFile = h5py.File(dataFileName, 'r')
        if name == 'velocity':
            vx = dataFile['VX'][:]
            vy = dataFile['VY'][:]
            data = sqrt(vx ** 2 + vy ** 2)
            dataFile.close()
            return data, vx, vy
        else:
            data = dataFile[name][:]
            dataFile.close()
            return data

    def getInterpolatedValue(self, xPosition, yPosition):
        return self.interp(xPosition, yPosition)
