from helperFiles.classes.MainWindow import *
from helperFiles.createDatasets import *

'''
Creator: James Stauder
Creation Date: 1/30/18
Last edit Date: 5/29/18
Purpose: GUI application that allows profile creation and model running on Greenland Glaciers.
'''


'''
Main function
Arguments list:
--help -prints out menu and exits
Purpose:
Return types, values:
Dependencies:
Creator: James Stauder
Date created: 1/31/18
Last edited: 2/02/18
'''


def main(argv):

    if len(argv) > 1:
        if "--help" in argv:
            printMainMenu()
            sys.exit()
    else:
        app = QApplication(sys.argv)
        mw = MainWindow()
        mw.show()
        datasetDict = createInitialDatasets()

        datasetDict['velocity'].imageItem.mouseClickEvent = mw.mouseClick
        datasetDict['velocity'].imageItem.hoverEvent = mw.mouseMove
        mw.datasetDict = datasetDict
        mw.createIntegrator()

        mw.addToImageItemContainer(datasetDict)

        sys.exit(app.exec_())

'''
Function: printMainMenu
Argument list: None
Purpose: Print the main menu
Return types, values: None
Dependencies: None
Creator: James Stauder
Date created: 1/31/18
Last edited: 1/31/18
'''


def printMainMenu():
    print("Greenland Ice Sheet modeling tool")
    print("usage: greenland.py [-h]")
    print("optional arguments:")
    print("     --help  show help message and exit")


if __name__ == '__main__':
    main(sys.argv)
