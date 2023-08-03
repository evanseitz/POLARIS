import sys, os, time
import psutil, logging
import multiprocessing

import numpy as np

from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton,\
     QGridLayout, QLabel, QFrame, QLineEdit, QFileDialog,\
     QDoubleSpinBox, QAbstractSpinBox, QMessageBox, QCheckBox

# below ordering matters:
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from pylab import imshow, show, loadtxt, axes
from math import factorial
import copy

# ignore matplotlib deprecation warnings:
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

import POLARIS_BE

################################################################################
# global assets:

pyDir = os.path.dirname(os.path.realpath(__file__)) #python file location
progname = 'POLARIS'
progversion = "0.1"
font_standard = QtGui.QFont('Arial', 12)


################################################################################
# plotting:

class FigCanvas(FigureCanvas):

    initiated = 0
    
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        FigCanvas.fig = Figure(figsize=(width, height), dpi=dpi)
        FigCanvas.axes = FigCanvas.fig.add_subplot(111)
        self.compute_initial_figure()
        FigureCanvas.__init__(self, FigCanvas.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                    QtWidgets.QSizePolicy.Expanding,
                                    QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass

class MplCanvas(FigCanvas):
    def __init__(self, *args, **kwargs):
        FigCanvas.__init__(self, *args, **kwargs)

    def compute_initial_figure(self):
        FigCanvas.fig.set_tight_layout(True)
        FigCanvas.fig.patch.set_facecolor('k')
        FigCanvas.axes.patch.set_facecolor('k')
        FigCanvas.axes.xaxis.set_ticklabels([])
        FigCanvas.axes.yaxis.set_ticklabels([])
        empty_plot = FigCanvas.axes.imshow(np.zeros(shape=(10,10)),cmap='gray', vmin=0, vmax=1)      
        FigCanvas.cbar = FigCanvas.fig.colorbar(empty_plot)
        FigCanvas.cbar.outline.set_facecolor('k')
        FigCanvas.text = FigCanvas.fig.text(.5,.5, 'No files detected...',
                                  ha='center', va='center',
                                  color='w', fontsize=12)

    def update_figure(self):
        if P1.df:
            if FigCanvas.initiated == 0:            
                FigCanvas.text.remove()
            
            FigCanvas.axes.cla()
            FigCanvas.cbar.remove()
            FigCanvas.fig.set_tight_layout(True)
            FigCanvas.fig.patch.set_facecolor('w')
            try:
                LS = loadtxt(P1.df,float,delimiter=',')
            except ValueError:
                LS = loadtxt(P1.df,float)

            for j,k in P1.points_xy:
                FigCanvas.axes.scatter([j],[k], c='k', s=5)
                #FigCanvas.axes.hold(True)
                FigCanvas.axes.scatter([j],[k], c='w', s=.5)
                #FigCanvas.axes.hold(True)
            im = FigCanvas.axes.imshow(LS, cmap='jet', origin='lower', )
            #FigCanvas.axes.grid(True)
            FigCanvas.axes.format_coord = lambda x, y: "x={0:.0f}, y={0:.0f}".format(x,y)
            
            FigCanvas.cbar = FigCanvas.fig.colorbar(im, ax=FigCanvas.axes,
                                        orientation='vertical',
                                        use_gridspec=True)
            #FigCanvas.axes.set_xlabel('Reaction Coordinate 1', size=6)
            #FigCanvas.axes.set_ylabel('Reaction Coordinate 2', size=6)
            for tick in FigCanvas.axes.xaxis.get_major_ticks():
                tick.label.set_fontsize(6)
            for tick in FigCanvas.axes.yaxis.get_major_ticks():
                tick.label.set_fontsize(6)
            FigCanvas.cbar.ax.tick_params(labelsize=6)
            FigCanvas.cbar.ax.set_title(label='Energy',size=6)
            
            self.show()
            self.draw()

            FigCanvas.initiated = 1
                                          
        else:
            pass

################################################################################
# coordinates (tab 1):         

class P1(QtWidgets.QWidget):
    df = False
    df_prev = False
    LS = False
    run = 0
    coords_x = {}
    coords_y = {}
    coords_e = {}
    points_x = np.zeros([10,1], dtype=int)
    points_y = np.zeros([10,1], dtype=int)
    points_xy = np.zeros([10,2], dtype=int)
    upper_bound = 1
    
    def __init__(self, parent=None):
        super(P1, self).__init__(parent)
        layout = QGridLayout(self)
        layout.setContentsMargins(20,20,20,20)
        layout.setSpacing(20)

        P1.points_xy = np.zeros([10,2], dtype=int)

######## left, top:
        self.label_edge1a = QLabel('') #body
        self.label_edge1a.setMargin(5)
        self.label_edge1a.setLineWidth(1)
        self.label_edge1a.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        layout.addWidget(self.label_edge1a, 0, 0, 14, 11)
        self.label_edge1a.show()

        self.label_edge1b = QLabel('') #head
        self.label_edge1b.setMargin(5)
        self.label_edge1b.setLineWidth(1)
        self.label_edge1b.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        layout.addWidget(self.label_edge1b, 0, 0, 1, 11)
        self.label_edge1b.show()

        self.label_LS = QLabel('Energy Landscape')
        self.label_LS.setFont(font_standard)
        self.label_LS.setMargin(5)
        self.label_LS.setFrameStyle(QFrame.Box | QFrame.Sunken)
        self.label_LS.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        layout.addWidget(self.label_LS, 0, 0, 1, 11)
        self.label_LS.show()

        self.entry_LS = QLineEdit('Filename')
        self.entry_LS.setDisabled(True)
        layout.addWidget(self.entry_LS, 1, 1, 1, 8)
        self.entry_LS.show()

        self.button_browse = QPushButton('Browse', self)
        self.button_browse.clicked.connect(self.browseBox)
        self.button_browse.setToolTip('Load a supported data file.')
        layout.addWidget(self.button_browse, 1, 9, 1, 1)
        self.button_browse.show()
        
        self.plot_canvas = MplCanvas(self, width=5, height=4, dpi=100)
        layout.addWidget(self.plot_canvas, 2, 1, 9, 9)
        
        self.navi_toolbar = NavigationToolbar(self.plot_canvas, self)
        layout.addWidget(self.navi_toolbar, 11, 1, 1, 9)
        self.navi_toolbar.setDisabled(True)

######## right, top:
        self.label_edge3a = QLabel('') #body
        self.label_edge3a.setMargin(5)
        self.label_edge3a.setLineWidth(1)
        self.label_edge3a.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        layout.addWidget(self.label_edge3a, 0, 11, 14, 6)
        self.label_edge3a.show()

        self.label_edge3b = QLabel('') #head
        self.label_edge3b.setMargin(5)
        self.label_edge3b.setLineWidth(1)
        self.label_edge3b.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        layout.addWidget(self.label_edge3b, 0, 11, 1, 6)
        self.label_edge3b.show()

        P1.label_coords = QLabel('Select Coordinates')
        P1.label_coords.setFont(font_standard)
        P1.label_coords.setMargin(5)
        P1.label_coords.setFrameStyle(QFrame.Box | QFrame.Sunken)
        P1.label_coords.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        P1.label_coords.setDisabled(True)
        layout.addWidget(P1.label_coords, 0, 11, 1, 6)
        P1.label_coords.show()

        P1.coords_x = {}
        P1.coords_y = {}
        P1.coords_e = {}
        P1.coords_check = {}
        ii = 0
        for i in range(10):
            ii += 1
            coord_x = QDoubleSpinBox(self)
            coord_y = QDoubleSpinBox(self)
            coord_e = QDoubleSpinBox(self)
            P1.coords_x[i] = coord_x
            P1.coords_y[i] = coord_y
            P1.coords_e[i] = coord_e
            coord_x.setDecimals(0)
            coord_y.setDecimals(0)
            coord_e.setDecimals(2)
            coord_e.setButtonSymbols(QAbstractSpinBox.NoButtons)
            coord_e.setSuffix('  kcal/mol')
            coord_x.setMinimum(0)
            coord_y.setMinimum(0)
            coord_x.setDisabled(True)
            coord_y.setDisabled(True)
            coord_e.setDisabled(True)
            layout.addWidget(coord_x, ii, 13, 1, 1)
            layout.addWidget(coord_y, ii, 14, 1, 1)
            layout.addWidget(coord_e, ii, 15, 1, 1)
            coord_check = QCheckBox('')
            coord_check.setDisabled(True)
            P1.coords_check[i] = coord_check
            layout.addWidget(coord_check, ii, 12, 1, 1)
            coord_check.show()
            
        P1.coords_x[0].valueChanged.connect(lambda: self.change_xy(0,'x'))
        P1.coords_x[1].valueChanged.connect(lambda: self.change_xy(1,'x'))
        P1.coords_x[2].valueChanged.connect(lambda: self.change_xy(2,'x'))
        P1.coords_x[3].valueChanged.connect(lambda: self.change_xy(3,'x'))
        P1.coords_x[4].valueChanged.connect(lambda: self.change_xy(4,'x'))
        P1.coords_x[5].valueChanged.connect(lambda: self.change_xy(5,'x'))
        P1.coords_x[6].valueChanged.connect(lambda: self.change_xy(6,'x'))
        P1.coords_x[7].valueChanged.connect(lambda: self.change_xy(7,'x'))
        P1.coords_x[8].valueChanged.connect(lambda: self.change_xy(8,'x'))
        P1.coords_x[9].valueChanged.connect(lambda: self.change_xy(9,'x'))
        P1.coords_y[0].valueChanged.connect(lambda: self.change_xy(0,'y'))
        P1.coords_y[1].valueChanged.connect(lambda: self.change_xy(1,'y'))
        P1.coords_y[2].valueChanged.connect(lambda: self.change_xy(2,'y'))
        P1.coords_y[3].valueChanged.connect(lambda: self.change_xy(3,'y'))
        P1.coords_y[4].valueChanged.connect(lambda: self.change_xy(4,'y'))
        P1.coords_y[5].valueChanged.connect(lambda: self.change_xy(5,'y'))
        P1.coords_y[6].valueChanged.connect(lambda: self.change_xy(6,'y'))
        P1.coords_y[7].valueChanged.connect(lambda: self.change_xy(7,'y'))
        P1.coords_y[8].valueChanged.connect(lambda: self.change_xy(8,'y'))
        P1.coords_y[9].valueChanged.connect(lambda: self.change_xy(9,'y'))
        P1.coords_check[0].stateChanged.connect(lambda: self.check_xy(0))
        P1.coords_check[1].stateChanged.connect(lambda: self.check_xy(1))
        P1.coords_check[2].stateChanged.connect(lambda: self.check_xy(2))
        P1.coords_check[3].stateChanged.connect(lambda: self.check_xy(3))
        P1.coords_check[4].stateChanged.connect(lambda: self.check_xy(4))
        P1.coords_check[5].stateChanged.connect(lambda: self.check_xy(5))
        P1.coords_check[6].stateChanged.connect(lambda: self.check_xy(6))
        P1.coords_check[7].stateChanged.connect(lambda: self.check_xy(7))
        P1.coords_check[8].stateChanged.connect(lambda: self.check_xy(8))
        P1.coords_check[9].stateChanged.connect(lambda: self.check_xy(9))

######## right, bot:

        P1.button_reset = QPushButton('   Reset Points   ', self)
        P1.button_reset.clicked.connect(self.coordsReset)
        P1.button_reset.setToolTip('Reset coordinates.')
        P1.button_reset.setDisabled(True)
        layout.addWidget(P1.button_reset, 11, 13, 1, 2)
        P1.button_reset.show()

        P1.button_calc = QPushButton('   Calculate Path   ', self)
        P1.button_calc.clicked.connect(self.leastAction)
        P1.button_calc.setToolTip('Calculate path of least action.')
        P1.button_calc.setDisabled(True)
        layout.addWidget(P1.button_calc, 11, 15, 1, 1)
        P1.button_calc.show()

    def check_xy(self, i):
        if P1.coords_check[i].isChecked():
            P1.coords_x[i].setDisabled(False)
            P1.coords_y[i].setDisabled(False)
        else:
            P1.coords_x[i].setDisabled(True)
            P1.coords_y[i].setDisabled(True)
            P1.coords_x[i].setValue(0)
            P1.coords_y[i].setValue(0)

    def change_xy(self, i, xy):
        if xy == 'x':
            P1.points_x[i] = int(P1.coords_x[i].value())
        if xy == 'y':
            P1.points_y[i] = int(P1.coords_y[i].value())
        P1.points_xy[i] = [P1.points_x[i], P1.points_y[i]]
        P1.coords_e[i].setValue(P1.LS[P1.points_y[i], P1.points_x[i]])
        self.plot_canvas.update_figure()

    def browseBox(self):
        P1.run += 1
        if P1.run > 1:
            msg = ("<span style='font-weight:normal;'>\
                    Performing this action will reset all information. \
                    <br /><br />\
                    Do you want to proceed?\
                    </span>")
            box = QMessageBox(self)
            box.setWindowTitle('%s Warning' % progname)
            box.setText('<b>Warning</b>')
            iconDir = os.path.join(pyDir, 'icons/70x70.png')
            if iconDir:
                box.setIconPixmap(QtGui.QPixmap(iconDir))
            box.setFont(font_standard)
            box.setIcon(QMessageBox.Warning)
            iconDir = os.path.join(pyDir, 'icons/70x70.png')
            if iconDir:
                box.setIconPixmap(QtGui.QPixmap(iconDir))
            box.setInformativeText(msg)
            box.setStandardButtons(QMessageBox.Yes|QMessageBox.No)
            reply = box.exec_()
            if reply == QMessageBox.Yes:
                self.browseFile()
            else:
                pass
        else:
            self.browseFile()

    def browseFile(self):
        if FigCanvas.initiated == 1:
            P1.df_prev = copy.copy(P1.df)
        P1.df = False
        P1.df = QtWidgets.QFileDialog.getOpenFileName(self, 'Choose Data File', '',
                                                     ('Data Files (*.csv *.txt)'))
        if P1.df[0]:
            try:
                tabs.setTabEnabled(1, True)
                P1.df = P1.df[0]
                try:
                    P1.LS = loadtxt(P1.df,float,delimiter=',')
                except ValueError:
                    P1.LS = loadtxt(P1.df,float)
                self.navi_toolbar.setDisabled(False)
                P2.p_list = np.zeros([5,1], dtype=int)
                for i in P2.nboxes:
                    if i == 4:
                        P2.nboxes[i].setValue(2)
                        P2.nPrboxes[i].setValue(nPr((4**(2)),i+1))
                    else:
                        P2.nboxes[i].setValue(1)
                        P2.nPrboxes[i].setValue(nPr((4**(1)),i+1))
                    P2.p_checks[i].setChecked(False)
                self.entry_LS.setDisabled(False)
                self.entry_LS.setText(P1.df)
                self.entry_LS.setDisabled(True)
                self.plot_canvas.update_figure()
                P1.label_coords.setDisabled(False)
                P1.button_reset.setDisabled(False)
                P1.button_calc.setDisabled(False)
                P2.processors.clearFocus()
                P2.processors.setValue(1)
                P2.TSW.clearFocus()
                #P2.TSW.setValue(1)
                P2.TSW.setChecked(False)

                for i in P1.coords_x:
                    P1.coords_x[i].setValue(0)
                    P1.coords_y[i].setValue(0)
                    P1.coords_check[i].setDisabled(False)
                    if i < 2:
                        P1.coords_x[i].setDisabled(False)
                        P1.coords_y[i].setDisabled(False)
                        P1.coords_check[i].setChecked(True)
                    else:
                        P1.coords_x[i].setDisabled(True)
                        P1.coords_y[i].setDisabled(True)
                        P1.coords_check[i].setChecked(False)

                for i in P1.coords_e:
                    P1.coords_e[i].setValue(np.amax(P1.LS[0,0]))
                df_dim = (int(np.shape(P1.LS)[0]), int(np.shape(P1.LS)[1]))
                for i in range(sys.maxsize):
                    if 2**i >= max(df_dim):
                        P1.upper_bound = 2**i
                        break
                for j in P1.coords_x:
                    P1.coords_x[j].setMaximum(df_dim[1]-1)
                for k in P1.coords_y:
                    P1.coords_y[k].setMaximum(df_dim[0]-1)
                for i in P2.nboxes:
                    P2.nboxes[i].setMaximum(np.log2(P1.upper_bound))
                    if i == 0:
                        P2.nboxes[i].setValue(np.log2(P1.upper_bound))
                        P2.nboxes[i].setDisabled(True)
                        P2.p_list[i] = np.log2(P1.upper_bound)
                        P2.p_checks[i].setChecked(True)

            except ValueError:
                box = QMessageBox(self)
                box.setWindowTitle('%s Error' % progname)
                box.setText('<b>Input Error</b>')
                iconDir = os.path.join(pyDir, 'icons/70x70.png')
                if iconDir:
                    box.setIconPixmap(QtGui.QPixmap(iconDir))
                box.setIcon(QMessageBox.Warning)
                box.setFont(font_standard)
                msg = "<span style='font-weight:normal;'>\
                        Incorrect file structure detected.\
                        </span>"
                box.setInformativeText(msg)
                box.setStandardButtons(QMessageBox.Ok)
                box.setDefaultButton(QMessageBox.Ok)
                ret = box.exec_()

                P1.df = P1.df_prev #bring back the previous

        else:
            P1.df = P1.df_prev #bring back the previous
            pass


            

    def coordsReset(self):
        for i in P1.coords_x:
            P1.coords_x[i].setValue(0)
            P1.coords_y[i].setValue(0)
            P1.coords_e[i].setValue(0)
            if i < 2:
                P1.coords_x[i].setDisabled(False)
                P1.coords_y[i].setDisabled(False)
                P1.coords_check[i].setChecked(True)
            else:
                P1.coords_x[i].setDisabled(True)
                P1.coords_y[i].setDisabled(True)
                P1.coords_check[i].setDisabled(False)
                P1.coords_check[i].setChecked(False)
                
        for i in P1.coords_e:
            P1.coords_e[i].setValue(np.amax(P1.LS[0,0]))

    def leastAction(self):
            
        final_xy = []
        flat_R = []
        final_N = []
        
        ii = 0
        for i in P1.points_xy:
            ii += 1
            if P1.coords_check[ii-1].isChecked():
                final_xy.append(i)
        flat_xy = ', '.join([str(i) for i in final_xy])

        ii = 0
        for i in P2.p_list:
            ii += 1
            if P2.p_checks[ii-1].isChecked():
                flat_R.append(ii)
                final_N.append(i)
        flat_N = [val for sublist in final_N for val in sublist]

        if (len(final_xy) > 1) and (len(final_N) > 0):
            
            final_xy = np.vstack(final_xy)
            
            box = QMessageBox(self)
            box.setWindowTitle('')
            box.setText('Confirm user inputs:<hr>')
            iconDir = os.path.join(pyDir, 'icons/70x70.png')
            if iconDir:
                box.setIconPixmap(QtGui.QPixmap(iconDir))
            box.setFont(font_standard)
            box.setInformativeText("<span style='font-weight:normal;'>\
                                    \
                                    <style>\
                                    p {\
                                        text-indent: 15px;\
                                    }\
                                    </style>\
                                    \
                                    <u>Coordinates:</u>\
                                    <p>\
                                        [x%s y%s] : <b>%s</b> </p>\
                                    <br /><br />\
                                    <u>Parameters:</u>\
                                    <p>\
                                        r :  <b>%s</b> </p>\
                                    <p>\
                                        n :  <b>%s</b> </p>\
                                    <br /><br />\
                                    <u>Constraints:</u>\
                                    <p>\
                                        Transition State Weighting: <b>%s</b> </p>\
                                    <p>\
                                        Number of Processors: <b>%s</b> </p>\
                                    <hr>\
                                    Would you like to proceed?\
                                    </span>\
                                    <br /><br />" % (u"\u1D62",u"\u1D62",
                                                     flat_xy,flat_R,flat_N,
                                                     P2.TSW.isChecked(),
                                                     int(P2.processors.value())))                              
            box.setStandardButtons(QMessageBox.Yes|QMessageBox.No)        
            reply = box.exec_()
            if reply == QMessageBox.Yes:
                
                POLARIS_BE.user_df = P1.df
                POLARIS_BE.user_pts = final_xy
                POLARIS_BE.user_R = flat_R
                POLARIS_BE.user_N = flat_N
                POLARIS_BE.user_proc = int(P2.processors.value())
                #POLARIS_BE.user_rate = int(P2.TSW.value())
                POLARIS_BE.user_rate = P2.TSW.isChecked()
                POLARIS_BE.init()

                # popup message upon completion:
                box = QMessageBox(self)
                box.setWindowTitle('')
                box.setText('<b>Path Complete</b>')
                iconDir = os.path.join(pyDir, 'icons/70x70.png')
                if iconDir:
                    box.setIconPixmap(QtGui.QPixmap(iconDir))
                box.setFont(font_standard)
                box.setInformativeText("<span style='font-weight:normal;'>\
                                        Please see project directory for a list\
                                        of coordinates and plots of the\
                                        trajectory through the energy landscape.\
                                        </span>")
                box.setStandardButtons(QMessageBox.Ok)        
                ret = box.exec_()

            else:
                pass
        elif len(final_xy) <= 1:
            box = QMessageBox(self)
            box.setWindowTitle('%s Error' % progname)
            box.setText('<b>Input Error</b>')
            iconDir = os.path.join(pyDir, 'icons/70x70.png')
            if iconDir:
                box.setIconPixmap(QtGui.QPixmap(iconDir))
            box.setFont(font_standard)
            box.setIcon(QMessageBox.Warning)
            box.setInformativeText("<span style='font-weight:normal;'>\
                                    At least two points must be selected\
                                    to perform this calculation.\
                                    <br /><br />\
                                    To proceed, first pick a minimum of two\
                                    unique points on the energy landscape\
                                    via the Coordinates tab.\
                                    </span>")
            box.setStandardButtons(QMessageBox.Ok)        
            ret = box.exec_()
        elif len(final_N) == 0:
            box = QMessageBox(self)
            box.setWindowTitle('%s Error' % progname)
            box.setText('<b>Input Error</b>')
            iconDir = os.path.join(pyDir, 'icons/70x70.png')
            if iconDir:
                box.setIconPixmap(QtGui.QPixmap(iconDir))
            box.setFont(font_standard)
            box.setIcon(QMessageBox.Warning)
            box.setInformativeText("<span style='font-weight:normal;'>\
                                    At least one parameter must be selected\
                                    to perform this calculation.\
                                    <br /><br />\
                                    To proceed, first pick one N value for\
                                    a minimum of one P value via the\
                                    Settings tab.\
                                    </span>")
            box.setStandardButtons(QMessageBox.Ok)        
            ret = box.exec_()


################################################################################
# settings (tab 2): 

class P2(QtWidgets.QWidget):

    nboxes = {}
    nPrboxes = {}
    p_checks = {}
    p_list = {}
    
    def __init__(self, parent=None):
        super(P2, self).__init__(parent)
        layout = QGridLayout(self)
        layout.setContentsMargins(20,20,20,20)
        layout.setSpacing(10)

        P2.p_list = np.zeros([5,1], dtype=int)

######## left:
        self.label_edge1a = QLabel('') #largest bin
        self.label_edge1a.setMargin(5)
        self.label_edge1a.setLineWidth(1)
        self.label_edge1a.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        layout.addWidget(self.label_edge1a, 0, 0, 11, 5)
        self.label_edge1a.show()

        self.label_edge1b = QLabel('') #header bin
        self.label_edge1b.setMargin(5)
        self.label_edge1b.setLineWidth(1)
        self.label_edge1b.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        layout.addWidget(self.label_edge1b, 0, 0, 1, 5)
        self.label_edge1b.show()

        self.label_settings = QLabel('Parameters')
        self.label_settings.setFont(font_standard)
        self.label_settings.setMargin(5)
        self.label_settings.setFrameStyle(QFrame.Box | QFrame.Sunken)
        self.label_settings.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        layout.addWidget(self.label_settings, 0, 0, 1, 5)
        self.label_settings.show()

        ii = 0
        for i in range(5):
            ii += 1

            rbox = QDoubleSpinBox(self)
            rbox.setButtonSymbols(QAbstractSpinBox.NoButtons)
            rbox.setPrefix('r: ')
            rbox.setDecimals(0)
            rbox.setDisabled(True)
            rbox.setValue(ii)
            layout.addWidget(rbox, 2*ii-1, 1, 1, 1)
            
            nbox = QDoubleSpinBox(self)
            P2.nboxes[i] = nbox
            nbox.setDecimals(0)
            if i == 4: #r5
                nbox.setMinimum(2) #r5, n1 -> error
            else:
                nbox.setMinimum(1)
            nbox.setPrefix('n: ')
            nbox.setDisabled(False)
            layout.addWidget(nbox, 2*ii-1, 2, 1, 1)

            nPrbox = QDoubleSpinBox(self)
            P2.nPrboxes[i] = nPrbox
            nPrbox.setButtonSymbols(QAbstractSpinBox.NoButtons)
            nPrbox.setPrefix('P(4%s, r): ' % (u"\u207F"))
            nPrbox.setDecimals(0)
            nPrbox.setMaximum(10**12)
            nPrbox.setDisabled(True)
            layout.addWidget(nPrbox, 2*ii-1, 3, 1, 1)

            p_check = QCheckBox('')
            p_check.setDisabled(False)
            P2.p_checks[i] = p_check
            P2.p_checks[i].setChecked(False)
            layout.addWidget(p_check, 2*ii-1, 4, 1, 1)
            p_check.show()

            hline = QLabel("")
            hline.setFrameStyle(QFrame.HLine | QFrame.Sunken)
            if i < 4:
                layout.addWidget(hline, 2*ii, 1, 1, 3, QtCore.Qt.AlignVCenter)

        P2.nboxes[0].valueChanged.connect(lambda: self.change_nPr(0))
        P2.nboxes[1].valueChanged.connect(lambda: self.change_nPr(1))
        P2.nboxes[2].valueChanged.connect(lambda: self.change_nPr(2))
        P2.nboxes[3].valueChanged.connect(lambda: self.change_nPr(3))
        P2.nboxes[4].valueChanged.connect(lambda: self.change_nPr(4))

        P2.p_checks[0].setChecked(True)
        P2.p_checks[0].stateChanged.connect(lambda: self.change_p(0))
        P2.p_checks[1].stateChanged.connect(lambda: self.change_p(1))
        P2.p_checks[2].stateChanged.connect(lambda: self.change_p(2))
        P2.p_checks[3].stateChanged.connect(lambda: self.change_p(3))
        P2.p_checks[4].stateChanged.connect(lambda: self.change_p(4))

######## right:
        
        self.label_edgeR = QLabel('') #largest bin
        self.label_edgeR.setMargin(5)
        self.label_edgeR.setLineWidth(1)
        self.label_edgeR.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        layout.addWidget(self.label_edgeR, 11, 0, 3, 5)
        self.label_edgeR.show()
        
        self.label_edge1b = QLabel('') #header bin
        self.label_edge1b.setMargin(5)
        self.label_edge1b.setLineWidth(1)
        self.label_edge1b.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        layout.addWidget(self.label_edge1b, 11, 0, 1, 5)
        self.label_edge1b.show()
        
        self.label_settings = QLabel('Constraints')
        self.label_settings.setFont(font_standard)
        self.label_settings.setMargin(5)
        self.label_settings.setFrameStyle(QFrame.Box | QFrame.Sunken)
        self.label_settings.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        layout.addWidget(self.label_settings, 11, 0, 1, 5)
        self.label_settings.show()

        P2.processors = QDoubleSpinBox(self)
        P2.processors.setFocus(False)
        P2.processors.setPrefix('Processors: ')
        P2.processors.setFont(font_standard)
        P2.processors.setDecimals(0)
        P2.processors.setMinimum(1)
        P2.processors.setMaximum(multiprocessing.cpu_count())
        layout.addWidget(P2.processors, 12, 3, 1, 1)
        P2.processors.show()
        '''
        P2.TSW = QDoubleSpinBox(self)
        P2.TSW.setFocus(False)
        P2.TSW.setPrefix('Transition State Weighting: ')
        P2.TSW.setFont(font_standard)
        P2.TSW.setDecimals(0)
        P2.TSW.setMinimum(1)
        P2.TSW.setMaximum(2)
        layout.addWidget(P2.TSW, 12, 1, 1, 1)
        P2.TSW.show()
        '''
        P2.TSW = QCheckBox('Transition State Weighting')
        P2.TSW.setFont(font_standard)
        P2.TSW.setChecked(False)
        #P2.TSW.setEnabled(False)
        layout.addWidget(P2.TSW, 12, 1, 1, 1)
        P2.TSW.show()

    def change_p(self, i):
        if P2.p_checks[i].isChecked():
            P2.nboxes[i].setDisabled(True)
            P2.p_list[i] = P2.nboxes[i].value()
        else:
            P2.nboxes[i].setDisabled(False)
            P2.p_list[i] = 0

    def change_nPr(self, i):
        P2.nPrboxes[i].setValue(nPr(4**(int(P2.nboxes[i].value())), i+1))

def nPr(n, r): #calculate total number of permutations
    total = int(factorial(n)/factorial(n-r))
    if total > 10**12:
        return 999999999999 #maxed out, i.e., never going to use
    else:
        return total
                        

################################################################################
# class hierarchy: 

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        #self.showMaximized()

        style = """QTabWidget::tab-bar{
            alignment: center;
        }"""

        tab1 = P1(self)
        tab2 = P2(self)
        global tabs
        tabs = QtWidgets.QTabWidget(self)
        tabs.addTab(tab1, '       Coordinates       ')
        tabs.addTab(tab2, '        Settings        ')
        tabs.setTabEnabled(1, False)
        self.setStyleSheet(style)
        self.setCentralWidget(tabs)
        self.showMaximized()

        mainMenu = self.menuBar()
        mainMenu.setNativeMenuBar(False)
        fileMenu = mainMenu.addMenu('&File')
        fileMenu.addAction('&About', self.about)
        fileMenu.addSeparator()
        fileMenu.addAction('&Restart', self.fileRestart)
        helpMenu = mainMenu.addMenu('&Help')
        helpMenu.addAction('&Coordinates', self.guide_coords)

        paraMenu = helpMenu.addMenu('&Settings')
        paraMenu.addAction('&Parameters', self.guide_params1)
        paraMenu.addAction('&Constraints', self.guide_params2)        

    def gotoP2(self):        
        tabs.setCurrentIndex(1)
        
    def closeEvent(self, ce): #safety message if user clicks to exit via window button
        msg = "<span style='font-weight:normal;'>\
               Performing this action will close the program.\
               All progress will be lost.\
               <br /><br />\
               Do you want to proceed?\
               </span>"
        box = QMessageBox(self)
        box.setWindowTitle('%s Warning' % progname)
        box.setText('<b>Exit Warning</b>')
        iconDir = os.path.join(pyDir, 'icons/70x70.png')
        if iconDir:
            box.setIconPixmap(QtGui.QPixmap(iconDir))
        box.setFont(font_standard)
        box.setIcon(QMessageBox.Warning)
        box.setInformativeText(msg)
        box.setStandardButtons(QMessageBox.Yes|QMessageBox.No)
        reply = box.exec_()
        if reply == QMessageBox.Yes:
            self.close()
        else:
            ce.ignore()

    def fileRestart(self):
        msg = "<span style='font-weight:normal;'>\
               Performing this action will restart the program,\
               resetting all user inputs.\
               <br /><br />\
               Do you want to proceed?\
               </span>"
        box = QMessageBox(self)
        box.setWindowTitle('%s Warning' % progname)
        box.setText('<b>Restart Warning</b>')
        iconDir = os.path.join(pyDir, 'icons/70x70.png')
        if iconDir:
            box.setIconPixmap(QtGui.QPixmap(iconDir))
        box.setFont(font_standard)
        box.setIcon(QMessageBox.Warning)
        box.setInformativeText(msg)
        box.setStandardButtons(QMessageBox.Yes|QMessageBox.No)
        reply = box.exec_()
        if reply == QMessageBox.Yes:
            try:
                p = psutil.Process(os.getpid())
                for handler in p.open_files() + p.connections():
                    os.close(handler.fd)
            except Exception as e:
                logging.error(e)

            python = sys.executable
            os.execl(python, python, * sys.argv)
        else:
            pass

    def guide_coords(self):
        box = QMessageBox(self)
        box.setWindowTitle('%s Guide' % progname)
        box.setText('<b>Navigating the coordinates tab</b>')
        iconDir = os.path.join(pyDir, 'icons/70x70.png')
        if iconDir:
            box.setIconPixmap(QtGui.QPixmap(iconDir))
        box.setFont(font_standard)
        box.setInformativeText("Energy Landscape:\
                                <span style='font-weight:normal;'>\
                                select a valid file to plot within the energy\
                                landscape window via the Browse button.\
                                Once an accepted file has been loaded,\
                                use the navigation toolbar below the plotted\
                                figure to zoom in and pan across the\
                                landscape (amongst other functions).\
                                <br /><br />\
                                <b>Select Coordinates:</b>\
                                Click the box next to each row to enable\
                                selection of individual points on the energy\
                                landscape. As integers are typed in or scrolled\
                                to via the up and down arrows, that row's respective\
                                point will update its location in the plotting\
                                window. A minimum of two points must be\
                                selected for computation to occur.\
                                <br /><br />\
                                POLARIS will calculate the path of\
                                least action between these chosen points in the\
                                order in which they appear, from the topmost row\
                                down to the bottom. Up to 10 points can be selected\
                                as intermediate transit locations, with unchecked\
                                rows ignored. If a complete cycle is preferred\
                                (closed loop), make sure to include the first point\
                                again as the last.\
                                </span>")
        box.setStandardButtons(QMessageBox.Ok)        
        ret = box.exec_()

    def guide_params1(self):
        box = QMessageBox(self)
        box.setWindowTitle('%s Guide' % progname)
        box.setText('<b>Navigating the settings tab</b>')
        iconDir = os.path.join(pyDir, 'icons/70x70.png')
        if iconDir:
            box.setIconPixmap(QtGui.QPixmap(iconDir))
        box.setFont(font_standard)
        box.setInformativeText('Permutational Order (r): \
                                <span style="font-weight:normal;">\
                                the number of intermediate minima nodes to traverse\
                                through between each pair of start and end coordinates\
                                during comparison of line approximations within a\
                                given permutation pool.\
                                <br /><br />\
                                <b>Segmentation Depth (n):</b>\
                                the number of minima to include as traversal options\
                                within the given permutation order, defined via the\
                                number of image subdivisions (4%s) performed.\
                                <br /><br />\
                                Select a segmentation level n for each permutational\
                                order r, with corresponding line permutations computed\
                                via P(4%s, r). The highest value of n is automatically\
                                calculated from the datafile and selected for the\
                                lowest r, such that for any two start and end point\
                                combinations, all points in the landscape will be\
                                tested as indepedent, singular midpoint options once,\
                                with the corresponding line approximation from each\
                                subsequently compared before moving on to the next\
                                highest value of r.\
                                </span>' % (u"\u207F", u"\u207F"))
        box.setStandardButtons(QMessageBox.Ok)        
        ret = box.exec_()

    def guide_params2(self):
        box = QMessageBox(self)
        box.setWindowTitle('%s Guide' % progname)
        box.setText('<b>Navigating the settings tab</b>')
        iconDir = os.path.join(pyDir, 'icons/70x70.png')
        if iconDir:
            box.setIconPixmap(QtGui.QPixmap(iconDir))
        box.setFont(font_standard)
        box.setInformativeText('<span style="font-weight:normal;">\
                                <b>Transition State Weighting:</b>\
                                option to compare competing lowest energy path\
                                bifurcations based on their rate limiting step\
                                (points of maximal energy that the object must\
                                pass through) instead of by their net integrated\
                                energies. In theory, the rate of any reaction\
                                is dependent on the highest energy transition state\
                                (i.e., the highest energy point along the lowest\
                                energy path), and not the net energy of the path.\
                                </span>')
        box.setStandardButtons(QMessageBox.Ok)        
        ret = box.exec_()

    def about(self):
        box = QMessageBox(self)
        box.setWindowTitle('%s About' % progname)
        box.setText('<b>POLARIS v.%s Python 3.x</b>' % (progversion))
        box.setFont(font_standard)
        iconDir = os.path.join(pyDir, 'icons/70x70.png')
        if iconDir:
            box.setIconPixmap(QtGui.QPixmap(iconDir))
        box.setInformativeText('<span style="font-weight:normal;">\
                                %s Evan Elliott Seitz, 2018-2020\
                                <br /><br />\
                                <b>LICENSE:</b>\
                                <br /><br />\
                                POLARIS is a free software: you can redistribute it and/or\
                                modify it under the terms of the GNU General Public License\
                                as published by the Free Software Foundation, either version\
                                3 of the License, or (at your option) any later version.\
                                <br /><br />\
                                This program is distributed in the hope that it will be useful\
                                but WITHOUT ANY WARRANTY; without even the implied warranty of\
                                MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\
                                See the GNU General Public License for more details.\
                                <br /><br />\
                                You should have received a copy of the GNU General Public\
                                License along with this program, which can be found at\
                                http://www.gnu.org/licenses.\
                                <br /><br />\
                                <b>Contact:</b>\
                                evan.e.seitz@gmail.com\
                                <br /><br />\
                                <b>DOI:</b>\
                                10.1021/acs.jcim.9b01108\
                                </span>' % (u"\u00A9"))
        box.setStandardButtons(QMessageBox.Ok)        
        ret = box.exec_()

if __name__ == '__main__':
    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)

    QtCore.QCoreApplication.setApplicationName(progname) #if non-macOS menu bar
        
    # set app icon for tray:
    iconDir = os.path.join(pyDir, 'icons')
    app_icon = QtGui.QIcon()
    app_icon.addFile(os.path.join(iconDir, '256x256.png'), QtCore.QSize(256,256))
    app.setWindowIcon(app_icon)
        
    w = MainWindow()
    w.setWindowTitle('%s' % progname)
    w.show()
    POLARIS_BE.user()
    sys.exit(app.exec_())
