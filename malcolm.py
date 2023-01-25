import sys
import signal
import random
import matplotlib
import numpy as np
import math
from pathlib import Path, PurePath
from time import time, sleep
import os
# to find the number of cores
import multiprocessing
import shutil
import pandas as pd

from nodes import node
from tasks import Task_status

matplotlib.use('Qt5Agg')

from PySide2.QtWidgets import QMainWindow, QApplication, QHBoxLayout, QVBoxLayout, QGridLayout, QWidget, QLabel, QPushButton, QLineEdit, QMessageBox, QDesktopWidget, QTabWidget, QDoubleSpinBox, QSpinBox, QComboBox, QCheckBox, QPlainTextEdit, QProgressBar, QFileDialog, QStyle, QSplashScreen


from PySide2.QtCore import Qt, QTimer

from PySide2.QtCore import QObject, QThread, Signal
from PySide2.QtGui import QGuiApplication, QPixmap, QIcon

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from itertools import cycle


import warnings
warnings.simplefilter('always', UserWarning)

import logging

from globvars import globvars

from rebars import rebar, hoop, spiral 
from cross_sections import cross_section
from elements import element
from materials import concrete_mat, rebar_mat, tendon_mat, window_CDPM2, window_Mises
from interaction_diagrams import id_ACI, id_MC2010, id_MSR, id_MSR_simple

from concrete import concrete

import oofem_input

# https://www.pythonguis.com/tutorials/pyside-creating-your-first-window/


class Worker(QObject):

    finished = Signal()
    progress = Signal(object, object)
    result = Signal(object, object, object)

    def __init__(self, task, oofem_folder, n_cores):
        super().__init__()

    
        self.task = task
        self.load = []
        self.dummy_eps = []
        self.dummy_kappa = []
        
        self.oofem_folder = oofem_folder
        self.n_cores = n_cores
 
    #@Slot()  # QtCore.Slot
    def run(self):


        current_path = os.getcwd()

        problem_path = self.task.file_path
        os.chdir(problem_path) 
        sys.path.append(self.oofem_folder)
        # the import needs to be in the "run" function, the imports from the above do not work
        import oofempy
        #oofempy.init(logLevel=3, numberOfThreads=self.n_cores)
        oofempy.init(logLevel=2, numberOfThreads=self.n_cores) 
        self.task.status = Task_status.PROGRESS
        # TODO: task should have its name
        #dr = oofempy.OOFEMTXTDataReader(problem_path)
        dr = oofempy.OOFEMTXTDataReader("malcolm_oofem.in")
        problem = oofempy.InstanciateProblem(dr, oofempy.problemMode.processor, False, None, False)
        domain = problem.giveDomain(1)

        problem.init()
        problem.checkProblemConsistency()

        problem.giveTimer().startTimer(oofempy.EngngModelTimerType.EMTT_AnalysisTimer)
        activeMStep = problem.giveMetaStep(1)
        problem.initMetaStepAttributes(activeMStep);


        load_level = 0.
        max_load_level = 0.

        n_steps = problem.giveNumberOfSteps()
        
        for timeStep in range(n_steps):
            if (load_level == max_load_level):

                problem.preInitializeNextStep()
                problem.giveNextStep()
                currentStep = problem.giveCurrentStep()
                problem.initializeYourself( currentStep )
                problem.solveYourselfAt( currentStep )
                problem.updateYourself( currentStep )

                w_eps = problem.giveUnknownComponent (oofempy.ValueModeType.VM_Total, problem.giveCurrentStep(False), domain, domain.giveGlobalDofManager(1000001).giveDofWithID(oofempy.DofIDItem.D_w))
                w_kappa = problem.giveUnknownComponent (oofempy.ValueModeType.VM_Total, problem.giveCurrentStep(False), domain, domain.giveGlobalDofManager(1000002).giveDofWithID(oofempy.DofIDItem.D_w))
                load_level = problem.giveLoadLevel()
                
                if (load_level > max_load_level):
                    max_load_level = load_level
                
                problem.terminate( currentStep )
                self.progress.emit(timeStep, load_level)
                self.load.append(load_level)
                self.dummy_eps.append(w_eps)
                self.dummy_kappa.append(w_kappa)

            else:
                break

           
        problem.terminateAnalysis()

        self.task.load_level = self.load
        self.task.max_load = max_load_level
        
        self.result.emit(self.task, self.dummy_eps, self.dummy_kappa)
        
        #print ("FINISHED TASK")
        self.task.status = Task_status.COMPLETED
        self.finished.emit()
        
        os.chdir(current_path)
        

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

        
class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.w_cdpm2 = None
        
        self.globvar = globvars()

        self.setWindowTitle(self.globvar.program_name)
        
        
        # all ids to be updated
        self.need_update_ids_flag = True

        self.diagram_ACI = id_ACI(self.globvar)
        self.diagram_MC2010 = id_MC2010(self.globvar)
        self.diagram_MSR = id_MSR(self.globvar)
        self.diagram_MSR_simple = id_MSR_simple(self.globvar)
        
        # init info labels
        self.label_info_topology = QLabel()
        
        # concrete entity
        self.concrete = concrete(globvar = self.globvar, cs = 1, mat = 1)

        # MATERIALS
        # material: CONCRETE 
        mat =  concrete_mat(globvar = self.globvar, nr = 1)
        self.globvar.materials.append(mat) 

        #  material: LATERAL REINFORCEMENT 
        mat =  rebar_mat(globvar = self.globvar, nr = 2, mat = self.globvar.rebarmat[0])
        self.globvar.materials.append(mat) 

        #  material: VERTICAL REINFORCEMENT 
        mat =  rebar_mat(globvar = self.globvar, nr = 3, mat = self.globvar.rebarmat[0])
        self.globvar.materials.append(mat) 
        
        self.cdpm2 = self.globvar.materials[0]
        self.mises_lat = self.globvar.materials[1]
        self.mises_vert = self.globvar.materials[2]

       
        ### LATERAL REINFORCEMENT - topology
        self.vec = [0., 0., 1.]

        rebar_DS = (self.globvar.rebars_CS[self.globvar.DS]).diam
        rebar_DL = (self.globvar.rebars_CS[self.globvar.DL]).diam
        
        # positions of small spirals to fit at the corners
        auxX = self.globvar.Bx/2. - self.globvar.cover - self.globvar.dS/2 - rebar_DS/2.
        auxY = self.globvar.By/2. - self.globvar.cover - self.globvar.dS/2 - rebar_DS/2.
        
        ### large spiral
        center  = [0., 0., 0.]
        start = [ center[0], center[1] - self.globvar.dL/2., center[2] ]

        spiral_L = spiral(globvar = self.globvar, XYZ = center, cs_nr = 2, cs_tag = self.globvar.DL, cs_mat = 2, vec = self.vec, pitch = self.globvar.H, axial_length = self.globvar.H, start = start, tag = 'L1')
        #spiral_L.check_topology_consistency(self.globvar)
        self.globvar.rebars.append(spiral_L)

        ### small spirals:
        # S1
        center = [auxX, -auxY, 0.]
        start = [ center[0], center[1] - self.globvar.dS/2., center[2] ]
            
        spiral_S = spiral(globvar = self.globvar, XYZ = center, cs_nr = 3, cs_tag = self.globvar.DS, cs_mat = 2, vec = self.vec, pitch = self.globvar.H, axial_length = self.globvar.H, start = start, tag = 'S1')
        #spiral_S.check_topology_consistency(self.globvar)
        self.globvar.rebars.append(spiral_S)

        # S2
        center = [auxX, auxY, 0.]
        start = [ center[0], center[1] - self.globvar.dS/2., center[2] ]
            
        spiral_S = spiral(globvar = self.globvar, XYZ = center, cs_nr = 3, cs_tag = self.globvar.DS, cs_mat = 2, vec = self.vec, pitch = self.globvar.H, axial_length = self.globvar.H, start = start, tag = 'S2')
        #spiral_S.check_topology_consistency(self.globvar)
        self.globvar.rebars.append(spiral_S)

        # S3
        center = [-auxX, auxY, 0.]
        start = [ center[0], center[1] - self.globvar.dS/2., center[2] ]
            
        spiral_S = spiral(globvar = self.globvar, XYZ = center, cs_nr = 3, cs_tag = self.globvar.DS, cs_mat = 2, vec = self.vec, pitch = self.globvar.H, axial_length = self.globvar.H, start = start, tag = 'S3')
        #spiral_S.check_topology_consistency(self.globvar)
        self.globvar.rebars.append(spiral_S)

        # S4
        center = [-auxX, -auxY, 0.]
        start = [ center[0], center[1] - self.globvar.dS/2., center[2] ]
            
        spiral_S = spiral(globvar = self.globvar, XYZ = center, cs_nr = 3, cs_tag = self.globvar.DS, cs_mat = 2, vec = self.vec, pitch = self.globvar.H, axial_length = self.globvar.H, start = start, tag = 'S4')
        #spiral_S.check_topology_consistency(self.globvar)
        self.globvar.rebars.append(spiral_S)
                      

        ## VERTICAL REINFORCEMENT
           
        # generation of vertical reinforcement for the first time
        self.update_vertical_rebars()
        # compute MSR intersection for the first time
        self.diagram_MSR.update_MSR_intersections(self.globvar)
        

        # in a new dialog windows for simplicity always define bc and related tfs with the same number
        # TODO
        # boundary conditions
        self.globvar.bcs.append("BoundaryCondition  1 loadTimeFunction 1 prescribedvalue 0.")
        self.globvar.bcs.append("NodalLoad 2 loadTimeFunction 2 Components 3 0. 0. -1.")

        # TODO
        # load time functions
        self.globvar.ltfs.append("ConstantFunction 1 f(t) 1.0")
        self.globvar.ltfs.append("ConstantFunction 2 f(t) 1.0")

        #########################
        # TOPOLOGY etc.
        #########################

        # column breadth
        self.widget_Bx = QDoubleSpinBox()
        self.widget_Bx.setMinimum(0.2)
        self.widget_Bx.setMaximum(1.0)
        self.widget_Bx.setDecimals(3)
        self.widget_Bx.setValue( self.globvar.Bx )
                
        self.widget_Bx.setSingleStep(50./1000.)  
        self.widget_Bx.valueChanged.connect(self.set_Bx)
        self.widget_Bx.valueChanged.connect(self.update_topology_plot)

        # spiral pitch
        self.widget_H = QDoubleSpinBox()
        self.widget_H.setMinimum(0.03)
        self.widget_H.setMaximum(0.2)
        self.widget_H.setDecimals(3)
        self.widget_H.setValue( self.globvar.H )
                
        self.widget_H.setSingleStep(10./1000.)  
        self.widget_H.valueChanged.connect(self.set_H)
        self.widget_H.valueChanged.connect(self.update_topology_plot)
        

        # axial diameter of small spiral
        self.widget_cover = QDoubleSpinBox()
        self.widget_cover.setMinimum(0.)
        self.widget_cover.setMaximum(0.1)
        self.widget_cover.setDecimals(3)
        self.widget_cover.setValue( self.globvar.cover )
                
        self.widget_cover.setSingleStep(5./1000.)  
        self.widget_cover.valueChanged.connect(self.set_cover)
        self.widget_cover.valueChanged.connect(self.update_topology_plot)

        
        # axial diameter of small spiral
        self.widget_dS = QDoubleSpinBox()

        current_DS = (self.globvar.rebars_CS[self.globvar.DS]).diam
        
        self.widget_dS.setMinimum(self.globvar.minimum_radius*2.)
        self.widget_dS.setMaximum( min(self.globvar.Bx, self.globvar.By )/2. - self.globvar.cover - current_DS )
        self.widget_dS.setValue( self.globvar.dS )
        
        self.widget_dS.setSingleStep(10./1000.) 
        self.widget_dS.valueChanged.connect(self.set_dS)
        self.widget_dS.valueChanged.connect(self.update_topology_plot)

        # rebar diameter of small spirals
        self.widget_DS = QComboBox()
        self.widget_DS.addItems(self.globvar.rebars_CS.keys())
        self.widget_DS.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DS))
        self.widget_DS.currentTextChanged.connect(self.set_DS)
        self.widget_DS.currentTextChanged.connect(self.update_topology_plot)          

        # rebar diameter of large spiral(s)
        self.widget_DL = QComboBox()
        self.widget_DL.addItems(self.globvar.rebars_CS.keys())
        self.widget_DL.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DL))
        self.widget_DL.currentTextChanged.connect(self.set_DL)
        self.widget_DL.currentTextChanged.connect(self.update_topology_plot)

        # rebar diameter of vertical rebars
        self.widget_DV = QComboBox()
        self.widget_DV.addItems(self.globvar.rebars_CS.keys())
        self.widget_DV.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DV))
        self.widget_DV.currentTextChanged.connect(self.set_DV)
        self.widget_DV.currentTextChanged.connect(self.update_topology_plot)   


        #########################
        # TOPOLOGY & DEFINITION CANVAS
        #########################
       
        self.topology_canvas = MplCanvas(self, width=5, height=4, dpi=100)
            
        topology_specs_layout = QGridLayout()

        row = 0
        topology_specs_layout.addWidget( QLabel("Geometry"), row, 0, 1, 3)
        
        row += 1
        label_Bx = QLabel("Bx = ")
        label_Bx_unit = QLabel("[m]")
        topology_specs_layout.addWidget(label_Bx, row, 0, 1, 1, Qt.AlignRight)
        topology_specs_layout.addWidget(self.widget_Bx, row, 1)
        topology_specs_layout.addWidget(label_Bx_unit, row, 2)

        row += 1
        label_cover = QLabel("c = ")
        label_cover_unit = QLabel("[m]")
        topology_specs_layout.addWidget(label_cover, row, 0, 1, 1, Qt.AlignRight)
        topology_specs_layout.addWidget(self.widget_cover, row, 1)
        topology_specs_layout.addWidget(label_cover_unit, row, 2)
        
        row += 1
        label_H = QLabel("H = ")
        label_H_unit = QLabel("[m]")
        topology_specs_layout.addWidget(label_H, row, 0, 1, 1, Qt.AlignRight)
        topology_specs_layout.addWidget(self.widget_H, row, 1)
        topology_specs_layout.addWidget(label_H_unit, row, 2)
     
        row += 1
        label_dS = QLabel("dS = ")
        label_dS_unit = QLabel("[m]")
        topology_specs_layout.addWidget(label_dS, row, 0, 1, 1, Qt.AlignRight)
        topology_specs_layout.addWidget(self.widget_dS, row, 1)
        topology_specs_layout.addWidget(label_dS_unit, row, 2)


        row = 0
        topology_specs_layout.addWidget( QLabel("Reinforcement"), row, 3, 1, 3)
        
        row += 1
        label_DS = QLabel("DS = ")
        label_DS_unit = QLabel("[m]")
        topology_specs_layout.addWidget(label_DS, row, 3, 1, 1, Qt.AlignRight)
        topology_specs_layout.addWidget(self.widget_DS, row, 4)
        topology_specs_layout.addWidget(label_DS_unit, row, 5)
                
        row += 1
        label_DL = QLabel("DL = ")
        label_DL_unit = QLabel("[m]")
        topology_specs_layout.addWidget(label_DL, row, 3, 1, 1, Qt.AlignRight)
        topology_specs_layout.addWidget(self.widget_DL, row, 4)
        topology_specs_layout.addWidget(label_DL_unit, row, 5)
            
        row += 1
        label_DV = QLabel("DV = ")
        label_DV_unit = QLabel("[m]")
        topology_specs_layout.addWidget(label_DV, row, 3, 1, 1, Qt.AlignRight)
        topology_specs_layout.addWidget(self.widget_DV, row, 4)
        topology_specs_layout.addWidget(label_DV_unit, row, 5)

        row = 5
        row += 1
        topology_specs_layout.addWidget(self.label_info_topology, row, 0, 1, 6, Qt.AlignLeft)


        row += 1
        lineEditLoading = QLineEdit()
        lineEditLoading.setMaxLength(1000)
        lineEditLoading.setPlaceholderText("specify loading file")
        lineEditLoading.textChanged.connect(self.set_loading_file)
        
        button_loading_file = QPushButton("Loading [M,N]")
        button_loading_file.clicked.connect(self.select_loading_file)
        button_loading_file.clicked.connect( lambda: lineEditLoading.setText(self.globvar.loading_file) )
        
        topology_specs_layout.addWidget(button_loading_file, row, 0, 1, 2)
        topology_specs_layout.addWidget(lineEditLoading, row, 2, 1, 4)

        #########################
        # INTERACTION DIAGRAMS etc.
        #########################
        
        ### basic material properties
        # concrete
        self.widget_fcm = QComboBox()
        self.widget_fcm.addItems(self.globvar.concretes.keys())
        self.widget_fcm.setCurrentIndex(list(self.globvar.concretes.keys()).index(self.globvar.fcm))
        self.widget_fcm.currentTextChanged.connect(self.set_fcm)
        self.widget_fcm.currentTextChanged.connect(self.update_diagram_plot)

        # steel - vertical reinforcement
        self.widget_fy_vert = QDoubleSpinBox()
        self.widget_fy_vert.setMinimum(200)
        self.widget_fy_vert.setMaximum(700)
        self.widget_fy_vert.setDecimals(0)
        self.widget_fy_vert.setValue( self.globvar.fy_vert )
                
        self.widget_fy_vert.setSingleStep(20.)  
        self.widget_fy_vert.valueChanged.connect(self.set_fy_vert)
        self.widget_fy_vert.valueChanged.connect(self.update_diagram_plot)
        
        # steel - lateral reinforcement
        self.widget_fy_lat = QDoubleSpinBox()
        self.widget_fy_lat.setMinimum(200)
        self.widget_fy_lat.setMaximum(700)
        self.widget_fy_lat.setDecimals(0)
        self.widget_fy_lat.setValue( self.globvar.fy_lat )
                
        self.widget_fy_lat.setSingleStep(20.)  
        self.widget_fy_lat.valueChanged.connect(self.set_fy_lat)
        self.widget_fy_lat.valueChanged.connect(self.update_diagram_plot)

        ### checkboxes for interaction diagrams
        checkbox_id_steel = QCheckBox()
        checkbox_id_steel.setCheckState(self.globvar.flag_id_steel)
        checkbox_id_steel.stateChanged.connect(self.set_id_steel)
        checkbox_id_steel.stateChanged.connect(self.update_diagram_plot)

        checkbox_id_ACI_concrete = QCheckBox()
        checkbox_id_ACI_concrete.setCheckState(self.globvar.flag_id_ACI_concrete)
        checkbox_id_ACI_concrete.stateChanged.connect(self.set_id_ACI_concrete)
        checkbox_id_ACI_concrete.stateChanged.connect(self.update_diagram_plot)
        
        checkbox_id_ACI_total = QCheckBox()
        checkbox_id_ACI_total.setCheckState(self.globvar.flag_id_ACI_total)
        checkbox_id_ACI_total.stateChanged.connect(self.set_id_ACI_total)
        checkbox_id_ACI_total.stateChanged.connect(self.update_diagram_plot)
        
        checkbox_id_MC2010_concrete = QCheckBox()
        checkbox_id_MC2010_concrete.setCheckState(self.globvar.flag_id_MC2010_concrete)
        checkbox_id_MC2010_concrete.stateChanged.connect(self.set_id_MC2010_concrete)
        checkbox_id_MC2010_concrete.stateChanged.connect(self.update_diagram_plot)

        checkbox_id_MC2010_total = QCheckBox()
        checkbox_id_MC2010_total.setCheckState(self.globvar.flag_id_MC2010_total)
        checkbox_id_MC2010_total.stateChanged.connect(self.set_id_MC2010_total)
        checkbox_id_MC2010_total.stateChanged.connect(self.update_diagram_plot)

        checkbox_id_MSR_concrete = QCheckBox()
        checkbox_id_MSR_concrete.setCheckState(self.globvar.flag_id_MSR_concrete)
        checkbox_id_MSR_concrete.stateChanged.connect(self.set_id_MSR_concrete)
        checkbox_id_MSR_concrete.stateChanged.connect(self.update_diagram_plot)
        
        checkbox_id_MSR_total = QCheckBox()
        checkbox_id_MSR_total.setCheckState(self.globvar.flag_id_MSR_total)
        checkbox_id_MSR_total.stateChanged.connect(self.set_id_MSR_total)
        checkbox_id_MSR_total.stateChanged.connect(self.update_diagram_plot)

        checkbox_id_FEM = QCheckBox()
        checkbox_id_FEM.setCheckState(self.globvar.flag_id_FEM)
        checkbox_id_FEM.stateChanged.connect(self.set_id_FEM)
        checkbox_id_FEM.stateChanged.connect(self.update_diagram_plot)

        self.checkbox_id_FEM_paths = QCheckBox()
        self.checkbox_id_FEM_paths.setCheckState(self.globvar.flag_id_FEM_paths)
        self.checkbox_id_FEM_paths.stateChanged.connect(self.set_id_FEM_paths)
        self.checkbox_id_FEM_paths.stateChanged.connect(self.update_diagram_plot)     

        # position of neutral axis
        widget_c_neutral_axis = QDoubleSpinBox()
        widget_c_neutral_axis.setMinimum(0.)
        widget_c_neutral_axis.setMaximum( 1.25 * self.globvar.By )
        widget_c_neutral_axis.setDecimals(3)
        widget_c_neutral_axis.setValue( self.globvar.c_neutral_axis )
                
        widget_c_neutral_axis.setSingleStep(0.05)  
        widget_c_neutral_axis.valueChanged.connect(self.set_c_neutral_axis)
        widget_c_neutral_axis.valueChanged.connect(self.update_topology_plot)
        widget_c_neutral_axis.valueChanged.connect(self.update_diagram_plot)

        # solution for neutral axis visible
        checkbox_active_neutral_axis = QCheckBox()
        checkbox_active_neutral_axis.setCheckState(self.globvar.flag_active_neutral_axis)
        checkbox_active_neutral_axis.stateChanged.connect(self.set_active_neutral_axis)
        checkbox_active_neutral_axis.stateChanged.connect(self.update_topology_plot)
        checkbox_active_neutral_axis.stateChanged.connect(self.update_diagram_plot)

        # display loading points
        self.checkbox_show_loading = QCheckBox()
        self.checkbox_show_loading.setCheckState(self.globvar.flag_loading_display)
        self.checkbox_show_loading.stateChanged.connect(self.set_loading_display)
        self.checkbox_show_loading.stateChanged.connect(self.update_diagram_plot) 

        #########################
        # INTERACTION DIAGRAM CANVAS
        #########################

        self.diagram_canvas = MplCanvas(self, width=6, height=4, dpi=100)  

        id_specs_layout = QGridLayout()

        row = 0

        row += 1
        label_fcm = QLabel("Concrete class (characteristic strength)")
        label_fcm_unit = QLabel("[MPa]")
        id_specs_layout.addWidget(label_fcm, row, 0, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(self.widget_fcm, row, 1, 1, 1)
        id_specs_layout.addWidget(label_fcm_unit, row, 2)

        row += 1
        label_fy_vert = QLabel("Yield stress of vertical reinforcement")
        label_fy_vert_unit = QLabel("[MPa]")
        id_specs_layout.addWidget(label_fy_vert, row, 0, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(self.widget_fy_vert, row, 1, 1, 1)
        id_specs_layout.addWidget(label_fy_vert_unit, row, 2)
         
        row += 1
        label_fy_lat = QLabel("Yield stress of lateral reinforcement")
        label_fy_lat_unit = QLabel("[MPa]")
        id_specs_layout.addWidget(label_fy_lat, row, 0, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(self.widget_fy_lat, row, 1, 1, 1)
        id_specs_layout.addWidget(label_fy_lat_unit, row, 2)
        
        row += 1
        label_id_steel = QLabel("ID longitudinal reinforcement")
        id_specs_layout.addWidget(label_id_steel, row, 0, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(checkbox_id_steel, row, 1)

        row += 1
        id_specs_layout.addWidget(QLabel("ID ACI concrete"), row, 0, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(checkbox_id_ACI_concrete, row, 1)
        id_specs_layout.addWidget(QLabel("total"), row, 2, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(checkbox_id_ACI_total, row, 3)

        row += 1
        id_specs_layout.addWidget(QLabel("ID MC2010 concrete"), row, 0, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(checkbox_id_MC2010_concrete, row, 1)
        id_specs_layout.addWidget(QLabel("total"), row, 2, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(checkbox_id_MC2010_total, row, 3)
        
        row += 1
        id_specs_layout.addWidget(QLabel("ID MSR concrete"), row, 0, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(checkbox_id_MSR_concrete, row, 1)
        id_specs_layout.addWidget(QLabel("total"), row, 2, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(checkbox_id_MSR_total, row, 3)
        
        row += 1
        id_specs_layout.addWidget(QLabel("ID FEM"), row, 0, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(checkbox_id_FEM, row, 1)
        id_specs_layout.addWidget(QLabel("loading path"), row, 2, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(self.checkbox_id_FEM_paths, row, 3)
        
        row += 1
        id_specs_layout.addWidget(QLabel("Draw neutral axis"), row, 0, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(checkbox_active_neutral_axis, row, 1)
        id_specs_layout.addWidget(QLabel("loading"), row, 2, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(self.checkbox_show_loading, row, 3)

        row += 1
        id_specs_layout.addWidget( QLabel("Neutral axis from compression face"), row, 0, 1, 1, Qt.AlignRight)
        id_specs_layout.addWidget(widget_c_neutral_axis, row, 1, 1, 1)
        id_specs_layout.addWidget(QLabel("[m]"), row, 2)
       
        
        #########################
        # FEM etc.
        #########################
        
        # finite element mesh
        self.widget_dx = QDoubleSpinBox()
        self.widget_dx.setMinimum( self.globvar.Bx/200. )
        self.widget_dx.setMaximum( self.globvar.Bx/8. )
        self.widget_dx.setDecimals(4)
        self.widget_dx.setSingleStep(5./1000.)
        self.widget_dx.setValue( self.globvar.elem_size_X )
        self.widget_dx.valueChanged.connect( self.set_dx )
        self.widget_dx.valueChanged.connect( self.update_topology_plot )
        

        self.widget_dy = QDoubleSpinBox()
        self.widget_dy.setMinimum( self.globvar.By/200. )
        self.widget_dy.setMaximum( self.globvar.By/8. )
        self.widget_dy.setDecimals(4)
        self.widget_dy.setSingleStep(5./1000.) 
        self.widget_dy.setValue( self.globvar.elem_size_Y )
        self.widget_dy.valueChanged.connect( self.set_dy )
        self.widget_dy.valueChanged.connect( self.update_topology_plot )

        self.widget_dz = QDoubleSpinBox()
        self.widget_dz.setMinimum( self.globvar.H/20. )
        self.widget_dz.setMaximum( self.globvar.H/1. )
        self.widget_dz.setDecimals(4)
        self.widget_dz.setSingleStep(1./1000.) 
        self.widget_dz.setValue( self.globvar.elem_size_Z )
        self.widget_dz.valueChanged.connect( self.set_dz )

        self.checkbox_ignore_cover = QCheckBox()
        self.checkbox_ignore_cover.setCheckState(self.globvar.flag_ignore_cover)
        self.checkbox_ignore_cover.stateChanged.connect(self.set_ignore_cover)
        self.checkbox_ignore_cover.stateChanged.connect(self.adjust_element_sizes)
        self.checkbox_ignore_cover.stateChanged.connect(self.update_topology_plot)
        

        button_adjust_element_sizes = QPushButton("Adjust elements")
        button_adjust_element_sizes.clicked.connect(self.adjust_element_sizes)
        button_adjust_element_sizes.clicked.connect(self.update_topology_plot)
        
        self.checkbox_show_mesh = QCheckBox()
        self.checkbox_show_mesh.setCheckState(self.globvar.flag_show_mesh)
        self.checkbox_show_mesh.stateChanged.connect(self.set_show_mesh)
        self.checkbox_show_mesh.stateChanged.connect(self.update_topology_plot)

            
        button_CDPM2 = QPushButton("Concrete")
        button_CDPM2.clicked.connect( self.set_CDPM2 )

        button_Mises_vert = QPushButton("Vert. reinf.")
        button_Mises_vert.clicked.connect( self.set_mises_vert ) 

        button_Mises_lat = QPushButton("Lat. reinf.")
        button_Mises_lat.clicked.connect( self.set_mises_lat )
        
        button_generate_oofem_input = QPushButton("Generate oofem input(s)")
        button_generate_oofem_input.clicked.connect(self.generate_oofem_input )

        button_run_oofem = QPushButton("Run FEM analysis")
        button_run_oofem.clicked.connect(self.run_oofem_problems )
      
        #########################
        # FEM MODEL LAYOUT
        #########################

        fem_model_layout = QGridLayout()

        row = 0

        row += 1
        label_dx = QLabel("dx = ")
        fem_model_layout.addWidget(label_dx, row, 0, 1, 1, Qt.AlignRight)
        fem_model_layout.addWidget(self.widget_dx, row, 1)
        self.label_Nx = QLabel("[m] (" + str(self.globvar.N_X) + " elements)")
        fem_model_layout.addWidget(self.label_Nx, row, 2, 1, 2)
        
        row += 1
        label_dy = QLabel("dy = ")
        fem_model_layout.addWidget(label_dy, row, 0, 1, 1, Qt.AlignRight)
        fem_model_layout.addWidget(self.widget_dy, row, 1)
        self.label_Ny = QLabel("[m] (" + str(self.globvar.N_Y) + " elements)")
        fem_model_layout.addWidget(self.label_Ny, row, 2, 1, 2)
        
        row += 1
        label_dz = QLabel("dz = ")
        fem_model_layout.addWidget(label_dz, row, 0, 1, 1, Qt.AlignRight)
        fem_model_layout.addWidget(self.widget_dz, row, 1)
        self.label_Nz = QLabel("[m] (" + str(self.globvar.N_Z) + " elements)")
        fem_model_layout.addWidget(self.label_Nz, row, 2, 1, 2)

        row += 1
        label_ignore_cover = QLabel("Ignore concrete cover")
        fem_model_layout.addWidget(label_ignore_cover, row, 0, 1, 2, Qt.AlignRight)
        fem_model_layout.addWidget(self.checkbox_ignore_cover, row, 2)
        
        row += 1
        fem_model_layout.addWidget(button_adjust_element_sizes, row, 0, 1, 2)
        label_show_mesh = QLabel("Show FE mesh")
        fem_model_layout.addWidget(label_show_mesh, row, 2, 1, 1, Qt.AlignRight)
        fem_model_layout.addWidget(self.checkbox_show_mesh, row, 3)
        
        row += 1
        label_fem_mat = QLabel("Materials:")
        fem_model_layout.addWidget(label_fem_mat, row, 0)
        fem_model_layout.addWidget(button_CDPM2, row, 1)
        fem_model_layout.addWidget(button_Mises_vert, row, 2)
        fem_model_layout.addWidget(button_Mises_lat, row, 3)      

        row += 1
        label_project_name = QLabel("Project name:")
        self.lineEditProject = QLineEdit()
        self.lineEditProject.setMaxLength(50)
        self.lineEditProject.setText(self.globvar.project_name)
        self.lineEditProject.textChanged.connect(self.set_project_name)
        
        fem_model_layout.addWidget(label_project_name, row, 0)
        fem_model_layout.addWidget(self.lineEditProject, row, 1, 1, 3)

        row += 1
        lineEditOOFEM = QLineEdit()
        lineEditOOFEM.setMaxLength(200)
        lineEditOOFEM.setPlaceholderText("specify OOFEM folder")
        lineEditOOFEM.textChanged.connect(self.set_oofem_folder)
        
        button_oofem_folder = QPushButton("OOFEM folder")
        button_oofem_folder.clicked.connect(self.select_oofem_folder)
        button_oofem_folder.clicked.connect( lambda: lineEditOOFEM.setText(self.globvar.oofem_folder) )
        
        fem_model_layout.addWidget(button_oofem_folder, row, 0)
        fem_model_layout.addWidget(lineEditOOFEM, row, 1, 1, 3)

        #########################
        # FEM CHECKBOXES LAYOUT
        #########################

        self.fem_checkboxes_layout = QGridLayout()
        
        self.fem_checkboxes = []
        
        row = 0
        label_checkboxes = QLabel("Normalized eccentricity selection:")
        self.fem_checkboxes_layout.addWidget(label_checkboxes, row, 0, 1, self.globvar.ecc_nr)
        

        row += 1
        for i in range(self.globvar.ecc_nr):
            label_i = QLabel(f"{( self.globvar.tasks[i].eccentricity_normalized[1] ):.1f}" )
            self.fem_checkboxes_layout.addWidget(label_i, row, i)

        row += 1
        for i in range(self.globvar.ecc_nr):

            if (self.globvar.tasks[i].status == Task_status.PREDEFINED ):
                continue
            
            checkbox = QCheckBox()
            if ( self.globvar.tasks[i].status == Task_status.UNSELECTED):
                checkbox.setCheckState( Qt.Unchecked)
            else:
                checkbox.setCheckState( Qt.Checked)

            checkbox.stateChanged.connect( self.update_ecc_selection )
            self.fem_checkboxes.append(checkbox)


        for checkbox in self.fem_checkboxes:
            index = self.fem_checkboxes.index(checkbox)
            #print(index)
            self.fem_checkboxes_layout.addWidget(checkbox, row, index)

        #########################
        # FEM ANALYSIS LAYOUT
        #########################
        
        self.fem_analysis_layout = QGridLayout()

        row = 0
        self.fem_analysis_layout.addWidget(button_generate_oofem_input, row, 0, 1, 4)

        row += 1
        self.fem_analysis_layout.addWidget(button_run_oofem, row, 0, 1, 1)
        
        cpu_label= QLabel("Threads # (less than " + str(multiprocessing.cpu_count()) + "):" )
        cpu_label.setAlignment(Qt.AlignVCenter)

        cpu_nr = QDoubleSpinBox()
        cpu_nr.setDecimals(0)
        recommended_cpu = math.ceil(multiprocessing.cpu_count()*3./4.)
        self.set_cpu_nr( recommended_cpu )
        cpu_nr.setValue( recommended_cpu )
        cpu_nr.setMinimum( 1 )
        cpu_nr.setMaximum( multiprocessing.cpu_count() )

        cpu_nr.valueChanged.connect(self.set_cpu_nr)
        
        self.fem_analysis_layout.addWidget(cpu_label, row, 1, 1, 2)
        self.fem_analysis_layout.addWidget(cpu_nr, row, 3)        
        
        row += 1

        label_console = QLabel("Analysis results:")
        self.fem_analysis_layout.addWidget(label_console, row, 0, 1, 4)
        
        row += 1
        self.console = QPlainTextEdit()
        self.console.setReadOnly(True)
        self.fem_analysis_layout.addWidget(self.console, row, 0, 1, 4)

        row += 1
        progress_label = QLabel("Load vs. estimate:")
        self.fem_analysis_layout.addWidget(progress_label, row, 0, 1, 1)
        self.progressBar = QProgressBar()
        self.progressBar.setMinimum(0)
        self.progressBar.setMaximum(100)
        self.fem_analysis_layout.addWidget(self.progressBar, row, 1, 1, 3)
        
        fem_group_layout = QVBoxLayout()
        fem_group_layout.addLayout(fem_model_layout)
        fem_group_layout.addLayout(self.fem_checkboxes_layout)
        fem_group_layout.addLayout(self.fem_analysis_layout)
        
        #########################
        # MALCOLM LAYOUT
        #########################
        
        malcolm_layout = QGridLayout()

        label_L = QLabel("Topology Specification")
        label_L.setAlignment(Qt.AlignCenter) 
        font_heading = label_L.font()
        font_heading.setPointSize(16)
        font_heading.setBold(True)
        label_L.setFont(font_heading)
        malcolm_layout.addWidget(label_L, 0, 0)
        malcolm_layout.addWidget(self.topology_canvas, 1, 0)
        malcolm_layout.addLayout(topology_specs_layout, 2, 0)

        label_C = QLabel("Interaction Diagram + Materials")
        label_C.setAlignment(Qt.AlignCenter) 
        label_C.setFont(font_heading)
        
        malcolm_layout.addWidget(label_C, 0, 1)
        malcolm_layout.addWidget(self.diagram_canvas, 1, 1)
        malcolm_layout.addLayout(id_specs_layout, 2, 1)

        label_R = QLabel("FEM Definition")
        label_R.setAlignment(Qt.AlignCenter) 
        label_R.setFont(font_heading)
        malcolm_layout.addWidget(label_R, 0, 2)
        malcolm_layout.addLayout(fem_group_layout, 1, 2, 2, 1)

        self.update_topology_plot()
        self.update_diagram_plot()

        # TODO: improve
        self.globvar.sets = max( len(self.globvar.materials), len(self.globvar.cross_sections) )

        widget = QWidget()
        widget.setLayout(malcolm_layout)
        
        self.setCentralWidget(widget)
        
               
        # prevent maximization and rescaling
        #self.setFixedSize(1280, 720)

        #https://stackoverflow.com/questions/60815433/deprecationwarning-qdesktopwidget-availablegeometryint-screen-const-is-deprec
        '''
        self.setGeometry(
            QStyle.alignedRect(
            Qt.LeftToRight,
            Qt.AlignCenter,
            self.size(),
            QGuiApplication.primaryScreen().availableGeometry(),
            ),
        )
        '''

        #self.show()

        
    # closes all windows when the main is closed
    def closeEvent(self, event):
        for window in QApplication.topLevelWidgets():
            window.close()
        #if self.w:
        #    self.w.close()
    
    def define_CDPM2(self, cdpm2):
        self.w_cdpm2 = window_CDPM2(cdpm2)
        self.w_cdpm2.show()
        self.w_cdpm2.setGeometry(
            QStyle.alignedRect(
            Qt.LeftToRight,
            Qt.AlignCenter,
            self.w_cdpm2.size(),
            QGuiApplication.primaryScreen().availableGeometry(),
            ),
        )


    def define_Mises(self, mises):
        self.w_mises = window_Mises(mises)
        self.w_mises.show()
        self.w_mises.setGeometry(
            QStyle.alignedRect(
            Qt.LeftToRight,
            Qt.AlignCenter,
            self.w_mises.size(),
            QGuiApplication.primaryScreen().availableGeometry(),
            ),
        )


        
    def warning_dialog_ok(self,message):
         button = QMessageBox.warning(self, "Warning", message, buttons=QMessageBox.Ok, defaultButton=QMessageBox.Ok)


    def info_dialog(self,message):
         button = QMessageBox.information(self, "Info", message, buttons=QMessageBox.Ok, defaultButton=QMessageBox.Ok)



    def update_topology_information(self):
        
        rho_vert = self.compute_vertical_reinforcement_ratio()
        self.globvar.rho_vert = rho_vert
        
        rho_lat = self.compute_lateral_reinforcement_ratio()
        self.globvar.rho_lat = rho_lat
        
        Aconf = self.compute_total_confined_area()
        A = self.globvar.Bx * self.globvar.By

        sigL_L = 0.
        sigL_S = 0.

        ke_L = 0.
        ke_S = 0.
        
        for rebar in self.globvar.rebars:
            if ( rebar.give_rebar_type() == 'spiral' and rebar.tag[0] == 'S' ):
                sigL_S = rebar.compute_confinement( self.globvar, self.globvar.fy_lat )
                ke_S = rebar.compute_confinement_effectiveness()
                break
            
        for rebar in self.globvar.rebars:
            if ( rebar.give_rebar_type() == 'spiral' and rebar.tag[0] == 'L' ):
                sigL_L = rebar.compute_confinement( self.globvar, self.globvar.fy_lat )
                ke_L = rebar.compute_confinement_effectiveness()
                break

        self.label_info_topology.setText(f"Reinforcement ratios: {(rho_vert*100.):.3f}% (vertical) {(rho_lat*100.):.3f}% (lateral) \nConfinement: sigL_L = {(sigL_L):.3f} MPa, sigL_S = {(sigL_S):.3f} MPa \nEffectiveness: ke_L = {(ke_L):.3f}, ke_S = {(ke_S):.3f}\nArea = {(A):.3f}m^2, A_eff = {(Aconf):.3f}m^2, A_eff/A = {(Aconf/A):.3f}")

                
    def update_topology_plot(self):

        self.update_topology_information()
        
        logger.info("Topology updated")

        # all ids to be updated      
        self.need_update_ids_flag = True
        
        self.topology_canvas.axes.cla()  # Clear the canvas.

        edge_color_opacity = 1 # 0<val<1
        face_color_opacity = 1. # 0<val<1

        face_color_opacity_V = 1.
        
        self.topology_canvas.axes.add_patch(matplotlib.patches.Rectangle((-self.globvar.Bx/2., -self.globvar.By/2.), self.globvar.Bx, self.globvar.By, edgecolor=(0, 0, 0, edge_color_opacity), facecolor=(192./255.,192./255.,192./255., face_color_opacity), linewidth=2))

        # draw FE mesh of concrete
        if (self.globvar.flag_show_mesh):

            Bx_net = self.globvar.Bx
            By_net = self.globvar.By

            # subtract concrete cover to produce mesh which approximately corresponds to the confined region
            if (self.globvar.flag_ignore_cover):
                Bx_net -= 2. * self.globvar.cover
                By_net -= 2. * self.globvar.cover
       
            x_coords = np.linspace(-Bx_net/2., Bx_net/2., num=self.globvar.N_X+1);
            y_coords = np.linspace(-By_net/2., By_net/2., num=self.globvar.N_Y+1);


            for x_coord in x_coords:
                line = matplotlib.lines.Line2D([x_coord, x_coord], [-By_net/2., By_net/2.], lw=0.5, color='blue', alpha=0.5)
                self.topology_canvas.axes.add_line(line)


            for y_coord in y_coords:
                line = matplotlib.lines.Line2D([-Bx_net/2., Bx_net/2.], [y_coord, y_coord], lw=0.5, color='blue', alpha=0.5)
                self.topology_canvas.axes.add_line(line)      

        for rebar in self.globvar.rebars:

            reb_type = rebar.give_rebar_type()
            if (reb_type == 'rebar'):
                # get first coordinate, expecting vertical orientation
                XYZ = rebar.XYZ[0]
                xy = [XYZ[0], XYZ[1]]
                radius = rebar.give_diameter(self.globvar)/2.
               
                circle_path = matplotlib.path.Path.circle(xy, radius, readonly=False)
                circle_patch = matplotlib.patches.PathPatch(circle_path, edgecolor=(0, 0, 0, edge_color_opacity), facecolor=(0.0, 1.0, 0.0, face_color_opacity_V), linewidth=1.)
                self.topology_canvas.axes.add_patch(circle_patch)
            
            elif (reb_type == 'hoop' or reb_type == 'spiral' and rebar.tag[0] == 'S' ):
            
                xy = [rebar.XYZ[0],rebar.XYZ[1]]
                annulus_patch = matplotlib.patches.Annulus(xy, rebar.radius + rebar.give_diameter(self.globvar)/2., rebar.give_diameter(self.globvar), edgecolor=(0, 0, 0, edge_color_opacity), facecolor=(1., 0.0, 0.0, face_color_opacity), linewidth=1)
                self.topology_canvas.axes.add_patch(annulus_patch)


            elif (reb_type == 'hoop' or reb_type == 'spiral' and rebar.tag[0] == 'L' ):
            
                xy = [rebar.XYZ[0],rebar.XYZ[1]]
                annulus_patch = matplotlib.patches.Annulus(xy, rebar.radius + rebar.give_diameter(self.globvar)/2., rebar.give_diameter(self.globvar), edgecolor=(0, 0, 0, edge_color_opacity), facecolor=(30./255.,144./255.,255./255., face_color_opacity), linewidth=1)
                self.topology_canvas.axes.add_patch(annulus_patch)
          
            else:
                raise ValueError


        #self.topology_canvas.axes.plot(self.xdata, self.ydata, color = 'chartreuse')
        # Trigger the topology_canvas to update and redraw.
        self.topology_canvas.axes.set_aspect('equal')
        self.topology_canvas.axes.autoscale_view()

        # remove border lines
        self.topology_canvas.axes.spines['top'].set_visible(False)
        self.topology_canvas.axes.spines['right'].set_visible(False)
        self.topology_canvas.axes.spines['bottom'].set_visible(False)
        self.topology_canvas.axes.spines['left'].set_visible(False)

        # DRAW NEUTRAL AXIS & CONFINEMENT FEATURES
        if ( self.globvar.flag_active_neutral_axis ):
            
            NO_y = self.globvar.By / 2. - self.globvar.c_neutral_axis
            self.topology_canvas.axes.plot( [-self.globvar.Bx/2., self.globvar.Bx/2.],  [NO_y, NO_y], linestyle='dashdot', color='black', linewidth=1)

            if (self.globvar.debug_flag):
                for single in self.diagram_MSR.single_confined_areas:
                    self.topology_canvas.axes.plot( single.CG[0], single.CG[1], color='red', marker='2', markeredgewidth=2, markersize=10)
            
                for double in self.diagram_MSR.double_confined_areas:
                    self.topology_canvas.axes.plot( double.CG[0], double.CG[1], color='magenta', marker='1', markeredgewidth=2, markersize=10)

                          
        self.topology_canvas.draw()
        self.update_diagram_plot()

        

    def update_diagram_plot(self):

        if ( self.need_update_ids_flag ):
            self.diagram_ACI.clear_results()
            self.diagram_MC2010.clear_results()
            self.diagram_MSR.clear_results()
            self.diagram_MSR_simple.clear_results()
        
        self.diagram_canvas.axes.cla()  # Clear the canvas.

        self.diagram_canvas.axes.set_xlabel('Bending moment [MNm]')
        self.diagram_canvas.axes.set_ylabel('Normal force [MN]')
        self.diagram_canvas.axes.invert_yaxis()

        self.diagram_canvas.axes.grid(True)

        self.diagram_canvas.axes.autoscale_view()
        self.diagram_canvas.axes.spines['top'].set_visible(False)
        self.diagram_canvas.axes.spines['right'].set_visible(False)


           
        ### FEM Results
        if ( self.globvar.flag_id_FEM ):


            # update value for pure tension:
            tension = self.globvar.tasks[-1]
            # check we have got the correct task
            if (tension.status == Task_status.PREDEFINED):
                tension.max_MN = self.diagram_MSR.compute_steel_MN(self.globvar, 0.)
                
            max_N = []
            max_M = []
            for task in self.globvar.tasks:
                #if ( not math.isnan(task.max_load) ):
                if ( task.status == Task_status.COMPLETED  or  task.status == Task_status.PREDEFINED):

                    max_M.append(task.max_MN[0])
                    max_N.append(task.max_MN[1])


            if(self.globvar.flag_id_FEM_paths):
                for task in self.globvar.tasks:
                    if ( task.status == Task_status.COMPLETED ):

                        self.annotate_task_results(task)

                        aux_M = task.M
                        aux_M.insert(0, 0.)
                        aux_N = task.N
                        aux_N.insert(0, 0.)
                        self.diagram_canvas.axes.plot(aux_M, aux_N, color = (192./255.,192./255.,192./255.), linestyle = 'solid', linewidth = 0.5, marker = 'x', markeredgecolor = 'black', markeredgewidth=0.5, markersize=2)


            if (len(max_M) >= 3):
                
                self.diagram_canvas.axes.plot(max_M, max_N, '-or', markeredgewidth=1.5, markersize=6, label = "FEM", markeredgecolor = 'black')

            # draw single points
            elif (len(max_M) > 1):
                self.diagram_canvas.axes.plot(max_M, max_N, 'or', markeredgewidth=1.5, markersize=6, label = "FEM", markeredgecolor = 'black')
                    


        # color tables: 
        # https://www.rapidtables.com/web/color/RGB_Color.html
        c_ACI = (153./255.,51./255.,155./255.)
        c_MC = (30./255.,144./255.,255./255.)
        c_MSR = (0./255.,204./255.,0./255.)
        c_MSR_simple = (255./255.,153./255.,153./255.)

        
        ### STEEL ONLY - is the same for all models, however we use the data from MC2010
        if ( self.globvar.flag_id_steel ):
            if not self.diagram_MC2010.MNc_steel:
                self.diagram_MC2010.compute_steel_id(self.globvar)
            M, N, c = zip(* self.diagram_MC2010.MNc_steel)
            self.diagram_canvas.axes.plot(M,N, label='steel (MC2010)', color='black', linestyle='dashed')
        

        ### ACI ###
       
        # compute interaction diagrams if missing
        if ( ( self.globvar.flag_id_ACI_concrete ) or ( self.globvar.flag_id_ACI_total ) ):
            if not self.diagram_ACI.MNc_concrete:
                self.diagram_ACI.compute_concrete_id(self.globvar)
            
            if not self.diagram_ACI.MNc_steel:
                self.diagram_ACI.compute_steel_id(self.globvar)

        if ( self.globvar.flag_id_ACI_concrete ):
            M, N, c = zip(* self.diagram_ACI.MNc_concrete)
            self.diagram_canvas.axes.plot(M,N, label='ACI concrete', color=c_ACI, linestyle ='dashdot')

        if ( self.globvar.flag_id_ACI_total ):
            Mc, Nc, c = zip(* self.diagram_ACI.MNc_concrete)
            Ms, Ns, c = zip(* self.diagram_ACI.MNc_steel)

            M = [sum(i) for i in zip(Mc, Ms)]
            N = [sum(i) for i in zip(Nc, Ns)]
            
            self.diagram_canvas.axes.plot(M,N, label='ACI', color=c_ACI, linestyle ='solid')


        ### MC2010 ###
        # compute interaction diagrams if missing
        if ( ( self.globvar.flag_id_MC2010_concrete ) or ( self.globvar.flag_id_MC2010_total ) ):
            if not self.diagram_MC2010.MNc_concrete:
                self.diagram_MC2010.compute_concrete_id(self.globvar)
            
            if not self.diagram_MC2010.MNc_steel:
                self.diagram_MC2010.compute_steel_id(self.globvar)

        if ( self.globvar.flag_id_MC2010_concrete ):
            M, N, c = zip(* self.diagram_MC2010.MNc_concrete)
            self.diagram_canvas.axes.plot(M,N, label='MC2010 concrete', color=c_MC, linestyle ='dashdot')

        if ( self.globvar.flag_id_MC2010_total ):
            Mc, Nc, c = zip(* self.diagram_MC2010.MNc_concrete)
            Ms, Ns, c = zip(* self.diagram_MC2010.MNc_steel)

            M = [sum(i) for i in zip(Mc, Ms)]
            N = [sum(i) for i in zip(Nc, Ns)]
            
            self.diagram_canvas.axes.plot(M,N, label='MC2010', color=c_MC, linestyle ='solid')

        ### MSR ###
        # compute interaction diagrams if missing
        if ( ( self.globvar.flag_id_MSR_concrete ) or ( self.globvar.flag_id_MSR_total ) ):
            if not self.diagram_MSR.MNc_concrete:
                self.diagram_MSR.compute_concrete_id(self.globvar)
            
            if not self.diagram_MSR.MNc_steel:
                self.diagram_MSR.compute_steel_id(self.globvar)

            if not self.diagram_MSR_simple.MNc_concrete:
                self.diagram_MSR_simple.compute_concrete_id(self.globvar)
            
            if not self.diagram_MSR_simple.MNc_steel:
                self.diagram_MSR_simple.compute_steel_id(self.globvar)

        if ( self.globvar.flag_id_MSR_concrete ):
            M, N, c = zip(* self.diagram_MSR.MNc_concrete)
            self.diagram_canvas.axes.plot(M,N, label='MSR concrete', color=c_MSR, linestyle ='dashdot', linewidth=2.5)
            M, N, c = zip(* self.diagram_MSR_simple.MNc_concrete)
            self.diagram_canvas.axes.plot(M,N, color=c_MSR_simple, linestyle ='dashdot', linewidth=1)
            
        if ( self.globvar.flag_id_MSR_total ):
            Mc, Nc, c = zip(* self.diagram_MSR.MNc_concrete)
            Ms, Ns, c = zip(* self.diagram_MSR.MNc_steel)

            M = [sum(i) for i in zip(Mc, Ms)]
            N = [sum(i) for i in zip(Nc, Ns)]
            
            self.diagram_canvas.axes.plot(M,N, label='MSR', color=c_MSR, linestyle ='solid', linewidth=2.5)

            Mc, Nc, c = zip(* self.diagram_MSR_simple.MNc_concrete)
            Ms, Ns, c = zip(* self.diagram_MSR_simple.MNc_steel)

            M = [sum(i) for i in zip(Mc, Ms)]
            N = [sum(i) for i in zip(Nc, Ns)]
            
            self.diagram_canvas.axes.plot(M,N, label='FEM estimate', color=c_MSR_simple, linestyle ='solid', linewidth=1)  


        # point plots for displayed neutral axis
        if ( self.globvar.flag_active_neutral_axis ):
            m_size = 5

            if ( self.globvar.flag_id_steel ):
                Ms, Ns = self.diagram_MC2010.compute_steel_MN(self.globvar, self.globvar.c_neutral_axis)
                self.diagram_canvas.axes.plot( Ms, Ns, color='red', marker='o', markeredgewidth=1, markersize=m_size, markeredgecolor = 'black')
            
            if (  self.globvar.flag_id_ACI_concrete ): 
                Mc, Nc = self.diagram_ACI.compute_concrete_MN(self.globvar, self.globvar.c_neutral_axis)
                self.diagram_canvas.axes.plot( Mc, Nc, color='red', marker='o', markeredgewidth=1, markersize=m_size, markeredgecolor = c_ACI)

            if ( self.globvar.flag_id_ACI_total ):
                Mc, Nc = self.diagram_ACI.compute_concrete_MN(self.globvar, self.globvar.c_neutral_axis)
                Ms, Ns = self.diagram_ACI.compute_steel_MN(self.globvar, self.globvar.c_neutral_axis)
                self.diagram_canvas.axes.plot( Mc+Ms, Nc+Ns, color='red', marker='o', markeredgewidth=1, markersize=m_size, markeredgecolor = c_ACI)    
                
            if (  self.globvar.flag_id_MC2010_concrete ): 
                Mc, Nc = self.diagram_MC2010.compute_concrete_MN(self.globvar, self.globvar.c_neutral_axis)
                self.diagram_canvas.axes.plot( Mc, Nc, color='red', marker='o', markeredgewidth=1, markersize=m_size, markeredgecolor = c_MC)

            if ( self.globvar.flag_id_MC2010_total ):
                Mc, Nc = self.diagram_MC2010.compute_concrete_MN(self.globvar, self.globvar.c_neutral_axis)
                Ms, Ns = self.diagram_MC2010.compute_steel_MN(self.globvar, self.globvar.c_neutral_axis)
                self.diagram_canvas.axes.plot( Mc+Ms, Nc+Ns, color='red', marker='o', markeredgewidth=1, markersize=m_size, markeredgecolor = c_MC)

            if (  self.globvar.flag_id_MSR_concrete ):
                Mc, Nc = self.diagram_MSR.interpolate_concrete_MN(self.globvar.c_neutral_axis)
                #Mc, Nc = self.diagram_MSR.compute_concrete_MN(self.globvar, self.globvar.c_neutral_axis)
                self.diagram_canvas.axes.plot( Mc, Nc, color='red', marker='o', markeredgewidth=1, markersize=m_size, markeredgecolor = c_MSR)

            if ( self.globvar.flag_id_MSR_total ):
                #Mc, Nc = self.diagram_MSR.compute_concrete_MN(self.globvar, self.globvar.c_neutral_axis)
                Mc, Nc = self.diagram_MSR.interpolate_concrete_MN(self.globvar.c_neutral_axis)
                Ms, Ns = self.diagram_MSR.compute_steel_MN(self.globvar, self.globvar.c_neutral_axis)
                self.diagram_canvas.axes.plot( Mc+Ms, Nc+Ns, markerfacecolor='red', marker='o', markeredgewidth=1, markersize=m_size, markeredgecolor = c_MSR)


        ### LOADING
        if (self.globvar.flag_loading_selected and  (self.globvar.flag_loading_display==Qt.Checked) ):

            for col in self.globvar.loading.columns:
                if 'N' in col:
                    N_key = col
                elif 'M' in col:
                    M_key = col

            #print("N_flag = " + str(N_key) )
            #print("M_flag = " + str(M_key) )

            c_loading = (255./255.,128./255.,0./255.)
            for row in self.globvar.loading.iterrows():
                M = row[1][M_key]
                N = row[1][N_key]
                self.diagram_canvas.axes.plot( M, N, markerfacecolor=c_loading, marker='x', markeredgewidth=1.5, markersize=6, markeredgecolor = c_loading)
                #print([M,N])
                

        self.diagram_canvas.axes.legend(loc='lower right')
                        
        self.diagram_canvas.draw()

        self.need_update_ids_flag = False

              

    def update_small_spirals(self):

        for reb in self.globvar.rebars[:]:
            if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'S') ):

                dS = reb.radius * 2.
                rebar_DS = reb.give_diameter(self.globvar)
                
                auxX = self.globvar.Bx/2. - self.globvar.cover - dS/2  - rebar_DS/2.
                auxY = self.globvar.By/2. - self.globvar.cover - dS/2. - rebar_DS/2.

                # position of center and start based on previous position
                center = [auxX * np.sign(reb.XYZ[0]), auxY * np.sign(reb.XYZ[1]), reb.XYZ[2]]
                start = [ center[0], center[1] - dS/2., center[2] ]

                spiral_S = spiral(globvar = self.globvar, XYZ = center, cs_nr = reb.cs_nr, cs_tag = reb.give_cs_tag(self.globvar), cs_mat = reb.give_mat_nr(self.globvar), vec = self.vec, pitch = self.globvar.H, axial_length = self.globvar.H, start = start, tag = reb.tag)
                
                #spiral_S.check_topology_consistency(self.globvar)
  
                # switch the items in the container
                self.globvar.rebars[ self.globvar.rebars.index(reb) ] = spiral_S


    def update_large_spirals(self):
     
        for reb in self.globvar.rebars[:]:
            if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'L') ):

                dL = min(self.globvar.Bx, self.globvar.By ) - 2.*self.globvar.cover
                dL -= reb.give_diameter(self.globvar)    
                
                # position of center and start based on previous position
                center = reb.XYZ
                start = [ center[0], center[1] - dL/2., center[2] ]

                spiral_L = spiral(globvar = self.globvar, XYZ = center, cs_nr = reb.cs_nr, cs_tag = reb.give_cs_tag(self.globvar), cs_mat = reb.give_mat_nr(self.globvar), vec = self.vec, pitch = self.globvar.H, axial_length = self.globvar.H, start = start, tag = reb.tag)
                #spiral_L.check_topology_consistency(self.globvar)
  
                # switch the items in the container
                self.globvar.rebars[ self.globvar.rebars.index(reb) ] = spiral_L

                
    def update_vertical_rebars(self):
        self.remove_vertical_rebars()
        self.create_vertical_rebars()


    ####
    # VARIABLES FOR TOPOLOGY
    ####        
                            
    # set cross-section width & update
    def set_Bx(self, s):

        try:
            Bx = float(s)
        except ValueError:
            Bx_default = self.globvar.Bx
            self.warning_dialog_ok("unsupported value, using Bx = " + str(Bx_default) )
            warnings.warn("unsupported value, using Bx = " + str(Bx_default))
            Bx = Bx_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.Bx == Bx ):
            return

        # change
        if ( self.change_model_definition() ):

            self.globvar.Bx = Bx
            self.globvar.By = Bx
            logger.info(f"Setting Bx = {(self.globvar.Bx):.3f}")
            logger.info(f"Setting By = {(self.globvar.By):.3f}")
            
            self.update_small_spirals()
            self.update_large_spirals()
            self.update_vertical_rebars()
            self.diagram_MSR.update_MSR_intersections(self.globvar)

            # determine minimum element size - default values for computational efficiency
            def_elem_size_horizontal = min(self.globvar.Bx, self.globvar.By)/15.
        
            # update bounds for FE mesh 
            self.widget_dx.setMinimum( self.globvar.Bx/200. )
            self.widget_dx.setMaximum( self.globvar.Bx/8. )
            self.widget_dx.setValue( def_elem_size_horizontal )
            
            self.widget_dy.setMinimum( self.globvar.By/200. )
            self.widget_dy.setMaximum( self.globvar.By/8. )
            self.widget_dy.setValue( def_elem_size_horizontal )

        else:
            self.widget_Bx.setValue( self.globvar.Bx )


    # set spiral pitch & update
    def set_H(self, s):

        try:
            H = float(s)
        except ValueError:
            H_default = self.globvar.H
            self.warning_dialog_ok("unsupported value, using H = " + str(H_default) )
            warnings.warn("unsupported value, using H = " + str(H_default))
            H = H_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.H == H ):
            return
        
        # change
        if ( self.change_model_definition() ):

            self.globvar.H = H
            logger.info(f"Setting H = {(self.globvar.H):.3f}")
            
            self.update_small_spirals()
            self.update_large_spirals()
            self.update_vertical_rebars()        
            self.diagram_MSR.update_MSR_intersections(self.globvar)

            # update bounds for FE mesh 
            self.widget_dz.setMinimum( self.globvar.H/20. )
            self.widget_dz.setMaximum( self.globvar.H/1. )
            self.widget_dz.setValue( self.globvar.H/4.)

        else:
            self.widget_H.setValue( self.globvar.H )
        
    # set concrete cover and update spirals topology
    def set_cover(self, s):

        try:
            cover = float(s)
        except ValueError:
            cover_default = self.globvar.cover
            self.warning_dialog_ok("unsupported value, using cover = " + str(cover_default) )
            warnings.warn("unsupported value, using cover = " + str(cover_default))
            cover = cover_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.cover == cover ):
            return
        
        # change
        if ( self.change_model_definition() ):
            
            self.globvar.cover = cover
            logger.info(f"Setting cover = {(self.globvar.cover):.3f}")

            self.update_small_spirals()
            self.update_large_spirals()
            self.update_vertical_rebars()
            self.diagram_MSR.update_MSR_intersections(self.globvar)

        else:
            self.widget_cover.setValue( self.globvar.cover )

    # set axial diameter to all small spirals
    def set_dS(self, s):

        try:
            dS = float(s)
        except ValueError:
            dS_default = self.globvar.rebars[0].radius * 0.3 * 2.
            self.warning_dialog_ok("unsupported value, using dS = " + str(dS_default) )
            warnings.warn("unsupported value, using dS = " + str(dS_default))
            dS = dS_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.dS == dS ):
            return
        
        # change
        if ( self.change_model_definition() ):

            self.globvar.dS = dS
            logger.info(f"Setting dS = {(self.globvar.dS):.3f}")
             
            for reb in self.globvar.rebars:
                if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'S') ):
                
                    reb.radius = dS/2.

            self.update_small_spirals()
            self.update_vertical_rebars()
            self.diagram_MSR.update_MSR_intersections(self.globvar)

        else:
            self.widget_dS.setValue( self.globvar.dS )

            

    # set rebar diameter to all small spirals
    def set_DS(self, s):

        try:
            rebar_DS = (self.globvar.rebars_CS[s]).diam

        except ValueError:
            rebar_DS_default = (self.globvar.rebars_CS[self.globvar.DS]).diam
            self.warning_dialog_ok("unsupported value, using DS = " + str(rebar_DS_default) )
            warnings.warn("unsupported value, using DS = " + str(rebar_DS_default))
            rebar_DS = rebar_DS_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.DS == s ):
            return

        # change
        if ( self.change_model_definition() ):

            self.globvar.DS = s
            logger.info("Setting DS = " + self.globvar.DS)
            area = self.globvar.rebars_CS[self.globvar.DS].area
                        
            # update cross-section first - otherwise a new one is created in "update_xxx_spirals"
            for reb in self.globvar.rebars:
                if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'S') ):
                    cs = reb.find_cross_section(self.globvar, reb.cs_nr)
                    cs.set_properties_from_cs_tag(self.globvar, self.globvar.DS)
                    break

            self.update_small_spirals()
            self.update_vertical_rebars()
            self.diagram_MSR.update_MSR_intersections(self.globvar)

        else:
            self.widget_DS.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DS))

        
    # set axial diameter to large spiral(s)
    def set_dL(self, s):

        try:
            dL = float(s)
        except ValueError:
            dL_default = min(self.globvar.Bx, self.globvar.By )/2. - 2.*self.globvar.cover 
            dL_default -= (self.globvar.rebars_CS[self.globvar.DL]).diam        
            
            self.warning_dialog_ok("unsupported value, using dL = " + str(dL_default) )
            warnings.warn("unsupported value, using dL = " + str(dL_default))
            dL = dL_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.dL == dL ):
            return            

        # change
        if ( self.change_model_definition() ):

            self.globvar.dL = dL
            logger.info(f"Setting dL = {(self.globvar.dL):.3f}")
            
            for reb in self.globvar.rebars:
                if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'L') ):
                
                    reb.radius = dL/2.

            self.update_large_spirals()
            self.update_vertical_rebars()
            self.diagram_MSR.update_MSR_intersections(self.globvar)

        else:
            self.widget_dL.setValue( self.globvar.dL )
        

    # set rebar diameter to all small spirals
    def set_DL(self, s):

        try:
            rebar_DL = (self.globvar.rebars_CS[s]).diam

        except ValueError:
            rebar_DL_default = (self.globvar.rebars_CS[self.globvar.DL]).diam
            self.warning_dialog_ok("unsupported value, using DL = " + str(rebar_DL_default) )
            warnings.warn("unsupported value, using DL = " + str(rebar_DL_default))
            rebar_DL = rebar_DL_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.DL == s ):
            return

        # change
        if ( self.change_model_definition() ):
            
            self.globvar.DL = s
            logger.info(f"Setting DL = " + self.globvar.DL)
            
            area = self.globvar.rebars_CS[self.globvar.DL].area

            for reb in self.globvar.rebars:
                if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'L') ):
                    cs = reb.find_cross_section(self.globvar, reb.cs_nr)
                    cs.set_properties_from_cs_tag(self.globvar, self.globvar.DL)
                    break
                    
            self.update_large_spirals()
            self.update_vertical_rebars()
            self.diagram_MSR.update_MSR_intersections(self.globvar)
       
        else:
            self.widget_DL.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DL))


    def set_DV(self, s):

        try:
            rebar_DV = (self.globvar.rebars_CS[s]).diam
        except ValueError:
            rebar_DV_default = (self.globvar.rebars_CS[self.globvar.DV]).diam
            self.warning_dialog_ok("unsupported value, using DV = " + str(rebar_DV_default) )
            warnings.warn("unsupported value, using DV = " + str(rebar_DV_default))
            rebar_DV = rebar_DV_default


        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.DV == s ):
            return

        # change
        if ( self.change_model_definition() ):
        
            self.globvar.DV = s
            logger.info(f"Setting DV = " + self.globvar.DV)

            for reb in self.globvar.rebars:
                if ( (reb.give_rebar_type() == 'rebar') and (reb.tag[0] == 'V') ):
                    cs = reb.find_cross_section(self.globvar, reb.cs_nr)
                    cs.set_properties_from_cs_tag(self.globvar, self.globvar.DV)
                    break
                    # cross-section needs to be updated!
                
            self.update_vertical_rebars()
            self.diagram_MSR.update_MSR_intersections(self.globvar)

            # update cross-section
            for reb in self.globvar.rebars:
                if ( (reb.give_rebar_type() == 'rebar') and (reb.tag[0] == 'V') ):
                    #cs = reb.find_cross_section(self.globvar, reb.cs)
                    #cs.area = reb.area
                    # TODO
                    break

        else:
            self.widget_DV.setCurrentIndex(list(self.globvar.rebars_CS.keys()).index(self.globvar.DV))


    ####
    # VARIABLE FOR INTERACTION DIAGRAM
    ####
        
    # set concrete grade
    def set_fcm(self, s):

        # selection from catalogue
        try:
            fcm_val = float(self.globvar.concretes[s] )

        except ValueError:
            fcm_default = (self.globvar.concretes[self.globvar.fcm])
            self.warning_dialog_ok("unsupported value, using fcm = " + str(fcm_default) )
            warnings.warn("unsupported value, using fcm = " + str(fcm_default))
            fcm_val = fcm_default
            
        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.fcm == s ):
            return
        
        # change
        if ( self.change_model_definition() ):
            self.globvar.fcm = s
            logger.info(f"Setting fcm = " + self.globvar.fcm)
            self.cdpm2.fcm = fcm_val
            self.cdpm2.predict_concrete_parameters()

        else:
            self.widget_fcm.setCurrentIndex(list(self.globvar.concretes.keys()).index(self.globvar.fcm))
        

    # set steel grade for vertical (longitudinal) reinforcement 
    def set_fy_vert(self, s):

        try:
            fy = float(s)
        except ValueError:
            fy_default = self.globvar.fy_vert
            self.warning_dialog_ok("unsupported value, using fy = " + str(fy_default) )
            warnings.warn("unsupported value, using fy = " + str(fy_default))
            fy = fy_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.fy_vert == fy ):
            return
        
        # change
        if ( self.change_model_definition() ):

            self.globvar.fy_vert = fy
            logger.info(f"Setting fy_vert = {(self.globvar.fy_vert):.3f}")
            self.mises_vert.sig_0 = fy

        else:
            self.widget_fy_vert.setValue( self.globvar.fy_vert )



    # set steel grade for lateral (transverse) reinforcement 
    def set_fy_lat(self, s):

        try:
            fy = float(s)
        except ValueError:
            fy_default = self.globvar.fy_lat
            self.warning_dialog_ok("unsupported value, using fy = " + str(fy_default) )
            warnings.warn("unsupported value, using fy = " + str(fy_default))
            fy = fy_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.fy_lat == fy ):
            return
        
        # change
        if ( self.change_model_definition() ):

            self.globvar.fy_lat = fy
            logger.info(f"Setting fy_lat = {(self.globvar.fy_lat):.3f}")
            self.mises_lat.sig_0 = fy

        else:
            self.widget_fy_lat.setValue( self.globvar.fy_lat )


    def set_dx(self,s):

        try:
            dx = float(s)
        except ValueError:
            dx_default = self.globvar.elem_size_X
            self.warning_dialog_ok("unsupported value, using = " + str(dx_default) )
            warnings.warn("unsupported value, using = " + str(dx_default))
            dx = dx_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.elem_size_X == dx ):
            return
        
        # change
        if ( self.change_model_definition() ):
            
            self.set_show_mesh(Qt.Unchecked) 
            self.checkbox_show_mesh.setCheckState(self.globvar.flag_show_mesh)    

        else:
            self.widget_dx.setValue( self.globvar.elem_size_X )
            logger.info(f"Setting elem_size_X = {(self.globvar.elem_size_X):.3f}")

    def set_dy(self,s):

        try:
            dy = float(s)
        except ValueError:
            dy_default = self.globvar.elem_size_Y
            self.warning_dialog_ok("unsupported value, using = " + str(dy_default) )
            warnings.warn("unsupported value, using = " + str(dy_default))
            dy = dy_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.elem_size_Y == dy ):
            return
        
        # change
        if ( self.change_model_definition() ):
            
            self.set_show_mesh(Qt.Unchecked) 
            self.checkbox_show_mesh.setCheckState(self.globvar.flag_show_mesh)    

        else:
            self.widget_dy.setValue( self.globvar.elem_size_Y )
            logger.info(f"Setting elem_size_Y = {(self.globvar.elem_size_Y):.3f}")


    def set_dz(self,s):

        try:
            dz = float(s)
        except ValueError:
            dz_default = self.globvar.elem_size_Z
            self.warning_dialog_ok("unsupported value, using = " + str(dz_default) )
            warnings.warn("unsupported value, using = " + str(dz_default))
            dz = dz_default

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.elem_size_Z == dz ):
            return
        
        # change
        if ( self.change_model_definition() ):
            pass
            #self.set_show_mesh(Qt.Unchecked) 
            #self.checkbox_show_mesh.setCheckState(self.globvar.flag_show_mesh)    

        else:
            self.widget_dz.setValue( self.globvar.elem_size_Z )
            logger.info(f"Setting elem_size_Z = {(self.globvar.elem_size_Z):.3f}")
            


    def set_CDPM2(self):
        # change
        if ( self.change_model_definition() ):
            self.define_CDPM2(self.cdpm2)


    def set_mises_vert(self):
        # change
        if ( self.change_model_definition() ):
            self.define_Mises(self.mises_vert)


    def set_mises_lat(self):
        # change
        if ( self.change_model_definition() ):
            self.define_Mises(self.mises_lat)

       

    def set_id_steel(self,s):
        self.globvar.flag_id_steel = s
        
    def set_id_ACI_concrete(self,s):
        self.globvar.flag_id_ACI_concrete = s

    def set_id_ACI_total(self,s):
        self.globvar.flag_id_ACI_total = s

    def set_id_MC2010_concrete(self,s):
        self.globvar.flag_id_MC2010_concrete = s

    def set_id_MC2010_total(self,s):
        self.globvar.flag_id_MC2010_total = s

    def set_id_MSR_concrete(self,s):
        self.globvar.flag_id_MSR_concrete = s

    def set_id_MSR_total(self,s):
        self.globvar.flag_id_MSR_total = s

    def set_id_FEM(self,s):
        self.globvar.flag_id_FEM = s

        if ( not self.globvar.flag_id_FEM ):
            self.checkbox_id_FEM_paths.setCheckState(Qt.Unchecked)
            
        
    def set_id_FEM_paths(self,s):
        self.globvar.flag_id_FEM_paths = s

    def set_active_neutral_axis(self,s):
        self.globvar.flag_active_neutral_axis = s

    def set_loading_display(self,s):
        self.globvar.flag_loading_display = s


    # set steel grade for lateral (transverse) reinforcement 
    def set_c_neutral_axis(self, s):

        try:
            c = float(s)
        except ValueError:
            c_default = self.globvar.c_neutral_axis
            self.warning_dialog_ok("unsupported value, using c = " + str(c_default) )
            warnings.warn("unsupported value, c = " + str(c_default))
            c = c_default

        self.globvar.c_neutral_axis = c        


    def set_show_mesh(self,s):
        self.globvar.flag_show_mesh = s

    def set_ignore_cover(self,s):

        # not to show the dialog again if the value has been set to its former setting
        if ( self.globvar.flag_ignore_cover == s ):
            return

        if ( self.change_model_definition() ):
            self.globvar.flag_ignore_cover = s
        else:
            self.checkbox_ignore_cover.setCheckState(self.globvar.flag_ignore_cover)
            
    
    def change_model_definition(self):
        # if the definition has been changed previously, do nothing
        if (self.globvar.flag_problem_changed):
            self.need_update_ids_flag = True
            return True

        # otherwise ask whether the content should be deleted
        else:
            button = QMessageBox.question(self, "Change in problem definition detected", " Delete computed results/generated inputs and proceed?")

        if button == QMessageBox.Yes:
            # change flags
            self.globvar.flag_problem_changed = True
            self.globvar.flag_output_generated = False
            
            # project directory will be cleared in a sequel upon input creation
            for task in self.globvar.tasks:
                task.reset()

            self.update_ecc_selection()

            self.console.clear()

            self.progressBar.setValue(0.)
            
            self.need_update_ids_flag = True
            return True
        else:
            return False
            
        

    def remove_vertical_rebars(self):

        # removing all vertical rebars, iterating over a copy of a list
        for reb in self.globvar.rebars[:]:
            if ( (reb.give_rebar_type() == 'rebar') and (reb.tag[0] == 'V') ):
                self.globvar.rebars.remove(reb)


    def create_vertical_rebars(self):

        DV = self.globvar.rebars_CS[self.globvar.DV].diam

        #cs = self.find_cross_section(globvar, self.cs)
        #cs.D =
        #cs.A = 
        
        # large spirals - vertical quarters - inner
        for reb in self.globvar.rebars:
            if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'L') ):

                # start with axial radius, subtract spiral rebar radius and radius of vertical rebar
                diam = reb.give_diameter(self.globvar)
                plan_radius = reb.radius - diam/2. - DV/2.

                for i in range(4):
                    if (i == 0): # X+
                        topology_origin = [reb.XYZ[0] + plan_radius, reb.XYZ[1], 0.]
                    elif (i == 1): # X-
                        topology_origin = [reb.XYZ[0] - plan_radius, reb.XYZ[1], 0.]
                    elif (i == 2): # Y+
                        topology_origin = [reb.XYZ[0], reb.XYZ[1] + plan_radius, 0.]
                    elif (i == 3): # Y-
                        topology_origin = [reb.XYZ[0], reb.XYZ[1] - plan_radius, 0.]
                        
                    topology_end = [ topology_origin[0], topology_origin[1], self.globvar.H ]
                    topology = [topology_origin, topology_end]

                    rebar_V = rebar(globvar = self.globvar, XYZ = topology, cs_nr = 4, cs_tag = self.globvar.DV, cs_mat = 3, tag = 'V_center_'+reb.tag)
                    #rebar_V.check_topology_consistency(self.globvar)
                    self.globvar.rebars.append(rebar_V)


        # large spirals - diagonals - inner
        for reb in self.globvar.rebars:
            if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'L') ):

                # start with axial radius, subtract spiral rebar radius and radius of vertical rebar
                diam = reb.give_diameter(self.globvar)
                plan_radius = reb.radius - diam/2. - DV/2.
                radius_projection = plan_radius/math.sqrt(2)

                for i in range(4):
                    if (i == 0): # X+, Y+
                        topology_origin = [reb.XYZ[0] + radius_projection, reb.XYZ[1] + radius_projection, 0.]
                    elif (i == 1): # X-, Y+
                        topology_origin = [reb.XYZ[0] - radius_projection, reb.XYZ[1] + radius_projection, 0.]
                    elif (i == 2): # X-, Y-
                        topology_origin = [reb.XYZ[0] - radius_projection, reb.XYZ[1] - radius_projection, 0.]
                    elif (i == 3): # X+, Y-
                        topology_origin = [reb.XYZ[0] + radius_projection, reb.XYZ[1] - radius_projection, 0.]
                        
                    topology_end = [ topology_origin[0], topology_origin[1], self.globvar.H ]
                    topology = [topology_origin, topology_end]

                    rebar_V = rebar(globvar = self.globvar, XYZ = topology, cs_nr = 4, cs_tag = self.globvar.DV, cs_mat = 3, tag = 'V_diag_'+reb.tag)
                    #rebar_V.check_topology_consistency(globvar = self.globvar)
                    self.globvar.rebars.append(rebar_V)
                
        
        # small spirals - diagonals - inner, outer
        # flag_inner_S = True
        flag_inner_S = False
        flag_outer_S = True
        
        for reb in self.globvar.rebars:
            if ( (reb.give_rebar_type() == 'spiral') and (reb.tag[0] == 'S') ):
                # start with axial radius
                plan_radius = reb.radius
                # math.sqrt( (reb.start[0]-reb.XYZ[0])**2 + (reb.start[1]-reb.XYZ[1])**2 )
                plan_radius -= reb.give_diameter(self.globvar)/2. # inner radius
                plan_radius -= DV/2. # to axis of vertical rebar

                radius_projection = plan_radius/math.sqrt(2)

                for i in range(2):
                    # spiral center X+ Y+

                    if ( reb.XYZ[0] > 0.):
                        aux_X =  1.
                    else:
                        aux_X = -1.

                    if ( reb.XYZ[1] > 0.):
                        aux_Y =  1.
                    else:
                        aux_Y = -1.
                    
                    if ( (i == 0) and flag_inner_S ):
                        topology_origin = [reb.XYZ[0] - aux_X*radius_projection, reb.XYZ[1] - aux_Y*radius_projection, 0.]

                    elif ( (i == 1) and flag_outer_S ): 
                        topology_origin = [reb.XYZ[0] + aux_X*radius_projection, reb.XYZ[1] + aux_Y*radius_projection, 0.]

                    topology_end = [ topology_origin[0], topology_origin[1], self.globvar.H ]
                    topology = [topology_origin, topology_end]

                    rebar_V = rebar(globvar = self.globvar, XYZ = topology, cs_nr = 4, cs_tag = self.globvar.DV, cs_mat = 3, tag = 'V_'+reb.tag)
                    #rebar_V.check_topology_consistency(self.globvar)
                    self.globvar.rebars.append(rebar_V)
                
        

        # large vs. small spirals - inner or outer intersections
        # eg. S1 vs. L1 etc.
        # flag_inner_L_vs_S = True
        flag_inner_L_vs_S = False
        
        for reb_L in self.globvar.rebars:
            # outer loop over all large spirals and find all potential intersecting small spirals
            if ( (reb_L.give_rebar_type() == 'spiral') and (reb_L.tag[0] == 'L') ):
                plan_radius_L = reb_L.radius
                #math.sqrt( (reb_L.start[0]-reb_L.XYZ[0])**2 + (reb_L.start[1]-reb_L.XYZ[1])**2 )

                if flag_inner_L_vs_S:
                    r_L_to_V = plan_radius_L - reb_L.give_diameter(self.globvar)/2. - DV/2.
                else:
                    r_L_to_V = plan_radius_L + reb_L.give_diameter(self.globvar)/2. + DV/2.

                # inner loop over inner spirals
                for reb_S in self.globvar.rebars:
                    # outer loop over all large spirals and find all potential intersecting small spirals
                    if ( (reb_S.give_rebar_type() == 'spiral') and (reb_S.tag[0] == 'S') ):
                        plan_radius_S = reb_S.radius
                        # plan_radius_S = math.sqrt( (reb_S.start[0]-reb_S.XYZ[0])**2 + (reb_S.start[1]-reb_S.XYZ[1])**2)
                        # calculate distance center_L to center_S
                        center_to_center = math.sqrt( (reb_L.XYZ[0]-reb_S.XYZ[0])**2 + (reb_L.XYZ[1]-reb_S.XYZ[1])**2 )

                        if (center_to_center < plan_radius_L + plan_radius_S - reb_L.give_diameter(self.globvar)/2. - reb_S.give_diameter(self.globvar)/2. - DV):
                            r_S_to_V = plan_radius_S - reb_S.give_diameter(self.globvar)/2. - DV/2.

                            DD = center_to_center/math.sqrt(2.)

                            alpha = math.acos ( (2*DD**2 + r_L_to_V**2 - r_S_to_V**2 )/ (2.*math.sqrt(2.) * DD * r_L_to_V) );
                            beta = math.pi/4.-alpha;

                            D_sin = r_L_to_V*math.sin(beta);
                            D_cos = r_L_to_V*math.cos(beta);

                            if ( reb_S.XYZ[0] > reb_L.XYZ[0] ):
                                aux_X = 1.
                            else:
                                aux_X = -1.

                            if ( reb_S.XYZ[1] > reb_L.XYZ[1] ):
                                aux_Y = 1.
                            else:
                                aux_Y = -1.
                                                  
                            ## nr 1
                            topology_origin = [reb_L.XYZ[0] + aux_X*D_sin, reb_L.XYZ[1] + aux_Y*D_cos, 0.]
                            topology_end = [ topology_origin[0], topology_origin[1], self.globvar.H ]
                            topology = [topology_origin, topology_end]

                            rebar_V = rebar(globvar = self.globvar, XYZ = topology, cs_nr = 4, cs_tag = self.globvar.DV, cs_mat = 3, tag = 'V_' + reb_L.tag + reb_S.tag)
                            #rebar_V.check_topology_consistency(self.globvar)
                            self.globvar.rebars.append(rebar_V)

                            ## nr 2
                            topology_origin = [reb_L.XYZ[0] + aux_X*D_cos, reb_L.XYZ[1] + aux_Y*D_sin, 0.]
                            topology_end = [ topology_origin[0], topology_origin[1], self.globvar.H ]
                            topology = [topology_origin, topology_end]
                            
                            rebar_V = rebar(globvar = self.globvar, XYZ = topology, cs_nr = 4, cs_tag = self.globvar.DV, cs_mat = 3, tag = 'V_' + reb_L.tag + reb_S.tag)
                            #rebar_V.check_topology_consistency(self.globvar)
                            self.globvar.rebars.append(rebar_V)
                




    def annotate_task_results(self,task):

        text = (f"e\u0302={(task.eccentricity_normalized[1]):.2f}\n[{(task.max_MN[0]):.2f}, {(task.max_MN[1]):.2f}]" )
        
        dN = -2.
        dM = -0.25

        self.diagram_canvas.axes.annotate(
        text,
        fontsize=8,
        xy=(task.max_MN[0], task.max_MN[1]), xycoords='data',
        xytext=(task.max_MN[0]+dM,task.max_MN[1]+dN),
        arrowprops=dict(arrowstyle="->",
        connectionstyle="arc3,rad=.2"))

     
        
    def create_finite_elements(self):

        self.globvar.ndofman = 0
        self.globvar.nelem = 0

        self.globvar.master_nodes.clear()
        self.globvar.brick_elements.clear()

        self.globvar.hanging_nodes.clear()
        self.globvar.truss_elements.clear()
        
        for rebar in self.globvar.rebars:
            rebar.rebar_points.clear()
            
        for i_cs in self.globvar.cross_sections:
            i_cs.elements.clear()
            

        self.concrete.create_concrete_mesh(self.globvar)
        
        for rebar in self.globvar.rebars:

            # create only vertices
            rebar.create_rebar_nodes()

            # create rebar elements
            rebar.create_rebar_elements(self.globvar)

        info = "generated FE mesh"
        logger.info(info)


    def compute_lateral_reinforcement_ratio(self):

        
        vol_concrete = self.globvar.Bx * self.globvar.By * self.globvar.H
        vol_steel = 0.

        for reb in self.globvar.rebars:
            if ( reb.give_rebar_type() == 'spiral') :

                vol_steel += reb.compute_rebar_volume(self.globvar)

        rho = vol_steel / vol_concrete
        return rho

                
    def compute_vertical_reinforcement_ratio(self):

        vol_concrete = self.globvar.Bx * self.globvar.By * self.globvar.H
        vol_steel = 0.
        
        for reb in self.globvar.rebars:
            if ( (reb.give_rebar_type() == 'rebar') and (reb.tag[0] == 'V') ):

                vol_steel += reb.compute_rebar_volume(self.globvar)

        rho = vol_steel / vol_concrete
        return rho

    
    def compute_total_confined_area(self):

        A_conf = 0.
        # summ all spiral areas
        for reb in self.globvar.rebars:
            if ( reb.give_rebar_type() == 'spiral') :
                A_conf += self.compute_circle_area(reb)

                # subtract all multiple-confined areas
                for reb_intersect in self.globvar.rebars:
                    if ( reb_intersect.give_rebar_type() == 'spiral'):
                        # to subtract the intersection only once
                        if (self.globvar.rebars.index(reb_intersect) > self.globvar.rebars.index(reb) ):

                            A_conf -= self.compute_intersecting_area_of_two_circles(reb,reb_intersect)

        return A_conf
                

    def compute_circle_area(self,circ):
        # compute area - axial-wise

        A_circle = math.pi * circ.radius**2

        return A_circle

                
    def compute_intersecting_area_of_two_circles(self,circ_1,circ_2):
        # assuming circles in horizontal plane
        # computed axial-wise
        # https://mathworld.wolfram.com/Circle-CircleIntersection.html

        # distance between centers
        d = math.sqrt( (circ_1.XYZ[0]-circ_2.XYZ[0])**2 + (circ_1.XYZ[1]-circ_2.XYZ[1])**2 )
        # radius of first and second circles
        R1 = circ_1.radius
        R2 = circ_2.radius

        A_intersect = 0.

        if ( d < R1+R2 ):
            d1 = (d**2 - R2**2 + R1**2) / (2.*d)
            d2 = d - d1

            A_intersect += R1**2 * math.acos(d1/R1) - d1 * math.sqrt( R1**2 - d1**2 )
            A_intersect += R2**2 * math.acos(d2/R2) - d2 * math.sqrt( R2**2 - d2**2 )

        return A_intersect
        


    def compute_unconfined_area(self):
        
        area_unconf = (self.globvar.Bx * self.globvar.By) - self.compute_confined_area()
        return area_unconf


    def adjust_element_sizes(self):

        Bx_net = self.globvar.Bx
        By_net = self.globvar.By

        # subtract concrete cover to produce mesh which approximately corresponds to the confined region
        if (self.globvar.flag_ignore_cover):
            Bx_net -= 2. * self.globvar.cover
            By_net -= 2. * self.globvar.cover
        
        self.globvar.N_X = round(Bx_net/ self.widget_dx.value() );
        self.globvar.N_Y = round(By_net / self.widget_dy.value() );
        self.globvar.N_Z = round(self.globvar.H / self.widget_dz.value() );

        self.label_Nx.setText("[m] (" + str(self.globvar.N_X) + " elements)")
        self.label_Ny.setText("[m] (" + str(self.globvar.N_Y) + " elements)")
        self.label_Nz.setText("[m] (" + str(self.globvar.N_Z) + " elements)")
        
        # final element size
        self.globvar.elem_size_X = Bx_net/self.globvar.N_X;
        self.globvar.elem_size_Y = By_net/self.globvar.N_Y;
        self.globvar.elem_size_Z = self.globvar.H/self.globvar.N_Z;

        self.widget_dx.setValue( self.globvar.elem_size_X )
        self.widget_dy.setValue( self.globvar.elem_size_Y )
        self.widget_dz.setValue( self.globvar.elem_size_Z )

    def set_project_name(self, s):

        try:
            name = str(s)
        except ValueError:
            warnings.warn("unsupported format of project name")
            name = "test"
        
        self.globvar.project_name = name
        logger.info(f"Setting project name to  = " + self.globvar.project_name)
        

    def select_oofem_folder(self):
        folder = "/home/pedro/Programs/oofem_official/python/"
        # temporary comment
        #folder = str(QFileDialog.getExistingDirectory(self, "Select OOFEM Directory"))
        self.set_oofem_folder(folder)

           
        #self.globvar.oofem_folder = 
        

    def set_oofem_folder(self, s):

        try:
            folder = str(s)
        except ValueError:
            warnings.warn("Unsupported format of oofem folder.")
            return

        self.globvar.oofem_folder = folder
        logger.info(f"Setting OOFEM folder to  = " + self.globvar.oofem_folder)


    def select_loading_file(self):
        #loading_file = "/home/pedro/Programs/oofem_official/python/"
        # temporary comment
        dialog = QFileDialog(self)
        dialog.setNameFilter(str("*.csv, tab-separated values [M,N] (*.csv)"))
        dialog.setViewMode(QFileDialog.Detail)
        dialog.setFileMode(QFileDialog.ExistingFile)
        if dialog.exec_():
            loading_file = dialog.selectedFiles()
        #loading_file = str(QFileDialog.getExistingDirectory(self, "Select *.csv file with [M,N] loading combinations"))
        #self.set_loading_file(loading_file[0])
        self.globvar.loading_file = loading_file[0]


    def set_loading_file(self, s):

        try:
            loading_file = str(s)
        except ValueError:
            warnings.warn("Unsupported format of loading file.")
            return
        
        if Path(loading_file).is_file():
        
            self.globvar.flag_loading_selected = True
            self.globvar.loading_file = loading_file
            logger.info(f"Setting loading file to  = " + self.globvar.loading_file)
            self.globvar.loading = pd.read_csv(loading_file, sep='\t', header=0, na_values=['nan'])
            self.checkbox_show_loading.setCheckState(Qt.Checked)
        else:

            self.globvar.flag_loading_selected = False
            self.globvar.loading_file = None
            self.globvar.loading = None

            self.checkbox_show_loading.setCheckState(Qt.Unchecked)
        
    def set_cpu_nr(self, s):

        try:
            cpu_nr = int(s)
        except ValueError:
            warnings.warn("unsupported format of cpu number")
            cpu_nr = 4
        
        self.globvar.cpu_nr = cpu_nr

    def update_ecc_selection(self):

    
        
        for checkbox in self.fem_checkboxes:
            index = self.fem_checkboxes.index(checkbox)

            checkbox_state = checkbox.checkState()
            current_task_status = self.globvar.tasks[index].status

            if (checkbox_state == Qt.Unchecked):
                if (current_task_status == Task_status.SELECTED):
                    self.globvar.tasks[index].status = Task_status.UNSELECTED
                else:
                    continue

            # Qt.Checked
            else:
                # was not selected before -> no results
                if (current_task_status == Task_status.UNSELECTED):
                    self.globvar.tasks[index].status = Task_status.SELECTED
                
    def generate_oofem_input(self):
       
        # disable changing name of the project
        self.lineEditProject.setEnabled(False)
        
        # adjusting element size (if the user forgets to do so)
        self.adjust_element_sizes()

        self.checkbox_show_mesh.setCheckState( Qt.Checked)

        self.create_finite_elements()
        
        project = Path(self.globvar.project_name) 

        # folder exists - make sure the definition has not been changed, otherwise it needs to be cleared
        if ( project.exists() ):
            if (self.globvar.flag_problem_changed):
                
                # remove entire tree
                shutil.rmtree(project)

                # create brand new project folder
                try:
                    project.mkdir(parents=False, exist_ok=False)
                except FileExistsError: 
                    self.warning_dialog_ok("folder exists and does not need to be created again" )
                
                # next change all tasks flags to correspond to the selection
                for task in self.globvar.tasks:

                    if (not task.status == Task_status.PREDEFINED):
                        task.status = Task_status.UNSELECTED
                    
                self.update_ecc_selection()
            
        # project folder does not exist and will be created
        else:
            project.mkdir(parents=False, exist_ok=False)
            
        Bx_net = self.globvar.Bx
        By_net = self.globvar.By

        # subtract concrete cover to produce mesh which approximately corresponds to the confined region
        if (self.globvar.flag_ignore_cover):
            Bx_net -= 2. * self.globvar.cover
            By_net -= 2. * self.globvar.cover
        
        for task in self.globvar.tasks:

            if ( not task.status == Task_status.PREDEFINED ):
                task_folder = project / str(task.eccentricity_normalized[1])

                if ( Path(task_folder).exists() ):
                    continue
            

            if ( task.status == Task_status.SELECTED ):
            
                # create folder
                try:
                    task_folder.mkdir(parents=False, exist_ok=False)
                except FileExistsError:
                    self.warning_dialog_ok("task folder already exists (" + task_folder + ")" )
                    return

                # modify current eccentricity
                task.eccentricity_actual = [task.eccentricity_normalized[0]*Bx_net/2., task.eccentricity_normalized[1]*By_net/2.]

                # modify analysis loading record

                # small eccentricity
                if ( task.eccentricity_normalized[1] <= 0.6 ):
                    engng_definition = "NonLinearStatic nsteps 200 stiffmode 2 rtolv 1e-4"
                    engng_definition += (f" stepLength {(1.e-4 * self.globvar.H):.6e} minStepLength {(0.125e-4 * self.globvar.H):.6e} hpc 2 {int(self.globvar.master_node_shortening)} 3 hpcw 1 -1.")
                else:
                    # significantly worse convergence due to tensile damage etc.
                    engng_definition = "NonLinearStatic nsteps 500 stiffmode 2 rtolv 5e-4"
                    engng_definition += (f" stepLength {(1.e-4 * self.globvar.H):.6e} minStepLength {(0.125e-4 * self.globvar.H):.6e} hpc 2 {int(self.globvar.master_node_bending_x)} 3 hpcw 1 -1.")
                    
                engng_definition += " Psi 0.0 MinIter 3 MaxIter 200 ReqIterations 20 lstype 4 smtype 8 initialguess 1 maxrestarts 1 renumber 1  nmodules 3\n"

                task.engng_definition = engng_definition

                # create oofem input file and write
                file_path = task_folder / self.globvar.filename_in            
 
                self.write_oofem_input(file_path, task)

                #task.file_path = file_path
                task.file_path = task_folder

        # CHANGE FLAG        
        self.globvar.flag_output_generated = True
        self.globvar.flag_problem_changed = False
  
    def write_oofem_input(self, file_path, task):

        # create project folder, name needs to be unique-> guarantee that folder is empty

        with file_path.open("w", encoding ="utf-8") as f:
            #f = open(self.globvar.filename_in, "w")
            oofem_input.write_oofem_input_header(self.globvar, task, f)
            oofem_input.write_oofem_engineering_model_modules_domain(self.globvar, task, f)
            oofem_input.write_oofem_input_mesh_nodes(self.globvar, task, f)
            oofem_input.write_oofem_input_hanging_nodes(self.globvar, task, f)
            oofem_input.write_oofem_input_dummy_nodes(self.globvar, task, f)
            oofem_input.write_oofem_input_brick_elements(self.globvar, task, f)
            oofem_input.write_oofem_input_truss_elements(self.globvar, task, f)
            oofem_input.write_oofem_input_cross_sections(self.globvar, task, f)
            oofem_input.write_oofem_input_materials(self.globvar, task, f)
            oofem_input.write_oofem_input_boundary_conditions(self.globvar, task, f)
            oofem_input.write_oofem_input_time_functions(self.globvar, task, f) 
            oofem_input.write_oofem_input_elements_sets(self.globvar, task, f)
            oofem_input.write_oofem_input_extractor_record(self.globvar, task, f)
            #f.close()
        

        
    def report_progress(self, step, load):
        self.console.appendPlainText(f"Finished time step {(step)} at load level = {(load):.2f} MN")
        self.progressBar.setValue(load)


    def set_FEM_results(self, task, dummy_eps, dummy_kappa):
       
        for w in dummy_eps:
            task.eps.append( w / self.globvar.H )

        for w in dummy_kappa:
            task.kappa.append( w / (self.globvar.H * self.globvar.By/2.) )

        for load in task.load_level:
            task.N.append(-load)
            task.M.append( load * task.eccentricity_actual[1] )

        task.max_MN = [task.max_load * task.eccentricity_actual[1], -task.max_load]

        self.update_diagram_plot()


    def execute_next_task(self,tasks):


       
        searching = True

        while (searching):
            task = next(tasks,"end")
            if (task == "end"):
                self.write_console_report()
                self.progressBar.setValue(0.)
                logger.info("Finished all selected tasks")
                return
            if (task.status == Task_status.SELECTED):
                searching = False
                      
        logger.info(f"Starting FEM analysis for eccentricity = {(task.eccentricity_normalized[1]):.2f} [-]")
        self.console.setPlainText(f"   Starting FEM analysis for eccentricity = {(task.eccentricity_normalized[1]):.2f} [-]")
        estimate = self.diagram_MSR_simple.find_load_for_eccentricity(task.eccentricity_actual[1])
        self.console.appendPlainText(f"   estimated load =  {(estimate):.3f} MN")
        self.console.appendPlainText("---------------------------------------------------------------------------")

        self.progressBar.setMaximum(estimate)
        
        self.thread = QThread()
        self.worker = Worker(task,self.globvar.oofem_folder,self.globvar.cpu_nr)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
                
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.worker.progress.connect(self.report_progress)

        self.worker.result.connect( self.set_FEM_results )
        self.thread.finished.connect( lambda: self.execute_next_task(tasks) )
        
        self.thread.start()
        
    def run_oofem_problems(self):
        # multiprocessing - works but window freezes 
        #https://superfastpython.com/run-function-in-new-process/
        #maybe in the future
        #https://stackoverflow.com/questions/26833093/how-to-terminate-qthread-in-python
        #https://realpython.com/python-pyqt-qthread/#using-qthread-to-prevent-freezing-guis

        logger.info("Started FEM analysis")
        
        if ( len(self.globvar.oofem_folder) ):
            oofem = Path(self.globvar.oofem_folder) 

            if ( oofem.exists() ):
                self.globvar.flag_oofem_selected = True
            else:
                self.globvar.flag_oofem_selected = False
        
        if ( not self.globvar.flag_oofem_selected ):
            message = "Invalid or not selected OOFEM folder."
            warnings.warn(message)
            self.warning_dialog_ok(message)
            return
        
        elif (not self.globvar.flag_output_generated):
            message = "FEM input files have not been generated."
            warnings.warn(message)
            self.warning_dialog_ok(message)
            return
        
        else:
            
            self.thread = QThread()
            tasks = iter(self.globvar.tasks)
            self.execute_next_task(tasks)
            self.globvar.flag_analyses_run = True       


    def change_status_to_progress(self,task):
        task.status = Task_status.PROGRESS
        
    def write_console_report(self):

        self.console.setPlainText("     SUMMARY OF FEM RESULTS:   ")
        self.console.appendPlainText("-----------------------------------------------------------------------")

        self.console.appendPlainText("e\u0302 [-]\tM [MNm]\tN [MN]\t% estimate")
        
        for task in self.globvar.tasks:
            if ( task.status == Task_status.COMPLETED ):
                estimate = self.diagram_MSR_simple.find_load_for_eccentricity(task.eccentricity_actual[1])
                self.console.appendPlainText(f"{(task.eccentricity_normalized[1]):.3f}\t{(task.max_MN[0]):.2f}\t{(task.max_MN[1]):.3f}\t{(100.*task.max_load/estimate):.1f}")
            
            elif ( task.status == Task_status.PREDEFINED):
                self.console.appendPlainText(f"tension\t{(task.max_MN[0]):.3f}\t{(task.max_MN[1]):.3f}")

        self.console.appendPlainText("-----------------------------------------------------------------------")               


if __name__ == "__main__":        
    # LOGGER INITIALIZATION    
    # logging warnings
    #  https://code-maven.com/python-warnings

 
    
    logger = logging.getLogger('py.warnings')
    logger.setLevel(logging.DEBUG)
    logging.captureWarnings(True)

    # console output
    sh = logging.StreamHandler()
    sh.setLevel(logging.DEBUG)
    sh.setFormatter(logging.Formatter('%(asctime)s %(levelname)-8s: %(message)s',
                                  "%Y-%m-%d %H:%M:%S") )
    # file output
    fh = logging.FileHandler(filename='malcolm.log')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s %(levelname)-8s: %(message)s',
                                  "%Y-%m-%d %H:%M:%S") )
    # add both outputs
    logger.addHandler(fh)
    logger.addHandler(sh)

    logger.info("Program started")

    
    
    app = QApplication(sys.argv)
    start = time()

    pixmap = QPixmap("malcolm_logo_text.png")

    
    splash = QSplashScreen(pixmap, Qt.WindowStaysOnTopHint)
    #splash.resize(1000,600)
    splash.resize(1280,758)
    splash.setGeometry(
        QStyle.alignedRect(
        Qt.LeftToRight,
        Qt.AlignCenter,
        splash.size(),
        QGuiApplication.primaryScreen().availableGeometry(),
        ),
    )
    splash.show()

    message = "(c) Petr Havlásek 2022, credits: Alena Plačková (Malcolm logo)\n"
    message += "Acknowledgment: TAČR project nr. TM01000059\n"
    message += "Malcolm v. 2.0\n"
    splash.showMessage(message, Qt.AlignLeft| Qt.AlignTop, Qt.black)
    
    while time() - start < 2.:
        sleep(0.001)
        app.processEvents()
    
    w = MainWindow()
    w.setFixedSize(1280, 720)
    w.setGeometry(
        QStyle.alignedRect(
        Qt.LeftToRight,
        Qt.AlignCenter,
        w.size(),
        QGuiApplication.primaryScreen().availableGeometry(),
        ),
    )
    w.show()
    
    splash.finish(w)

    app.setWindowIcon(QIcon('malcolm_icon.png'))

    #splash.raise_()
    sys.exit(app.exec_())

    #logger.info("Program finished")
    #app.closeAllWindows()
    #app.exec_()
    
 
