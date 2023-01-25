import math

import globvars
import warnings

from PySide2.QtWidgets import QWidget, QHBoxLayout, QVBoxLayout, QGridLayout, QLabel, QDesktopWidget, QLineEdit, QDialogButtonBox, QPushButton, QMessageBox
from PySide2.QtGui import QGuiApplication
from PySide2.QtCore import Qt

class material:
    def __init__(self, globvar, nr):

        self.nr = nr
        self.material_type = None

    def give_material_type(self):
        return self.material_type
       


class concrete_mat(material):
    
    def __init__(self, globvar, nr, fcm = None):
        super().__init__(globvar, nr)

        self.fcm = globvar.concretes[globvar.fcm]
                
        self.material_type = globvar.material_type[0]
        self.predict_concrete_parameters()

        

    def predict_concrete_parameters(self):

        # fcm = mean compressive strength on cylinders [MPa]
           
        # characteristic value of compressive strength [MPa]
        fck = self.fcm-8.;
        # mean tensile strength [MPa], math.log = natural logarithm
        if (fck>50.):
            self.fctm = 2.12*math.log(1.+self.fcm/10.);
        else:
            self.fctm = 0.3*fck**(2./3.);

        # fracture energy [N/m]
        self.Gf = 73.*self.fcm**0.18;
        # initial modulus [MPa]
        self.Eci = 21.5e3*(self.fcm/10.)**(1./3.);
        # critical crack opening [m]
        self.wf = 4.444*self.Gf/self.fctm*1.e-6;
        # Poisson's ratio
        self.nu = 0.2;
        # eccentricity
        self.ecc = 0.525;
        # initial value of hardening variable
        self.kinit = 0.3;
        # hardening modulus
        self.Hp = 0.01;
        # dilation factor
        self.dilation = 0.85;
        # hardening parameter Ahard
        self.Ahard = 0.08;
        # hardening parameter Ahard
        self.Bhard = 0.003;
        # hardening parameter Ahard
        self.Chard = 2.;
        # hardening parameter Ahard
        self.Dhard = 1.e-6;
        # softening parameter
        self.Asoft = 15;
        # softening parameter for compression
        self.efc = 1.e-4;
        # type of softening (1=bilinear)
        # stype = 1;
        # softening parameter for tension
        self.wf1 = 0.15;
        # softening parameter for tension
        self.ft1 = 0.3;


    def update_concrete_strength(self, s):

      
        try:
            strength = float(s)
        except ValueError:
                       
            warnings.warn("unsupported value, using fcm = 30 MPa")
            strength = 30.
            
        self.fcm = strength


class window_CDPM2(QWidget):

    def __init__(self, cdpm2):
        super().__init__()

        self.setWindowTitle("CDPM2 material parameters")

        param_length = 8
        
        label_fcm = QLabel("fcm = ")
        self.label_edit_fcm = QLabel(str(cdpm2.fcm))
        label_fcm_unit = QLabel("[MPa]")

        label_number = QLabel("Definition of material number: " + str(cdpm2.nr) )
        
        label_fctm = QLabel("fctm = ")
        self.line_edit_fctm = QLineEdit()
        self.line_edit_fctm.setMaxLength(param_length)
        self.line_edit_fctm.setText(str(cdpm2.fctm))
        label_fctm_unit = QLabel("[MPa]")

        label_Gf = QLabel("Gf = ")
        self.line_edit_Gf = QLineEdit()
        self.line_edit_Gf.setMaxLength(param_length)
        self.line_edit_Gf.setText(str(cdpm2.Gf))
        label_Gf_unit = QLabel("[N/m]")

        label_Eci = QLabel("Eci = ")
        self.line_edit_Eci = QLineEdit()
        self.line_edit_Eci.setMaxLength(param_length)
        self.line_edit_Eci.setText(str(cdpm2.Eci))
        label_Eci_unit = QLabel("[MPa]")

        label_nu = QLabel("nu = ")
        self.line_edit_nu = QLineEdit()
        self.line_edit_nu.setMaxLength(param_length)
        self.line_edit_nu.setText(str(cdpm2.nu))
        label_nu_unit = QLabel("[-]")

        label_wf = QLabel("wf = ")
        self.line_edit_wf = QLineEdit()
        self.line_edit_wf.setMaxLength(param_length)
        self.line_edit_wf.setText(str(cdpm2.wf))
        label_wf_unit = QLabel("[m]")
        
        label_ecc = QLabel("ecc = ")
        self.line_edit_ecc = QLineEdit()
        self.line_edit_ecc.setMaxLength(param_length)
        self.line_edit_ecc.setText(str(cdpm2.ecc))
        label_ecc_unit = QLabel("[m]")

        label_kinit = QLabel("kinit = ")
        self.line_edit_kinit = QLineEdit()
        self.line_edit_kinit.setMaxLength(param_length)
        self.line_edit_kinit.setText(str(cdpm2.kinit))
        label_kinit_unit = QLabel("[-]")

        label_Hp = QLabel("Hp = ")
        self.line_edit_Hp = QLineEdit()
        self.line_edit_Hp.setMaxLength(param_length)
        self.line_edit_Hp.setText(str(cdpm2.Hp))
        label_Hp_unit = QLabel("[-]")

        label_dilation = QLabel("dilation = ")
        self.line_edit_dilation = QLineEdit()
        self.line_edit_dilation.setMaxLength(param_length)
        self.line_edit_dilation.setText(str(cdpm2.dilation))
        label_dilation_unit = QLabel("[-]")

        label_Ahard = QLabel("Ahard = ")
        self.line_edit_Ahard = QLineEdit()
        self.line_edit_Ahard.setMaxLength(param_length)
        self.line_edit_Ahard.setText(str(cdpm2.Ahard))
        label_Ahard_unit = QLabel("[-]")

        label_Bhard = QLabel("Bhard = ")
        self.line_edit_Bhard = QLineEdit()
        self.line_edit_Bhard.setMaxLength(param_length)
        self.line_edit_Bhard.setText(str(cdpm2.Bhard))
        label_Bhard_unit = QLabel("[-]")
        
        label_Chard = QLabel("Chard = ")
        self.line_edit_Chard = QLineEdit()
        self.line_edit_Chard.setMaxLength(param_length)
        self.line_edit_Chard.setText(str(cdpm2.Chard))
        label_Chard_unit = QLabel("[-]")

        label_Dhard = QLabel("Dhard = ")
        self.line_edit_Dhard = QLineEdit()
        self.line_edit_Dhard.setMaxLength(param_length)
        self.line_edit_Dhard.setText(str(cdpm2.Dhard))
        label_Dhard_unit = QLabel("[-]")

        label_Asoft = QLabel("Asoft = ")
        self.line_edit_Asoft = QLineEdit()
        self.line_edit_Asoft.setMaxLength(param_length)
        self.line_edit_Asoft.setText(str(cdpm2.Asoft))
        label_Asoft_unit = QLabel("[-]")
        
        label_efc = QLabel("efc = ")
        self.line_edit_efc = QLineEdit()
        self.line_edit_efc.setMaxLength(param_length)
        self.line_edit_efc.setText(str(cdpm2.efc))
        label_efc_unit = QLabel("[-]")

        label_wf1 = QLabel("wf1 = ")
        self.line_edit_wf1 = QLineEdit()
        self.line_edit_wf1.setMaxLength(param_length)
        self.line_edit_wf1.setText(str(cdpm2.wf1))
        label_wf1_unit = QLabel("[-]")

        label_ft1 = QLabel("ft1 = ")
        self.line_edit_ft1 = QLineEdit()
        self.line_edit_ft1.setMaxLength(param_length)
        self.line_edit_ft1.setText(str(cdpm2.ft1))
        label_ft1_unit = QLabel("[-]")

        
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Close)
        self.buttonBox.accepted.connect(lambda: self.update_cdpm2_user_parameters(cdpm2) )
        self.buttonBox.accepted.connect(self.close)
        self.buttonBox.rejected.connect(self.close)
        
        cdpm2_layout = QGridLayout()

        row = 0
        cdpm2_layout.addWidget(label_number, row, 0, 1, 3)

        row += 1
        cdpm2_layout.addWidget(label_fcm, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.label_edit_fcm, row, 1)
        cdpm2_layout.addWidget(label_fcm_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_fctm, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_fctm, row, 1)
        cdpm2_layout.addWidget(label_fctm_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_Gf, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_Gf, row, 1)
        cdpm2_layout.addWidget(label_Gf_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_Eci, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_Eci, row, 1)
        cdpm2_layout.addWidget(label_Eci_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_nu, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_nu, row, 1)
        cdpm2_layout.addWidget(label_nu_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_wf, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_wf, row, 1)
        cdpm2_layout.addWidget(label_wf_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_ecc, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_ecc, row, 1)
        cdpm2_layout.addWidget(label_ecc_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_kinit, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_kinit, row, 1)
        cdpm2_layout.addWidget(label_kinit_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_Hp, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_Hp, row, 1)
        cdpm2_layout.addWidget(label_Hp_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_dilation, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_dilation, row, 1)
        cdpm2_layout.addWidget(label_dilation_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_Ahard, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_Ahard, row, 1)
        cdpm2_layout.addWidget(label_Ahard_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_Bhard, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_Bhard, row, 1)
        cdpm2_layout.addWidget(label_Bhard_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_Chard, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_Chard, row, 1)
        cdpm2_layout.addWidget(label_Chard_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_Dhard, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_Dhard, row, 1)
        cdpm2_layout.addWidget(label_Dhard_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_Asoft, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_Asoft, row, 1)
        cdpm2_layout.addWidget(label_Asoft_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_efc, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_efc, row, 1)
        cdpm2_layout.addWidget(label_efc_unit, row, 2)
        
        row += 1
        cdpm2_layout.addWidget(label_wf1, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_wf1, row, 1)
        cdpm2_layout.addWidget(label_wf1_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(label_ft1, row, 0, 1, 1, Qt.AlignRight)
        cdpm2_layout.addWidget(self.line_edit_ft1, row, 1)
        cdpm2_layout.addWidget(label_ft1_unit, row, 2)

        row += 1
        cdpm2_layout.addWidget(self.buttonBox, row, 0, 1, 3)

        self.setLayout(cdpm2_layout)

       
        
    def check_float_value(self, value):
        try:
            test = float(value)
            return True
        except:
            QMessageBox.warning(self, "Warning", "unsupported value:" + str(value), buttons=QMessageBox.Ok, defaultButton=QMessageBox.Ok)
            return False


    def update_text_fields(self, cdpm2):
        #self.line_edit_fctm.setText(str(cdpm2.fctm))
        self.line_edit_Gf.setText(str(cdpm2.Gf))
        self.line_edit_Eci.setText(str(cdpm2.Eci))
        self.line_edit_wf.setText(str(cdpm2.wf))
        self.line_edit_nu.setText(str(cdpm2.nu))
        self.line_edit_ecc.setText(str(cdpm2.ecc))
        self.line_edit_kinit.setText(str(cdpm2.kinit))
        self.line_edit_Hp.setText(str(cdpm2.Hp))
        self.line_edit_dilation.setText(str(cdpm2.dilation))
        self.line_edit_Ahard.setText(str(cdpm2.Ahard))
        self.line_edit_Bhard.setText(str(cdpm2.Bhard))
        self.line_edit_Chard.setText(str(cdpm2.Chard))
        self.line_edit_Dhard.setText(str(cdpm2.Dhard))
        self.line_edit_Asoft.setText(str(cdpm2.Asoft))
        self.line_edit_efc.setText(str(cdpm2.efc))
        self.line_edit_wf1.setText(str(cdpm2.wf1))
        self.line_edit_ft1.setText(str(cdpm2.ft1))


        
    def update_cdpm2_user_parameters (self, cdpm2):

        try:
            cdpm2.fctm = float( self.line_edit_fctm.text() )
            cdpm2.Gf = float( self.line_edit_Gf.text() )
            cdpm2.Eci = float( self.line_edit_Eci.text() )
            cdpm2.wf = float( self.line_edit_wf.text() )
            cdpm2.nu = float( self.line_edit_nu.text() )
            cdpm2.ecc = float( self.line_edit_ecc.text() )
            cdpm2.kinit = float( self.line_edit_kinit.text() )
            cdpm2.Hp = float( self.line_edit_Hp.text() )
            cdpm2.dilation = float( self.line_edit_dilation.text() )
            cdpm2.Ahard= float( self.line_edit_Ahard.text() )
            cdpm2.Bhard = float( self.line_edit_Bhard.text() )
            cdpm2.Chard = float( self.line_edit_Chard.text() )
            cdpm2.Dhard = float( self.line_edit_Dhard.text() )
            cdpm2.Asoft = float( self.line_edit_Asoft.text() )
            cdpm2.efc = float( self.line_edit_efc.text() )
            cdpm2.wf1 = float( self.line_edit_wf1.text() )
            cdpm2.ft1 = float( self.line_edit_ft1.text() )

        except ValueError:
            warnings.warn("incorrectly defined user value" )
            return




class rebar_mat(material):
    
    def __init__(self, globvar, nr, mat):
        super().__init__(globvar, nr)

        self.material_type = globvar.material_type[1]
        self.predict_steel_parameters(globvar,mat)
        self.sig_0 = 500;

        
    def predict_steel_parameters(self, globvar, mat = None):
        # TODO

        if ( mat == globvar.rebarmat[0]):
            self.E = 210.e3
            self.H = 0.;
            self.omega_c = 0.;
            self.a = 0.;
        
        elif ( mat == globvar.rebarmat[1]):
            self.E = 210.e3
            self.H = 7.5e3;
            self.omega_c = 1.;
            self.a = 6.;

        elif ( mat == globvar.rebarmat[2]):
            self.E = 210.e3
            self.H = 8.2e3;
            self.omega_c = 1.;
            self.a = 6.;

        elif ( mat == globvar.rebarmat[3]):
            self.sig_0 = 500;
            self.H = 9.e3;
            self.omega_c = 1.;
            self.a = 13.;

        elif ( mat == globvar.rebarmat[4]):
            self.E = 210.e3
            self.H = 5.5e3;
            self.omega_c = 1.;
            self.a = 7.25;

        elif ( mat == globvar.rebarmat[5]):
            self.E = 210.e3
            self.H = 5.5e3;
            self.omega_c = 1.;
            self.a = 6.15;

        elif ( mat == globvar.rebarmat[6]):
            self.E = 210.e3
            self.H = 0.;
            self.omega_c = 0.;
            self.a = 0.;

        else:
            warnings.warn("unknown material definition, using fy = 500 MPa")
            self.E = 210.e3
            self.sig_0 = 500;
            self.H = 0.;
            self.omega_c = 0.;
            self.a = 0.;
            

class tendon_mat(material):
    # TODO
    
    def __init__(self, globvar, nr, mat):
        super().__init__(globvar, nr)

        self.material_type = globvar.material_type[2]
        self.predict_steel_parameters(mat)

    def predict_steel_parameters(self, mat = None):

        if ( mat == globvar.tendonmat[0]):
            self.E = 200.e3
            self.sig_0 = 1800.;
            self.H = 0.;
            self.omega_c = 0.;
            self.a = 0.;

        else:
            warnings.warn("unknown material definition, using fy = 1800 MPa")
            self.E = 200.e3
            self.sig_0 = 1800.;
            self.H = 0.;
            self.omega_c = 0.;
            self.a = 0.;



class window_Mises(QWidget):

    def __init__(self, mises):
        super().__init__()

        self.setWindowTitle("Mises material parameters")

        param_length = 8

        label_number = QLabel("Definition of material number: " + str(mises.nr) )

        label_sig_0 = QLabel("sig_0 = ")
        self.label_edit_sig_0 = QLabel(str(mises.sig_0))
        label_sig_0_unit = QLabel("[MPa]")
        
        label_E = QLabel("E = ")
        self.line_edit_E = QLineEdit()
        self.line_edit_E.setMaxLength(param_length)
        self.line_edit_E.setText(str(mises.E))
        label_E_unit = QLabel("[MPa]")

        label_H = QLabel("H = ")
        self.line_edit_H = QLineEdit()
        self.line_edit_H.setMaxLength(param_length)
        self.line_edit_H.setText(str(mises.H))
        label_H_unit = QLabel("[MPa]")

        label_omega_c = QLabel("omega_c = ")
        self.line_edit_omega_c = QLineEdit()
        self.line_edit_omega_c.setMaxLength(param_length)
        self.line_edit_omega_c.setText(str(mises.omega_c))
        label_omega_c_unit = QLabel("[-]")

        label_a = QLabel("a = ")
        self.line_edit_a = QLineEdit()
        self.line_edit_a.setMaxLength(param_length)
        self.line_edit_a.setText(str(mises.a))
        label_a_unit = QLabel("[-]")

        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Close)
        self.buttonBox.accepted.connect(lambda: self.update_mises_user_parameters(mises) )
        self.buttonBox.accepted.connect(self.close)
        self.buttonBox.rejected.connect(self.close)
        
        mises_layout = QGridLayout()

        row = 0
        mises_layout.addWidget(label_number, row, 0, 1, 3)
        
        row += 1
        mises_layout.addWidget(label_sig_0, row, 0, 1, 1, Qt.AlignRight)
        mises_layout.addWidget(self.label_edit_sig_0, row, 1)
        mises_layout.addWidget(label_sig_0_unit, row, 2)
        
        row += 1
        mises_layout.addWidget(label_E, row, 0, 1, 1, Qt.AlignRight)
        mises_layout.addWidget(self.line_edit_E, row, 1)
        mises_layout.addWidget(label_E_unit, row, 2)

        row += 1
        mises_layout.addWidget(label_H, row, 0, 1, 1, Qt.AlignRight)
        mises_layout.addWidget(self.line_edit_H, row, 1)
        mises_layout.addWidget(label_H_unit, row, 2)

        row += 1
        mises_layout.addWidget(label_omega_c, row, 0, 1, 1, Qt.AlignRight)
        mises_layout.addWidget(self.line_edit_omega_c, row, 1)
        mises_layout.addWidget(label_omega_c_unit, row, 2)

        row += 1
        mises_layout.addWidget(label_a, row, 0, 1, 1, Qt.AlignRight)
        mises_layout.addWidget(self.line_edit_a, row, 1)
        mises_layout.addWidget(label_a_unit, row, 2)

        row += 1
        mises_layout.addWidget(self.buttonBox, row, 0, 1, 3)

        self.setLayout(mises_layout)


    def update_mises_user_parameters (self, mises):

        try:
            mises.E = float( self.line_edit_E.text() )
            mises.H = float( self.line_edit_H.text() )
            mises.omega_c = float( self.line_edit_omega_c.text() )
            mises.a = float( self.line_edit_a.text() )

        except ValueError:
            warnings.warn("incorrectly defined user value" )
            return

        
