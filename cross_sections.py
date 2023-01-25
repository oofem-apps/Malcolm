import math

class cross_section:
    def __init__(self, cs, mat):
        self.nr = cs
        self.mat = mat
        # list of element numbers
        self.elements = []
        # number of element set coincides with cross-section number
        self.elem_set = cs

        self.cs_types = ("bulk", "rebar")
        self.cs_type = self.cs_types[0]
       

class rebar_cross_section(cross_section):
    def __init__(self, cs, mat, diam, area, cs_tag = None):
        super().__init__(cs, mat)
        
        self.diam = diam
        self.area = area
        self.cs_tag = cs_tag
        
        self.cs_type = self.cs_types[1]
        
        if (area is None):
            self.area = math.pi*(self.diam/2.)**2

    def give_diameter(self):
        return self.diam

    def give_area(self):
        return self.area

    def set_properties_from_cs_tag(self, globvar, cs_tag):
        self.cs_tag = cs_tag
        self.diam = (globvar.rebars_CS[self.cs_tag]).diam
        self.area = (globvar.rebars_CS[self.cs_tag]).area
        if (globvar.debug_flag):
            print("updated diameter D = " + str(self.diam) + " A = " +str(self.area) )
        
class rebar_CS:
    # rebar-specific cross-section - just a container
    def __init__(self, diam, area = None):
        
        self.diam = diam
        self.area = area
                
        if (area is None):
            self.area = math.pi*(self.diam/2.)**2
        
