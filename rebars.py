import math
import numpy as np

import warnings

from cross_sections import rebar_cross_section
from nodes import node
from elements import truss_elem

class rebar:
    def __init__(self, globvar, XYZ, cs_nr, cs_tag, cs_mat, div_length = None, tag = None):

        self.rebar_type = ("rebar", "hoop", "spiral")

        self.reb_type = self.rebar_type[0]
        
        # control points coordinates (or center/origin for a spiral and hoop reinforcement)
        self.XYZ = XYZ
       
        # number of cross-section
        self.cs_nr = cs_nr

        # label, tag
        self.tag = tag

        # division length
        if (div_length is None):
            self.div_length = globvar.global_div_length
            
        else:
            self.div_length = div_length

        # empty list of point objects
        self.rebar_points = []

        # add cross-sections
        do_create_cs = True

      
        if globvar.cross_sections:
            for i_cs in globvar.cross_sections:

                # rebar cross-section tag (i.e. #4)
                if (i_cs.cs_type == "rebar"):
                    if ( (i_cs.nr == self.cs_nr) and (i_cs.mat == cs_mat) and (i_cs.cs_tag == cs_tag) ):
                        do_create_cs = False
                        break

        if (do_create_cs):
            diam =  (globvar.rebars_CS[cs_tag]).diam
            area =  (globvar.rebars_CS[cs_tag]).area
            truss_cs = rebar_cross_section(cs_nr, cs_mat, diam, area, cs_tag ) 
            globvar.cross_sections.append(truss_cs)

    def give_rebar_type(self):
        return self.reb_type

    def find_cross_section(self, globvar, nr):

        for i_cs in globvar.cross_sections:
            if (i_cs.nr == nr):
                break

        return i_cs


    def give_area(self, globvar):

        cs = self.find_cross_section(globvar, self.cs_nr)
        area = cs.give_area()
        return area

    def give_diameter(self, globvar):
        cs = self.find_cross_section(globvar, self.cs_nr)
        diam = cs.give_diameter()
        return diam

    def give_mat_nr(self, globvar):
        cs = self.find_cross_section(globvar, self.cs_nr)
        mat_nr = cs.mat
        return mat_nr
    
    def give_cs_tag(self, globvar):
        cs = self.find_cross_section(globvar, self.cs_nr)
        cs_tag = cs.cs_tag
        return cs_tag    

    
    def compute_rebar_volume(self, globvar):

        reinf_length = self.compute_rebar_length()
        area = self.give_area(globvar)
        
        vol_rebar = reinf_length * area
        return vol_rebar
        

    def compute_rebar_length(self):

        reinf_length = 0
        prev_xyz = None
        for xyz in self.XYZ:
            if (prev_xyz == None):
                prev_xyz = xyz
                continue

            else:
                reinf_length += math.sqrt ( (xyz[0]-prev_xyz[0])**2 +  (xyz[1]-prev_xyz[1])**2 +  (xyz[2]-prev_xyz[2])**2 )
                prev_xyz = xyz

        return reinf_length
    
    
    def create_rebar_nodes(self):

        for seg in range( len(self.XYZ)-1 ):

            seg_length = math.sqrt( (self.XYZ[seg+1][0] - self.XYZ[seg][0])**2 + (self.XYZ[seg+1][1] - self.XYZ[seg][1])**2 + (self.XYZ[seg+1][2] - self.XYZ[seg][2])**2 )

            div_nr = round(seg_length / self.div_length)
            
            # unit vector
            vec = [ (self.XYZ[seg+1][0] - self.XYZ[seg][0]) / seg_length,
                    (self.XYZ[seg+1][1] - self.XYZ[seg][1]) / seg_length,
                    (self.XYZ[seg+1][2] - self.XYZ[seg][2]) / seg_length ]

       
            # number of points = number of elements + 1
            for i in range(div_nr):

                # new_point = self.XYZ0 + [x*i/self.div_nr for x in vec]
                new_point =  [ origin_i + vec_i for origin_i, vec_i in zip( self.XYZ[seg], [x*seg_length*i/div_nr for x in vec] ) ]
                # print(new_point)
                self.rebar_points.append(new_point)

        # add last ending point
        new_point = self.XYZ[ len(self.XYZ)-1 ]
        self.rebar_points.append(new_point)  


    #    def create_rebar_elements(self, last_node_nr, last_elem_nr):
    def create_rebar_elements(self, globvar):

        node_nr = self.find_last_node_nr(globvar)
        elem_nr = self.find_last_elem_nr(globvar)
        
        # create first "solo" node
        node_nr += 1
        new_node = node( nr = node_nr, XYZ = self.rebar_points[0], tag = self.tag )
        globvar.hanging_nodes.append(new_node)
        # global counter
        globvar.ndofman += 1
        
        # create elements + remaining nodes
        for point in self.rebar_points[1:]:
            elem_nr += 1
            new_truss = truss_elem ( globvar = globvar, nr = elem_nr, nodes = [node_nr, node_nr+1], cs = self.cs_nr, tag = self.tag )
            
            globvar.truss_elements.append(new_truss)
            # global counter
            globvar.nelem += 1
                                     
            node_nr += 1
            new_node = node( nr = node_nr, XYZ = point, tag = self.tag )
            globvar.hanging_nodes.append(new_node)
            # global counter
            globvar.ndofman += 1


    def find_last_node_nr(self, globvar):
        
        if (globvar.hanging_nodes == []):
            last_node_nr = globvar.ndofman
        else:
            last_node_nr = globvar.hanging_nodes[-1].nr

        return last_node_nr

    def find_last_elem_nr(self, globvar):
        
        if (globvar.hanging_nodes == []):
            last_elem_nr = globvar.nelem
        else:
            last_elem_nr = globvar.truss_elements[-1].nr

        return last_elem_nr

    def check_topology_consistency(self, globvar):
        for xyz in self.XYZ:
            diam = self.give_diameter(globvar)
            limit_x = globvar.Bx/2. - globvar.cover - diam/2.
            limit_y = globvar.By/2. - globvar.cover - diam/2.

            if xyz[0] < -limit_x:
                warnings.warn("specified topology beyond allowed bounds, x = " + str(xyz[0]))
                xyz[0] = -limit_x
            elif xyz[0] > limit_x:
                warnings.warn("specified topology beyond allowed bounds, x = " + str(xyz[0]))
                xyz[0] = limit_x
            else:
                pass

            if xyz[1] < -limit_y:
                warnings.warn("specified topology beyond allowed bounds, y = " + str(xyz[1]))
                xyz[1] = -limit_y
            elif xyz[1] > limit_y:
                warnings.warn("specified topology beyond allowed bounds, y = " + str(xyz[1]))
                xyz[1] = limit_y
            else:
                pass
           
            if xyz[2] < 0.:
                warnings.warn("specified topology beyond allowed bounds, z = " + str(xyz[2]))
                xyz[2] = 0.
            elif xyz[2] > globvar.H:
                warnings.warn("specified topology beyond allowed bounds, z = " + str(xyz[2]))
                xyz[2] = globvar.H
            else:
                pass

    
# hoops reinforcement - can be connected or even partial
class hoop(rebar):
    def __init__(self, globvar, XYZ, cs_nr, cs_tag, cs_mat, vec, start = None, angle = None, radius = None, div_length = None, tag = None, connect = False):
        super().__init__(globvar, XYZ, cs_nr, cs_tag, cs_mat, div_length, tag)

        self.reb_type = self.rebar_type[1]
        
        # start = starting point of the hoop reinforcement
        # angle = specifies length of the arch, in degrees
        # radius = axial radius of the curved reinforcement
        # connect = deside whether the ends of the reinforcement should join or not
        
        # vector normal to the plane where hoop reinforcement is defined
        # in the case of spiral reinforcement it is the direction of the spiral

        # reb_type = super.rebar_type[1]
                
        norm_vec = 0
        for vec_i in vec:
            norm_vec += vec_i**2 
        norm_vec = math.sqrt(norm_vec)

        # normalize normal vector 
        self.vec = [ vec_i/norm_vec for vec_i in vec ]
        
        # starting point for not entire circles
        if (start is not None):
            self.start = start

            start_vec = [self.start[0] - self.XYZ[0],
                         self.start[1] - self.XYZ[1],
                         self.start[2] - self.XYZ[2] ]


            norm_start = 0
            for vec_i in start_vec:
                norm_start += vec_i**2 
            norm_start = math.sqrt(norm_start)

            # vec should be already normalized
            start_dot_vec = 0.
            for start_vec_i, vec_i in zip(start_vec, self.vec):
                start_dot_vec += start_vec_i * vec_i 

            # angle between vec (plane normal) and starting point of the reinforcement
            theta = math.acos (  start_dot_vec  / (norm_start * 1.) )

            self.radius = norm_start * math.sin(theta)
            
        else:

            self.radius = radius

            if ( self.vec[2] == 0. ):

                start_vec = [ self.radius * self.vec[1],
                              -self.radius * self.vec[0],
                              0. ]

            else: 
                norm = math.sqrt(1 + self.vec[1]**2 / self.vec[2]**2 )              
                
                start_vec = [ 0.,
                          1. * self.radius/norm,
                          -self.vec[1]/self.vec[2] * self.radius/norm ]

            self.start = [self.XYZ[0] + start_vec[0],
                          self.XYZ[1] + start_vec[1],
                          self.XYZ[2] + start_vec[2] ]
            
            
        if (angle is not None):
            # angle to be specified in degrees and is converted to radians
            self.angle = angle/180. * math.pi
        else:
            # full circle assumed
            self.angle = 2. * math.pi


            
        self.connect = connect


    # https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    def rotation_matrix(self, axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis)
        axis = axis / math.sqrt(np.dot(axis, axis))
        a = math.cos(theta / 2.0)
        b, c, d = -axis * math.sin(theta / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                         [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                         [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

    def create_rebar_nodes(self):
 
        start_vec = [self.start[0] - self.XYZ[0],
                     self.start[1] - self.XYZ[1],
                     self.start[2] - self.XYZ[2] ]


        reinf_length = self.angle * self.radius

        div_nr = round(reinf_length / self.div_length)

        div_angle = self.angle / div_nr

        # number of points = number of elements + 1
        number_of_nodes = div_nr + 1

        # subtract one node in the case of a full hoop
        if ( ( self.angle == 2. * math.pi ) or ( self.connect ) ):
            number_of_nodes -= 1
        
        for i in range(number_of_nodes):

            alpha = i/div_nr * self.angle
            rotated = np.dot( self.rotation_matrix(self.vec, alpha), start_vec )
                             
            new_point = [ origin_i + rot_i for origin_i, rot_i in zip (self.XYZ, rotated) ]

            self.rebar_points.append(new_point)
            

    def compute_rebar_length(self):

        reinf_length = 2. * math.pi * self.radius

        return reinf_length
    

    def create_rebar_elements(self, globvar):
        super().create_rebar_elements(globvar)

        if ( self.connect ):

            elem_nr = globvar.truss_elements[-1].nr
            # last node of the rebar
            end_node_nr = globvar.hanging_nodes[-1].nr
            # first node of the rebar
            start_node_nr = last_node_nr + 1

            elem_nr += 1
            new_truss = truss_elem ( globvar = globvar, nr = elem_nr, nodes = [end_node_nr, start_node_nr], cs = self.cs_nr, tag = self.tag)
            globvar.truss_elements.append(new_truss)
            # global counter
            globvar.nelem += 1


    def check_topology_consistency(self, globvar):

        pass


            
class spiral(hoop):
    def __init__(self, globvar, XYZ, cs_nr, cs_tag, cs_mat, vec, pitch, axial_length, start = None, radius = None, div_length = None, tag = None, spin = True):
        super().__init__(globvar, XYZ, cs_nr, cs_tag, cs_mat, vec, start, None, radius, div_length, tag, connect = False)

        self.reb_type = self.rebar_type[2]
        
        # axial pitch
        self.pitch = pitch
        # length of the spiral in the axial direction
        self.axial_length = axial_length
        # positive or negative spin
        self.spin = spin
        
    def create_rebar_nodes(self):

        reinf_length = self.compute_rebar_length()
        
        nr_loops = self.axial_length / self.pitch
        self.angle = nr_loops * 2.*math.pi
        
        # number of elements of the entire spiral
        div_nr = round(reinf_length / self.div_length)

        # corresponding angle to one element
        div_angle = self.angle / div_nr

        # axial distance corresponding to one element
        div_axial = div_angle * self.pitch / (2. * math.pi)

        # vector from the center to the starting point
        start_vec = [self.start[0] - self.XYZ[0],
                     self.start[1] - self.XYZ[1],
                     self.start[2] - self.XYZ[2] ]
                                  
        # number of points = number of elements + 1
        number_of_nodes = div_nr + 1

        
        
        for i in range(number_of_nodes):

            alpha = i/div_nr * self.angle

            # change orientation of the spin to negative
            # (orientation of vec needs to be preserved because we move in that direction)
            if (self.spin == False):
                alpha *= -1.

            # start_vec shifted in the direction of the normal vector
            current_vec = [ start_vec[0] + self.vec[0] * i * div_axial,
                            start_vec[1] + self.vec[1] * i * div_axial,
                            start_vec[2] + self.vec[2] * i * div_axial ]

            rotated = np.dot( self.rotation_matrix(self.vec, alpha), current_vec )
                             
            new_point = [ origin_i + rot_i for origin_i, rot_i in zip (self.XYZ, rotated) ]

            self.rebar_points.append(new_point)


    def compute_rebar_length(self):

        nr_loops = self.axial_length / self.pitch
        loop_length = math.sqrt( (2. * math.pi * self.radius)**2 + self.pitch**2 )
        reinf_length = loop_length * nr_loops

        return reinf_length
            

    def compute_confinement_effectiveness(self):
        # formula for hoop reinforcement
        ke = ( 1.- self.pitch / (4.*self.radius) )**2.
        #return 1.0
        return ke

    def compute_effective_radius(self):
        # assuming 2nd order parabola and 45-angle intercept
        r_eff = self.radius - self.pitch/4.
        
        return r_eff

    def compute_confinement(self, globvar, fy):

        area = self.give_area(globvar)
        sig_L = 2. * area * fy / (2.*self.radius*self.pitch)
        return sig_L


    def compute_effective_confinement(self, globvar, fy):

        sig_L = self.compute_confinement(globvar, fy)
        ke = self.compute_confinement_effectiveness()
        
        sig_L_prime = sig_L * ke
        return sig_L_prime

    
    def check_topology_consistency(self, globvar):

        # first requirements on radius
        diam = self.give_diameter(globvar)
        max_radius = ( min(globvar.Bx, globvar.By) - 2.* globvar.cover - diam ) / 2.

        if (self.radius < globvar.minimum_radius):
            warnings.warn("specified radius beyond allowed bounds, d = " + str(self.radius) + " < " + str(globvar.minimum_radius) )
            self.radius = globvar.minimum_radius

        elif (self.radius > max_radius):
            warnings.warn("specified radius beyond allowed bounds, d = " + str(self.radius) + " > " + str(max_radius) )
            self.radius = max_radius
            

        # next requirements on spiral pitch
        if( self.pitch > globvar.H ):
            warnings.warn("specified topology beyond allowed bounds, pitch = " + str(self.pitch))
            self.pitch = globvar.H
        
        # final requirements on position    
        limit_x = globvar.Bx/2. - globvar.cover - diam/2  - self.radius
        limit_y = globvar.By/2. - globvar.cover - diam/2. - self.radius

        if (self.XYZ[0] < -limit_x):
            warnings.warn("specified topology beyond allowed bounds, x = " + str(self.XYZ[0]))
            self.XYZ[0] = -limit_x
        elif self.XYZ[0] > limit_x:
            warnings.warn("specified topology beyond allowed bounds, x = " + str(self.XYZ[0]))
            self.XYZ[0] = limit_x
        else:
            pass

        if self.XYZ[1] < -limit_y:
            warnings.warn("specified topology beyond allowed bounds, y = " + str(self.XYZ[1]))
            self.XYZ[1] = -limit_y
        elif self.XYZ[1] > limit_y:
            warnings.warn("specified topology beyond allowed bounds, y = " + str(self.XYZ[1]))
            self.XYZ[1] = limit_y
        else:
            pass

        if self.XYZ[2] < 0:
            warnings.warn("specified topology beyond allowed bounds, z = " + str(self.XYZ[2]))
            self.XYZ[2] = 0.
        elif self.XYZ[2] > globvar.H - self.pitch:
            warnings.warn("specified topology beyond allowed bounds, y = " + str(self.XYZ[2]))
            self.XYZ[2] = globvar.H - self.pitch
        else:
            pass

        if ( self.XYZ[2] + self.axial_length > globvar.H ):
            warnings.warn("specified topology beyond allowed bounds, axial length = " + str(self.axial_length) )
            self.axial_length = globvar.H - self.XYZ[2]
            
